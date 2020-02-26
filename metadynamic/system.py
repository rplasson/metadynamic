#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2019 by Raphaël Plasson
#
# This file is part of metadynamic
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.

import gc
from os import path
from math import ceil
from numpy import nan, seterr
from itertools import repeat
from os import getpid
from typing import Dict, Any, List, Union
from psutil import Process

from json import load, dump, JSONEncoder

from metadynamic.ends import (
    Finished,
    TimesUp,
    RuntimeLim,
    NotFound,
    NoMore,
    HappyEnding,
    BadEnding,
    InitError,
    Interrupted,
    SignalCatcher,
    FileNotFound,
    InternalError,
)
from metadynamic.logger import Logged
from metadynamic.proba import Probalistic
from metadynamic.ruleset import Ruled
from metadynamic.collector import Collectable
from metadynamic.chemical import Collected, trigger_changes
from metadynamic.inputs import Param, LockedError
from metadynamic.inval import invalidstr
from metadynamic.json2dot import Json2dot
from metadynamic.hdf5 import Saver
from metadynamic.mpi import Parallel


class Encoder(JSONEncoder):
    def default(self, obj: Any) -> Any:
        if isinstance(obj, Collectable):
            return obj.serialize()
        return super().default(obj)


class RunStatus(Logged):
    def __init__(self, param: Param):
        self.param = param

    def initialize(self, num: int) -> None:
        self.num: int = num
        self.tnext: float = 0.0
        self.dstep: int = 0
        self.step: int = 0
        self.time: float = 0.0
        self.infonames = ["#", "thread", "ptime", "memuse", "step", "dstep", "time"]

    @property
    def memuse(self) -> float:
        """Return the memory used by the process in Mb"""
        return float(Process(getpid()).memory_info().rss) / 1024 / 1024

    @property
    def info(self) -> List[Union[float, int]]:
        return [
            1,
            self.num,
            self.log.runtime(),
            self.memuse,
            self.step,
            self.dstep,
            self.time,
        ]

    def logstat(self) -> None:
        self.log.info(f"#{self.step}'{self.dstep}'': {self.time} -> {self.tnext}")

    def inc(self, dt: float) -> None:
        self.time += dt
        self.dstep += 1

    def next_step(self) -> None:
        if self.finished:
            self.tnext += self.param.tstep
        self.step += 1

    def checkend(self) -> None:
        if self.time >= self.param.tend:
            raise TimesUp(f"t={self.time}")
        if self.log.runtime() >= self.param.rtlim:
            raise RuntimeLim(f"t={self.time}")

    @property
    def finished(self) -> bool:
        return self.time >= self.tnext


class Statistic(Collected, Saver, Parallel):
    def __init__(self, param: Param, status: RunStatus, comment: str):
        self.param = param
        self.status = status
        self.statnames = list(self.param.statparam.keys())
        self.lines = self.status.infonames + self.param.save + self.statnames
        self.mapnames = list(self.param.mapsparam.keys())
        self.mapdict: Dict[str, Dict[float, List[float]]] = {
            name: {} for name in self.mapnames
        }
        self.tsnapshot = (
            self.param.sstep if self.param.sstep >= 0 else 2 * self.param.tend
        )
        self._snapfilenames: List[str] = []
        self._snaptimes: List[float] = []
        self._nbcomp: int = 0
        self._nbreac: int = 0
        self._nbsnap: int = 0
        self.comment = comment

    def conc_of(self, compound: str) -> float:
        try:
            return self.comp_collect.active[compound].pop / self.param.vol
        except KeyError:  # compounds doesn't exist => conc=0
            return 0.0

    @property
    def concentrations(self) -> List[float]:
        return [self.conc_of(comp) for comp in self.param.save]

    def startwriter(self) -> None:
        if self.param.hdf5 == "":
            self.log.error("No hdf5 filename given, can't save")
            raise FileNotFound("No hdf5 filename given, can't save")
        else:
            self.writer.init_stat(
                datanames=self.lines,
                mapnames=self.mapnames,
                params=self.param,
                comment=self.comment,
                nbcol=ceil(self.param.tend / self.param.tstep) + 1,
            )

    def writestat(self) -> None:
        res = (
            self.status.info
            + self.concentrations
            + [
                self.collstat(
                    collection=stats.collection,
                    prop=stats.prop,
                    weight=stats.weight,
                    method=stats.method,
                    full=stats.full,
                )
                for stats in self.param.statparam.values()
            ]
        )
        self.writer.add_data(res)
        self.log.debug(str(res))

    def calcmap(self) -> None:
        for name, maps in self.param.mapsparam.items():
            collmap = self.collmap(
                collection=maps.collection,
                prop=maps.prop,
                weight=maps.weight,
                sort=maps.sort,
                method=maps.method,
                full=maps.full,
            )
            for cat in collmap.keys() | self.mapdict[name].keys():
                if cat in self.mapdict[name]:
                    if cat in collmap:
                        self.mapdict[name][cat].append(collmap[cat])
                    else:
                        self.mapdict[name][cat].append(nan)
                else:
                    self.mapdict[name][cat] = [nan] * self.status.step + [collmap[cat]]

    def writemap(self) -> None:
        for name in self.mapnames:
            # correct sizes
            categories = self.mpi.sortlist(list(self.mapdict[name].keys()))
            self.writer.mapsize(name, categories)
            # write maps
            self.writer.add_map(name, self.mapdict[name])

    def calcsnapshot(self, final: bool = False) -> None:
        if self.param.snapshot:
            if final:
                timestr = "end"
            else:
                if self.status.tnext < self.tsnapshot:
                    return None
                timestr = str(self.tsnapshot)
                self.tsnapshot += self.param.sstep
            basename, ext = path.splitext(self.param.snapshot)
            filled = "{base}-{n}_{t}" if self.status.num >= 0 else "{base}-{t}"
            basename = filled.format(base=basename, n=self.status.num, t=timestr)
            filename = basename + ext
            with open(filename, "w") as outfile:
                comp = self.comp_collect.save()
                nbcomp = len(comp)
                reac = self.reac_collect.save()
                nbreac = len(reac)
                dump(
                    {"Compounds": comp, "Reactions": reac},
                    outfile,
                    cls=Encoder,
                    indent=4,
                )
            self._snapfilenames.append(filename)
            self._snaptimes.append(self.param.tend if final else self.tsnapshot)
            self._nbsnap += 1
            self._nbcomp = max(self._nbcomp, nbcomp)
            self._nbreac = max(self._nbreac, nbreac)
            if self.param.printsnap:
                if not Json2dot(filename).write(f"{basename}.{self.param.printsnap}"):
                    self.log.error(
                        f"Couldn't create snapshot {basename}.{self.param.printsnap} from {filename}"
                    )

    def writesnap(self) -> None:
        # Correct snapshot sizes
        nbsnap = self.mpi.max(self._nbsnap)
        nbcomp = self.mpi.max(self._nbcomp)
        nbreac = self.mpi.max(self._nbreac)
        self.log.debug(f"resize snapshot with {nbsnap}-{nbcomp}-{nbreac}")
        self.writer.snapsize(nbcomp, nbreac, nbsnap)
        # Write snapshots
        col = 0
        for time, filename in zip(self._snaptimes, self._snapfilenames):
            with open(filename, "r") as reader:
                data = load(reader)
            self.writer.add_snapshot(data["Compounds"], data["Reactions"]  if self.param.store_snapreac else {}, col, time)
            col += 1

    def end(self, the_end: Finished) -> None:
        self.writer.add_end(the_end, self.log.runtime())
        if isinstance(the_end, HappyEnding):
            self.log.info(str(the_end))
        elif isinstance(the_end, BadEnding):
            self.log.error(str(the_end))
        else:
            self.log.warning(str(the_end))

    def close_log(self) -> None:
        cutline = self.mpi.max(self.writer.logcount[self.mpi.rank])
        self.writer.close_log(cutline)

    def data_recut(self) -> None:
        maxcol = self.mpi.max(self.writer.currentcol)
        self.writer.data_resize(maxcol)

    def close(self) -> None:
        self.log.info(f"File {self.writer.filename} to be written and closed...")
        self.data_recut()
        self.close_log()
        self.writer.close()
        self.log.info(f"...written and closed, done.")


class System(Probalistic, Collected, Saver, Parallel):
    def __init__(
        self,
        filename: str,
        logfile: str = invalidstr,
        loglevel: str = "INFO",
        comment: str = "",
    ):
        seterr(divide="ignore", invalid="ignore")
        self.initialized = False
        self.param: Param = Param.readfile(filename)
        Parallel.setmpi(taginit=100)
        Saver.setsaver(self.param.hdf5, self.param.maxstrlen, self.param.lengrow)
        Logged.setlogger(logfile, loglevel)
        self.writer.init_log(self.param.maxlog)
        self.log.info("Parameter files loaded.")
        self.signcatch = SignalCatcher()
        self.status = RunStatus(self.param)
        self.comment = comment

    def initialize(self) -> None:
        if self.initialized:
            raise InitError("Double Initialization")
        if self.param.gcperio:
            gc.disable()
        Probalistic.setprobalist(vol=self.param.vol)
        # Add all options for collections
        Collected.setcollections(dropmode_reac=self.param.dropmode)
        Ruled.setrules(self.param.rulemodel, self.param.consts)
        self.log.info("System created")
        for compound, pop in self.param.init.items():
            self.comp_collect[compound].change_pop(pop)
        trigger_changes()
        self.log.info(f"Initialized with {self.param}")
        self.initialized = True
        self.param.lock()

    def _process(self) -> None:
        # Check if a cleanup should be done
        if self.param.autoclean:
            self.probalist.clean()
        # Then process self.maxsteps times
        for _ in repeat(None, self.param.maxsteps):
            # ... but stop is step ends
            if self.status.finished:
                break
            # choose a random event
            chosen, dt = self.probalist.choose()
            # check if there even was an event to choose
            if chosen is None:
                raise NotFound(f"t={self.status.time}")
            # perform the (chosen one) event
            chosen.process()
            if self.probalist.probtot == 0:
                raise NoMore(f"t={self.status.time}")
            if not self.signcatch.alive:
                raise Interrupted(
                    f" by {self.signcatch.signal} at t={self.status.time}"
                )
            # update time for next step
            self.status.inc(dt)
        if not self.status.finished:
            self.log.warning(f"maxsteps per process (={self.param.maxsteps}) too low")

    def _run(self, num: int = -1) -> str:
        self.status.initialize(num)
        statistic = Statistic(self.param, self.status, self.comment)
        if num >= 0:
            self.log.connect(f"Reconnected from thread {num+1}", num + 1)
        if not self.initialized:
            self.log.info("Will initialize")
            self.initialize()
        statistic.startwriter()
        self.log.info(f"Run #{num}={getpid()} launched")
        self.signcatch.listen()
        # Process(getpid()).cpu_affinity([num % cpu_count()])
        with self.gate.context() as gate:
            while gate.cont.ok:
                try:
                    self.status.logstat()
                    self._process()
                    if self.status.memuse > self.param.maxmem/gate.mem_divide:
                        self.log.warning(f"{self.status.memuse}>{self.param.maxmem/gate.mem_divide}, call oom")
                        gate.close("oom")
                    if self.param.gcperio:
                        gc.collect()
                    statistic.writestat()
                    statistic.calcmap()
                    statistic.calcsnapshot()
                    self.status.next_step()
                    self.status.checkend()
                    gate.checkpoint()
                except Finished as the_end:
                    end = f"{the_end} ({self.log.runtime()} s)"
                    statistic.end(the_end)
                    break
            else:
                gate_end = Interrupted("by kind request of OOM killer")
                end = f"{gate_end} ({self.log.runtime()} s)"
                statistic.end(gate_end)
        self.log.info(f"Run #{num}={getpid()} finished")
        statistic.calcsnapshot(final=True)
        statistic.writesnap()
        statistic.writemap()
        statistic.close()
        if num >= 0:
            self.log.disconnect(f"Disconnected from #{num}")
        return end

    def purge(self, num: int) -> None:
        self.comp_collect.purge()
        self.reac_collect.purge()
        Probalistic.setprobalist(vol=self.param.vol)
        trigger_changes()
        gc.collect()
        self.log.debug(f"Collection purged for #{num}")

    def run(self) -> List[str]:
        if self.mpi.ismpi:
            self.log.info(f"Launching MPI run from thread #{self.mpi.rank}")
            self.log.disconnect(reason="Launching MPI....")
            res = [self._run(self.mpi.rank)]
            self.signcatch.reset()
            return res
        self.log.info("Launching single run.")
        res = [self._run()]
        self.signcatch.reset()
        return res

    def set_param(self, **kwd) -> None:
        try:
            self.param.set_param(**kwd)
        except LockedError:
            raise InternalError("Parameter set after initialization")
        self.log.info(f"Parameters changed: {self.param}")
