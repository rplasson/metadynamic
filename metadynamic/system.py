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
    HappyEnding,
    BadEnding,
    Interrupted,
    SignalCatcher,
    FileNotFound,
    InternalError,
    OOMError,
)
from metadynamic.logger import LOGGER
from metadynamic.mpi import MPI_GATE, MPI_STATUS
from metadynamic.collector import Collectable
from metadynamic.chemical import CRN
from metadynamic.inputs import Param, LockedError
from metadynamic.inval import invalidstr
from metadynamic.json2dot import Json2dot
from metadynamic.hdf5 import ResultWriter


class Encoder(JSONEncoder):
    def default(self, obj: Any) -> Any:
        if isinstance(obj, Collectable):
            return obj.serialize()
        return super().default(obj)


class RunStatus:
    infonames = ["#", "thread", "ptime", "memuse", "step", "dstep", "time"]

    def __init__(self) -> None:
        self.tnext: float = 0.0
        self.dstep: int = 0
        self.step: int = 0
        self.time: float = 0.0
        self.tstep: float = 1.0
        self.tend: float = 0.0
        self.rtlim: float = 0.0
        self.maxmem: float = 0.0
        self.gcperio: bool = False

    def initialize(self, param: Param) -> None:
        self.tnext = 0.0
        self.dstep = 0
        self.step = 0
        self.time = 0.0
        self.tstep = param.tstep
        self.tend = param.tend
        self.rtlim = param.rtlim
        self.maxmem = param.maxmem
        self.gcperio = param.gcperio
        if self.gcperio:
            gc.disable()

    @property
    def memuse(self) -> float:
        """Return the memory used by the process in Mb"""
        return float(Process(getpid()).memory_info().rss) / 1024 / 1024

    @property
    def info(self) -> List[Union[float, int]]:
        return [
            1,
            MPI_STATUS.rank,
            LOGGER.runtime,
            self.memuse,
            self.step,
            self.dstep,
            self.time,
        ]

    def logstat(self) -> None:
        LOGGER.info(f"#{self.step}'{self.dstep}'': {self.time} -> {self.tnext}")

    def inc(self, dt: float) -> None:
        self.time += dt
        self.dstep += 1

    def next_step(self) -> None:
        if self.finished:
            self.tnext += self.tstep
        self.step += 1

    def checkend(self) -> None:
        if self.time >= self.tend:
            raise TimesUp(f"t={self.time}")
        if LOGGER.runtime >= self.rtlim:
            raise RuntimeLim(f"t={self.time}")

    def checkmem(self) -> None:
        if self.gcperio:
            gc.collect()
        if self.memuse > self.maxmem / MPI_GATE.mem_divide:
            LOGGER.warning(f"{self.memuse}>{self.maxmem/MPI_GATE.mem_divide}, call oom")
            MPI_GATE.close("oom")

    def finalize(self) -> None:
        if self.gcperio:
            gc.enable()

    @property
    def finished(self) -> bool:
        return self.time >= self.tnext


class Statistic:
    def __init__(
            self, crn: CRN, writer: ResultWriter, param: Param, status: RunStatus, comment: str
    ):
        self.crn: CRN = crn
        self.writer: ResultWriter = writer
        self.param: Param = param
        self.status: RunStatus = status
        self.statnames: List[str] = list(self.param.statparam.keys())
        self.lines: List[str] = RunStatus.infonames + self.param.save + self.statnames
        self.mapnames: List[str] = list(self.param.mapsparam.keys())
        self.mapdict: Dict[str, Dict[float, List[float]]] = {
            name: {} for name in self.mapnames
        }
        self.tsnapshot: float = (
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
            return self.crn.comp_collect.active[compound].pop / self.param.vol
        except KeyError:  # compounds doesn't exist => conc=0
            return 0.0

    @property
    def concentrations(self) -> List[float]:
        return [self.conc_of(comp) for comp in self.param.save]

    def startwriter(self) -> None:
        if self.param.hdf5 == "":
            LOGGER.error("No hdf5 filename given, can't save")
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
                self.crn.collstat(
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
        LOGGER.debug(str(res))

    def calcmap(self) -> None:
        for name, maps in self.param.mapsparam.items():
            collmap = self.crn.collmap(
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
            categories = MPI_STATUS.sortlist(list(self.mapdict[name].keys()))
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
            basename = f"{basename}-{MPI_STATUS.rank}_{timestr}"
            filename = basename + ext
            with open(filename, "w") as outfile:
                comp = self.crn.comp_collect.save()
                nbcomp = len(comp)
                reac = self.crn.reac_collect.save()
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
                    LOGGER.error(
                        f"Couldn't create snapshot {basename}.{self.param.printsnap} from {filename}"
                    )

    def writesnap(self) -> None:
        # Correct snapshot sizes
        nbsnap = MPI_STATUS.max(self._nbsnap)
        nbcomp = MPI_STATUS.max(self._nbcomp)
        nbreac = MPI_STATUS.max(self._nbreac)
        LOGGER.debug(f"resize snapshot with {nbsnap}-{nbcomp}-{nbreac}")
        self.writer.snapsize(nbcomp, nbreac, nbsnap)
        # Write snapshots
        col = 0
        for time, filename in zip(self._snaptimes, self._snapfilenames):
            with open(filename, "r") as reader:
                data = load(reader)
            self.writer.add_snapshot(
                data["Compounds"],
                data["Reactions"] if self.param.store_snapreac else {},
                col,
                time,
            )
            col += 1

    def end(self, the_end: Finished) -> None:
        self.writer.add_end(the_end, LOGGER.runtime)
        if isinstance(the_end, HappyEnding):
            LOGGER.info(str(the_end))
        elif isinstance(the_end, BadEnding):
            LOGGER.error(str(the_end))
        else:
            LOGGER.warning(str(the_end))

    def close_log(self) -> None:
        cutline = MPI_STATUS.max(self.writer.logcount[MPI_STATUS.rank])
        self.writer.close_log(cutline)

    def data_recut(self) -> None:
        maxcol = MPI_STATUS.max(self.writer.currentcol)
        self.writer.data_resize(maxcol)

    def close(self) -> None:
        LOGGER.info(f"File {self.writer.filename} to be written and closed...")
        self.data_recut()
        self.close_log()
        self.writer.close()
        LOGGER.info(f"...written and closed, done.")


class System:
    def __init__(
        self,
        filename: str,
        logfile: str = invalidstr,
        loglevel: str = "INFO",
        comment: str = "",
    ):
        LOGGER.info("Creating the system.")
        seterr(divide="ignore", invalid="ignore")
        self.initialized = False
        self.param: Param = Param.readfile(filename)
        self.writer = ResultWriter(
            self.param.hdf5, self.param.maxstrlen, self.param.lengrow
        )
        self.writer.init_log(self.param.maxlog)
        LOGGER.setlevel(loglevel)
        LOGGER.settxt(logfile)
        LOGGER.setsaver(self.writer)
        LOGGER.timeformat = self.param.timeformat
        self.signcatch = SignalCatcher()
        self.status = RunStatus()
        self.comment = comment
        MPI_GATE.init(taginit=100)
        LOGGER.info("System created")

    def initialize(self) -> None:
        if not self.initialized:
            self.crn = CRN(self.param)
            self.param.lock()
            LOGGER.info("System fully initialized")
            self.initialized = True

    def _process(self) -> None:
        # Check if a cleanup should be done
        if self.param.autoclean:
            self.crn.clean()
        # Process self.maxsteps times
        for _ in repeat(None, self.param.maxsteps):
            # ... but stop if step ends
            if self.status.finished:
                break
            # ... or if process is stopped by a signal
            if not self.signcatch.alive:
                raise Interrupted(
                    f" by {self.signcatch.signal} at t={self.status.time}"
                )
            # perform a random event
            self.status.inc(self.crn.stepping())
        else:
            #  had to performed loop until maxsteps
            LOGGER.warning(f"maxsteps per process (={self.param.maxsteps}) too low")

    def run(self) -> str:
        LOGGER.reset_timer()
        if MPI_STATUS.ismpi:
            LOGGER.info(f"Launching MPI run from thread #{MPI_STATUS.rank}")
        else:
            LOGGER.info("Launching single run.")
        # Setup working environment
        self.signcatch.listen()
        self.status.initialize(self.param)
        self.initialize()
        statistic = Statistic(self.crn, self.writer, self.param, self.status, self.comment)
        statistic.startwriter()
        # Ready, let's run
        LOGGER.info(f"Run #{MPI_STATUS.rank}={getpid()} launched")
        with MPI_GATE.launch():
            # Do not stop except if asked to exit the gate
            while MPI_GATE.cont:
                try:
                    # Log stat info
                    self.status.logstat()
                    # perform a step
                    self._process()
                    # Write data statistics
                    statistic.writestat()
                    statistic.calcmap()
                    statistic.calcsnapshot()
                    # Check memory
                    self.status.checkmem()
                    self.status.next_step()
                    self.status.checkend()
                    # Pass gate checkpoint
                    MPI_GATE.checkpoint()
                # exit the gate because the work is finished.
                except Finished as the_end:
                    end = f"{the_end} ({LOGGER.runtime} s)"
                    statistic.end(the_end)
                    break
            # exit the gate because asked by collective decision.
            else:
                gate_end = OOMError("Stopped by kind request of OOM killer")
                end = f"{gate_end} ({LOGGER.runtime} s)"
                statistic.end(gate_end)
        # Write and log final data
        LOGGER.info(f"Run #{MPI_STATUS.rank}={getpid()} finished")
        statistic.calcsnapshot(final=True)
        statistic.writesnap()
        statistic.writemap()
        # Then cleanly exit
        statistic.close()
        self.status.finalize()
        self.signcatch.reset()
        return end

    def set_param(self, **kwd) -> None:
        try:
            self.param.set_param(**kwd)
        except LockedError:
            raise InternalError("Parameter set after initialization")
        LOGGER.info(f"Parameters changed: {self.param}")
