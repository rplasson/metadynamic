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
from math import ceil
from multiprocessing import get_context
from itertools import repeat
from os import getpid
from typing import Dict, Tuple, Any
from psutil import Process
from mpi4py import MPI

from pandas import DataFrame
from json import dump, JSONEncoder

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
)
from metadynamic.logger import Logged
from metadynamic.proba import Probalistic
from metadynamic.processing import Result
from metadynamic.ruleset import Ruled
from metadynamic.collector import Collectable
from metadynamic.chemical import Collected, trigger_changes
from metadynamic.inputs import Param
from metadynamic.inval import invalidstr, invalidfloat, isvalid
from metadynamic.json2dot import Json2dot
from metadynamic.hdf5 import ResultWriter


class Encoder(JSONEncoder):
    def default(self, obj: Any) -> Any:
        if isinstance(obj, Collectable):
            return obj.serialize()
        return super().default(obj)


class System(Probalistic, Collected):
    def __init__(
        self, filename: str, logfile: str = invalidstr, loglevel: str = "INFO"
    ):
        self.initialized = False
        Logged.setlogger(logfile, loglevel)
        self.param: Param = Param.readfile(filename)
        self.log.info("Parameter files loaded.")
        self.signcatch = SignalCatcher()

    @property
    def memuse(self) -> float:
        return float(Process(getpid()).memory_info().rss) / 1024 / 1024

    def conc_of(self, compound: str) -> float:
        try:
            return self.comp_collect.active[compound].pop / self.param.vol
        except KeyError:  # compounds doesn't exist => conc=0
            return 0.0

    @property
    def poplist(self) -> Dict[str, int]:
        return {
            comp.description: comp.pop for comp in self.comp_collect.active.values()
        }

    @property
    def lendist(self) -> Dict[int, int]:
        res = self.comp_collect.dist("length", lenweight=True, full=False)
        res[-1] = 1
        return res

    @property
    def pooldist(self) -> Dict[int, int]:
        res = self.comp_collect.dist("length", lenweight=False, full=True)
        res[-1] = 1
        return res

    def statlist(self) -> Tuple[Dict[str, int], Dict[str, Dict[int, int]]]:
        stat = {}
        dist = {}
        dist["lendist"] = self.lendist
        dist["pooldist"] = self.pooldist
        stat["nbcomp"] = len(self.comp_collect.active)
        stat["nbreac"] = len(self.reac_collect.active)
        stat["poolsize"] = len(self.comp_collect.pool)
        stat["poolreac"] = len(self.reac_collect.pool)
        stat["maxlength"] = max(dist["lendist"])
        return stat, dist

    def initialize(self) -> None:
        if self.initialized:
            raise InitError("Double Initialization")
        if self.param.gcperio:
            gc.disable()
        self.time = 0.0
        # If seed set, always restart from the same seed. For timing/debugging purpose
        self.step = 0
        Probalistic.setprobalist(vol=self.param.vol)
        self.probalist.seed(self.param.seed)
        # Add all options for collections
        Collected.setcollections(dropmode_reac=self.param.dropmode)
        Ruled.setrules(self.param.rulemodel, self.param.consts)
        self.log.info("System created")
        for compound, pop in self.param.init.items():
            self.comp_collect[compound].change_pop(pop)
        trigger_changes()
        self.log.info(f"Initialized with {self.param}")
        self.initialized = True

    def _process(self, tstop: float) -> bool:
        # Check if a cleanup should be done
        if self.param.autoclean:
            self.probalist.clean()
        # Check if end of time is nigh
        if self.time >= self.param.tend:
            raise TimesUp(f"t={self.time}")
        if self.log.runtime() >= self.param.rtlim:
            raise RuntimeLim(f"t={self.time}")
        # Then process self.maxsteps times
        for _ in repeat(None, self.param.maxsteps):
            # ... but stop is step ends
            if self.time >= tstop:
                return True
            # choose a random event
            chosen, dt = self.probalist.choose()
            # check if there even was an event to choose
            if chosen is None:
                raise NotFound(f"t={self.time}")
            # perform the (chosen one) event
            chosen.process()
            if self.probalist.probtot == 0:
                raise NoMore(f"t={self.time}")
            if not self.signcatch.alive:
                raise Interrupted(f" by {self.signcatch.signal} at t={self.time}")
            # update time for next step
            self.time += dt
            self.step += 1
        self.log.warning(f"maxsteps per process (={self.param.maxsteps}) too low")
        return False

    def _run(self, num: int = -1, ismpi=False) -> Tuple[DataFrame, int, int, str]:
        lines = (
            ["#", "thread", "ptime", "memuse", "step", "time"]
            + self.param.save
            + ["maxlength", "nbcomp", "poolsize", "nbreac", "poolreac"]
        )
        if num >= 0:
            self.log.connect(f"Reconnected from thread {num+1}", num + 1)
        if not self.initialized:
            self.log.info("Will initialize")
            self.initialize()
        if ismpi:
            writer = ResultWriter(
                self.param.hdf5, lines, ceil(self.param.tend / self.param.tstep) + 1
            )
            writer.add_parameter(self.paramasdict())
        table = DataFrame(index=lines)
        lendist = DataFrame()
        pooldist = DataFrame()
        tnext = 0.0
        tsnapshot = self.param.sstep if self.param.sstep >= 0 else 2 * self.param.tend
        step = 0
        self.log.info(f"Run {num}={getpid()} launched")
        self.signcatch.listen()
        # Process(getpid()).cpu_affinity([num % cpu_count()])
        while True:
            try:
                self.log.info(f"#{step}: {self.time} -> {tnext}")
                finished = self._process(tnext)
                stat, dist = self.statlist()
                res = (
                    [1, num, self.log.runtime(), self.memuse, self.step, self.time]
                    + [self.conc_of(comp) for comp in self.param.save]
                    + [
                        stat["maxlength"],
                        stat["nbcomp"],
                        stat["poolsize"],
                        stat["nbreac"],
                        stat["poolreac"],
                    ]
                )
                table[step] = res
                if ismpi:
                    writer.add_data(res)
                self.log.debug(str(res))
                lendist = lendist.join(
                    DataFrame.from_dict({step: dist["lendist"]}), how="outer",
                ).fillna(0)
                pooldist = pooldist.join(
                    DataFrame.from_dict({step: dist["pooldist"]}), how="outer",
                ).fillna(0)
                if self.param.gcperio:
                    gc.collect()
                if finished:
                    if tnext > tsnapshot:
                        self.make_snapshot(num, tsnapshot)
                        tsnapshot += self.param.sstep
                    tnext += self.param.tstep
                step += 1
            except Finished as the_end:
                end = f"{the_end} ({self.log.runtime()} s)"
                if isinstance(the_end, HappyEnding):
                    self.log.info(str(the_end))
                elif isinstance(the_end, BadEnding):
                    self.log.error(str(the_end))
                else:
                    self.log.warning(str(the_end))
                break
        retval: Tuple[DataFrame, int, int, str] = (
            table,
            lendist.astype(int),
            pooldist.astype(int),
            end,
        )
        self.log.debug(f"Run {num}={getpid()} finished")
        self.make_snapshot(num)
        if num >= 0:
            # Clean memory as much as possible to leave room to still alive threads
            self.comp_collect.purge()
            self.reac_collect.purge()
            Probalistic.setprobalist(vol=self.param.vol)
            trigger_changes()
            gc.collect()
            self.log.debug(f"Collection purged for {num}")
            self.log.disconnect(f"Disconnected from thread {num}")
        if ismpi:
            writer.close()
        return retval

    def make_snapshot(self, num: int, time: float = invalidfloat) -> None:
        timestr = str(time) if isvalid(time) else "end"
        if self.param.snapshot:
            basename, ext = self.param.snapshot.split(".")
            filled = "{base}-{n}_{t}" if num >= 0 else "{base}-{t}"
            basename = filled.format(base=basename, n=num, t=timestr)
            filename = f"{basename}.{ext}"
            with open(filename, "w") as outfile:
                dump(
                    {
                        "Compounds": self.comp_collect.save(),
                        "Reactions": self.reac_collect.save(),
                    },
                    outfile,
                    cls=Encoder,
                    indent=4,
                )
            if self.param.printsnap:
                Json2dot(filename).write(f"{basename}.{self.param.printsnap}")

    def run(self) -> Result:
        if MPI.COMM_WORLD.size > 1:
            rank = MPI.COMM_WORLD.rank
            self.log.info(f"Launching MPI run from thread #{rank}")
            self.log.disconnect(reason="Launching MPI....")
            res = [self._run(rank, ismpi=True)]
            self.signcatch.reset()
            return Result(res)
        if self.param.nbthread == 1:
            self.log.info("Launching single run.")
            res = [self._run()]
            self.signcatch.reset()
            return Result(res)
        self.log.info("Launching threaded run.")
        return self.multirun(self.param.nbthread)

    def multirun(self, nbthread: int = -1) -> Result:
        ctx = get_context(self.param.context)
        self.signcatch.ignore()
        if nbthread == -1:
            nbthread = ctx.cpu_count()
        self.log.info(f"Launching multithread from {getpid()}")
        self.log.disconnect(reason="Launching multithreading....")
        with ctx.Pool(nbthread) as pool:
            running = pool.map_async(self._run, range(nbthread))
            self.log.connect(reason="...Multithreading launched")
            res = running.get()
            self.log.info("Multirun finished!")
            self.signcatch.reset()
        return Result(res)

    def set_param(self, **kwd) -> None:
        if self.initialized:
            raise InitError("Parameter set after initialization")
        self.param.set_param(**kwd)
        self.log.info(f"Parameters changed: {self.param}")
