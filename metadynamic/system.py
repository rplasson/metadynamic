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


"""Simulation interface.

It provides L{RunStatus}, a class for controlling the status of a simulation (time, step, memory),
L{Statistic}: a class for computing, getting and saving statistics on a simulation and L{System}, a
class for creating, running, and directly controlling a simulation.

"""


import gc

from math import ceil
from itertools import repeat
from os import getpid
from typing import Dict, Any, List, Union
from psutil import Process

import numpy as np


from metadynamic.ends import (
    Finished,
    TimesUp,
    RuntimeLim,
    HappyEnding,
    BadEnding,
    Interrupted,
    SignalCatcher,
    InternalError,
    OOMError,
)
from metadynamic.logger import LOGGER
from metadynamic.mpi import MPI_GATE, MPI_STATUS
from metadynamic.outputs import Output
from metadynamic.chemical import Crn
from metadynamic.inputs import Param, LockedError
from metadynamic.hdf5 import ResultWriter
from metadynamic.result import ResultReader
from metadynamic.inval import isvalid, invalidint


class RunStatus:
    """Control the status of a simulation (time, step, memory)."""

    infonames = ["thread", "ptime", "memuse", "step", "dstep", "time"]
    """List of tracked information"""

    def __init__(self) -> None:
        """Create an empty object, to be initialized later."""
        self.tnext: float = 0.0
        """next time to be reached"""

        self.dstep: int = 0
        """stochatic step number"""

        self.step: int = 0
        """timestep number"""

        self.time: float = 0.0
        """simulated time"""

        self.tstep: float = 1.0
        """timestep"""

        self.tend: float = 0.0
        """final simulated time"""

        self.rtlim: float = 0.0
        """maximum run time"""

        self.maxmem: float = 0.0
        """maximum memory use"""

        self.gcperio: bool = False
        """periodic garbage collection flag
        (True: only GC after each time step,
        False: regular GC process) """

    def initialize(self, param: Param) -> None:
        """First initialization of the run from the set of parameters 'param'.

        @param param: parameters
        @type param: Param

        """
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
        """Memory used by the process in Mb."""
        return float(Process(getpid()).memory_info().rss) / 1024 / 1024

    @property
    def info(self) -> List[Union[float, int]]:
        """Information summary.

        thread number, runtime, memory use, simulation step, stochastic step, time

        """
        return [
            MPI_STATUS.rank,
            LOGGER.runtime,
            self.memuse,
            self.step,
            self.dstep,
            self.time,
        ]

    def logstat(self) -> None:
        """Log (at INFO level) run information."""
        LOGGER.info(
            f"#{self.step}'{self.dstep}'': {self.time} -> {self.tnext}  <{int(self.memuse)}Mb>"
        )

    def inc(self, delta_t: float) -> None:
        """
        Increment time by 'delta_t' and dstep by 1.

        @param delta_t: time increment
        @type delta_t: float
        """
        self.time += delta_t
        self.dstep += 1

    def next_step(self) -> None:
        """Increment tnext by tstep if tnext was reached, and step by 1."""
        if self.finished:
            self.tnext += self.tstep
        self.step += 1

    def checkend(self) -> None:
        """
        Check if the run reached its final limit.

        @raise TimesUp: if time reached tend
        @raise RuntimeLim: if runtime reached rtlim
        """
        if self.time >= self.tend:
            raise TimesUp(f"t={self.time}")
        if LOGGER.runtime >= self.rtlim:
            raise RuntimeLim(f"t={self.time}")

    def checkmem(self) -> None:
        """
        Check if the process reached its memory limit.

        If reached, close the MPI gate for 'oom' reason
        for selectively stopping some threads at the next
        gathering at the MPI gate.
        """
        if self.gcperio:
            gc.collect()
        if self.memuse > self.maxmem / MPI_GATE.mem_divide:
            LOGGER.warning(f"{self.memuse}>{self.maxmem/MPI_GATE.mem_divide}, call oom")
            MPI_GATE.close("oom")

    def checkstep(self) -> None:
        """
        Follow a checklist before starting next step.

         1. check the memory,
         2. update tnext and step,
         3. check if run ends.
        """
        self.checkmem()
        self.next_step()
        self.checkend()

    def gcclean(self) -> None:
        """Call garbage collector for memory cleaning."""
        gc.collect()
        if self.gcperio:
            gc.enable()

    @property
    def finished(self) -> bool:
        """Return  True if present time reached the tnext limit."""
        return self.time >= self.tnext


class Statistic:
    """Compute, get and save statustics on a simulation run."""

    def __init__(
        self,
        crn: Crn,
        writer: ResultWriter,
        param: Param,
        status: RunStatus,
        comment: str,
    ):
        """
        Create the statistic object.

        @param crn: chemical reaction network of the simulation
        @type crn: Crn
        @param writer: writer to the .hdf5 result file
        @type writer: ResultWriter
        @param param: run parameters
        @type param: Param
        @param status: run status
        @type status: RunStatus
        @param comment: Comment string describing the run
        @type comment: str
        """
        self.crn: Crn = crn
        """Chemical Reaction Network"""
        self.writer: ResultWriter = writer
        """Writer to HDF5 file"""
        self.param: Param = param
        """Run  parameters"""
        self.status: RunStatus = status
        """Run status"""
        self.statnames: List[str] = list(self.param.statparam.keys())
        """list of statistics names"""
        self.lines: List[str] = RunStatus.infonames + self.param.save + self.statnames
        """list of fields to be saved in HDF5 file"""
        self.mapnames: List[str] = list(self.param.mapsparam.keys())
        """list of statistics map names"""
        self._mapdict: Dict[str, Dict[float, List[float]]] = {
            name: {} for name in self.mapnames
        }
        self.tsnapshot: float = self.param.sstep if self.param.sstep >= 0 else 2 * self.param.tend
        """time of next snapshot"""
        self._snaptimes: List[float] = []
        self._snapcomp: List[Dict[str, Any]] = []
        self._snapreac: List[Dict[str, Any]] = []
        self._nbcomp: int = 0
        self._nbreac: int = 0
        self._nbsnap: int = 0
        self.comment = comment
        """Run comment"""

    def conc_of(self, compound: str) -> float:
        """Get the concentration of a given compound.

        @param compound: compound name
        @type compound: str
        @return: concentration value
        @rtype: float
        """
        try:
            return self.crn.comp_collect.active[compound].pop / self.param.vol
        except KeyError:  # compounds doesn't exist => conc=0
            return 0.0

    @property
    def concentrations(self) -> List[float]:
        """List of concentrations of compounds to be saved."""
        return [self.conc_of(comp) for comp in self.param.save]

    def startwriter(self) -> None:
        """Initialize the hdf5 writer."""
        self.writer.init_stat(
            datanames=self.lines,
            mapnames=self.mapnames,
            params=self.param,
            ruleparam=self.crn.model.param.asdict(),
            comment=self.comment,
            nbcol=ceil(self.param.tend / self.param.tstep) + 1,
        )

    def writestat(self) -> None:
        """Write general statistics at present time to hdf5."""
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
        """Calculate map statistics at present time, saved for later writing."""
        for name, maps in self.param.mapsparam.items():
            collmap = self.crn.collmap(
                collection=maps.collection,
                prop=maps.prop,
                weight=maps.weight,
                sort=maps.sort,
                method=maps.method,
                full=maps.full,
            )
            for cat in collmap.keys() | self._mapdict[name].keys():
                if cat in self._mapdict[name]:
                    if cat in collmap:
                        self._mapdict[name][cat].append(collmap[cat])
                    else:
                        self._mapdict[name][cat].append(np.nan)
                else:
                    self._mapdict[name][cat] = [np.nan] * self.status.step + [
                        collmap[cat]
                    ]

    def writemap(self) -> None:
        """Write all previously calculated map statistics to hdf5."""
        for name in self.mapnames:
            # write maps
            categories = MPI_STATUS.sortlist(list(self._mapdict[name].keys()))
            self.writer.add_map(name, categories, self._mapdict[name])

    def calcsnapshot(self, final: bool = False) -> None:
        """Snapshot the CRN at present time, if time reached tsnapshot.

        Also saves if it is the final snapshot.

        Saved for later writing in hdf5.

        @param final: True if this is the last snapshot
        @type final: bool
        """
        if not final:
            if self.status.tnext < self.tsnapshot:
                return None
            self.tsnapshot += self.param.sstep
        comp = self.crn.comp_collect.asdict()
        nbcomp = len(comp)
        reac = self.crn.reac_collect.asdict() if self.param.store_snapreac else {}
        nbreac = len(reac)
        self._snapcomp.append(comp)
        self._snapreac.append(reac)
        self._snaptimes.append(self.param.tend if final else self.tsnapshot)
        self._nbsnap += 1
        self._nbcomp = max(self._nbcomp, nbcomp)
        self._nbreac = max(self._nbreac, nbreac)
        return None

    def writesnap(self) -> None:
        """Write all previous snapshots to hdf5."""
        # Correct snapshot sizes
        nbsnap = MPI_STATUS.max(self._nbsnap)
        nbcomp = MPI_STATUS.max(self._nbcomp)
        nbreac = MPI_STATUS.max(self._nbreac)
        LOGGER.debug(f"resize snapshot with {nbsnap}-{nbcomp}-{nbreac}")
        self.writer.snapsize(nbcomp, nbreac, nbsnap)
        # Write snapshots
        for col, (time, comp, reac) in enumerate(
            zip(self._snaptimes, self._snapcomp, self._snapreac)
        ):
            self.writer.add_snapshot(comp, reac, col, time)

    def calc(self) -> None:
        """Calculate all statistics: stat, map and snapshots."""
        self.writestat()
        self.calcmap()
        self.calcsnapshot()

    def end(self, the_end: Finished) -> None:
        """Write ending message to hdf5.

        @param the_end: ending message
        @type the_end: Finished
        """
        self.writer.add_end(the_end, LOGGER.runtime)
        if isinstance(the_end, HappyEnding):
            LOGGER.info(str(the_end))
        elif isinstance(the_end, BadEnding):
            LOGGER.error(str(the_end))
        else:
            LOGGER.warning(str(the_end))

    def close(self) -> None:
        """Write all remaining data  in hdf5 file, clean it, then close it."""
        LOGGER.info(f"File {self.writer.filename} to be written and closed...")
        self.writesnap()
        self.writemap()
        self.writer.close()
        LOGGER.info(f"...written and closed, done.")


class System:
    """Create, run and control a simulation."""

    @classmethod
    def fromjson(cls, filename: str, **kwd: Any) -> "System":
        """
        Create a System object from a Param json file.

        Additional parameters can be passed,
        they will override the one defined in the json file.

        @param filename: Name of the Param json file
        @type filename: str
        @param kwd: Additional parameters
        @return: newly created System object
        @rtype: System
        """
        param = Param.readfile(filename)
        param.set_param(**kwd)
        return cls(param)

    @classmethod
    def fromhdf5(
        cls, filename: str, snapnum: int = invalidint, snapstep: int = -1, **kwd: Any,
    ) -> "System":
        """
        Create a System object from a result hdf5 file.

        Additional parameters can be passed,
        they will override the one defined in the hdf5 file.

        if a 'snapnum' is provided, initial conditions will be set from the recorded snapshot as
        saved by the corresponding thread number at the corresponding 'snapstep' (defaulting to the
        final snapshot)

        @param filename: Name of the Param json file
        @type filename: str
        @param snapnum: snapshot number to start from (Default value = invalindint)
        @type snapnum: int
        @param snapstep: snapshot step to start from (Default value = -1)
        @type snapstep: int
        @param kwd: Additional parameters
        @return: newly created System object
        @rtype: System
        """
        res = ResultReader(filename)
        param = res.parameters
        param.set_param(**kwd)
        if isvalid(snapnum):
            if snapnum < 0:
                snapnum = MPI_STATUS.rank % res.size
            param.set_param(init=res.getsnap_comp(snapnum, snapstep))
        return cls(param)

    def __init__(
        self, param: Param,
    ):
        """Create a System from parameters 'param'.

        @param param: parameters
        @type param: Param
        """
        LOGGER.debug("Creating the system.")
        np.seterr(divide="ignore", invalid="ignore")
        self.initialized = False
        """Initialized flag"""
        self.param: Param = param
        """run parameters"""
        MPI_STATUS.init(self.param.timeformat)
        self.output: Output = Output(self.param)
        """Access to output files and folders"""
        self.writer: ResultWriter = ResultWriter(
            self.output.h5file, self.param.maxstrlen, self.param.lengrow
        )
        """Writer to HDF5 file."""
        self.writer.init_log(self.param.maxlog)
        LOGGER.setlevel(self.param.loglevel)
        LOGGER.connect(self.output.logfile)
        LOGGER.setsaver(self.writer)
        LOGGER.timeformat = self.param.timeformat
        self.signcatch: SignalCatcher = SignalCatcher()
        """Signal catcher for dealing with system interruptions"""
        self.status: RunStatus = RunStatus()
        """Interface to run status"""
        self.comment = self.param.comment
        """run comment"""
        MPI_GATE.init(taginit=100, sleeptime=param.sleeptime)
        LOGGER.info("System created")
        self.crn: Crn
        """Chemical Reaction Network"""
        self.statistic: Statistic
        """interface to run statistics"""

    def _initialize(self) -> None:
        """Initialize the system before starting the run.

        Intended to be launched from self.run, not manually

        @raise InternalError: in case of double initialization.

        """
        if not self.initialized:
            self.signcatch.listen()
            self.status.initialize(self.param)
            self.crn = Crn(self.param)
            self.param.lock()
            self.statistic = Statistic(
                self.crn, self.writer, self.param, self.status, self.comment
            )
            self.statistic.startwriter()
            self.initialized = True
            LOGGER.info("System fully initialized")
        else:
            raise InternalError("Attempt to initialize already initialized system!")

    def _release(self, end: Finished, fullrelease: bool = True) -> str:
        """
        Clean data before run ending.

        A double release will raise an InternalError.
        Intended to be launched from self.run, not manually

        @param end: ending message
        @type end: Finished
        @param fullrelease: if true, also clean memory (Default value = True)
        @type fullrelease: bool
        @return: end message
        @rtype: str
        """
        if self.initialized:
            self.statistic.end(end)
            self.statistic.calcsnapshot(final=True)
            if fullrelease:
                self.crn.close()
                self.param.unlock()
                self.initialized = False
            self.status.gcclean()
            self.signcatch.reset()
            LOGGER.info(f"System {'fully ' if fullrelease else ''}cleaned")
            return f"{end} ({LOGGER.runtime} s)"
        raise InternalError("Attempt to release non-initialized system!")

    def _process(self) -> None:
        """Perform a run step."""
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

    def run(self, release: bool = True) -> str:
        """
        Launch a run.

        The final state can be kept by setting 'release' to false (useful for run directly launched
        from the command line for further direct data manipulation). It is by default released in
        order to release memory as soon as a run is finished during batch runs.

        @param release: Shall the final state be released at run end
        @type release: bool
        @return: end message
        @rtype: str
        """
        LOGGER.reset_timer()
        if MPI_STATUS.ismpi:
            LOGGER.info(f"Launching MPI run from thread #{MPI_STATUS.rank}")
        else:
            LOGGER.info("Launching single run.")
        # Setup working environment
        self._initialize()
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
                    self.statistic.calc()
                    # Check status before next step
                    self.status.checkstep()
                    # Pass gate checkpoint
                    MPI_GATE.checkpoint()
                # exit the gate because the work is finished.
                except Finished as the_end:
                    end = self._release(the_end, release)
                    break
            # exit the gate because asked by collective decision.
            else:
                end = self._release(
                    OOMError("Stopped by kind request of OOM killer"), release
                )
        # Write and log final data, then exit.
        LOGGER.info(f"Run #{MPI_STATUS.rank}={getpid()} finished")
        self.statistic.close()
        return end

    def set_param(self, **kwd: Any) -> None:
        """Set run parameters, overriding the previous ones."""
        try:
            self.param.set_param(**kwd)
        except LockedError:
            raise InternalError("Parameter set after initialization")
        LOGGER.info(f"Parameters changed: {self.param}")
