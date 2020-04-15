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

"""
metadynamic.mpi
===============

Control the coordination between thread in MPI runs

Provides
--------

 - L{Cont}: switchable "continue" flag

 - L{MpiGate}: MPI barrier that can be opened/closed by any thread, for on-demand requests of
   synchronization

 - L{MpiStatus}: Interface to general MPI operations.

 - L{MPI_STATUS}: global MpiStatus object

 - L{MPI_GATE}: global MpiGate object

"""

from datetime import datetime
from socket import gethostname
from time import sleep
from typing import List, Dict, Any, Callable, Iterable, Optional, Type
from types import TracebackType
from mpi4py import MPI


import numpy as np


def nop() -> None:
    """'no operation function' (surprisingly useful)"""


class Cont:
    """switchable "continue" flag"""

    def __init__(self) -> None:
        """Create an object in initial 'continue' state"""
        self._cont = True
        """'continue' state flag"""

    def stop(self) -> None:
        """set object to 'stop' state"""
        self._cont = False

    def reset(self) -> None:
        """reset object to initial 'continue' state"""
        self._cont = True

    def __str__(self) -> str:
        """convert to a string representing the object state"""
        return "continue" if self._cont else "stop"

    def __bool__(self) -> bool:
        """True in 'continue' state, False in 'stop' state"""
        return self._cont


class MpiGate:
    """
    MPI barrier that can be opened/closed by any thread, for on-demand requests of synchronization
    """

    def __init__(
        self,
        taginit: int = 1,
        operations: Optional[Dict[str, Callable[[], None]]] = None,
        sleeptime: float = 0.1,
    ) -> None:
        """
        create a gate

        @param taginit: initial tag number. each request will be tagged from tag increments
            starting from this value
        @type taginit: int
        @param operations: Dictionary of synchronisation operations {name: operation}
        @type operation: Dict[str, Callable[[], None]]
        @param sleeptime: time spent (in seconds) in idle loops when waiting before the exit gate
        """
        self.comm: MPI.Intracomm = MPI.COMM_WORLD
        """MPI communicator"""
        self.size: int = int(self.comm.size)
        """total number of MPI threads"""
        self.rank: int = int(self.comm.rank)
        """thread number"""
        self.procs: Iterable[int] = range(self.size)
        """list of thread numbers"""
        self.snd_state: List[Optional[MPI.Request]] = [None] * (self.size)
        """list of send states (to which a message was sent?)"""
        self.rcv_state: List[bool] = [False] * (self.size)
        """list of received states (from which a message was received?)"""
        self.nb_running: int = self.size
        """number of threads still runnning"""
        self.mem_divide: int = self.size
        """number of thread actaually using memory"""
        self.gatenum: int
        """Present gate operation n umber"""
        self.sleeptime: float
        """time spent (in seconds) in idle loops when waiting before the exit gate"""
        self.init(taginit, sleeptime)
        self.running: bool = False
        """running flag state"""
        self.launched: bool = False
        """launched flag state"""
        self.cont = Cont()
        """coninue toggable state"""
        self.msg: List[int] = [0] * (self.size)
        """list of MPI message box"""
        self._op: Dict[int, Callable[[], None]] = {
            1: nop,
            2: self.check_all_out,
            3: self.cont.stop,
            4: self.oomkill,
        }
        """dictionary of synchronisation operation functions"""
        self._opnum: Dict[str, int] = {"nop": 1, "final": 2, "exit": 3, "oom": 4}
        """dictionary of synchronisation operation numbers"""
        if operations:
            for name, oper in operations.items():
                self.register_function(name, oper)

    def init(self, taginit: int = -1, sleeptime: float = -1) -> None:
        """
        Initialize a gate to its starting position

        @param taginit: initial tag number (if <0 value given, stays untouched)
        @type taginit: int
        @param sleeptime: sleep time  (if <0 value given, stays untouched)
        @type sleeptime: float
        """
        self.nb_running = self.size
        self.mem_divide = self.size
        if taginit > 0:
            self.gatenum = taginit
        if sleeptime > 0:
            self.sleeptime = sleeptime

    def __enter__(self) -> "MpiGate":
        """Operations performed when entering the gate (via context manager)"""
        self.init()
        self.cont.reset()
        return self

    def __exit__(
        self,
        exctype: Optional[Type[BaseException]],
        excinst: Optional[BaseException],
        exctb: Optional[TracebackType],
    ) -> None:
        """
        Operations performed at gate exit.

        The parameters (automatically sent by the context manager) are ignored.
        """
        self.running = False
        self.close("final")
        while self.nb_running > 0:
            sleep(self.sleeptime)
            self.checkpoint()
        self.launched = False

    def register_function(self, opname: str, func: Callable[[], None]) -> None:
        """
        register a new synchronisation operation function

        @param opname: name of the operation (will override a previously registered one)
        @type opname: str
        @param fun: synchronisation operation function
        @type func: Callable[[], None]
        """
        funcnum = len(self._op) + 1
        self._op[funcnum] = func
        self._opnum[opname] = funcnum

    def operate(self, funcnum: int) -> None:
        """
        Perform synchronisation operation function num 'funcnum'

        @param funcnum: synchronisation operation function number
        @type funcnum: int
        """
        self._op[funcnum]()

    def oomkill(self) -> None:
        """
        Synchronisation operation to be called when too much memory is used.

        Half of still running processes will be kindly asked to stop.
        """
        self.mem_divide = self.nb_running
        runningpos = np.where(self.comm.allgather(self.running))
        for running in runningpos[0][::2]:
            if running == self.rank:
                self.cont.stop()
            self.mem_divide -= 1
        if self.mem_divide == 0:
            self.mem_divide = 1

    def check_all_out(self) -> None:
        """
        Synchronisation operation function to be called when a thread enters the exit.

        The number of still running threads is updated.
        """
        self.nb_running = self.comm.allreduce(self.running, op=MPI.SUM)

    def check_msg(self) -> None:
        """
        """
        for src in self.procs:
            received = self.comm.Iprobe(source=src, tag=self.gatenum)
            self.rcv_state[src] = received

    def read_msg(self) -> List[int]:
        self.check_msg()
        for src in self.procs:
            if self.rcv_state[src]:
                msg = self.comm.recv(source=src, tag=self.gatenum)
                self.msg[src] = msg
                self.rcv_state[src] = False
        res = self.msg.copy()
        self.msg = [0] * self.size
        return res

    def wait_sent(self) -> None:
        for dest in self.procs:
            req = self.snd_state[dest]
            if req is not None:
                req.Wait()
                self.snd_state[dest] = None

    def send_msg(self, msg: str) -> None:
        for dest in self.procs:
            if self.snd_state[dest] is None:
                msgnum = self._opnum[msg]
                self.snd_state[dest] = self.comm.isend(msgnum, dest=dest, tag=self.gatenum)

    @property
    def closed(self) -> bool:
        self.check_msg()
        return sum(self.rcv_state) > 0

    def close(self, msg: str = "nop") -> None:
        if not self.launched:
            raise ValueError("Cannot close a gate that has not been launched.")
        self.send_msg(msg)

    def open(self) -> None:
        # Wait for everyone for exchanging messages
        if not self.launched:
            raise ValueError("Cannot open a gate that has not been launched.")
        self.comm.Barrier()
        # get operation list, sorted to be processed in order.
        # 0 is removed as it corresponds to 'no message'
        operations = list(set(self.read_msg()) - {0})
        operations.sort()
        self.wait_sent()
        self.gatenum += 1
        # Wait for everyone for processing operations
        self.comm.Barrier()
        for oper in operations:
            self.operate(oper)

    def checkpoint(self) -> None:
        if self.closed:
            self.open()

    def launch(self) -> "MpiGate":
        self.running = True
        self.launched = True
        return self


class MpiStatus:
    def __init__(self, rootnum: int = 0) -> None:
        self.comm = MPI.COMM_WORLD
        self.hostname = gethostname()
        self.size: int = int(self.comm.size)
        self.rank: int = int(self.comm.rank)
        self.ismpi: bool = self.size > 1
        self.rootnum: int = rootnum
        self._starttime: str = "[not-started]"

    def init(self, timeformat: str) -> None:
        self._starttime = self.bcast(datetime.now().strftime(timeformat))

    @property
    def starttime(self) -> str:
        return self._starttime

    @property
    def root(self) -> bool:
        return self.rank == self.rootnum

    def max(self, val: Any) -> Any:
        if self.ismpi:
            return self.comm.allreduce(val, op=MPI.MAX)
        return val

    def bcast(self, val: Any) -> Any:
        if self.root:
            res = val
        else:
            res = None
        res = self.comm.bcast(res, root=self.rootnum)
        return res

    def sortlist(self, data: List[float]) -> List[float]:
        if not self.ismpi:
            data.sort()
            return data
        sendbuf = np.array(data)
        sendcounts = np.array(self.comm.gather(len(sendbuf), self.rootnum))
        if self.root:
            recvbuf = np.empty(sum(sendcounts), dtype=float)
        else:
            recvbuf = None
        self.comm.Gatherv(
            sendbuf=sendbuf, recvbuf=(recvbuf, sendcounts), root=self.rootnum
        )
        if self.root:
            start = 0
            gathered: List[np.ndarray] = []
            for i in sendcounts:
                gathered.append(recvbuf[start : start + i])
                start = start + i
            fused: List[float] = list(set().union(*[set(i) for i in gathered]))
            fused.sort()
        else:
            fused = []
        fused = self.comm.bcast(fused, root=self.rootnum)
        return fused


MPI_STATUS = MpiStatus()
"""global MpiStatus object"""

MPI_GATE = MpiGate()
"""global MpiGate object"""
