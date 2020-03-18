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

import numpy as np

from mpi4py import MPI
from time import sleep
from typing import List, Dict, Any, Callable, Iterable, Optional, Type
from types import TracebackType


def nop() -> None:
    pass


class Cont:
    def __init__(self) -> None:
        self._cont = True

    def stop(self) -> None:
        self._cont = False

    def reset(self) -> None:
        self._cont = True

    def __bool__(self) -> bool:
        return self._cont


class MpiGate:
    def __init__(
        self,
        taginit: int = 1,
        operations: Optional[Dict[str, Callable[[], None]]] = None,
    ) -> None:
        self.comm: MPI.Intracomm = MPI.COMM_WORLD
        self.size: int = int(self.comm.size)
        self.rank: int = int(self.comm.rank)
        self.procs: Iterable[int] = range(self.size)
        self.snd_state: List[Optional[MPI.Request]] = [None] * (self.size)
        self.rcv_state: List[bool] = [False] * (self.size)
        self.gatenum: int
        self.init(taginit)
        self.running: bool = False
        self.launched: bool = False
        self.nb_running: int = self.size
        self.mem_divide: int = self.size
        self.cont = Cont()
        self.msg: List[int] = [0] * (self.size)
        self._op: Dict[int, Callable[[], None]] = {
            1: nop,
            2: self.check_all_out,
            3: self.cont.stop,
            4: self.oomkill,
        }
        self._opnum: Dict[str, int] = {"nop": 1, "final": 2, "exit": 3, "oom": 4}
        if operations:
            for name, op in operations.items():
                self.register_function(name, op)

    def init(self, taginit: int = 1) -> None:
        self.gatenum = taginit

    def __enter__(self) -> "MpiGate":
        self.cont.reset()
        return self

    def __exit__(
        self,
        exctype: Optional[Type[BaseException]],
        excinst: Optional[BaseException],
        exctb: Optional[TracebackType],
    ) -> None:
        self.exit()

    def register_function(self, opname: str, func: Callable[[], None]) -> None:
        funcnum = len(self._op) + 1
        self._op[funcnum] = func
        self._opnum[opname] = funcnum

    def operate(self, funcnum: int) -> None:
        self._op[funcnum]()

    def oomkill(self) -> None:
        self.mem_divide = self.nb_running
        runningpos = np.where(self.comm.allgather(self.running))
        for running in runningpos[0][::2]:
            if running == self.rank:
                self.cont.stop()
            self.mem_divide -= 1

    def check_all_out(self) -> None:
        self.nb_running = self.comm.allreduce(self.running, op=MPI.SUM)

    def check_msg(self, tag: int) -> None:
        for src in self.procs:
            received = self.comm.Iprobe(source=src, tag=tag)
            self.rcv_state[src] = received

    def read_msg(self, tag: int) -> List[int]:
        self.check_msg(tag)
        for src in self.procs:
            if self.rcv_state[src]:
                msg = self.comm.recv(source=src, tag=tag)
                self.msg[src] = msg
                self.rcv_state[src] = False
        res = self.msg.copy()
        self.msg = [0] * self.size
        return res

    def wait_sent(self, tag: int) -> None:
        for dest in self.procs:
            req = self.snd_state[dest]
            if req is not None:
                req.Wait()
                self.snd_state[dest] = None

    def send_msg(self, tag: int, msg: str) -> None:
        for dest in self.procs:
            if self.snd_state[dest] is None:
                msgnum = self._opnum[msg]
                self.snd_state[dest] = self.comm.isend(msgnum, dest=dest, tag=tag)

    @property
    def closed(self) -> bool:
        self.check_msg(tag=self.gatenum)
        return sum(self.rcv_state) > 0

    def close(self, msg: str = "nop") -> None:
        if not self.launched:
            raise ValueError("Cannot close a gate that has not been launched.")
        self.send_msg(tag=self.gatenum, msg=msg)

    def open(self) -> None:
        # Wait for everyone for exchanging messages
        if not self.launched:
            raise ValueError("Cannot open a gate that has not been launched.")
        self.comm.Barrier()
        # get operation list, sorted to be processed in order.
        # 0 is removed as it corresponds to 'no message'
        operations = list(set(self.read_msg(tag=self.gatenum)) - {0})
        operations.sort()
        self.wait_sent(self.gatenum)
        self.gatenum += 1
        # Wait for everyone for processing operations
        self.comm.Barrier()
        for op in operations:
            self.operate(op)

    def checkpoint(self) -> None:
        if self.closed:
            self.open()

    def exit(self, sleeptime: float = 0.1) -> None:
        self.running = False
        self.close("final")
        while self.nb_running > 0:
            self.checkpoint()
            sleep(sleeptime)
        self.launched = False

    def launch(self) -> "MpiGate":
        self.running = True
        self.launched = True
        return self


class MpiStatus:
    def __init__(self, rootnum: int = 0) -> None:
        self.comm = MPI.COMM_WORLD
        self.size: int = int(self.comm.size)
        self.rank: int = int(self.comm.rank)
        self.ismpi: bool = self.size > 1
        self.rootnum: int = rootnum

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
MPI_GATE = MpiGate()
