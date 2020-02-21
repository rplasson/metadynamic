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

from mpi4py import MPI
from time import sleep
from numpy import array, ndarray, empty
from typing import List, Any, Optional


class MpiGate:
    def __init__(self) -> None:
        self.comm: MPI.Intracomm = MPI.COMM_WORLD
        self.size: int = int(self.comm.size)
        self.rank: int = int(self.comm.rank)
        self.other: List[int] = list(range(self.size))
        self.other.pop(self.rank)
        self.snd_state: List[Optional[MPI.Request]] = [None] * (self.size)
        self.sender: bool = False
        self.rcv_state: List[bool] = [False] * (self.size)
        self.gatenum: int = 1
        self.out: List[bool] = [False] * (self.size)
        self.msg: List[int] = [0] * (self.size)

    def checktag(self, src: int, tag: int = 0) -> None:
        if src == self.rank:
            raise ValueError(f"Don't talk to yourself")
        received = self.comm.Iprobe(source=src, tag=tag)
        self.rcv_state[src] = received

    def readtag(self, src: int, tag: int = 0) -> None:
        if src == self.rank:
            raise ValueError(f"Don't talk to yourself")
        if self.rcv_state[src]:
            msg = self.comm.recv(source=src, tag=tag)
            if msg == 2:
                self.out[src] = True
            self.msg[src] = msg
            self.rcv_state[src] = False

    def sendtag(self, dest: int, tag: int = 0, final: bool = False) -> None:
        if dest == self.rank:
            raise ValueError(f"Don't talk to yourself")
        if self.snd_state[dest] is None:
            msg = 2 if final else 1
            self.snd_state[dest] = self.comm.isend(msg, dest=dest, tag=tag)
            self.sender = True

    def wait_sent_tag(self, dest: int, tag: int = 0) -> None:
        if dest == self.rank:
            raise ValueError(f"Don't talk to yourself")
        req = self.snd_state[dest]
        if req is not None:
            req.Wait()
            self.snd_state[dest] = None
        self.sender = False

    def check_mail(self, tag: int = 0) -> int:
        for src in self.other:
            self.checktag(src, tag)
        return sum(self.rcv_state)

    def read_mail(self, tag: int = 0) -> List[int]:
        self.check_mail(tag)
        for src in self.other:
            self.readtag(src, tag)
        res = self.msg.copy()
        self.msg = [0] * self.size
        return res

    def wait_sent(self, tag: int = 0) -> None:
        for dest in self.other:
            self.wait_sent_tag(dest, tag)

    def send_mails(self, tag: int = 0, final: bool = False) -> None:
        for dest in self.other:
            self.sendtag(dest, tag, final)

    def closed(self) -> bool:
        return self.sender or (self.check_mail(self.gatenum) > 0)

    def close(self, final: bool = False) -> None:
        self.send_mails(self.gatenum, final)

    def open(self) -> None:
        self.comm.Barrier()
        self.read_mail(self.gatenum)
        self.wait_sent(self.gatenum)
        self.gatenum += 1
        self.comm.Barrier()

    def exit(self, sleeptime: float = 0.1) -> None:
        self.out[self.rank] = True
        self.close(final=True)
        while True:
            if self.closed():
                self.open()
                if sum(self.out) == self.size:
                    break
            sleep(sleeptime)


class MpiStatus:
    def __init__(self, rootnum: int = 0) -> None:
        self.comm = MPI.COMM_WORLD
        self.size: int = int(self.comm.size)
        self.rank: int = int(self.comm.rank)
        self.ismpi: bool = self.size > 1
        self.rootnum: int = rootnum
        self.gate = MpiGate()

    @property
    def root(self) -> bool:
        return self.rank == self.rootnum

    def max(self, val: Any) -> Any:
        if self.ismpi:
            return self.comm.allreduce(val, op=MPI.MAX)
        return val

    def sortlist(self, data: List[float]) -> List[float]:
        if not self.ismpi:
            data.sort()
            return data
        sendbuf = array(data)
        sendcounts = array(self.comm.gather(len(sendbuf), self.rootnum))
        if self.root:
            recvbuf = empty(sum(sendcounts), dtype=float)
        else:
            recvbuf = None
        self.comm.Gatherv(
            sendbuf=sendbuf, recvbuf=(recvbuf, sendcounts), root=self.rootnum
        )
        if self.root:
            start = 0
            gathered: List[ndarray] = []
            for i in sendcounts:
                gathered.append(recvbuf[start : start + i])
                start = start + i
            fused: List[float] = list(set().union(*[set(i) for i in gathered]))
            fused.sort()
        else:
            fused = []
        fused = self.comm.bcast(fused, root=self.rootnum)
        return fused
