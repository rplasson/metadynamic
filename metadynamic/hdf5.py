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

from time import sleep
from typing import Dict, Any, List, Tuple
from h5py import File, Group, Dataset, string_dtype
from mpi4py import MPI
from numpy import nan, nanmean, nanstd, ndarray, array, convolve, ones, empty
from pandas import DataFrame

from metadynamic.ends import InitError, Finished
from metadynamic.inval import invalidstr, invalidint, isvalid


class MpiStatus:
    def __init__(self, rootnum: int = 0) -> None:
        self.comm = MPI.COMM_WORLD
        self.size: int = int(self.comm.size)
        self.rank: int = int(self.comm.rank)
        self.ismpi: bool = self.size > 1
        self.rootnum: int = rootnum
        if not self.ismpi:
            self.comm = None
            # for mulithread pickling in non -mpi multithread
            # temporary hack, non-mpi multithread to be removed at term.

    @property
    def root(self) -> bool:
        return self.rank == self.rootnum

    def barrier(self, sleeptime: float = 0.01, tag: int = 0) -> None:
        # From  https://goo.gl/NofOO9
        if self.ismpi:
            mask = 1
            while mask < self.size:
                dst = (self.rank + mask) % self.size
                src = (self.rank - mask + self.size) % self.size
                req = self.comm.isend(None, dst, tag)
                while not self.comm.Iprobe(src, tag):
                    sleep(sleeptime)
                self.comm.recv(None, src, tag)
                req.Wait()
                mask <<= 1

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


class ResultWriter:
    def __init__(
        self,
        filename: str,
        datanames: List[str],
        mapnames: List[str],
        nbcol: int,
        maxstrlen: int = 256,
    ) -> None:
        self.filename = filename
        self.mpi = MpiStatus()
        size = self.mpi.size
        if self.mpi.ismpi:
            self.h5file = File(filename, "w", driver="mpio", comm=self.mpi.comm)
        else:
            self.h5file = File(filename, "w")
        self.nbcol = nbcol
        self.params: Group = self.h5file.create_group("Parameters")
        self.dataset: Group = self.h5file.create_group("Dataset")
        self.dataset.attrs["datanames"] = datanames
        self.data: Dataset = self.dataset.create_dataset(
            "results", (size, len(datanames), nbcol), fillvalue=nan
        )
        self.end: Dataset = self.dataset.create_dataset(
            "end",
            (size,),
            dtype=[
                ("num", "int32"),
                ("message", string_dtype(length=maxstrlen)),
                ("runtime", "float32"),
            ],
        )
        self.snapshots: Group = self.h5file.create_group("Snapshots")
        self.timesnap: Dataset = self.snapshots.create_dataset(
            "time", (size, 1), maxshape=(size, None), dtype="float32",
        )
        self.compsnap: Dataset = self.snapshots.create_dataset(
            "compounds",
            (size, 1, 1),
            maxshape=(size, None, None),
            dtype=[("name", string_dtype(length=maxstrlen)), ("pop", "int32")],
        )
        self.reacsnap: Dataset = self.snapshots.create_dataset(
            "reactions",
            (size, 1, 1),
            maxshape=(size, None, None),
            dtype=[
                ("name", string_dtype(length=maxstrlen)),
                ("const", "float32"),
                ("rate", "float32"),
            ],
        )
        self._snapsized: bool = False
        self.maps: Group = self.h5file.create_group("Maps")
        for name in mapnames:
            self.maps.create_dataset(
                name,
                (size, 1, nbcol + 1),
                maxshape=(size, None, nbcol + 1),
                fillvalue=nan,
                dtype="float32",
            )
        self.map_cat: Dict[str, List[float]] = {}
        self._currentcol = 0

    def mapsize(self, name: str, categories: List[float]) -> None:
        self.map_cat[name] = categories
        mapsize = len(categories)
        self.maps[name].resize((self.mpi.size, mapsize, self.nbcol + 1))
        self.maps[name][self.mpi.rank, :, 0] = categories

    def add_map(self, name: str, data: Dict[float, List[float]]) -> None:
        for catnum, cat in enumerate(self.map_cat[name]):
            try:
                self.maps[name][self.mpi.rank, catnum, 1:] = data[cat]
            except KeyError:
                pass  # No problem, some categories may have been reached by only some processes

    def snapsize(self, maxcomp: int, maxreac: int, maxsnap: int) -> None:
        self.timesnap.resize((self.mpi.size, maxsnap))
        self.compsnap.resize((self.mpi.size, maxcomp, maxsnap))
        self.reacsnap.resize((self.mpi.size, maxreac, maxsnap))
        self._snapsized = True

    def close(self) -> None:
        self.h5file.close()

    def add_parameter(self, params: Dict[str, Any], name: str = "") -> None:
        self.dict_as_attr(self.params, params, name)

    def add_data(self, result: List[float]) -> None:
        try:
            self.data[self.mpi.rank, :, self._currentcol] = result
            self._currentcol += 1
        except ValueError:
            raise ValueError(
                f"No more space in file for #{self.mpi.rank} at column {self._currentcol}"
            )

    def add_end(self, ending: Finished, time: float) -> None:
        self.end[self.mpi.rank] = (ending.num, ending.message.encode(), time)

    def add_snapshot(
        self,
        complist: Dict[str, int],
        reaclist: Dict[str, Tuple[float, float]],
        col: int,
        time: float,
    ) -> None:
        if self._snapsized:
            self.timesnap[self.mpi.rank, col] = time
            for line, data in enumerate(complist.items()):
                self.compsnap[self.mpi.rank, line, col] = data
            for line, (name, (const, rate)) in enumerate(reaclist.items()):
                self.reacsnap[self.mpi.rank, line, col] = (name.encode(), const, rate)
        else:
            raise InitError(f"Snapshots data in hdf5 file wasn't properly sized")

    def dict_as_attr(self, group: Group, datas: Dict[str, Any], name: str = "") -> None:
        for key, val in datas.items():
            if name:
                key = f"{name}->{key}"
            try:
                group.attrs[key] = val
            except TypeError:
                self.dict_as_attr(group, val, name=key)


class ResultReader:
    def __init__(self, filename: str) -> None:
        self.h5file: File = File(filename, "r")
        self.params: Group = self.h5file["Parameters"]
        self.dataset: Group = self.h5file["Dataset"]
        self.datanames: List[str] = list(self.dataset.attrs["datanames"])
        self.data: Dataset = self.dataset["results"]
        self.end: Dataset = self.dataset["end"]
        self.snapshots: Group = self.h5file["Snapshots"]
        self.timesnap: Dataset = self.snapshots["time"]
        self.compsnap: Dataset = self.snapshots["compounds"]
        self.reacsnap: Dataset = self.snapshots["reactions"]
        self.maps: Group = self.h5file["Maps"]

    def _loc(self, field: str) -> int:
        return self.datanames.index(field)

    def get(
        self,
        field: str = invalidstr,
        procnum: str = invalidstr,
        meanlength: int = invalidint,
    ) -> ndarray:
        loc = self._loc(field) if isvalid(field) else slice(None, None, None)
        if isvalid(procnum):
            if procnum.isnumeric():
                res = self.data[int(procnum), loc]
            elif procnum[0] in ["m", "s", "+", "-"]:
                mean = nanmean(self.data[:, loc, :], axis=0)
                if procnum == "m":
                    res = mean
                else:
                    std = nanstd(self.data[:, loc, :], axis=0)
                    if procnum == "s":
                        res = std
                    else:
                        res = mean + float(procnum) * std
            else:
                raise ValueError(f"'procnum'={procnum} is invalid")
            return (
                # Running mean from https://stackoverflow.com/questions/13728392/moving-average-or-running-mean
                convolve(res, ones((meanlength,)) / meanlength, mode="valid")
                if isvalid(meanlength)
                else res
            )
        return self.data[:, loc]

    def getmap(
        self, field: str, procnum: str = invalidstr, meanlength: int = invalidint
    ) -> ndarray:
        try:
            data = self.maps[field][:, :, 1:]
        except KeyError:
            raise ValueError(f"{field} is not a recorded map name")
        if isvalid(procnum):
            if procnum.isnumeric():
                res = data[int(procnum)]
            elif procnum[0] in ["m", "s", "+", "-"]:
                mean = nanmean(data, axis=0)
                if procnum == "m":
                    res = mean
                else:
                    std = nanstd(data, axis=0)
                    if procnum == "s":
                        res = std
                    else:
                        res = mean + float(procnum) * std
            else:
                raise ValueError(f"'procnum'={procnum} is invalid")
            return (
                # Running mean from https://stackoverflow.com/questions/13728392/moving-average-or-running-mean
                convolve(res, ones((meanlength,)) / meanlength, mode="valid")
                if isvalid(meanlength)
                else res
            )
        return data

    def ending(self, num: int) -> Tuple[int, str, float]:
        endnum, message, time = self.end[num]
        return endnum, message.decode(), time

    def table(
        self, maps=invalidstr, procnum: str = "m", meanlength: int = invalidint
    ) -> DataFrame:
        if isvalid(maps):
            data = self.getmap(field=maps, procnum=procnum, meanlength=meanlength)
            index = self.maps[maps][0, :, 0]
        else:
            data = self.get(procnum=procnum, meanlength=meanlength)
            index = self.datanames
        return DataFrame(data, index=index)

    def xy(
        self,
        x: str = "time",
        y: str = "ptime",
        procnum: str = "m",
        meanlength: int = invalidint,
    ) -> Tuple[ndarray, ndarray]:
        x = self.get(field=x, procnum=procnum, meanlength=meanlength)
        y = self.get(field=y, procnum=procnum, meanlength=meanlength)
        return x, y
