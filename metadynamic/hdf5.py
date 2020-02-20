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

from datetime import datetime
from socket import gethostname
from time import sleep
from typing import Dict, Any, List, Tuple, Mapping
from h5py import File, Group, Dataset, string_dtype
from mpi4py import MPI
from graphviz import Digraph
from numpy import (
    nan,
    nanmean,
    nanstd,
    ndarray,
    array,
    convolve,
    ones,
    empty,
    nan_to_num,
    meshgrid,
)
from pandas import DataFrame

from metadynamic.ends import Finished, FileCreationError, InternalError
from metadynamic.inval import invalidstr, invalidint, invalidfloat, isvalid
from metadynamic.inputs import Readerclass, StatParam, MapParam, Param
from metadynamic.caster import Caster
from metadynamic.json2dot import Data2dot
from metadynamic import __version__


comp_cast = Caster(Dict[str, int])
reac_cast = Caster(Dict[str, List[float]])


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
    def __init__(self, filename: str, maxstrlen: int = 256) -> None:
        if not isvalid(filename) or filename == "":
            raise FileNotFoundError(f"Plese enter a valid output file name")
        self.filename = filename
        self.maxstrlen = maxstrlen
        self.mpi = MpiStatus()
        try:
            if self.mpi.ismpi:
                self.h5file = File(filename, "w", driver="mpio", comm=self.mpi.comm)
            else:
                self.h5file = File(filename, "w")
        except OSError as err:
            raise FileCreationError(f"'{filename}': {err}")
        except ValueError as err:
            raise FileNotFoundError(f"Couldn't find file {filename} : {err}")
        self._init_stat: bool = False
        self._init_log: bool = False

    def init_log(self, maxlog: int) -> None:
        self.maxlog = maxlog
        size = self.mpi.size
        self.logging: Group = self.h5file.create_group("Logging")
        self.logcount: Dataset = self.logging.create_dataset(
            "count", (size,), fillvalue=0, dtype="int32"
        )
        self.logs: Dataset = self.logging.create_dataset(
            "logs",
            (size, maxlog),
            dtype=[
                ("level", "int32"),
                ("time", string_dtype(length=18)),
                ("runtime", "float32"),
                ("message", string_dtype(length=self.maxstrlen)),
            ],
        )
        self._init_log = True

    def write_log(self, level: int, time: str, runtime: float, msg: str) -> None:
        if self._init_log:
            rank = self.mpi.rank
            col = self.logcount[rank]
            self.logs[rank, col] = (level, time, runtime, msg[: self.maxstrlen])
            self.logcount[rank] = col + 1
            if col == self.maxlog:
                self._init_log = False

    def init_stat(
        self,
        datanames: List[str],
        mapnames: List[str],
        params: Dict[str, Any],
        statparam: Dict[str, StatParam],
        mapparam: Dict[str, MapParam],
        comment: str,
        nbcol: int,
    ) -> None:
        size = self.mpi.size
        self.nbcol = nbcol
        self.run: Group = self.h5file.create_group("Run")
        self.run.attrs["version"] = __version__
        self.run.attrs["hostname"] = gethostname()
        self.run.attrs["date"] = datetime.now().strftime("%H:%M:%S, %d/%m/%Y")
        self.run.attrs["threads"] = self.mpi.size
        self.run.attrs["comment"] = comment
        self.params: Group = self.h5file.create_group("Parameters")
        self.dict_as_attr(self.params, params)
        self.statparam: Group = self.params.create_group("Stats")
        self.multiread_as_attr(self.statparam, statparam)
        self.mapparam: Group = self.params.create_group("Maps")
        self.multiread_as_attr(self.mapparam, mapparam)
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
                ("message", string_dtype(length=self.maxstrlen)),
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
            dtype=[("name", string_dtype(length=self.maxstrlen)), ("pop", "int32")],
        )
        self.reacsnap: Dataset = self.snapshots.create_dataset(
            "reactions",
            (size, 1, 1),
            maxshape=(size, None, None),
            dtype=[
                ("name", string_dtype(length=self.maxstrlen)),
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
        self._init_stat = True

    def test_initialized(self) -> None:
        if not self._init_stat:
            raise InternalError("Attempt to write in HDF5 file before intialization")

    def mapsize(self, name: str, categories: List[float]) -> None:
        self.test_initialized()
        self.map_cat[name] = categories
        mapsize = len(categories)
        self.maps[name].resize((self.mpi.size, mapsize, self.nbcol + 1))
        self.maps[name][self.mpi.rank, :, 0] = categories

    def add_map(self, name: str, data: Dict[float, List[float]]) -> None:
        self.test_initialized()
        for catnum, cat in enumerate(self.map_cat[name]):
            try:
                self.maps[name][self.mpi.rank, catnum, 1:] = data[cat]
            except KeyError:
                pass  # No problem, some categories may have been reached by only some processes

    def snapsize(self, maxcomp: int, maxreac: int, maxsnap: int) -> None:
        self.test_initialized()
        self.timesnap.resize((self.mpi.size, maxsnap))
        self.compsnap.resize((self.mpi.size, maxsnap, maxcomp))
        self.reacsnap.resize((self.mpi.size, maxsnap, maxreac))
        self._snapsized = True

    def close(self) -> None:
        self.run.attrs["end"] = datetime.now().strftime("%H:%M:%S, %d/%m/%Y")
        self._init_log = False
        self._init_stat = False
        self.h5file.close()

    def add_data(self, result: List[float]) -> None:
        self.test_initialized()
        try:
            self.data[self.mpi.rank, :, self._currentcol] = result
            self._currentcol += 1
        except ValueError:
            raise InternalError(
                f"No more space in file for #{self.mpi.rank} at column {self._currentcol}"
            )

    def add_end(self, ending: Finished, time: float) -> None:
        self.test_initialized()
        self.end[self.mpi.rank] = (
            ending.num,
            ending.message.encode()[: self.maxstrlen],
            time,
        )

    def add_snapshot(
        self,
        complist: Dict[str, int],
        reaclist: Dict[str, Tuple[float, float]],
        col: int,
        time: float,
    ) -> None:
        self.test_initialized()
        if self._snapsized:
            self.timesnap[self.mpi.rank, col] = time
            for line, data in enumerate(complist.items()):
                self.compsnap[self.mpi.rank, col, line] = (
                    data[0][: self.maxstrlen],
                    data[1],
                )
            for line, (name, (const, rate)) in enumerate(reaclist.items()):
                self.reacsnap[self.mpi.rank, col, line] = (
                    name.encode()[: self.maxstrlen],
                    const,
                    rate,
                )
        else:
            raise InternalError(f"Snapshots data in hdf5 file wasn't properly sized")

    def dict_as_attr(self, group: Group, datas: Dict[str, Any], name: str = "") -> None:
        for key, val in datas.items():
            if name:
                key = f"{name}->{key}"
            try:
                group.attrs[key] = val
            except TypeError:
                self.dict_as_attr(group, val, name=key)

    def multiread_as_attr(self, group: Group, datas: Mapping[str, Readerclass]) -> None:
        for key, val in datas.items():
            subgroup = group.create_group(key)
            self.dict_as_attr(subgroup, val.asdict())


class ResultReader:
    def __init__(self, filename: str) -> None:
        self.filename = filename
        self.h5file: File = File(filename, "r")
        self.run: Group = self.h5file["Run"]
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
        self.mapnames = self.maps.keys()
        self.logging: Group = self.h5file["Logging"]
        self.logcount: Dataset = self.logging["count"]
        self.logs: Dataset = self.logging["logs"]

    def __getitem__(self, field: str) -> ndarray:
        if field in self.datanames:
            return self.get(field)
        if field in self.mapnames:
            return self.getmap(field)
        raise KeyError

    def _loc(self, field: str) -> int:
        try:
            return self.datanames.index(field)
        except ValueError:
            raise ValueError(f"{field} is not a recorded stat name")

    def get(
        self, field: str = invalidstr, method: str = "m", meanlength: int = invalidint,
    ) -> ndarray:
        loc = self._loc(field) if isvalid(field) else slice(None, None, None)
        if method == "*":
            return self.data[:, loc]
        if method[0] == "p":
            res = self.data[int(method[1:]), loc]
        elif method[0] in ["m", "s", "+", "-"]:
            mean = nanmean(self.data[:, loc, :], axis=0)
            if method == "m":
                res = mean
            else:
                std = nanstd(self.data[:, loc, :], axis=0)
                if method == "s":
                    res = std
                else:
                    res = mean + float(method) * std
        else:
            raise ValueError(f"'method'={method} is invalid")
        return (
            # Running mean from https://stackoverflow.com/questions/13728392/moving-average-or-running-mean
            convolve(res, ones((meanlength,)) / meanlength, mode="valid")
            if isvalid(meanlength)
            else res
        )

    def getmap(
        self, field: str, method: str = "m", meanlength: int = invalidint
    ) -> ndarray:
        try:
            data = self.maps[field][:, :, 1:]
        except KeyError:
            raise ValueError(f"{field} is not a recorded map name")
        if method == "*":
            return data
        if method[0] == "p":
            res = data[int(method[1:])]
        elif method[0] in ["m", "s", "+", "-"]:
            mean = nanmean(data, axis=0)
            if method == "m":
                res = mean
            else:
                std = nanstd(data, axis=0)
                if method == "s":
                    res = std
                else:
                    res = mean + float(method) * std
        else:
            raise ValueError(f"'method'={method} is invalid")
        return (
            # Running mean from https://stackoverflow.com/questions/13728392/moving-average-or-running-mean
            convolve(res, ones((meanlength,)) / meanlength, mode="valid")
            if isvalid(meanlength)
            else res
        )

    def getsnap(self, num: int, step: int, parameterfile: str = "") -> Digraph:
        comp = comp_cast(self.snapshots["compounds"][num, step])
        if "" in comp:
            comp.pop("")
        reac = reac_cast(
            {i: (j, k) for i, j, k in self.snapshots["reactions"][num, step]}
        )
        if "" in reac:
            reac.pop("")
        return Data2dot(comp, reac, parameterfile).crn.dot

    def categories(self, field: str) -> ndarray:
        return self.maps[field][0, :, 0]

    def ending(self, num: int) -> Tuple[int, str, float]:
        endnum, message, time = self.end[num]
        return endnum, message.decode(), time

    def endmsg(self, num: int) -> str:
        endnum, message, time = self.ending(num)
        return f"#{num}: ending n°{endnum} at runtime t={time}s; {message}"

    def table(
        self, maps: str = invalidstr, method: str = "m", meanlength: int = invalidint
    ) -> DataFrame:
        if isvalid(maps):
            data = self.getmap(field=maps, method=method, meanlength=meanlength)
            index = self.categories(maps)
        else:
            data = self.get(method=method, meanlength=meanlength)
            index = self.datanames
        return DataFrame(data, index=index)

    def xy(
        self,
        y: str = "ptime",
        x: str = "time",
        method: str = "m",
        xmethod: str = "m",
        meanlength: int = invalidint,
    ) -> Tuple[ndarray, ndarray]:
        x = self.get(field=x, method=xmethod, meanlength=meanlength)
        y = self.get(field=y, method=method, meanlength=meanlength)
        return x, y

    def xypm(
        self,
        y: str = "ptime",
        x: str = "time",
        delta: float = 1.0,
        meanlength: int = invalidint,
    ) -> Tuple[ndarray, ndarray, ndarray]:
        x = self.get(field=x, method="m", meanlength=meanlength)
        yp = self.get(field=y, method=f"+{delta}", meanlength=meanlength)
        ym = self.get(field=y, method=f"-{delta}", meanlength=meanlength)
        return x, yp, ym

    def xyerr(
        self,
        y: str = "ptime",
        x: str = "time",
        delta: float = 1.0,
        meanlength: int = invalidint,
    ) -> Tuple[ndarray, ndarray, ndarray]:
        x = self.get(field=x, method="m", meanlength=meanlength)
        y = self.get(field=y, method=f"m", meanlength=meanlength)
        err = self.get(field=y, method=f"s", meanlength=meanlength)*delta
        return x, y, err

    def xyz(
        self,
        field: str,
        method: str = "m",
        tmethod: str = "m",
        meanlength: int = invalidint,
        nanval: float = 0.0,
        posinf: float = invalidfloat,
        neginf: float = invalidfloat,
    ) -> Tuple[ndarray, ndarray, ndarray]:
        time = self.get(field="time", method=tmethod, meanlength=meanlength)
        categories = self.categories(field)
        x, y = meshgrid(time, categories)
        z = nan_to_num(
            self.getmap(field=field, method=method, meanlength=meanlength),
            nan=nanval,
            posinf=posinf if isvalid(posinf) else None,
            neginf=neginf if isvalid(neginf) else None,
        )
        return x, y, z

    @property
    def printinfo(self) -> str:
        version = self.run.attrs["version"]
        hostname = self.run.attrs["hostname"]
        start = self.run.attrs["date"]
        threads = self.run.attrs["threads"]
        comment = self.run.attrs["comment"]
        end = self.run.attrs["end"]
        endline = "\n"
        return (
            f"----------------\n"
            f"{comment}\n"
            f"----------------\n"
            f"metadynamic version {version}, "
            f"ran on {threads} threads on {hostname}\n"
            f"from {start} to {end}\n"
            f"results saved in '{self.filename}'\n"
            f"----------------\n"
            f"{endline.join([self.endmsg(i) for i in range(threads)])}\n"
            f"----------------"
        )

    @property
    def runinfo(self) -> Dict[str, Any]:
        return dict(self.run.attrs)

    @property
    def parameters(self) -> Param:
        params = dict(self.params.attrs)
        res = params.copy()
        for key, val in params.items():
            if "->" in key:
                prekey, postkey = key.split("->")
                try:
                    res[prekey][postkey] = val
                except KeyError:
                    res[prekey] = {postkey: val}
                res.pop(key)
        return Param.readdict(res)

    @property
    def statparam(self) -> Dict[str, StatParam]:
        return StatParam.multipledict(
            {key: dict(val.attrs) for key, val in self.params["Stats"].items()}
        )

    @property
    def mapparam(self) -> Dict[str, MapParam]:
        return MapParam.multipledict(
            {key: dict(val.attrs) for key, val in self.params["Maps"].items()}
        )

    @property
    def fulllog(self) -> ndarray:
        # add simple system for log display/search (depending on level etc) ?
        maxcol = max(self.logcount)
        return self.logs[:, :maxcol]


class Saver:
    writer: ResultWriter

    @classmethod
    def setsaver(cls, filename: str, maxstrlen: int = 256) -> None:
        cls.writer = ResultWriter(filename, maxstrlen)
