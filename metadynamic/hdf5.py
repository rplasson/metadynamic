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
from typing import Dict, Any, List, Tuple, Mapping
from h5py import File, Group, Dataset, string_dtype
from numpy import nan

from metadynamic.ends import Finished, FileCreationError, InternalError
from metadynamic.inval import isvalid
from metadynamic.inputs import Readerclass, Param, RulesetParam
from metadynamic.caster import Caster
from metadynamic.mpi import Parallel
from metadynamic import __version__


comp_cast = Caster(Dict[str, int])
reac_cast = Caster(Dict[str, List[float]])


class ResultWriter(Parallel):
    def __init__(self, filename: str, maxstrlen: int = 256, lengrow: int = 10) -> None:
        if not isvalid(filename) or filename == "":
            raise FileNotFoundError(f"Plese enter a valid output file name")
        self.filename = filename
        self.maxstrlen = maxstrlen
        self.lengrow = lengrow
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
        self.dlog = maxlog
        size = self.mpi.size
        self.logging: Group = self.h5file.create_group("Logging")
        self.logcount: Dataset = self.logging.create_dataset(
            "count", (size,), fillvalue=0, dtype="int32"
        )
        self.logs: Dataset = self.logging.create_dataset(
            "logs",
            (size, maxlog),
            maxshape=(size, None),
            dtype=[
                ("level", "int32"),
                ("time", string_dtype(length=18)),
                ("runtime", "float32"),
                ("message", string_dtype(length=self.maxstrlen)),
            ],
        )
        self.gate.register_function("addlog", self.add_log_line)
        self._init_log = True

    def add_log_line(self) -> None:
        self.maxlog = self.maxlog + self.dlog
        self.logs.resize(self.maxlog, axis=1)

    def write_log(self, level: int, time: str, runtime: float, msg: str) -> None:
        if self._init_log:
            rank = self.mpi.rank
            col = self.logcount[rank]
            self.logs[rank, col] = (level, time, runtime, msg[: self.maxstrlen])
            self.logcount[rank] = col + 1
            if (self.maxlog - col) < self.lengrow:
                self.gate.close("addlog")

    def close_log(self, cutline: int) -> None:
        self.logs.resize(cutline, axis=1)
        self._init_log = False

    def init_stat(
        self,
        datanames: List[str],
        mapnames: List[str],
        params: Param,
        comment: str,
        nbcol: int,
    ) -> None:
        size = self.mpi.size
        self.nbcol = nbcol + self.lengrow
        self.dcol = nbcol
        self.run: Group = self.h5file.create_group("Run")
        self.run.attrs["version"] = __version__
        self.run.attrs["hostname"] = gethostname()
        self.run.attrs["date"] = datetime.now().strftime("%H:%M:%S, %d/%m/%Y")
        self.run.attrs["threads"] = self.mpi.size
        self.run.attrs["comment"] = comment
        self.params: Group = self.h5file.create_group("Parameters")
        self.dict_as_attr(self.params, params.asdict())
        self.statparam: Group = self.params.create_group("Stats")
        self.multiread_as_attr(self.statparam, params.statparam)
        self.mapparam: Group = self.params.create_group("Maps")
        self.multiread_as_attr(self.mapparam, params.mapsparam)
        self.ruleparam: Group = self.params.create_group("Rules")
        self.dict_as_attr(
            self.ruleparam, RulesetParam.readfile(params.rulemodel).asdict()
        )
        self.dataset: Group = self.h5file.create_group("Dataset")
        self.dataset.attrs["datanames"] = datanames
        self.data: Dataset = self.dataset.create_dataset(
            "results",
            (size, len(datanames), self.nbcol),
            maxshape=(size, len(datanames), None),
            fillvalue=nan,
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
        self.reacsnapsaved: Dataset = self.snapshots.create_dataset(
            "reactions_saved", (size, 1), maxshape=(size, None), dtype=bool,
        )
        self._snapsized: bool = False
        self.maps: Group = self.h5file.create_group("Maps")
        for name in mapnames:
            self.maps.create_dataset(
                name,
                (size, 1, self.nbcol + 1),
                maxshape=(size, None, None),
                fillvalue=nan,
                dtype="float32",
            )
        self.map_cat: Dict[str, List[float]] = {}
        self.currentcol = 0
        self.gate.register_function("addcol", self.add_col)
        self._init_stat = True

    def test_initialized(self) -> None:
        if not self._init_stat:
            raise InternalError("Attempt to write in HDF5 file before intialization")

    def add_col(self) -> None:
        self.nbcol = self.nbcol + self.dcol
        self.data_resize(self.nbcol)

    def data_resize(self, maxcol: int) -> None:
        self.data.resize(maxcol, axis=2)
        for datamap in self.maps.values():
            datamap.resize(maxcol + 1, axis=2)

    def mapsize(self, name: str, categories: List[float]) -> None:
        self.test_initialized()
        self.map_cat[name] = categories
        mapsize = len(categories)
        self.maps[name].resize(mapsize, axis=1)
        self.maps[name][self.mpi.rank, :, 0] = categories

    def add_map(self, name: str, data: Dict[float, List[float]]) -> None:
        self.test_initialized()
        for catnum, cat in enumerate(self.map_cat[name]):
            try:
                length = len(data[cat])
                self.maps[name][self.mpi.rank, catnum, 1 : length + 1] = data[cat]
            except KeyError:
                pass  # No problem, some categories may have been reached by only some processes

    def snapsize(self, maxcomp: int, maxreac: int, maxsnap: int) -> None:
        self.test_initialized()
        self.timesnap.resize((self.mpi.size, maxsnap))
        self.compsnap.resize((self.mpi.size, maxsnap, maxcomp))
        self.reacsnap.resize((self.mpi.size, maxsnap, maxreac))
        self.reacsnapsaved.resize((self.mpi.size, maxsnap))
        self._snapsized = True

    def close(self) -> None:
        self.run.attrs["end"] = datetime.now().strftime("%H:%M:%S, %d/%m/%Y")
        self._init_log = False
        self._init_stat = False
        self.h5file.close()

    def add_data(self, result: List[float]) -> None:
        self.test_initialized()
        try:
            self.data[self.mpi.rank, :, self.currentcol] = result
            self.currentcol += 1
            if (self.nbcol - self.currentcol) < self.lengrow:
                self.gate.close("addcol")
        except ValueError:
            raise InternalError(
                f"No more space in file for #{self.mpi.rank} at column {self.currentcol}"
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
            self.reacsnapsaved[self.mpi.rank, col] = len(reaclist) > 0
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
