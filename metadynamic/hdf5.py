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
metadynamic.hdf5
================

interface to hdf5 result file


Provides:
---------

 - L{ResultWriter}: storage for all the simulation data and results
"""

from datetime import datetime
from typing import Dict, Any, List, Tuple, Mapping, Callable
from h5py import File, Group, Dataset, string_dtype

import numpy as np

from metadynamic.ends import Finished, FileCreationError, InternalError
from metadynamic.inval import isvalid, invalidfloat
from metadynamic.inputs import Readerclass, Param
from metadynamic.caster import Caster
from metadynamic.mpi import MPI_GATE, MPI_STATUS
from metadynamic.version import __version__


comp_cast: Callable[[Any], Dict[str, int]] = Caster(Dict[str, int])
"""Caster to a compound field"""
reac_cast: Callable[[Any], Dict[str, List[float]]] = Caster(Dict[str, List[float]])
"""Caster to a reaction field"""


class ResultWriter:
    """Storage for all the simulation data and results"""

    def __init__(
        self,
        filename: str,
        maxstrlen: int = 256,
        lengrow: int = 10,
        timeformat: str = "%H:%M:%S, %d/%m/%y",
    ) -> None:
        """
        Open the hdf5 result file

        @param filename: name of the hdf5 file
        @type filename: str
        @param maxstrlen: maximum length of strings stored in the file (Default value = 256)
        @type maxlength: int
        @param lengrow: maximal remaining data space left empty before adding more space (default value = 10)
        @type lengrow: int
        @param timeformat: time/date formatting string
        @type timeformat: str
        """
        if not isvalid(filename) or filename == "":
            raise FileNotFoundError(f"Plese enter a valid output file name")
        self.filename: str = MPI_STATUS.bcast(filename)
        """name of hdf5 file"""
        self.maxstrlen: int = maxstrlen
        """maximum length of strings stored in the file"""
        self.lengrow: int = lengrow
        """maximal remaining data space left empty before adding more space"""
        self.timeformat: str = timeformat
        """time/date formatting string"""
        self.h5file: File
        """hdf5 file object"""
        try:
            if MPI_STATUS.ismpi:
                self.h5file = File(filename, "w", driver="mpio", comm=MPI_STATUS.comm)
            else:
                self.h5file = File(filename, "w")
        except OSError as err:
            raise FileCreationError(f"'{filename}': {err}")
        except ValueError as err:
            raise FileNotFoundError(f"Couldn't find file {filename} : {err}")
        # result data
        self._init_stat: bool = False
        """flag indicating if the writer is initialized"""
        self.nbcol: int
        """Available number of columns for writing data"""
        self.dcol: int
        """Number of column increment to add when space goes missing"""
        self.run: Group
        """hdf5 Group 'Run' for storing generic run informations"""
        self.params: Group
        """hdf5 Group 'Parameters' for storing run parameters"""
        self.statparam: Group
        """hdf5 Group 'Stats' for storing statistics parameters"""
        self.mapparam: Group
        """hdf5 Group 'Maps' for storing maps parameters"""
        self.ruleparam: Group
        """hdf5 Group 'Rules' for storing ruleset parameters"""
        self.dataset: Group
        """hdf5 Group 'Dataset' for storing result data"""
        self.data: Dataset
        """hdf5 Dataset 'Dataset/results' (recording results)"""
        self.end: Dataset
        """hdf5 Dataset 'Dataset/end' (recording end messages)"""
        self.snapshots: Group
        """hdf5 Group 'Snapshots' for storing snapshots"""
        self.timesnap: Dataset
        """hdf5 Dataset 'Snapshots/time' (times of snapshots)"""
        self.compsnap: Dataset
        """hdf5 Dataset 'Snapshots/compounds' (compounds snapshots)"""
        self.reacsnap: Dataset
        """hdf5 Dataset 'Snapshots/reactions' (reactions snapshots)"""
        self.reacsnapsaved: Dataset
        """hdf5 Dataset 'Snapshots/reactions_saved' (were reactions snapshotted?)"""
        self._snapsized: bool
        """If True, file space have been correctly sized for storing snapshots"""
        self.maps: Group
        """hdf5 Group 'Maps' for storing maps statistics"""
        self.currentcol: int
        """Cuurent column to save data"""
        # Logging data
        self._init_log: bool = False
        """flag indicating if logging into the writer is initialized"""
        self.maxlog: int
        """maximum log line to be recorded in file"""
        self.dlog: int
        """log line increment to add in file when space goes missing"""
        self.logging: Group
        """hdf5 Group 'Logging' for storing log data"""
        self.logcount: Dataset
        """hdf5 Dataset 'Logging/count' (number of recorded log lines)"""
        self.logs: Dataset
        """hdf5 Dataset 'Logging/logs' (recorded log lines)"""

    def init_log(self, maxlog: int) -> None:
        """
        Init logging interface to hdf5 file

        @param maxlog: (initial) maximum log line to reserve in hdf5 file
            if more space is needed, maxlog more lines will be reserved.
            unused lines will be removed at run end.
        @type maxlog: int
        """
        self.maxlog = maxlog
        self.dlog = maxlog
        size = MPI_STATUS.size
        self.logging = self.h5file.create_group("Logging")
        self.logcount = self.logging.create_dataset(
            "count", (size,), fillvalue=0, dtype="int32"
        )
        self.logs = self.logging.create_dataset(
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
        MPI_GATE.register_function("addlog", self.add_log_line)
        self._init_log = True

    def add_log_line(self) -> None:
        """Reserve more lines for logging
        (function intended to be called as a sync. op. of MPI gate)"""
        self.maxlog = self.maxlog + self.dlog
        try:
            self.logs.resize(self.maxlog, axis=1)
        except ValueError:
            self.maxlog = self.maxlog - self.dlog
            self._init_log = False

    def write_log(self, level: int, time: str, runtime: float, msg: str) -> None:
        """
        Write a log line if the file

        @param level: logging level number
        @type level: int
        @param time: time at logging event
        @type time: str
        @param runtime: runtime at logging event
        @type runtime: float
        @param msg: logged message
        @type msg: str
        """
        try:
            if self._init_log:
                rank = MPI_STATUS.rank
                col = self.logcount[rank]
                try:
                    self.logs[rank, col] = (level, time, runtime, msg[: self.maxstrlen])
                except ValueError:
                    # No more room in log, stop logging
                    self._init_log = False
                self.logcount[rank] = col + 1
                if (self.maxlog - col) < self.lengrow:
                    try:
                        MPI_GATE.close("addlog")
                    except ValueError:
                        # run out of room for log outside the gate, cannot sync withn other threads
                        self._init_log = False
        except OSError:
            # Big problem.... stop logging... probably an overflow of error messages...
            self._init_log = False

    def close_log(self) -> None:
        """Stop logging in file"""
        cutline = MPI_STATUS.max(self.logcount[MPI_STATUS.rank])
        self.logs.resize(cutline, axis=1)
        self._init_log = False

    def init_stat(
        self,
        datanames: List[str],
        mapnames: List[str],
        params: Param,
        ruleparam: Dict[str, Any],
        comment: str,
        nbcol: int,
    ) -> None:
        """
        Initialize statistics recording

        @param datanames: list of data field names
        @type datanames: List[str]
        @param mapnames: list of map field names
        @type mapnames: List[str]
        @param params: run parameters
        @type params: Param
        @param ruleparam: rules parameters
        @type ruleparam: Dict[str, Any]
        @param comment: run comment
        @type comment: str
        @param nbcol: initial number of columns to reserve in file
        @type nbcol: int
        """
        size = MPI_STATUS.size
        self.nbcol = nbcol + self.lengrow
        self.dcol = nbcol
        self.run = self.h5file.create_group("Run")
        self.run.attrs["version"] = __version__
        self.run.attrs["hostname"] = MPI_STATUS.hostname
        self.run.attrs["date"] = MPI_STATUS.starttime
        self.run.attrs["threads"] = MPI_STATUS.size
        self.run.attrs["comment"] = comment
        self.params = self.h5file.create_group("Parameters")
        self.dict_as_attr(self.params, params.asdict())
        self.statparam = self.params.create_group("Stats")
        self.multiread_as_attr(self.statparam, params.statparam)
        self.mapparam = self.params.create_group("Maps")
        self.multiread_as_attr(self.mapparam, params.mapsparam)
        self.ruleparam = self.params.create_group("Rules")
        self.dict_as_attr(self.ruleparam, ruleparam)
        self.dataset = self.h5file.create_group("Dataset")
        self.dataset.attrs["datanames"] = datanames
        self.data = self.dataset.create_dataset(
            "results",
            (size, len(datanames), self.nbcol),
            maxshape=(size, len(datanames), None),
            fillvalue=np.nan,
        )
        self.end = self.dataset.create_dataset(
            "end",
            (size,),
            dtype=[
                ("num", "int32"),
                ("message", string_dtype(length=self.maxstrlen)),
                ("runtime", "float32"),
            ],
        )
        self.snapshots = self.h5file.create_group("Snapshots")
        self.timesnap = self.snapshots.create_dataset(
            "time", (size, 1), maxshape=(size, None), dtype="float32",
        )
        self.compsnap = self.snapshots.create_dataset(
            "compounds",
            (size, 1, 1),
            maxshape=(size, None, None),
            dtype=[("name", string_dtype(length=self.maxstrlen)), ("pop", "int32")],
        )
        self.reacsnap = self.snapshots.create_dataset(
            "reactions",
            (size, 1, 1),
            maxshape=(size, None, None),
            dtype=[
                ("name", string_dtype(length=self.maxstrlen)),
                ("const", "float32"),
                ("rate", "float32"),
            ],
        )
        self.reacsnapsaved = self.snapshots.create_dataset(
            "reactions_saved", (size, 1), maxshape=(size, None), dtype=bool,
        )
        self._snapsized = False
        self.maps = self.h5file.create_group("Maps")
        for name in mapnames:
            self.maps.create_dataset(
                name,
                (size, 1, self.nbcol + 1),
                maxshape=(size, None, None),
                fillvalue=np.nan,
                dtype="float32",
            )
        self.currentcol = 0
        MPI_GATE.register_function("addcol", self.add_col)
        self._init_stat = True

    def test_initialized(self) -> None:
        """
        Test if the file was intialized for storing statistics

        @raise: InternalError if not initialized
        """
        if not self._init_stat:
            raise InternalError("Attempt to write in HDF5 file before intialization")

    def add_col(self) -> None:
        """Reserve additional columns for storing results
        (function intended to be called as a sync. op. of MPI gate)"""
        self.nbcol = self.nbcol + self.dcol
        self.data_resize(self.nbcol)

    def data_resize(self, nbcol: float = invalidfloat) -> None:
        """
        Resize data size of hdf5 datasets

        @param nbcol: number of column to resize to (if invalid, cutout empty lines)
            (Default value = invalidfloat)
        @type nbcol: float
        """
        if not isvalid(nbcol):
            nbcol = MPI_STATUS.max(self.currentcol)
        self.data.resize(nbcol, axis=2)
        for datamap in self.maps.values():
            datamap.resize(nbcol + 1, axis=2)

    def add_map(self, name: str, categories: List[float], data: Dict[float, List[float]]) -> None:
        self.test_initialized()
        mapsize = len(categories)
        self.maps[name].resize(mapsize, axis=1)
        self.maps[name][MPI_STATUS.rank, :, 0] = categories
        for catnum, cat in enumerate(categories):
            try:
                length = len(data[cat])
                self.maps[name][MPI_STATUS.rank, catnum, 1 : length + 1] = data[cat]
            except KeyError:
                pass  # No problem, some categories may have been reached by only some processes

    def snapsize(self, maxcomp: int, maxreac: int, maxsnap: int) -> None:
        self.test_initialized()
        self.timesnap.resize((MPI_STATUS.size, maxsnap))
        self.compsnap.resize((MPI_STATUS.size, maxsnap, maxcomp))
        self.reacsnap.resize((MPI_STATUS.size, maxsnap, maxreac))
        self.reacsnapsaved.resize((MPI_STATUS.size, maxsnap))
        self._snapsized = True

    def close(self) -> None:
        self.data_resize()
        self.close_log()
        self.run.attrs["end"] = datetime.now().strftime(self.timeformat)
        self._init_log = False
        self._init_stat = False
        self.h5file.close()

    def add_data(self, result: List[float]) -> None:
        self.test_initialized()
        try:
            self.data[MPI_STATUS.rank, :, self.currentcol] = result
            self.currentcol += 1
            if (self.nbcol - self.currentcol) < self.lengrow:
                MPI_GATE.close("addcol")
        except ValueError:
            raise InternalError(
                f"No more space in file for #{MPI_STATUS.rank} at column {self.currentcol}"
            )

    def add_end(self, ending: Finished, time: float) -> None:
        self.test_initialized()
        self.end[MPI_STATUS.rank] = (
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
            self.timesnap[MPI_STATUS.rank, col] = time
            self.reacsnapsaved[MPI_STATUS.rank, col] = len(reaclist) > 0
            for line, data in enumerate(complist.items()):
                self.compsnap[MPI_STATUS.rank, col, line] = (
                    data[0][: self.maxstrlen],
                    data[1],
                )
            for line, (name, (const, rate)) in enumerate(reaclist.items()):
                self.reacsnap[MPI_STATUS.rank, col, line] = (
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
