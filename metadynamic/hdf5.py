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

from typing import Dict, Any, List, Tuple
from h5py import File, Group, Dataset
from mpi4py import MPI


class ResultWriter:
    def __init__(self, filename: str, datanames: List[str], nbcol: int) -> None:
        self.h5file: File = File(filename, "w", driver="mpio", comm=MPI.COMM_WORLD)
        self.size = MPI.COMM_WORLD.size
        self.params: Group = self.h5file.create_group("Parameters")
        self.dataset: Group = self.h5file.create_group("Dataset")
        self.dataset.attrs["datanames"] = datanames
        self.data: Dataset = self.dataset.create_dataset(
            "Results", (self.size, len(datanames), nbcol)
        )
        self.h5file.create_group("Snapshots")
        self.snapshots: List[Group] = [
            self.h5file["Snapshots"].create_group(f"#{num}") for num in range(self.size)
        ]
        self._currentcol = 0

    def close(self) -> None:
        self.h5file.close()

    @property
    def rank(self) -> int:
        return int(MPI.COMM_WORLD.rank)  # cast is for mypy...

    def add_parameter(self, params: Dict[str, Any], name: str = "") -> None:
        if self.rank == 0:  # Avoid multiple write (?)
            self.dict_as_attr(self.params, params, name)

    def add_data(self, result: List[float]) -> None:
        try:
            self.data[self.rank, :, self._currentcol] = result
            self._currentcol += 1
        except ValueError:
            raise ValueError(
                f"No more space in file for #{self.rank} at column {self._currentcol}"
            )

    def add_snapshot(
        self,
        complist: Dict[str, int],
        reaclist: Dict[str, Tuple[float, float]],
        time: float,
    ) -> None:
        group: Group = self.snapshots[self.rank].create_group(f"t={time}")
        compgroup: Group = group.create_group("Compounds")
        reacgroup: Group = group.create_group("Reactions")
        self.dict_as_attr(compgroup, complist)
        self.dict_as_attr(reacgroup, reaclist)

    def dict_as_attr(self, group: Group, datas: Dict[str, Any], name: str = "") -> None:
        for key, val in datas.items():
            if name:
                key = f"{name}->{key}"
            try:
                group.attrs[key] = val
            except TypeError:
                self.dict_as_attr(group, val, name=key)
