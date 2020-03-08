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


from numpy import array, ndarray, nanmean, nanstd, convolve, ones, meshgrid, nan_to_num
from pandas import DataFrame
from h5py import File, Group, Dataset
from graphviz import Digraph
from typing import List, Dict, Tuple, Any

from metadynamic.inval import invalidstr, invalidint, isvalid, invalidfloat
from metadynamic.inputs import StatParam, MapParam, Param, RulesetParam
from metadynamic.caster import Caster
from metadynamic.json2dot import Data2dot
from metadynamic.system import CRN


comp_cast = Caster(Dict[str, int])
reac_cast = Caster(Dict[str, List[float]])


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
        self.reacsnapsave: Dataset = self.snapshots["reactions_saved"]
        self.maps: Group = self.h5file["Maps"]
        self.mapnames = self.maps.keys()
        self.logging: Group = self.h5file["Logging"]
        self.logcount: Dataset = self.logging["count"]
        self.logs: Dataset = self.logging["logs"]
        self.size = self.run.attrs["threads"]

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

    def getcrn(self, comp: Dict[str, int]) -> CRN:
        param = self.parameters
        param.set_param(init=comp)
        return CRN(param)

    def getsnap(self, num: int, step: int, parameterfile: str = "") -> Digraph:
        comp = comp_cast(self.compsnap[num, step])
        if "" in comp:
            comp.pop("")
        if self.reacsnapsave[num, step]:
            reac = reac_cast(
                {i: (j, k) for i, j, k in self.snapshots["reactions"][num, step]}
            )
        else:
            reac = reac_cast(self.getcrn(comp).reac_collect.asdict())
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

    def xyproc(
        self,
        y: str = "ptime",
        x: str = "time",
        method: str = "m",
        xmethod: str = "m",
        meanlength: int = invalidint,
    ) -> Tuple[ndarray, ndarray]:
        x = self.get(field=x, method=xmethod, meanlength=meanlength)
        y = array(
            [
                self.get(field=y, method=f"p{proc}", meanlength=meanlength)
                for proc in range(self.size)
            ]
        ).T
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
        err = self.get(field=y, method=f"s", meanlength=meanlength) * delta
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
        try:
            z = nan_to_num(
                self.getmap(field=field, method=method, meanlength=meanlength),
                nan=nanval,
                posinf=posinf if isvalid(posinf) else None,
                neginf=neginf if isvalid(neginf) else None,
            )
        except TypeError:
            # older numpy version cannot set replacement values for nan/posinf/neginf
            z = nan_to_num(
                self.getmap(field=field, method=method, meanlength=meanlength)
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
    def ruleset(self) -> RulesetParam:
        ruleset = dict(self.params["Rules"].attrs)
        res = ruleset.copy()
        for key, val in ruleset.items():
            if "->" in key:
                key1, key2, key3 = key.split("->")
                if key1 not in res:
                    res[key1] = {key2: {key3: val}}
                elif key2 not in res[key1]:
                    res[key1][key2] = {key3: val}
                else:
                    res[key1][key2][key3] = val
                res.pop(key)
        return RulesetParam.readdict(res)

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
