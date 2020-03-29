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

from pandas import DataFrame
from h5py import File, Group, Dataset
from graphviz import Digraph
from typing import List, Dict, Tuple, Any, Callable

from metadynamic.inval import invalidstr, invalidint, isvalid, invalidfloat
from metadynamic.inputs import StatParam, MapParam, Param, RulesetParam
from metadynamic.caster import Caster
from metadynamic.network import Data2dot
from metadynamic.chemical import Crn


comp_cast: Callable[[Any], Dict[str, int]] = Caster(Dict[str, int])
reac_cast: Callable[[Any], Dict[str, List[float]]] = Caster(Dict[str, List[float]])


class ResultReader:
    """Interface for a hdf5 result file"""

    def __init__(self, filename: str) -> None:
        """Connect the reader to a .hdf5 result file

        @param filename:str name of the .hdf5 file
        """
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

    def __getitem__(self, field: str) -> np.ndarray:
        if field in self.datanames:
            return self.get(field)
        if field in self.mapnames:
            return self.getmap(field)
        raise KeyError

    def _loc(self, field: str) -> int:
        """Returns the index number of the given field name.

        @param field: Name of the field
        @type field: str
        @param field: str: 
        @return: field index number
        @rtype: int

        """
        try:
            return self.datanames.index(field)
        except ValueError:
            raise ValueError(f"{field} is not a recorded stat name")

    def get(
        self, field: str = invalidstr, method: str = "m", meanlength: int = invalidint,
    ) -> np.ndarray:
        """Returns the saved result data 'field'. If no field is provided, the full
        dataset will be returned.
        
        If 'method' is set to:
            - 'm', the average values over all processes is returned (default).
            - 's', the standard deviation  of the values ove all processes is returned.
            - '+X', where X is a numerical value, the returned value is 'mean + X·std'
            - '-X', where X is a numerical value, the returned value is 'mean - X·std'
            - '*', all individual values from each process is returned
            - 'pX', where X is an integer, the values from process number X is returned.
            - 'sum', the values of each process are summed together.
        
        If meanlength is set, a running mean of the corresponding length is returned.

        @param field: data field  (Default value = invalidstr)
        @param method: method for processing data over processes  (Default value = 'm')
        @param meanlength: running mean length  (Default value = invalidint)
        @return: the processed set of data
        @rtype: ndarray

        """
        loc = self._loc(field) if isvalid(field) else slice(None, None, None)
        if method == "*":
            return self.data[:, loc]
        if method[0] == "p":
            res = self.data[int(method[1:]), loc]
        elif method[0] == "sum":
            res = np.nansum(self.data[:, loc, :], axis=0)
        elif method[0] in ["m", "s", "+", "-"]:
            mean = np.nanmean(self.data[:, loc, :], axis=0)
            if method == "m":
                res = mean
            else:
                std = np.nanstd(self.data[:, loc, :], axis=0)
                if method == "s":
                    res = std
                else:
                    res = mean + float(method) * std
        else:
            raise ValueError(f"'method'={method} is invalid")
        return (
            # Running mean from https://stackoverflow.com/questions/13728392/moving-average-or-running-mean
            np.convolve(res, np.ones((meanlength,)) / meanlength, mode="valid")
            if isvalid(meanlength)
            else res
        )

    def getmap(
        self, field: str, method: str = "m", meanlength: int = invalidint
    ) -> np.ndarray:
        """Returns the saved property map data 'field'. If no field is provided, the full
        dataset will be returned.
        
        If 'method' is set to:
            - 'm', the average values over all processes is returned (default).
            - 's', the standard deviation  of the values ove all processes is returned.
            - '+X', where X is a numerical value, the returned value is 'mean + X·std'
            - '-X', where X is a numerical value, the returned value is 'mean - X·std'
            - '*', all individual values from each process is returned
            - 'pX', where X is an integer, the values from process number X is returned.
        
        If meanlength is set, a running mean of the corresponding length is returned.

        @param field: data field  (Default value = invalidstr)
        @param method: method for processing data over processes  (Default value = 'm')
        @param meanlength: running mean length  (Default value = invalidint)
        @param field: str: 
        @param method: str:  (Default value = "m")
        @param meanlength: int:  (Default value = invalidint)
        @return: the processed set of data
        @rtype: ndarray

        """
        try:
            data = self.maps[field][:, :, 1:]
        except KeyError:
            raise ValueError(f"{field} is not a recorded map name")
        if method == "*":
            return data
        if method[0] == "p":
            res = data[int(method[1:])]
        elif method[0] in ["m", "s", "+", "-"]:
            mean = np.nanmean(data, axis=0)
            if method == "m":
                res = mean
            else:
                std = np.nanstd(data, axis=0)
                if method == "s":
                    res = std
                else:
                    res = mean + float(method) * std
        else:
            raise ValueError(f"'method'={method} is invalid")
        return (
            # Running mean from https://stackoverflow.com/questions/13728392/moving-average-or-running-mean
            np.convolve(res, np.ones((meanlength,)) / meanlength, mode="valid")
            if isvalid(meanlength)
            else res
        )

    def getcrn(self, comp: Dict[str, int]) -> Crn:
        """Generate a Chemical Reaction Network from a dictionnary
        as {compound_name : population}

        @param comp: compounds and population
        @type comp: Dict[str: int]
        @param comp: Dict[str: 
        @param int]: 

        """
        param = self.parameters
        param.set_param(init=comp)
        return Crn(param)

    def getsnap_comp(self, num: int, step: int) -> Dict[str, int]:
        """

        @param num: int:
        @param step: int:

        """
        comp = comp_cast(self.compsnap[num, step])
        if "" in comp:
            comp.pop("")
        return comp

    def getsnap(self, num: int, step: int, parameterfile: str = "") -> Digraph:
        """

        @param num: int:
        @param step: int:
        @param parameterfile: str:  (Default value = "")

        """
        comp = self.getsnap_comp(num, step)
        if self.reacsnapsave[num, step]:
            reac = reac_cast(
                {i: (j, k) for i, j, k in self.snapshots["reactions"][num, step]}
            )
        else:
            reac = reac_cast(self.getcrn(comp).reac_collect.asdict())
        if "" in reac:
            reac.pop("")
        return Data2dot(comp, reac, parameterfile).crn.dot

    def categories(self, field: str) -> np.ndarray:
        """

        @param field: str:

        """
        return self.maps[field][0, :, 0]

    def ending(self, num: int) -> Tuple[int, str, float]:
        """

        @param num: int:

        """
        endnum, message, time = self.end[num]
        return endnum, message.decode(), time

    def endmsg(self, num: int) -> str:
        """

        @param num: int:

        """
        endnum, message, time = self.ending(num)
        return f"#{num}: ending n°{endnum} at runtime t={time}s; {message}"

    def table(
        self, maps: str = invalidstr, method: str = "m", meanlength: int = invalidint
    ) -> DataFrame:
        """

        @param maps: str:  (Default value = invalidstr)
        @param method: str:  (Default value = "m")
        @param meanlength: int:  (Default value = invalidint)

        """
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
    ) -> Tuple[np.ndarray, np.ndarray]:
        """

        @param y: str:  (Default value = "ptime")
        @param x: str:  (Default value = "time")
        @param method: str:  (Default value = "m")
        @param xmethod: str:  (Default value = "m")
        @param meanlength: int:  (Default value = invalidint)

        """
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
    ) -> Tuple[np.ndarray, np.ndarray]:
        """

        @param y: str:  (Default value = "ptime")
        @param x: str:  (Default value = "time")
        @param method: str:  (Default value = "m")
        @param xmethod: str:  (Default value = "m")
        @param meanlength: int:  (Default value = invalidint)

        """
        x = self.get(field=x, method=xmethod, meanlength=meanlength)
        y = np.array(
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
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """

        @param y: str:  (Default value = "ptime")
        @param x: str:  (Default value = "time")
        @param delta: float:  (Default value = 1.0)
        @param meanlength: int:  (Default value = invalidint)

        """
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
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """

        @param y: str:  (Default value = "ptime")
        @param x: str:  (Default value = "time")
        @param delta: float:  (Default value = 1.0)
        @param meanlength: int:  (Default value = invalidint)

        """
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
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """

        @param field: str:
        @param method: str:  (Default value = "m")
        @param tmethod: str:  (Default value = "m")
        @param meanlength: int:  (Default value = invalidint)
        @param nanval: float:  (Default value = 0.0)
        @param posinf: float:  (Default value = invalidfloat)
        @param neginf: float:  (Default value = invalidfloat)

        """
        time = self.get(field="time", method=tmethod, meanlength=meanlength)
        categories = self.categories(field)
        x, y = np.meshgrid(time, categories)
        try:
            z = np.nan_to_num(
                self.getmap(field=field, method=method, meanlength=meanlength),
                nan=nanval,
                posinf=posinf if isvalid(posinf) else None,
                neginf=neginf if isvalid(neginf) else None,
            )
        except TypeError:
            # older numpy version cannot set replacement values for nan/posinf/neginf
            z = np.nan_to_num(
                self.getmap(field=field, method=method, meanlength=meanlength)
            )
        return x, y, z

    @property
    def printinfo(self) -> str:
        """ """
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
        """ """
        return dict(self.run.attrs)

    @property
    def parameters(self) -> Param:
        """ """
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
        """ """
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
        """ """
        return StatParam.multipledict(
            {key: dict(val.attrs) for key, val in self.params["Stats"].items()}
        )

    @property
    def mapparam(self) -> Dict[str, MapParam]:
        """ """
        return MapParam.multipledict(
            {key: dict(val.attrs) for key, val in self.params["Maps"].items()}
        )

    @property
    def fulllog(self) -> np.ndarray:
        """ """
        # add simple system for log display/search (depending on level etc) ?
        maxcol = max(self.logcount)
        return self.logs[:, :maxcol]
