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

"""Interface for processing data from finished runs.

It provides L{ResultReader}, a class for reading and extracting data from a .hdf5 result file.

"""

from typing import List, Dict, Tuple, Any, Callable, Optional
from tempfile import NamedTemporaryFile
from pandas import DataFrame
from h5py import File, Group, Dataset
from graphviz import Digraph

import numpy as np

from metadynamic.inval import invalidstr, invalidint, isvalid
from metadynamic.inputs import StatParam, MapParam, Param, RulesetParam, multiple_tojson
from metadynamic.caster import Caster
from metadynamic.network import Data2dot
from metadynamic.chemical import Crn

comp_cast: Callable[[Any], Dict[str, int]] = Caster(Dict[str, int])
"""Caster to a compound field"""
reac_cast: Callable[[Any], Dict[str, List[float]]] = Caster(Dict[str, List[float]])
"""Caster to a reaction field"""


class ResultReader:
    """Interface for a hdf5 result file."""

    def __init__(self, filename: str) -> None:
        """
        Create a reader connected to a .hdf5 result file.

        @param filename: name of the .hdf5 file
        @type filename: str
        """
        self.filename: str = filename
        """name of HDF5 file"""
        self.h5file: File = File(filename, "r")
        """HDF5 file"""
        self.run: Group = self.h5file["Run"]
        """'Run' group"""
        self.params: Group = self.h5file["Parameters"]
        """'Parameters' group"""
        self.dataset: Group = self.h5file["Dataset"]
        """'Dataset' group"""
        self.datanames: List[str] = list(self.dataset.attrs["datanames"])
        """list of data names"""
        self.data: Dataset = self.dataset["results"]
        """'Dataset/results' dataset"""
        self.end: Dataset = self.dataset["end"]
        """'Dataset/end' dataset"""
        self.snapshots: Group = self.h5file["Snapshots"]
        """'Snapshots' group"""
        self.timesnap: Dataset = self.snapshots["time"]
        """'Dataset/time' dataset"""
        self.compsnap: Dataset = self.snapshots["compounds"]
        """'Dataset/compounds' dataset"""
        self.reacsnap: Dataset = self.snapshots["reactions"]
        """'Dataset/reactions' dataset"""
        self.reacsnapsave: Dataset = self.snapshots["reactions_saved"]
        """'Dataset/reactions_saved' dataset"""
        self.maps: Group = self.h5file["Maps"]
        """'Maps' group"""
        self.mapnames: List[str] = list(self.maps.keys())
        """list of map names"""
        self.logging: Group = self.h5file["Logging"]
        """'Logging' group"""
        self.logcount: Dataset = self.logging["count"]
        """'Logging/count' dataset"""
        self.logs: Dataset = self.logging["logs"]
        """'Logging/logs' dataset"""
        self.size: int = self.run.attrs["threads"]
        """number of threads"""

    def __getitem__(self, field: str) -> np.ndarray:
        """
        Return either a data or a map array corresponding to 'field'.

        Valid for field in self.datanames or self.mapnames
        """
        if field in self.datanames:
            return self.get(field)
        if field in self.mapnames:
            return self.getmap(field)
        raise KeyError

    def _loc(self, field: str) -> int:
        """
        Return the index number of the given field name.

        Valid fields are listed in self.datanames

        @param field: Name of the field
        @type field: str
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
        """Return the saved result data named 'field'.

        Valid fields are listed in self.datanames If no field is provided, the full dataset will be
        returned.

        'method' can be set to:
            - 'm', the average values over all processes is returned (default).
            - 's', the standard deviation  of the values ove all processes is returned.
            - '+X', where X is a numerical value, the returned value is 'mean + X·std'
            - '-X', where X is a numerical value, the returned value is 'mean - X·std'
            - '*', all individual values from each process is returned
            - 'pX', where X is an integer, the values from process number X is returned.
            - 'sum', the values of each process are summed together.

        If meanlength is set, a running mean of the corresponding length is returned.

        @param field: data field  (Default value = invalidstr)
        @type field: str
        @param method: method for processing data over processes  (Default value = 'm')
        @type method: str
        @param meanlength: running mean length  (Default value = invalidint)
        @type method: int
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
            # Running mean taken from
            # https://stackoverflow.com/questions/13728392/moving-average-or-running-mean
            np.convolve(res, np.ones((meanlength,)) / meanlength, mode="valid")
            if isvalid(meanlength)
            else res
        )

    def getmap(
        self, field: str, method: str = "m", meanlength: int = invalidint
    ) -> np.ndarray:
        """
        Return the saved property map data named 'field'.

        Valid fields are listed in self.mapnames

        'method' can be set to:
            - 'm', the average values over all processes is returned (default).
            - 's', the standard deviation  of the values ove all processes is returned.
            - '+X', where X is a numerical value, the returned value is 'mean + X·std'
            - '-X', where X is a numerical value, the returned value is 'mean - X·std'
            - '*', all individual values from each process is returned
            - 'pX', where X is an integer, the values from process number X is returned.

        If meanlength is set, a running mean of the corresponding length is returned.

        @param field: data field  (Default value = invalidstr)
        @type field: str
        @param method: method for processing data over processes  (Default value = 'm')
        @type method: str
        @param meanlength: running mean length  (Default value = invalidint)
        @type meanlength: int
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
            # Running mean taken from
            # https://stackoverflow.com/questions/13728392/moving-average-or-running-mean
            np.convolve(res, np.ones((meanlength,)) / meanlength, mode="valid")
            if isvalid(meanlength)
            else res
        )

    def getcrn(self, comp: Dict[str, int]) -> Crn:
        """Generate a Chemical Reaction Network from a dictionnary as {compound_name : population}.

        @param comp: compounds and population
        @type comp: Dict[str: int]
        @return: the corresponding CRN
        @rtype: Crn

        """
        param = self.parameters
        param.set_param(init=comp)
        return Crn(param)

    def getsnap_comp(self, num: int, step: int) -> Dict[str, int]:
        """Return the compounds snapshots saved by thread 'num' at step 'step'.

        @param num: thread number
        @type num: int
        @param step: step number
        @type step: int
        @return: the snapshot as a dictionnary {compound name: population}
        @rtype: Dict[str, int]

        """
        comp = comp_cast(self.compsnap[num, step])
        if "" in comp:
            comp.pop("")
        return comp

    def getsnap(self, num: int, step: int, parameterfile: str = "") -> Digraph:
        """Return the snapshot saved by thread 'num' at step 'step'.

        It follows the format as described in 'parameterfile' that should point to a file readable
        as a DotParam json file.

        @param num: thread number
        @type num: int
        @param step: step number
        @type step: int
        @param parameterfile: filename of a Dotparam json file  (Default value = "")
        @return: the full snapshot as a graphviz object
        @rtype: Digraph

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
        """Get the list of available categories for th given 'field'.

        Valid fields are listed in self.mapnames

        @param field: Name of the map field.
        @type field: str
        @return: array of field categories
        @rtype: ndarray

        """
        return self.maps[field][0, :, 0]

    def ending(self, num: int) -> Tuple[int, str, float]:
        """Return the details of ending information sent by thread 'num'.

        @param num: thread number
        @type num: int
        @return: ending number, ending message, ending time
        @rtype: int, str, float

        """
        endnum, message, time = self.end[num]
        return endnum, message.decode(), time

    def endmsg(self, num: int) -> str:
        """Return formatted ending information sent by thread 'num'.

        @param num: thread number
        @type num: int
        @return: formatted ending information
        @rtype: str

        """
        endnum, message, time = self.ending(num)
        return f"#{num}: ending n°{endnum} at runtime t={time}s; {message}"

    def table(
        self, maps: str = invalidstr, method: str = "m", meanlength: int = invalidint
    ) -> DataFrame:
        """Return the requested data formatted in a pandas DataFrame.

        If maps is set, a dataframe of the corresponding map field will be returned
        Valid values are listed self.mapnames.
        else, a dataframe with all datanames will be returned

        'method' and 'meanlength' parameters are as defined in getmap and get methods.

        @param maps: map field name  (Default value = invalidstr)
        @type maps:
        @param method: method for processing data over processes  (Default value = 'm')
        @type method: str
        @param meanlength: running mean length  (Default value = invalidint)
        @type meanlength: int
        @return: the processed set of data
        @rtype: DataFrame

        """
        if isvalid(maps):
            data = self.getmap(field=maps, method=method, meanlength=meanlength)
            index = self.categories(maps)
        else:
            data = self.get(method=method, meanlength=meanlength)
            index = self.datanames
        return DataFrame(data, index=index)

    def x_y(
        self,
        y: str = "ptime",
        x: str = "time",
        method: str = "m",
        xmethod: str = "m",
        meanlength: int = invalidint,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Return two data fields as named in 'x' and 'y'.

        The returned value can be directly plotted with matplotlib.

        Valid field names for x and y are listed in self.datanames

        'method', 'methodx' and 'meanlength' parameters
        are as defined in get method.

        @param y: field name for y  (Default value = "ptime")
        @type y: str
        @param x: field name for x  (Default value = "time")
        @type x: str
        @param method: method for processing y  (Default value = 'm')
        @type method: str
        @param xmethod: method for processing x  (Default value = 'm')
        @type xmethod: str
        @param meanlength: running mean length  (Default value = invalidint)
        @type meanlength: int
        @return: two 1D data arrays for x and y
        @rtype: ndarray, ndarray

        """
        x = self.get(field=x, method=xmethod, meanlength=meanlength)
        y = self.get(field=y, method=method, meanlength=meanlength)
        return x, y

    def x_y_proc(
        self,
        y: str = "ptime",
        x: str = "time",
        xmethod: str = "m",
        meanlength: int = invalidint,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Return two data fields as named in 'x' and 'y'.

        The returned values can be directly plotted with matplotlib.

        'y' contains the individual datas obtained from each thread.

        Valid field names for x and y are listed in self.datanames

        'methodx' and 'meanlength' parameters are as defined in get method.

        @param y: field name for y  (Default value = "ptime")
        @type y: str
        @param x: field name for x  (Default value = "time")
        @type x: str
        @param xmethod: method for processing x  (Default value = 'm')
        @type xmethod: str
        @param meanlength: running mean length  (Default value = invalidint)
        @type meanlength: int
        @return: 1D data arrays for x and 2D data arrays for y
        @rtype: ndarray, ndarray

        """
        x = self.get(field=x, method=xmethod, meanlength=meanlength)
        y = np.array(
            [
                self.get(field=y, method=f"p{proc}", meanlength=meanlength)
                for proc in range(self.size)
            ]
        ).T
        return x, y

    def x_y_pm(
        self,
        y: str = "ptime",
        x: str = "time",
        delta: float = 1.0,
        meanlength: int = invalidint,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return three data fields as named in 'x' and 'y', as x, y+err, y-err.

        The returned values can be plotted with matplotlib with, e.g., fill_between.

        Valid field names for x and y are listed in self.datanames

        the error is equal to 'delta' times the standard deviation

        'meanlength' parameter is as defined in get method.

        @param y: field name for y  (Default value = "ptime")
        @type y: str
        @param x: field name for x  (Default value = "time")
        @type x: str
        @param delta: error factor  (Default value = 1.0)
        @type delta: float
        @param meanlength: running mean length  (Default value = invalidint)
        @type meanlength: int
        @return: three 1D data arrays for x, y+err, y-err
        @rtype: ndarray, ndarray, ndarray

        """
        return (
            self.get(field=x, method="m", meanlength=meanlength),
            self.get(field=y, method=f"+{delta}", meanlength=meanlength),
            self.get(field=y, method=f"-{delta}", meanlength=meanlength),
        )

    def x_y_err(
        self,
        y: str = "ptime",
        x: str = "time",
        delta: float = 1.0,
        meanlength: int = invalidint,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return three data fields as named in 'x' and 'y', as x, y, yerr.

        The returned values can be plotted with matplotlib with, e.g., errorbar.

        Valid field names for x and y are listed in self.datanames

        the error is equal to 'delta' times the standard deviation

        'meanlength' parameter is as defined in get method.

        @param y: field name for y  (Default value = "ptime")
        @type y: str
        @param x: field name for x  (Default value = "time")
        @type x: str
        @param delta: error factor  (Default value = 1.0)
        @type delta: float
        @param meanlength: running mean length  (Default value = invalidint)
        @type meanlength: int
        @return: three 1D data arrays for x, y, yerr
        @rtype: ndarray, ndarray, ndarray

        """
        x = self.get(field=x, method="m", meanlength=meanlength)
        y = self.get(field=y, method=f"m", meanlength=meanlength)
        err = self.get(field=y, method=f"s", meanlength=meanlength) * delta
        return x, y, err

    def x_y_z(
        self,
        field: str,
        method: str = "m",
        tmethod: str = "m",
        meanlength: int = invalidint,
        nanval: Optional[float] = 0.0,
        posinf: Optional[float] = None,
        neginf: Optional[float] = None,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return three 2D arrays.

        The retyurned values can be used for plotting a map with, e.g., matplotlib plot_surface, as
        X the time, Y the field category, and Z the data from 'field'

        Valid field names are listed in self.mapnames

        'method', 'tmethod' and 'meanlength' parameters
        are as defined in getmap method.

        nanval, posinf and neginf can be set for replacing
        NaN or infinites by custom values (e.g. 0.0 or None)
        (only works with recent versions of numpy)

        @param field: map field name
        @type param: str
        @param method: method for processing field  (Default value = 'm')
        @type method: str
        @param tmethod: method for processing time  (Default value = 'm')
        @type tmethod: str
        @param meanlength: running mean length  (Default value = invalidint)
        @type meanlength: int
        @param nanval: replacement value for NaN  (Default value = 0.0)
        @type nanval: float or None
        @param posinf: replacement value for +inf  (Default value = None)
        @type posinf: float or None
        @param neginf: replacement value for -inf  (Default value = None)
        @type neginf: float or None
        @return: three 2D mesh data arrays for X, Y, Z
        @rtype: ndarray, ndarray, ndarray

        """
        time = self.get(field="time", method=tmethod, meanlength=meanlength)
        categories = self.categories(field)
        x, y = np.meshgrid(time, categories)
        try:
            z = np.nan_to_num(
                self.getmap(field=field, method=method, meanlength=meanlength),
                nan=nanval,
                posinf=posinf,
                neginf=neginf,
            )
        except TypeError:
            # older numpy version cannot set replacement values for nan/posinf/neginf
            z = np.nan_to_num(
                self.getmap(field=field, method=method, meanlength=meanlength)
            )
        return x, y, z

    @property
    def printinfo(self) -> str:
        """Return formated information of the full run.

        @rtype: str

        """
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
        """Structured information of the full run.

        @rtype: dict

        """
        return dict(self.run.attrs)

    @property
    def parameters(self) -> Param:
        """Parameters used for the run.

        @rtype: Param

        """
        # Read parameters saved as HDF5 attributes
        params = dict(self.params.attrs)
        res = params.copy()
        # "unflatten" enbedded as dictionary saved as "'prekey->postkey' : val" entries
        for key, val in params.items():
            if "->" in key:
                prekey, postkey = key.split("->")
                try:
                    res[prekey][postkey] = val
                except KeyError:
                    res[prekey] = {postkey: val}
                res.pop(key)
        # recover stat and map parameters saved in self.statparam/self.mapsparam
        # into temporary files
        with NamedTemporaryFile() as tempstat, NamedTemporaryFile() as tempmap:
            # write temporary stat and map json parameter files
            multiple_tojson(tempstat.name, self.statparam)
            multiple_tojson(tempmap.name, self.mapparam)
            res["stat"] = tempstat.name
            res["maps"] = tempmap.name
            param = Param.readdict(res)
        return param

    @property
    def ruleset(self) -> RulesetParam:
        """Parameters used for setting the ruleset.

        @rtype: RulesetParam

        """
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
        """Parameters used for setting the statistics.

        @rtype: Dict[str, StatParam]

        """
        return StatParam.multipledict(
            {key: dict(val.attrs) for key, val in self.params["Stats"].items()}
        )

    @property
    def mapparam(self) -> Dict[str, MapParam]:
        """Parameters used for setting the maps.

        @rtype: Dict[str, MapParam]

        """
        return MapParam.multipledict(
            {key: dict(val.attrs) for key, val in self.params["Maps"].items()}
        )

    @property
    def fulllog(self) -> np.ndarray:
        """Array containing the full log of the run processes.

        fulllog[nb_thread, nb_log_message] gives a given message as message_level, time, runtime,
        message

        @rtype: ndarray

        """
        # add simple system for log display/search (depending on level etc) ?
        maxcol = max(self.logcount)
        return self.logs[:, :maxcol]
