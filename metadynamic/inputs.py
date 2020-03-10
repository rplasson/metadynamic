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

from json import load, dump, JSONDecodeError
from typing import List, Dict, TypeVar, Type, Any
from dataclasses import dataclass, field
from metadynamic.caster import Caster

from metadynamic.ends import BadFile, FileNotFound, BadJSON

R = TypeVar("R", bound="Readerclass")


class LockedError(Exception):
    pass


class Castreader(Caster):
    def __call__(self, value: Any) -> Any:
        if issubclass(self.dest, Readerclass):
            return self.dest.readdict(value)
        if self.dest is dict and issubclass(self.args[1].dest, Readerclass):
            return self.args[1].dest.multipledict(value)
        return super().__call__(value)


@dataclass
class Readerclass:
    _list_param: Dict[str, Any] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        pass

    @staticmethod
    def _fromfile(filename: str) -> Dict[str, Any]:
        if filename == "":
            return {}
        try:
            with open(filename) as json_data:
                parameters: Dict[str, Any] = load(json_data)
        except FileNotFoundError:
            raise FileNotFound(f"Unknown file {filename}")
        except JSONDecodeError as jerr:
            raise BadJSON(f"({jerr})")
        return parameters

    @classmethod
    def readdict(
        cls: Type[R],
        parameters: Dict[str, Any],
        checktype: bool = True,
        autocast: bool = True,
    ) -> R:
        """Return a Readerclass object, updated by the data from parameters dict"""
        new = cls()
        if checktype:
            new.set_checktype()
        else:
            new.unset_checktype()
        if autocast:
            new.set_autocast()
        else:
            new.unset_autocast()
        new.set_param(**parameters)
        return new

    @classmethod
    def multipledict(
        cls: Type[R],
        parameters: Dict[str, Dict[str, Any]],
        checktype: bool = True,
        autocast: bool = True,
    ) -> Dict[str, R]:
        res: Dict[str, R] = {}
        for name, params in parameters.items():
            res[name] = cls.readdict(
                parameters=params, checktype=checktype, autocast=autocast
            )
        return res

    @classmethod
    def readfile(
        cls: Type[R], filename: str, checktype: bool = True, autocast: bool = True,
    ) -> R:
        """Return a Readerclass object, updated by the data from filename"""
        return cls.readdict(
            parameters=cls._fromfile(filename), checktype=checktype, autocast=autocast
        )

    @classmethod
    def readmultiple(
        cls: Type[R], filename: str, checktype: bool = True, autocast: bool = True,
    ) -> Dict[str, R]:
        """Return a dictionnary of Readerclass object,
           each updated by specific entry from filename"""
        return cls.multipledict(
            parameters=cls._fromfile(filename), checktype=checktype, autocast=autocast
        )

    @classmethod
    def list_param(cls) -> Dict[str, Castreader]:
        if not hasattr(cls, "_list_param"):
            cls._list_param = {
                key: Castreader(val.type)
                for key, val in cls.__dataclass_fields__.items()
                if val.init
            }
        return cls._list_param

    @classmethod
    def conv_param(cls, param: str) -> Castreader:
        return cls.list_param()[param]

    def checked_items(self, key: str, val: Any) -> Any:
        err = ""
        if key not in self.list_param().keys():
            err += f"'{key}' parameter unknown.\n"
        else:
            if self.autocast:
                try:
                    val = self.conv_param(key)(val)
                except ValueError:
                    err += f"Couldn't cast '{val}' into {self.conv_param(key).dest}. "
            if self.checktype:
                if not isinstance(val, self.conv_param(key).dest):
                    err += f"{key} parameter should be of type {self.conv_param(key).dest}, not {type(val)}\n"
        if err != "":
            raise BadFile(err)
        return val

    def set_param(self, **kwd: Any) -> None:
        if self.locked:=
            raise LockedError
        for key, val in kwd.items():
            val = self.checked_items(key, val)
            setattr(self, key, val)
        self.__post_init__()

    def asdict(self) -> Dict[str, Any]:
        res = {}
        for key in self.list_param().keys():
            val = getattr(self, key)
            if isinstance(val, Readerclass):
                val = val.asdict()
            elif isinstance(val, dict):
                val = {
                    subkey: subval.asdict()
                    if isinstance(subval, Readerclass)
                    else subval
                    for subkey, subval in val.items()
                }
            res[key] = val
        return res

    def tojson(self, filename: str) -> None:
        with open(filename, "w") as out:
            dump(self.asdict(), out, indent=4)

    def lock(self) -> None:
        self._locked = True

    def unlock(self) -> None:
        self._locked = False

    @property
    def locked(self) -> bool:
        if not hasattr(self, "_locked"):
            self._locked = False
        return self._locked

    def set_autocast(self) -> None:
        self._autocast = True

    def unset_autocast(self) -> None:
        self._autocast = False

    @property
    def autocast(self) -> bool:
        if not hasattr(self, "_autocast"):
            self._autocast = True
        return self._autocast

    def set_checktype(self) -> None:
        self._checktype = True

    def unset_checktype(self) -> None:
        self._checktype = False

    @property
    def checktype(self) -> bool:
        if not hasattr(self, "_checktype"):
            self._checktype = True
        return self._checktype


@dataclass
class CategoryParam(Readerclass):
    func: str = ""
    descr: str = ""


@dataclass
class PropertyParam(Readerclass):
    func: str = ""
    descr: str = ""


@dataclass
class RuleParam(Readerclass):
    reactants: List[str] = field(default_factory=list)
    builder_func: str = ""
    builder_const: str = "kinvar"
    builder_variant: str = "novariant"
    descr: str = ""


@dataclass
class RulesetParam(Readerclass):
    rulemodel: str = "metadynamic.polymers"  # rulemodel to be used
    categories: Dict[str, CategoryParam] = field(default_factory=dict)  # categories
    properties: Dict[str, PropertyParam] = field(default_factory=dict)  # properties
    rules: Dict[str, RuleParam] = field(default_factory=dict)  # rules


@dataclass
class StatParam(Readerclass):
    prop: str = "count"
    weight: str = "count"
    method: str = "m"
    full: bool = False
    collection: str = "compounds"


@dataclass
class MapParam(Readerclass):
    prop: str = "count"
    weight: str = "count"
    sort: str = "count"
    method: str = "+"
    full: bool = False
    collection: str = "compounds"


@dataclass
class Param(Readerclass):
    # chemical
    conc: float = 0.1  # Concentration
    ptot: int = field(init=False)
    vol: float = field(init=False)
    init: Dict[str, int] = field(default_factory=dict)  # initial concentrations
    rulemodel: str = "polymers-ruleset.json"  # rule model to be used
    consts: Dict[str, List[float]] = field(default_factory=dict)  # kinetic constants
    # simulation
    tend: float = 1.0  # final simulation time
    tstep: float = 0.01  # timestep
    sstep: float = -1  # snapshot step (if <0: only at end, if = 0: at each step)
    rtlim: float = 900.0  # Limit runtime
    maxsteps: int = 10000  # maximum time steps
    # System
    autoclean: bool = True  # Perform a periodic cleaning of probabilities if True
    dropmode: str = ""  # drop mode (can be 'keep', 'drop' or 'soft')
    gcperio: bool = True  # if True, only call garbage collector at each timestep.
    endbarrier: float = 0.01  # If non zero, final threads will wait in idle loops of corresponding values
    maxmem: int = 2048  # Max memory (in Mb).
    # IO
    save: List[str] = field(
        default_factory=list
    )  # list of compounds to be saved at each time step
    stat: str = ""  # json filename describing statistics
    statparam: Dict[str, StatParam] = field(init=False)
    maps: str = ""  # json filename describing stat maps
    mapsparam: Dict[str, MapParam] = field(init=False)
    snapshot: str = ""  # filename for final snapshot
    printsnap: str = "pdf"  # filetype of snapshots
    hdf5: str = ""  # filename for hdf5 file
    store_snapreac: bool = False  # Store reaction snapshots? (can take lots of time)
    maxstrlen: int = 256  # max string length to be stored in hdf5
    lengrow: int = 10  # number of length left before requesting a resize
    maxlog: int = 100  # max log lines per process to be saved
    timeformat: str = "%H:%M:%S, %d/%m/%y"  # timeformat used in log files

    def __post_init__(self) -> None:
        self.ptot = sum([pop * len(comp) for comp, pop in self.init.items()])
        self.vol = self.ptot / self.conc
        self.statparam = StatParam.readmultiple(self.stat)
        self.mapsparam = MapParam.readmultiple(self.maps)


@dataclass
class DotParam(Readerclass):
    # type
    binode: bool = False
    # Graph rendering
    margin: float = 0.0
    concentrate: bool = False
    maxsize: float = 10.0
    # compounds
    min_fontsize: float = 1.0
    max_fontsize: float = 20.0
    min_c_width: float = 0.05
    max_c_width: float = 0.5
    c_penwidth: float = 0.5
    c_margin: float = 0.01
    c_color: str = "blue"
    c_powerscale: float = 0.5
    font_powerscale: float = 1.0
    # reactions
    min_r_width: float = 0.01
    max_r_width: float = 0.1
    r_color: str = "red"
    r_powerscale: float = 1.0
    # flows
    min_f_width: float = 0.001
    max_f_width: float = 5.0
    f_color: str = "black"
    cutoff: float = 0.05
    f_powerscale: float = 1.0
