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
from dataclasses import dataclass, field, asdict

from metadynamic.ends import BadFile, FileNotFound, BadJSON
from metadynamic.logger import Logged

R = TypeVar("R", bound="Readerclass")


class LockedError(Exception):
    pass


@dataclass
class Readerclass(Logged):
    _default_section: str = field(repr=False, default="")
    _autocast: bool = field(repr=False, default=True)
    _checktype: bool = field(repr=False, default=True)

    def __post_init__(self) -> None:
        pass

    @classmethod
    def readfile(
        cls: Type[R],
        filename: str,
        section: str = "",
        checktype: bool = True,
        autocast: bool = True,
    ) -> R:
        """Return a SysParam object, updated by the data from filename"""
        if section == "":
            section = cls._default_section
        try:
            with open(filename) as json_data:
                parameters = load(json_data)
                if section:
                    parameters = parameters[section]
        except FileNotFoundError:
            raise FileNotFound(f"Unknown file {filename}")
        except JSONDecodeError as jerr:
            raise BadJSON(f"({jerr})")
        # class memders???
        cls._checktype = checktype
        cls._autocast = autocast
        # Validate file entries
        for key, val in parameters.items():
            val = cls.checked_items(key, val)
        # OK, initialize data
        return cls(**parameters)

    @classmethod
    def checked_items(cls, key: str, val: Any) -> Any:
        err = ""
        if key not in cls.list_param().keys():
            err += f"'{key}' parameter unknown.\n"
        else:
            if cls._autocast:
                try:
                    val = cls.list_param()[key](val)
                except ValueError:
                    err += f"Couldn't cast {val} into {cls.list_param()[key]}"
            if cls._checktype:
                if not isinstance(val, cls.list_param()[key]):
                    err += f"{key} parameter should be of type {cls.list_param()[key]}, not {type(val)}\n"
        if err != "":
            raise BadFile(err)
        return val

    @classmethod
    def list_param(cls) -> Dict[str, type]:
        if not hasattr(cls, "_list_param"):
            cls._list_param = {
                key: val.__origin__ if hasattr(val, "__origin__") else val
                for key, val in cls.__annotations__.items()
                if key[0] != "_"
            }
        return cls._list_param

    def set_param(self, **kwd) -> None:
        if self.locked:
            raise LockedError
        for key, val in kwd.items():
            val = self.checked_items(key, val)
            setattr(self, key, val)
        self.__post_init__()

    def asdict(self) -> Dict[str, Any]:
        return {key: getattr(self, key) for key in self.list_param().keys()}

    def tojson(self, filename: str) -> None:
        with open(filename, "w") as out:
            dump(self.asdict(), out)

    def lock(self) -> None:
        self._locked = True

    def unlock(self) -> None:
        self._locked = False

    @property
    def locked(self):
        if not hasattr(self, "_locked"):
            self._locked = False
        return self._locked


@dataclass
class Param(Readerclass):
    # chemical
    conc: float = 0.1  # Concentration
    ptot: int = field(init=False)
    vol: float = field(init=False)
    init: Dict[str, int] = field(default_factory=dict)  # initial concentrations
    rulemodel: str = "metadynamic.polymers"  # rule model to be used
    consts: Dict[str, List[float]] = field(default_factory=dict)  # kinetic constants
    # simulation
    tend: float = 1.0  # final simulation time
    tstep: float = 0.01  # timestep
    sstep: float = -1  # snapshot step (if <0: only at end, if = 0: at each step)
    rtlim: float = 900.0  # Limit runtime
    maxsteps: int = 10000  # maximum time steps
    seed: int = 0  # random initial seed (now ignored)
    # System
    nbthread: int = 1  # number of thread (-1 is to use as many threads as detected cores)
    autoclean: bool = True  # Perform a periodic cleaning of probabilities if Treu
    dropmode: str = ""  # drop mode (can be 'keep', 'drop' or 'soft')
    gcperio: bool = True  # if True, only call garbage collector at each timestep.
    context: str = "fork"  # thread context to be used (now, only "fork" is implemented)
    # IO
    save: List[str] = field(
        default_factory=list
    )  # list of compounds to be saved at each time step
    snapshot: str = ""  # filename for final snapshot
    printsnap: str = "pdf"  # filetype of snapshots
    hdf5: str = ""  # filename for hdf5 file

    def __post_init__(self) -> None:
        self.ptot = sum([pop * len(comp) for comp, pop in self.init.items()])
        self.vol = self.ptot / self.conc


@dataclass
class DotParam(Readerclass):
    # type
    binode: bool = False
    # compounds
    min_fontsize: int = 20
    max_fontsize: int = 200
    min_c_width: float = 0.5
    max_c_width: float = 5.0
    c_color: str = "blue"
    c_powerscale: float = 0.5
    font_powerscale: float = 1.0
    # reactions
    min_r_width: float = 0.02
    max_r_width: float = 0.2
    r_color: str = "red"
    r_powerscale: float = 1.0
    # flows
    min_f_width: float = 0.1
    max_f_width: float = 20.0
    f_color: str = "black"
    cutoff: float = 0.05
    f_powerscale: float = 1.0
