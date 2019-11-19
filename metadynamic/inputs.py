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

from json import load, JSONDecodeError
from typing import List, Dict
from dataclasses import dataclass, field

from metadynamic.ends import BadFile, FileNotFound, BadJSON


@dataclass
class Readerclass:
    _default_section: str = "Parameter"

    def __post_init__(self) -> None:
        pass

    @classmethod
    def readfile(
        cls, filename: str, section: str = "", checktype: bool = True
    ) -> "Readerclass":
        """Return a SysParam object, updated by the data from filename"""
        if section == "":
            section = cls._default_section
        try:
            with open(filename) as json_data:
                parameters = load(json_data)[section]
        except FileNotFoundError:
            raise FileNotFound(f"Unknown file {filename}")
        except JSONDecodeError as jerr:
            raise BadJSON(f"({jerr})")
        # Validate file entries
        err = ""
        list_param = cls.list_param()
        for key, val in parameters.items():
            if key not in list_param.keys():
                err += f"{key} parameter unknown.\n"
            elif checktype:
                if not isinstance(val, list_param[key]):
                    err += f"{key} parameter should be of type {list_param[key]}, not {type(val)}\n"
        if err != "":
            raise BadFile(err)
        # OK, initialize data
        return cls(**parameters)

    @classmethod
    def list_param(cls) -> Dict[str, type]:
        return {
            key: val.__origin__ if hasattr(val, "__origin__") else val
            for key, val in cls.__annotations__.items()
            if key[0] != "_"
        }

    def set_param(self, **kwd) -> None:
        list_param = self.list_param()
        for key, val in kwd.items():
            setattr(self, key, list_param[key](val))
        self.__post_init__()


@dataclass
class SysParam(Readerclass):
    _default_section = "System"
    conc: float = 0.1
    tend: float = 1.0
    tstep: float = 0.01
    rtlim: float = 900.0  # Limit to 15min runtime
    maxsteps: int = 10000
    ptot: int = field(init=False)
    vol: float = field(init=False)
    seed: int = 0
    save: List[str] = field(default_factory=list)
    init: Dict[str, int] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.ptot = sum([pop * len(comp) for comp, pop in self.init.items()])
        self.vol = self.ptot / self.conc


@dataclass
class RunParam(Readerclass):
    _default_section = "Run"
    autoclean: bool = True
    minprob: float = 1e-10
    dropmode: str = ""
    gcperio: bool = True
    nbthread: int = 1
    context: str = "fork"
    consts: Dict[str, float] = field(default_factory=dict)
    altconsts: Dict[str, float] = field(default_factory=dict)
    catconsts: Dict[str, float] = field(default_factory=dict)
