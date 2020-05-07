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

"""Definition of structured json files for parameter inputs.

This module provides a generic L{Readerclass} for reading parameters by defining dataclasses that
will mirror the json file structure, as well as all the json files required for defining matadynamic
simulations.

"""

from json import load, dump, JSONDecodeError
from typing import List, Dict, TypeVar, Type, Any
from dataclasses import dataclass, field
from psutil import virtual_memory

from metadynamic.caster import Caster
from metadynamic.ends import BadFile, FileNotFound, BadJSON

R = TypeVar("R", bound="Readerclass")
"""Generic type for L{Readerclass} and  its subclasses"""


class LockedError(Exception):
    """Exception raised when attempting to modify a locked Readerclass."""


class Castreader(Caster):
    """Extend a L{Caster} for dealing with L{Readerclass}, converting them as dictionary."""

    def __call__(self, value: Any) -> Any:
        """Convert 'value' to self.dest type."""
        if issubclass(self.dest, Readerclass):
            return self.dest.readdict(value)
        if self.dest is dict and issubclass(self.args[1].dest, Readerclass):
            return self.args[1].dest.multipledict(value)
        return super().__call__(value)


@dataclass
class Readerclass:
    """Dataclass with interface for reading its data from json files."""

    _list_param: Dict[str, Castreader] = field(init=False, repr=False)
    """all class parameters as a dictionary {parameter name: parameter caster}"""
    __dataclass_fields__: Dict[str, Any] = field(init=False, repr=False)
    """list of fields types"""

    def __post_init__(self) -> None:
        """(Re)Calculate data after fields init or change."""
        self._locked: bool
        """lock flag"""
        self._autocast: bool
        """autocast flag"""
        self._checktype: bool
        """checktype flag"""

    @staticmethod
    def _fromfile(filename: str) -> Dict[str, Any]:
        """Read a json file named 'filename' as a dict.

        @param filename: name of json file
        @type filename: str
        @return: json data as dict
        @rtype: Dict[str, Any]

        """
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
        """Return a Readerclass object, updated by the data from parameters dict.

        @param parameters: parameters to be set (override the default values)
        @type parameters: Dict[str, Any]
        @param checktype: if True, an error will be raised if the parameters are given
            in a faulty type (Default value = True)
        @type checktype: bool
        @param autocast: if True, given parameters will be catsed to the correct type
            (Default value = True)
        @type autocast: bool
        @return: new object
        @rtype: Readerclass

        """
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
        """Return a dictionary of Readerclass object.

        Each one is updated by the data from a dictionary of parameters dict

        @param parameters: dictionary of parameters to be set (override the default values)
        @type parameters: Dict[str, Dict[str, Any]]
        @param checktype: if True, an error will be raised if the parameters are given
            in a faulty type (Default value = True)
        @type checktype: bool
        @param autocast: if True, given parameters will be catsed to the correct type
            (Default value = True)
        @type autocast: bool
        @return: dictionary of new objects
        @rtype: Dict[str, Readerclass]

        """
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
        """Return a Readerclass object, updated by the data from a json Param file.

        @param filename: name of json file
        @type filename: str
        @param checktype: if True, an error will be raised if the parameters are given
            in a faulty type (Default value = True)
        @type checktype: bool
        @param autocast: if True, given parameters will be catsed to the correct type
            (Default value = True)
        @type autocast: bool
        @return: new object
        @rtype: Readerclass

        """
        return cls.readdict(
            parameters=cls._fromfile(filename), checktype=checktype, autocast=autocast
        )

    @classmethod
    def readmultiple(
        cls: Type[R], filename: str, checktype: bool = True, autocast: bool = True,
    ) -> Dict[str, R]:
        """Return a dictionnary of Readerclass object.

        Each one is updated by specific entry from a json file

        @param filename: name of json file
        @type filename: str
        @param checktype: if True, an error will be raised if the parameters are given
            in a faulty type (Default value = True)
        @type checktype: bool
        @param autocast: if True, given parameters will be catsed to the correct type
            (Default value = True)
        @type autocast: bool
        @return: dictionary of new objects
        @rtype: Dict[str, Readerclass]

        """
        return cls.multipledict(
            parameters=cls._fromfile(filename), checktype=checktype, autocast=autocast
        )

    @classmethod
    def list_param(cls) -> Dict[str, Castreader]:
        """List all the parameters of the class, with a caster to their type.

        @return: dictionary {parameter name: parameter caster}
        @rtype: Dict[str, Castreader]

        """
        if not hasattr(cls, "_list_param"):
            cls._list_param = {
                key: Castreader(val.type)
                for key, val in cls.__dataclass_fields__.items()
                if val.init
            }
        return cls._list_param

    @classmethod
    def conv_param(cls, param: str) -> Castreader:
        """Return the caster to the type defined for 'param'.

        @param param: parameter name
        @type param: str
        @return caster to the parameter type
        @rtype: Castreader

        """
        return cls.list_param()[param]

    def checked_items(self, key: str, val: Any) -> Any:
        """Autocast and check 'val' to the consistent type defined for 'key'.

        The operations performed will depend on self.autocast and self.checktype flags

        @param key: name of the parameter
        @type key: str
        @param val: value to be checked
        @return: converted value consistent with 'key' type
        @raise BadFile: raised if 'val' cannot be casted, or if 'val' if not of the correct type
        (in case of an autocast set to False, or of a faulty Caster

        """
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
                    err += f"{key} parameter should be of type {self.conv_param(key).dest}, "
                    err += f"not {type(val)}\n"
        if err != "":
            raise BadFile(err)
        return val

    def set_param(self, **kwd: Any) -> None:
        """Set the object parameters.

        This function must be used, instead of directly setting parameters,
        so that the check/autocast/lock features can be correctly used.

        @raise LockedError: raised if attempted on a locked object
        @raise BadFile: raised if parameters value are of uncorrect type

        """
        if self.locked:
            raise LockedError
        for key, val in kwd.items():
            val = self.checked_items(key, val)
            setattr(self, key, val)
        self.__post_init__()

    def asdict(self) -> Dict[str, Any]:
        """Convert the object data to a dictionary.

        @return:  full set of parameters as a dictionary
        @rtype: Dict[str, Any]

        """
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
        """Write all parameters in a json file.

        @param filename: name of json file
        @type filename: str

        """
        with open(filename, "w") as out:
            dump(self.asdict(), out, indent=4)

    def lock(self) -> None:
        """Lock the object, preventing parameter changes."""
        self._locked = True

    def unlock(self) -> None:
        """Unlock the object, enabling parameter changes."""
        self._locked = False

    @property
    def locked(self) -> bool:
        """Lock state.

        @return: lock state value @rtype: bool

        """
        if not hasattr(self, "_locked"):
            self._locked = False
        return self._locked

    def set_autocast(self) -> None:
        """Set the autocast feature."""
        self._autocast = True

    def unset_autocast(self) -> None:
        """Unset the autocast feature."""
        self._autocast = False

    @property
    def autocast(self) -> bool:
        """Autocast state.

        @return: autocast flag value
        @rtype: bool

        """
        if not hasattr(self, "_autocast"):
            self._autocast = True
        return self._autocast

    def set_checktype(self) -> None:
        """Set the checktype feature."""
        self._checktype = True

    def unset_checktype(self) -> None:
        """Unset the checktype feature."""
        self._checktype = False

    @property
    def checktype(self) -> bool:
        """Checktype state.

        @return: checktype flag value
        @rtype: bool

        """
        if not hasattr(self, "_checktype"):
            self._checktype = True
        return self._checktype


@dataclass
class RuleParam(Readerclass):
    """Parameters for a reaction rule."""

    reactants: List[str] = field(default_factory=list)
    """list of reactants categories"""
    builder_func: str = ""
    """Name of the reaction builder function"""
    builder_const: str = "kinvar"
    """Name of the reaction constant builder function"""
    builder_variant: str = "novariant"
    """Name of the reaction variant builder function"""
    descr: str = ""
    """Description of the reaction"""
    robust: bool = True
    """Reaction robustness flag (products are recomputed
    at each reaction proceesing step when non-robust)"""


@dataclass
class RulesetParam(Readerclass):
    """Parameters for a reaction ruleset."""

    rulemodel: str = "metadynamic.models.polymers"
    """rulemodel module to be used"""
    categories: List[str] = field(default_factory=list)
    """list of categories to be imported from the rulemodel"""
    relations: List[str] = field(default_factory=list)
    """list of relations to be imported from the rulemodel"""
    properties: List[str] = field(default_factory=list)
    """list of properties to be imported from the rulemodel"""
    rules: Dict[str, RuleParam] = field(default_factory=dict)
    """list of rules to be imported from the rulemodel, followinfg L{RuleParam} format"""


@dataclass
class StatParam(Readerclass):
    """Parameters for statistics."""

    prop: str = "count"
    """property used for the statistic calculation"""
    weight: str = "count"
    """property used as weight (e.g. for weighted average)"""
    method: str = "m"
    """method; '+' returns the weighted sum, 'm' the weighted average,
       'max' the weighted max value, 'min' the weighted min value"""
    full: bool = False
    """statistics to be performed on active (False) are the full pool (True)"""
    collection: str = "compounds"
    """Statisitic to be performed on 'compounds' or 'reactions'"""


@dataclass
class MapParam(Readerclass):
    """Parameters for map statistics."""

    prop: str = "count"
    """property used for the statistic calculation"""
    weight: str = "count"
    """property used as weight (e.g. for weighted average)"""
    sort: str = "count"
    """property used for sorting the statistics in categories"""
    method: str = "+"
    """method; '+' returns the weighted sum, 'm' the weighted average"""
    full: bool = False
    """statistics to be performed on active (False) are the full pool (True)"""
    collection: str = "compounds"
    """Statisitic to be performed on 'compounds' or 'reactions'"""


@dataclass
class Param(Readerclass):
    """Run parameters."""

    # Description
    name: str = "run"
    """Run name"""
    comment: str = ""
    """Run comment"""
    # Save
    savedir: str = ""
    """Where the final .hdf5 file will be saved. Defaults to working dir"""
    logdir: str = ""
    """Where the logs will be saved. If empty, log to standard output"""
    loglevel: str = "INFO"
    """logging level (CRITICAL, ERROR, WARNING, INFO, or DEBUG)"""
    # chemical
    conc: float = 0.1
    """Concentration"""
    ptot: int = field(init=False)
    """Total population"""
    vol: float = field(init=False)
    """System volume"""
    init: Dict[str, int] = field(default_factory=dict)
    """initial concentrations"""
    rulemodel: str = "metadynamic.models.polymers"
    """name of the module containing the rule model code"""
    reactions: List[str] = field(default_factory=list)
    """list of reaction types to load from rulemodel (if empty, load all)"""
    parameters: Dict[str, float] = field(default_factory=dict)
    """kinetic constants"""
    # simulation
    tend: float = 1.0
    """final simulation time"""
    tstep: float = 0.01
    """timestep"""
    sstep: float = -1
    """snapshot step (if <0: only at end, if = 0: at each step)"""
    rtlim: float = 900.0
    """Limit runtime"""
    maxsteps: int = 10000
    """maximum time steps"""
    # System
    autoclean: bool = True
    """"Perform a periodic cleaning of probabilities if True"""
    dropmode: str = ""
    """drop mode (can be 'keep', 'drop' or 'soft')"""
    gcperio: bool = True
    """if True, only call garbage collector at each timestep."""
    sleeptime: float = 0.1
    """finished threads will wait in idle loops of corresponding values the end of others"""
    maxmem: int = 0
    """Max memory (in Mb). If 0, set to maxmem_percent/100 * total physical memory"""
    maxmem_percent: int = 95
    """Maximum percentage of memory to be used (used if maxmem set to 0)"""
    # IO
    save: List[str] = field(default_factory=list)
    """list of compounds to be saved at each time step"""
    stat: str = ""
    """json filename describing statistics"""
    statparam: Dict[str, StatParam] = field(init=False)
    """List of statistic parameters, with fields as described in L{StatParam}"""
    maps: str = ""
    """json filename describing stat maps"""
    mapsparam: Dict[str, MapParam] = field(init=False)
    """List of statistic parameters, with fields as described in L{MapParam}"""
    store_snapreac: bool = False
    """Store reaction snapshots? (can take lots of time for large CRNs)"""
    maxstrlen: int = 256
    """max string length to be stored in hdf5"""
    lengrow: int = 20
    """number of length left before requesting a resize"""
    maxlog: int = 100
    """max log lines per process to be saved"""
    timeformat: str = "[%d.%m.%Y-%H:%M:%S]"
    """timeformat used in log files"""

    def __post_init__(self) -> None:
        self.ptot = sum([pop * len(comp) for comp, pop in self.init.items()])
        self.vol = self.ptot / self.conc
        self.statparam = StatParam.readmultiple(self.stat)
        self.mapsparam = MapParam.readmultiple(self.maps)
        if self.maxmem == 0:
            self.maxmem = int(
                self.maxmem_percent * virtual_memory().total / 1024 / 1024 / 100
            )

    def set_param(self, **kwd: Any) -> None:
        if "parameters" in kwd:
            oldparams = self.parameters.copy()
            for key, val in kwd["parameters"].items():
                oldparams[key] = val
            kwd["parameters"] = oldparams
        super().set_param(**kwd)


@dataclass
class DotParam(Readerclass):
    """Parameters for graphviz CRN representation."""

    # type
    binode: bool = False
    """reaction graph with only compounds as single nodes (False)
    or both compounds and reaction as dual nodes (True)"""
    # Graph rendering
    margin: float = 0.0
    """graph margin"""
    concentrate: bool = False
    """if true, try to reduce graph size"""
    maxsize: float = 10.0
    """maximum size"""
    # compounds
    min_fontsize: float = 1.0
    """minimum compound font size"""
    max_fontsize: float = 20.0
    """maximum compound font size"""
    min_c_width: float = 0.05
    """minimal compounds width"""
    max_c_width: float = 0.5
    """maximal compound width"""
    c_penwidth: float = 0.5
    """compound pen width"""
    c_margin: float = 0.01
    """compound margin"""
    c_color: str = "blue"
    """compound color"""
    c_powerscale: float = 0.5
    """compound size scale as pop**c_powerscale"""
    font_powerscale: float = 1.0
    """compound font scale as pop**font_powerscale"""
    # reactions
    min_r_width: float = 0.01
    """minimum reaction width"""
    max_r_width: float = 0.1
    """maximum reaction width"""
    r_color: str = "red"
    """reaction color"""
    r_powerscale: float = 1.0
    """reaction scale as const**r_powerscale"""
    # flows
    min_f_width: float = 0.001
    """minimum flow width"""
    max_f_width: float = 5.0
    """maximum flow width"""
    f_color: str = "black"
    """flow color"""
    cutoff: float = 0.05
    """cut data below cutoff fraction"""
    f_powerscale: float = 1.0
    """scale flow as prob**f+powerscale"""
