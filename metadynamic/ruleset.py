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
metadynamic.ruleset
===================

All elements for designing a set of rule,
i.e. the Core for model design !


Provides:
---------

 - classes for the description of a rule model:
    - L{Parameters}: Maintain a set of parameters {name:values}
      that can be linked with each other.
    - L{Descriptor}: Tool for describing a compound from its name
      (categories and properties)
    - L{Rule}: Describes a given rule for building a reaction
    - L{Model}: Full description of a rule model
 - generators to be used for building models:
    - L{Paramrel} generators:
      - L{parmul}
      - L{arrhenius}
      - L{linproc}
    - L{ConstBuilder} generators:
      - L{kinvar}
      - L{kalternate}
      - L{kdualchoice}
    - L{VariantBuilder} generators:
      - L{novariant_gen}
      - L{singlevariant}
      - L{rangevariant}
    - L{ProdBuilder} generators:
      - L{joiner}
      - L{splitter}

"""

from types import ModuleType
from typing import Callable, Dict, KeysView, Tuple, Set, Iterable, List
from itertools import product
from importlib import import_module
from dataclasses import dataclass

import numpy as np

# from functools import partial, cached_property from 3.8 only
# from functools import reduce


from metadynamic.ends import InitError
from metadynamic.inval import invalidint
from metadynamic.logger import LOGGER
from metadynamic.mpi import MPI_STATUS
from metadynamic.inputs import RulesetParam


# Type alias (~~ data struct)

Compset = Tuple[str, ...]
"""Set of componds"""

Paramdict = Dict[str, float]
"""Parameters, as {name:value}"""

Paramrel = Callable[[Paramdict], float]
"""Function computing a value from a set of parameters
(i.e. calculating a value from a given Paramdict)"""

Stoechio = Iterable[Tuple[str, int]]
"""Stoechimetry, as a list of (compound name,stechiometry number)"""

Categorizer = Callable[[str], bool]
"""Function identifying if a name belongs to a specific category?'
(i.e. returning a bool from a given name string)"""

Propertizer = Callable[[str], float]
"""Function computing a property from a name
(i.e. returning a float from a given name string)"""

ProdBuilder = Callable[[Compset, int], Compset]
"""Function computing products from (reactants, variant)
(i.e. computing a new Compset from a Compset and a reaction variant)"""


class Parameters:
    """
    Maintain a set of parameters {name:values}
    that can be linked with each other.
    """

    def __init__(self, paramdict: Paramdict) -> None:
        """
        Create a Parameters object from a dictionary

        Defaults parameters will be created:
         - T = 300.0 (temperature)
         - pH = 7.0
         - num = MPI rank
         - ntot = MPI size

        @param paramdict: dictionary of parameters as {name:value}
        @type object: Paramdict
        """
        self._paramdict: Paramdict = {}
        """internal parameter values storage"""
        self._relation: Dict[str, Paramrel] = {}
        """relation functions between parameters"""
        self._updating: List[str] = []
        """list of parameters being updated (avoids circular updates)"""
        self.add_set_param("T", 300.0)
        self.add_set_param("pH", 7.0)
        self.add_set_param("num", MPI_STATUS.rank)
        self.add_set_param("ntot", MPI_STATUS.size)
        for key, val in paramdict.items():
            self.add_set_param(key, val)

    def add_set_param(self, key: str, val: float) -> None:
        """
        Create a new parameter and set it.

        @param key: parameter name
        @type key: str
        @param val: parameter value
        @type val: float

        @raise KeyError: if the parameter already exists
        """
        self.add_param(key)
        self.set_param(key, val)

    def add_param(self, key: str) -> None:
        """
        Create a new parameter (set to 0.0)

        @param key: parameter name
        @type key: str

        @raise KeyError: if the parameter already exists
        """
        if key not in self._paramdict:
            self._paramdict[key] = 0.0
        else:
            raise KeyError(f"Key {key} already registered")

    def set_param(self, key: str, val: float, terminate: bool = True) -> None:
        """
        Set a parameter value.

        This will automatically update parameters defined by a relation

        If an attempt is done to fix the value of a parameter defined by a
        relation, the request will be ignored *except* if the 'terminate'
        flag is set to False (do not use this flag manually, this is intended to be used
        by _param_init; direct use may result in an inconsistent set of parameters)

        @param key: parameter name
        @type key: str
        @param val: parameter value
        @type val: float
        @param terminate: do not use manually (see above)
        @type terminate: bool

        @raise KeyError: if the parameter does not exists
        """
        if terminate and key in self._relation:
            LOGGER.warning(
                f"'{key}' was set by a relation, direct setting request to '{val}' ignored"
            )
        elif key in self._paramdict:
            if key not in self._updating:
                self._paramdict[key] = val
                self._updating.append(key)
                self._param_init()
        else:
            raise KeyError(f"Key {key} not registered")
        if terminate:
            self._updating.clear()

    def __getitem__(self, key: str) -> float:
        """Return the value of key. If the key is non existent, it is created and set to 0"""
        try:
            return self._paramdict[key]
        except KeyError:  # better dealing of missing keys ???
            self.add_param(key)
            return 0.0

    def add_relation(self, key: str, relation: Paramrel) -> None:
        """
        Add a new parameter and set it by a relation.

        If the key already exists, it will be overriden.

        @param key: parameter name
        @type key: str
        @param relation: relation function
        @type relation: Paramrel
        """
        try:
            self.add_param(key)
        except KeyError:
            LOGGER.warning(
                f"Setting '{key}' relation will override already set value ({self[key]})"
            )
        self._relation[key] = relation
        self._param_init()

    def _param_init(self) -> None:
        """Initialize all relation-defined parameters"""
        for param, func in self._relation.items():
            self.set_param(param, func(self._paramdict), terminate=False)

    def __repr__(self) -> str:
        return dict.__repr__(self._paramdict)

    def __str__(self) -> str:
        return str(self._paramdict)


ConstBuilder = Callable[[Compset, Parameters, int], float]
"""Function computing a kinetic constant from (reactants, parameters, variant)
(i.e. computing a float from a Compset, a Parameters and a reaction variant)"""

VariantBuilder = Callable[[Compset], Iterable[int]]
"""Function computing a list of possible reaction variants from reactants
(i.e. computing an iterable of int from a Compset)"""

Builder = Tuple[ProdBuilder, ConstBuilder, VariantBuilder]
"""Full reaction builder, with all elements for describing a reaction
simply a tuple (prodbuilder,constbuilder,variantbuilder)"""

ReacDescr = Tuple[str, Compset, int]
"""reaction description as (rule name, reactants, variant)"""

ReacProp = Tuple[Stoechio, Stoechio, float, bool]
"""reaction property, as
(reactants stoechiometry, products stoechiometry, constant, robustness)"""


class Descriptor:
    """tool for describing a compound from its name (categories and properties)"""

    def __init__(self) -> None:
        self.cat_dict: Dict[str, Categorizer] = {}
        """Collection of Categorizers"""
        self.prop_dict: Dict[str, Propertizer] = {}
        """Collection of Propertizers"""

    @property
    def catlist(self) -> KeysView[str]:
        """
        List of possible category names

        @rtype: KeysView[str]
        """
        return self.cat_dict.keys()

    def prop(self, propname: str, name: str) -> float:
        """
        Return the property from a name

        @param propname: name of the property to be computed
        @type propname: str
        @param name: name of the compound to be evaluated
        @type name: str
        @return: computed property value:
        @rtype: float
        """
        try:
            return self.prop_dict[propname](name)
        except KeyError:
            return float(self.cat_dict[propname](name))

    def categories(self, name: str) -> Set[str]:
        """
        Return the set of categories a compound belongs to

        @param name: name of the compound to be evaluated
        @type name: str
        @return: set of categories
        @rtype: Set[str]
        """
        return {catname for catname, rule in self.cat_dict.items() if rule(name)}

    def __repr__(self) -> str:
        return f"Descriptor: {self.cat_dict.keys()}"

    def add_cat(self, catname: str, rule: Categorizer) -> None:
        """
        Add a new Categorizer to the Descriptor

        @param catname: name of the categorizer
        @type catname: str
        @param rule: Categorizer function
        @type rule: Categorizer
        """
        if catname in self.cat_dict:
            raise KeyError(f"Category {catname} already defined")
        self.cat_dict[catname] = rule

    def add_prop(self, propname: str, func: Propertizer) -> None:
        """
        Add a new Properizer to the Descriptor

        @param propname: name of the propertizer
        @type propname: str
        @param func: Propertizer function
        @type func: Propertizer
        """
        if propname in self.prop_dict:
            raise KeyError(f"Propert {propname} already defined")
        self.prop_dict[propname] = func


@dataclass
class Rule:
    """Describes a given rule for building a reaction"""

    name: str
    """rule name"""
    reactants: Compset
    """set of reactant category names"""
    builder: Builder
    """Builder function"""
    descr: str
    """reaction rule description"""
    parameters: Parameters
    """parameter set"""
    robust: bool
    """robustness"""

    def _build_products(self, reactants: Compset, variant: int) -> Compset:
        """
        Build the set of products from a set of reactants

        @param reactants: set of reactants
        @type reactants: Compset
        @param variant: reaction variant
        @type variant: int
        @return: set of products
        @rtype: Compset
        """
        products: Compset = self.builder[0](reactants, variant)
        if "" in products:
            raise InitError(
                f"Reaction from {reactants} lead to null compound: {products}"
            )
        return products

    def _build_constant(self, reactants: Compset, variant: int) -> float:
        """
        Build the chemical constant from a set of reactants

        @param reactants: set of reactants
        @type reactants: Compset
        @param variant: reaction variant
        @type variant: int
        @return: constant value
        @rtype: float
        """
        return self.builder[1](reactants, self.parameters, variant)

    def build(self, description: ReacDescr) -> ReacProp:
        """
        Build a reaction from its description

        @param description: reaction description
        @type description: ReacDescr
        @return: reaction properties
        @rtype: ReacProp
        """
        _, reactants, variant = description
        products: Compset = self._build_products(reactants, variant)
        constant: float = self._build_constant(reactants, variant)
        return (
            self.getstoechio(reactants),
            self.getstoechio(products),
            constant,
            self.robust,
        )

    def rebuild_prod(self, description: ReacDescr) -> Stoechio:
        """
        (re)build the list of products generated by the reaction

        @param description: reaction description
        @type description: ReacDescr
        @return: reactants stoechiometry
        @rtype: Stoechio
        """
        _, reactants, variant = description
        return self.getstoechio(self._build_products(reactants, variant))

    @staticmethod
    def getstoechio(compounds: Compset) -> Stoechio:
        """
        Return the stoechiometry of a compound set

        @param compounds: set of compounds
        @type compounds: Compset
        @return: compounds stoechiometry
        @rtype: Stoechio
        """
        #  Extend common cases for faster computations
        length = len(compounds)
        if length == 1:
            return ((compounds[0], 1),)
        if length == 2:
            c_0 = compounds[0]
            c_1 = compounds[1]
            if c_0 == c_1:
                return ((c_0, 2),)
            return ((c_0, 1), (c_1, 1))
        # General computation from order 3
        return ((reac, compounds.count(reac)) for reac in set(compounds))

    def __str__(self) -> str:
        return self.descr


class Model:
    """full description of a rule model"""

    def __init__(
        self, modelparam: str, reactions: List[str], paramdict: Paramdict
    ) -> None:
        """
        Create a model

        @param modelparam: name of the module containing the rule model code
        @type modelparam: str
        @param reactions: list of reaction types to load from rulemodel (if empty, load all)
        @type reactions: List[str]
        @param paramdict:parameters values
        @type paramdict: Paramdict
        """
        self.modelparam: str = modelparam
        """name of the module containing the rule model code"""
        self.rulepath: ModuleType
        """Direct access to model code"""
        self.param: RulesetParam
        """ruleset parameters"""
        self._load_param()
        self.reactions: Iterable[
            str
        ] = reactions if reactions else self.param.rules.keys()
        """list of reactions types"""
        self.parameters = Parameters(paramdict)
        """set of model parameters (constants, temperature, etc.)"""
        self.descriptor: Descriptor = Descriptor()
        """set of categorizers/propertizers"""
        self._create_descriptors()
        # create rules
        self.categories: Dict[str, Set[str]] = {}
        """Dictionary of categories"""
        for catname in self.descriptor.catlist:
            self.categories[catname] = set()
        self.rules: Dict[str, Rule] = {}
        """Dictionary of rules"""
        self._create_rules()

    def _load_param(self) -> None:
        """load parameters of the rule module"""
        try:
            # model given as a module, e.g. metadynamic.models.polymers
            self.rulepath = import_module(self.modelparam)
            try:
                ruleset = getattr(self.rulepath, "default_ruleset")
            except AttributeError:
                InitError(
                    f"'{self.modelparam}' module didn't define a 'default_ruleset'"
                )
            self.param = RulesetParam.readdict(ruleset)
        except AttributeError:
            raise InitError(
                f"Model {self.modelparam} does not provides a 'default_ruleset' dict"
            )
        except ModuleNotFoundError:
            # model given as a json file
            self.param = RulesetParam.readfile(self.modelparam)
            self.rulepath = import_module(self.param.rulemodel)

    def _create_descriptors(self) -> None:
        """create the descriptors"""
        for rel in self.param.relations:
            try:
                self.parameters.add_relation(rel, getattr(self.rulepath, rel))
            except AttributeError:
                LOGGER.error(f"relation '{rel}' not found in {self.modelparam}")
        for catname in self.param.categories:
            try:
                self.descriptor.add_cat(catname, getattr(self.rulepath, catname))
            except AttributeError:
                LOGGER.error(f"category '{catname}' not found in {self.modelparam}")
        for propname in self.param.properties:
            try:
                self.descriptor.add_prop(propname, getattr(self.rulepath, propname))
            except AttributeError:
                LOGGER.error(f"'[property '{propname}' not found in {self.modelparam}")

    def _create_rules(self) -> None:
        """Read and create rules"""
        for rulename in self.reactions:
            # Get rule from parameter file
            if rulename in self.param.rules:
                ruleparam = self.param.rules[rulename]
                try:
                    # then create it from rule module
                    rule = Rule(
                        name=rulename,
                        reactants=tuple(ruleparam.reactants),
                        builder=(
                            getattr(self.rulepath, ruleparam.builder_func),
                            getattr(self.rulepath, ruleparam.builder_const),
                            getattr(self.rulepath, ruleparam.builder_variant),
                        ),
                        descr=ruleparam.descr,
                        parameters=self.parameters,
                        robust=ruleparam.robust,
                    )
                except AttributeError:
                    # raise an error if the rule from file is not in the module
                    raise InitError(
                        f"The rule '{rulename}' from '{self.modelparam}'"
                        f"is not defined in '{self.rulepath}'"
                    )
                # Register the created rule
                self._add_rule(rulename, rule)
            else:
                LOGGER.warning(
                    f"Reaction '{rulename}' not found in {self.rulepath}, ignored."
                )

    def _add_rule(self, rulename: str, rule: Rule) -> None:
        """
        Add a new reaction rule

        @param rulename: name of the rule
        @type param: str
        @param rule: rule object
        @type rule: Rule
        """
        if rulename in self.rules:
            raise KeyError(f"Rule {rulename} already defined")
        self.rules[rulename] = rule
        for reac in rule.reactants:
            try:
                self.categories[reac].add(rulename)
            except KeyError:
                raise ValueError(f"Unrecognize category {reac}")

    def get_related(
        self, comp_name: str, coll_cat: Dict[str, Set[str]]
    ) -> Set[ReacDescr]:
        """
        Build the description of all reactions that can be performed
        from a given compound

        @param comp_name: name of the compound to react
        @type comp_name: str
        @param coll_cat: dictionary of compounds to react with, sorted in categories
        @type coll_cat: Dict[str, Set[str]]
        @return: a set of reaction description
        @rtype: Set[ReacDescr]
        """
        # get the categories to which belongs comp_name
        comp_categories = self.descriptor.categories(comp_name)
        res: Set[ReacDescr] = set()
        # Scan all registered rules
        # NOTE: dense and critical code. to be checked again (performance are important!)
        for rulename, rule in self.rules.items():
            # Check if the compound is concerned by the rule
            # necessary? avoid to scan further, but => overhead...
            if comp_categories & set(rule.reactants):
                # scan all the rules reactant
                for pos, reactant_category in enumerate(rule.reactants):
                    # Check if the compound can be in this pos
                    if reactant_category in comp_categories:
                        # OK => scan all possible combinations, with fixing comp_name in this pos
                        combinations = product(
                            *[
                                # fixed position => original compound
                                [comp_name] if pos2 == pos
                                # other positions => scam all other concerned compounds
                                else coll_cat[other_category]
                                for pos2, other_category in enumerate(rule.reactants)
                            ]
                        )
                        # Then for each combination...
                        for reactants in combinations:
                            # ...Compute all variants...
                            for variant in rule.builder[2](reactants):
                                # ...and add them to the result...
                                res.add((rulename, reactants, variant))
        return res

    def buildreac(self, reacdescr: ReacDescr) -> ReacProp:
        """
        Build a reaction from its description

        @param reacdescr: reaction description
        @type reacdescr: ReacDescr
        @return: reaction properties
        @rtype: ReacProp
        """
        return self.rules[reacdescr[0]].build(reacdescr)

    def rebuild_prod(self, reacdescr: ReacDescr) -> Stoechio:
        """
        (re)build the list of products generated by the reaction

        @param reacdescr: reaction description
        @type reacdescr: ReacDescr
        @return: reactants stoechiometry
        @rtype: Stoechio
        """
        return self.rules[reacdescr[0]].rebuild_prod(reacdescr)


# - Generic elements - #

# Paramrel Generators #


def parmul(name: str, factor: str) -> Paramrel:
    """
    Generate a parameter relation →'name'×'factor'

    @param name: name of the base parameter.
    @type name: str
    @param factor: name of the factor parameter.
    @type factor: str
    @return: Paramrel function
    @rtype: Paramrel
    """
    return lambda k: k[name] * k[factor]


def arrhenius(k_0: str, eact: str) -> Paramrel:
    """
    Generate a parameter relation →'k_0'·e^(-'eact'/(R·'T'))
    With R=8.314

    @param k_0: name of the Arrhenius pre-exponential factor
    @type k0: str
    @param eact: name of the activation energy
    @type eact: str
    @return: Paramrel function
    @rtype: Paramrel
    """
    return lambda k: k[k_0] * float(np.exp(-k[eact] / 8.314 / k["T"]))


def linproc(start: str, end: str) -> Paramrel:
    """
    Generate a parameter relation returning a value proportional to
    the thread number, ranging from 'start' to 'end'

    @param start: name of the parameter given to the first thread
    @type start: str
    @param end: name of the parameter given to the last thread
    @type end: str
    @return: Paramrel function
    @rtype: Paramrel
    """
    return lambda k: float(np.linspace(k[start], k[end], k["ntot"])[k["num"]])


# ConstBuilder Generators #


def kinvar(name: str) -> ConstBuilder:
    """
    Generates a Constbuilder giving the same constant 'name' in all situations

    @param name: name of the constant to return
    @type name: str
    @return: ConstBuilder function
    @rtype: ConstBuilder
    """
    return lambda names, k, variant: k[name]


def kalternate(
    condition: Callable[[Compset, int], bool], name_t: str, name_f: str
) -> ConstBuilder:
    """
    Generates a Constbuilder giving constant 'name_t' when 'condition' is True,
    and 'name_f' when 'condition' is False.

    @param condition: function that takes as arguments a set of compound names
        (a tuple of strings) and a variant (integer)
    @type condition: Callable[[Compset, int], bool]
    @param name_t: name of the constant to return when condition is True
    @type name_t: str
    @param name_f: name of the constant to return when condition is False
    @type name_f: str
    @return: ConstBuilder function
    @rtype: ConstBuilder
    """
    return (
        lambda names, k, variant: k[name_t] if condition(names, variant) else k[name_f]
    )


def kdualchoice(
    condition_1: Callable[[Compset, int], bool],
    condition_2: Callable[[Compset, int], bool],
    name_tt: str,
    name_ff: str,
    name_tf: str,
    name_ft: str = "",
) -> ConstBuilder:
    """
    Generates a Constbuilder giving constant 'name_xy' where x and y are either
    t or f (True or False), depending on the respective booleam results
    of 'condition_1' and 'condition_2'.

    @param condition_1: first condition, as a function that takes as arguments
        a set of compound names (a tuple of strings) and a variant (integer)
    @type condition_1: Callable[[Compset, int], bool]
    @param condition_2: second condition
    @type condition_2: Callable[[Compset, int], bool]
    @param name_tt: name of the constant to return when conditions are (True, True)
    @type name_tt: str
    @param name_ff: name of the constant to return when conditions are (False, False)
    @type name_ff: str
    @param name_tf: name of the constant to return when conditions are (True, False)
    @type name_tf: str
    @param name_ft: name of the constant to return when conditions are (False, True)
        if not set (default), it will be set equal to 'name_tf'.
    @type name_ft: str
    @return: ConstBuilder function
    @rtype: ConstBuilder
    """
    if not name_ft:
        name_ft = name_tf

    def kdual(names: Compset, k: Parameters, variant: int) -> float:  # ConstBuilder
        c_1 = condition_1(names, variant)
        c_2 = condition_2(names, variant)
        if c_1:
            if c_2:
                return k[name_tt]
            return k[name_tf]
        if c_2:
            return k[name_ft]
        return k[name_ff]

    return kdual


# VariantBuilder generatore #


def novariant_gen() -> VariantBuilder:
    """
    Generate a VariantBuilder for reactions with no variant
    (i.e. only one possible outcome).

    The returned VariantBuilder will return a single invalid variant value,
    thus is intended to be used with ProdBuilder and ConstBuilder
    that do not use variant for computing the reaction property.

    @return: VariantBuilder function
    @rtype: VariantBuilder
    """
    return lambda reactants: (invalidint,)


def singlevariant(num: int) -> VariantBuilder:
    """
    Generate a VariantBuilder for reactions with a single variant
    (i.e. only one possible outcome).

    The returned VariantBuilder will return a single valid variant value 'num',
    thus is intended to be used with ProdBuilder and ConstBuilder
    that do use variant for computing the reaction property, but in cases where
    only one of all possible variants is expected to be used.

    @param num: variant number to return
    @type num: int
    @return: VariantBuilder function
    @rtype: VariantBuilder
    """
    return lambda reactants: (num,)


def rangevariant(
    reacnum: int, first_offset: int = 0, last_offset: int = 0
) -> VariantBuilder:
    """
    Generate a VariantBuilder for reactions with as many variants
    as the length of the reactant number 'reacnum',
    ranging from 0 to len(reactants[reacnum]).

    A first_offset can be added to change the first variant.
    A last_offset can be added to change the last variant.

    @param reacnum: number of the reactant to be used for computing the variant
    @type reacnum: int
    @param first_offset: offset to be added to the first variant
    @type first_offset: int
    @param last_offset: offset to be added to the last variant
    @type last_offset: int
    @return: VariantBuilder function
    @rtype: VariantBuilder
    """
    return lambda reactants: range(first_offset, len(reactants[reacnum]) + last_offset)


# ProdBuilder generatore #


def joiner(sep: str) -> ProdBuilder:
    """
    Generate a ProdBuilder that will fuse all reactants into a single
    compounds joining all the names in one using a 'sep' as a separator string.

        >>> chainer=joiner("-")
        >>> chainer(["A","B","C"], invalidint)
        "A-B-C"

    @param sep: string separator
    @type str: str
    @return: ProdBuilder function
    @rtype: ProdBuilder
    """
    return lambda names, variant: (sep.join(names),)


def splitter(sep: str) -> ProdBuilder:
    """
    Generate a ProdBuilder that will cut the first reactant into individual compounds
    by splitting at each 'sep' occurence:

        >>> cutter=splitter("-")
        >>> cutter(["A-B-C"], invalidint)
        ["A", "B", "C"]

    @param sep: string separator
    @type str: str
    @return: ProdBuilder function
    @rtype: ProdBuilder
    """
    return lambda names, variant: tuple(names[0].split(sep))
