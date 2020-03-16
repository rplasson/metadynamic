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

from typing import Callable, Dict, KeysView, Tuple, Set, Iterable, List
from itertools import product
from importlib import import_module
from dataclasses import dataclass, field

# from functools import partial, cached_property from 3.8 only
# from functools import reduce


from metadynamic.ends import InitError
from metadynamic.inval import invalidint

# from metadynamic.logger import LOGGER
from metadynamic.inputs import RulesetParam


# Type alias (~~ data struct)
Compset = Tuple[str, ...]
Paramset = List[float]
Paramdict = Dict[str, float]
Paramrel = Callable[[Paramdict], float]
Stoechio = Iterable[Tuple[str, int]]
Categorizer = Callable[[str], bool]
Propertizer = Callable[[str], float]
# reactants, variant -> products
ProdBuilder = Callable[[Compset, int], Compset]
# reactants, parameters, variant -> constant
ConstBuilder = Callable[[Compset, Paramset, int], float]
# reactants -> variants
VariantBuilder = Callable[[Compset], Iterable[int]]
Builder = Tuple[ProdBuilder, ConstBuilder, VariantBuilder]
# rule, reactants, variant
ReacDescr = Tuple[str, Compset, int]
# products, constant, stoechiometry
ReacProp = Tuple[Stoechio, Stoechio, float]
ChemDescr = str


class Parameters:
    def __init__(self):
        self._paramdict: Paramdict = {}
        self._relation: Dict[str, Paramrel] = {}
        self._updating: List[str] = []

    def add_param(self, key: str) -> None:
        if key not in self._paramdict:
            self._paramdict[key] = 0.0
        else:
            raise KeyError(f"Key {key} already registered")

    def set_param(self, key: str, val: float, terminate=True) -> None:
        if key in self._paramdict:
            if key not in self._updating:
                self._paramdict[key] = val
                self._updating.append(key)
                self.param_init()
        else:
            raise KeyError(f"Key {key} not registered")
        if terminate:
            self._updating.clear()

    def __getitem__(self, key: str) -> float:
        return self._paramdict[key]

    def add_relation(self, key: str, relation: Paramrel) -> None:
        self._relation[key] = relation
        self.param_init()

    def param_init(self):
        for param, func in self._relation.items():
            self.set_param(param, func(self._paramdict), terminate=False)

    def __repr__(self) -> str:
        return dict.__repr__(self._paramdict)

    def __str__(self) -> str:
        return str(self._paramdict)


class Descriptor:
    def __init__(self) -> None:
        self.cat_dict: Dict[str, Categorizer] = {}
        self.prop_dict: Dict[str, Propertizer] = {}

    @property
    def catlist(self) -> KeysView[str]:
        return self.cat_dict.keys()

    def prop(self, propname: str, name: str) -> float:
        try:
            return self.prop_dict[propname](name)
        except KeyError:
            return float(self.cat_dict[propname](name))

    def categories(self, name: str) -> Set[str]:
        return {catname for catname, rule in self.cat_dict.items() if rule(name)}

    def __repr__(self) -> str:
        return f"Descriptor: {self.cat_dict.keys()}"

    def add_cat(self, catname: str, rule: Categorizer) -> None:
        if catname in self.cat_dict:
            raise KeyError(f"Category {catname} already defined")
        self.cat_dict[catname] = rule

    def add_prop(self, propname: str, func: Propertizer) -> None:
        if propname in self.prop_dict:
            raise KeyError(f"Propert {propname} already defined")
        self.prop_dict[propname] = func


@dataclass
class Rule:
    name: str
    reactants: Compset
    builder: Builder
    descr: str
    constants: Paramset = field(default_factory=list)
    initialized: bool = False

    def set_constants(self, const_list: Paramset) -> None:
        self.constants = const_list
        self.initialized = True

    def build(self, description: ReacDescr) -> ReacProp:
        _, reactants, variant = description
        if not self.initialized:
            raise InitError("Rule {self} used before constant initialization")
        products: Compset = self.builder[0](reactants, variant)
        if "" in products:
            raise InitError(f"Reaction {description} lead to null compund: {products}")
        constant: float = self.builder[1](reactants, self.constants, variant)
        return self.getstoechio(reactants), self.getstoechio(products), constant

    @staticmethod
    def getstoechio(compounds: Compset) -> Stoechio:
        #  Extend common cases for faster computations
        length = len(compounds)
        if length == 1:
            return ((compounds[0], 1),)
        if length == 2:
            c0 = compounds[0]
            c1 = compounds[1]
            if c0 == c1:
                return ((c0, 2),)
            return ((c0, 1), (c1, 1))
        # General computation from order 3
        return ((reac, compounds.count(reac)) for reac in set(compounds))

    def __str__(self) -> str:
        return self.descr


class Ruleset:
    def __init__(self, descriptor: Descriptor):
        self.descriptor: Descriptor = descriptor
        self.categories: Dict[str, Set[str]] = {}
        for catname in descriptor.catlist:
            self.categories[catname] = set()
        self.rules: Dict[str, Rule] = {}

    def add_rule(self, rulename: str, rule: Rule) -> None:
        if rulename in self.rules:
            raise KeyError(f"Rule {rulename} already defined")
        self.rules[rulename] = rule
        for reac in rule.reactants:
            try:
                self.categories[reac].add(rulename)
            except KeyError:
                raise ValueError(f"Unrecognize category {reac}")

    def initialize(self, paramdict: Dict[str, Paramset]) -> None:
        for rulename, parameters in paramdict.items():
            self.rules[rulename].set_constants(parameters)

    def get_related(
        self, comp_name: str, coll_cat: Dict[str, Set[str]]
    ) -> Set[ReacDescr]:
        # get the categories to which belongs comp_name
        comp_categories = self.descriptor.categories(comp_name)
        res: Set[ReacDescr] = set()
        # Scan all registered rules
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
        return self.rules[reacdescr[0]].build(reacdescr)


class Model:
    def __init__(self, modelparam: str, paramdict: Dict[str, Paramset]) -> None:
        # load parameters an rule module
        self.param = RulesetParam.readfile(modelparam)
        self.rulepath = import_module(self.param.rulemodel)
        # create descriptors
        self.descriptor: Descriptor = Descriptor()
        for catname, catparam in self.param.categories.items():
            self.descriptor.add_cat(catname, getattr(self.rulepath, catparam.func))
        for propname, propparam in self.param.properties.items():
            self.descriptor.add_prop(propname, getattr(self.rulepath, propparam.func))
        # create rules
        self.ruleset: Ruleset = Ruleset(self.descriptor)
        for rulename in paramdict:
            # Gate rule from parameter file
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
                    )
                except KeyError:
                    # raise an error if the rule from paramdict is not in parameter file
                    raise InitError(
                        f"The rule '{rulename}' from '{paramdict}' is not defined in '{modelparam}'"
                    )
                except AttributeError:
                    # raise an error if the rule from file is not in the module
                    raise InitError(
                        f"The rule '{rulename}' from '{modelparam}' is not defined in '{self.rulepath}'"
                    )
                # Register the created rule
                self.ruleset.add_rule(rulename, rule)
        # intialize the rules parameters
        self.ruleset.initialize(paramdict)


# Generic elements

# Invariant constant
kinvar: ConstBuilder = lambda names, k, variant: k[0]


# Reaction with no variants
novariant: VariantBuilder = lambda reactants: (invalidint,)


# Reactins with a single variant
def singlevariant(num: int) -> VariantBuilder:
    return lambda reactants: (num,)
