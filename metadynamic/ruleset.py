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
from dataclasses import dataclass

# from functools import partial, cached_property from 3.8 only
# from functools import reduce


from metadynamic.ends import InitError
from metadynamic.inval import invalidint
from metadynamic.logger import LOGGER
from metadynamic.inputs import RulesetParam


# Type alias (~~ data struct)
Compset = Tuple[str, ...]
Paramdict = Dict[str, float]
Paramrel = Callable[[Paramdict], float]
Stoechio = Iterable[Tuple[str, int]]
Categorizer = Callable[[str], bool]
Propertizer = Callable[[str], float]
# reactants, variant -> products
ProdBuilder = Callable[[Compset, int], Compset]
# reactants, parameters, variant -> constant


class Parameters:
    def __init__(self, paramdict: Paramdict) -> None:
        self._paramdict: Paramdict = {}
        self._relation: Dict[str, Paramrel] = {}
        self._updating: List[str] = []
        for key, val in paramdict.items():
            self.add_param(key)
            self.set_param(key, val)

    def add_param(self, key: str) -> None:
        if key not in self._paramdict:
            self._paramdict[key] = 0.0
        else:
            raise KeyError(f"Key {key} already registered")

    def set_param(self, key: str, val: float, terminate: bool = True) -> None:
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
        try:
            return self._paramdict[key]
        except KeyError:  # better dealing of missing keys ???
            self.add_param(key)
            return 0.0

    def add_relation(self, key: str, relation: Paramrel) -> None:
        self._relation[key] = relation
        self.param_init()

    def param_init(self) -> None:
        for param, func in self._relation.items():
            self.set_param(param, func(self._paramdict), terminate=False)

    def __repr__(self) -> str:
        return dict.__repr__(self._paramdict)

    def __str__(self) -> str:
        return str(self._paramdict)


ConstBuilder = Callable[[Compset, Parameters, int], float]
# reactants -> variants
VariantBuilder = Callable[[Compset], Iterable[int]]
Builder = Tuple[ProdBuilder, ConstBuilder, VariantBuilder]
# rule, reactants, variant
ReacDescr = Tuple[str, Compset, int]
# products, constant, stoechiometry
ReacProp = Tuple[Stoechio, Stoechio, float]
ChemDescr = str


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
    parameters: Parameters
    #    initialized: bool = False

    #    def set_param(self, paramdict: Paramdict) -> None:
    #        self.paramdict = paramdict
    #        self.initialized = True      ### Still useful???

    def build(self, description: ReacDescr) -> ReacProp:
        _, reactants, variant = description
        #        if not self.initialized:
        #            raise InitError("Rule {self} used before constant initialization")
        products: Compset = self.builder[0](reactants, variant)
        if "" in products:
            raise InitError(f"Reaction {description} lead to null compund: {products}")
        constant: float = self.builder[1](reactants, self.parameters, variant)
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

    #    def initialize(self, paramdict: Paramdict) -> None:
    #        for rules in self.rules.values():
    #            rules.set_param(paramdict)

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
    def __init__(
        self, modelparam: str, reactions: List[str], paramdict: Paramdict
    ) -> None:
        # load parameters an rule module
        self.param = RulesetParam.readfile(modelparam)
        self.rulepath = import_module(self.param.rulemodel)
        # create descriptors
        self.descriptor: Descriptor = Descriptor()
        self.parameters = Parameters(paramdict)
        for catname, catparam in self.param.categories.items():
            self.descriptor.add_cat(catname, getattr(self.rulepath, catparam.func))
        for propname, propparam in self.param.properties.items():
            self.descriptor.add_prop(propname, getattr(self.rulepath, propparam.func))
        # create rules
        self.ruleset: Ruleset = Ruleset(self.descriptor)
        # read all rules from reactions if not empty, else read them from [aram.rules
        for rulename in reactions if reactions else self.param.rules:
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
                    )
                except AttributeError:
                    # raise an error if the rule from file is not in the module
                    raise InitError(
                        f"The rule '{rulename}' from '{modelparam}' is not defined in '{self.rulepath}'"
                    )
                # Register the created rule
                self.ruleset.add_rule(rulename, rule)
            else:
                LOGGER.warning(
                    f"Reaction '{rulename}' not found in {self.rulepath}, ignored."
                )


# Generic elements

# Invariant constant
def kinvar(name: str) -> ConstBuilder:
    """Build an invariable kinetic constant named 'name'"""
    return lambda names, k, variant: k[name]


def kalternate(
    condition: Callable[[Compset, int], bool], name_t: str, name_f: str
) -> ConstBuilder:
    """Build an kinetic constant 'name_t' when 'condition' is True,
       and 'name_f' when 'condition' is False.
       'condition' is a function that takes as arguments a set of compound names
       (a tuple of strings) and a variant (integer)"""
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
    """Build an kinetic constant 'name_xy' where x and y are either
       t or f (True or False), depending on the respective booleam results
       of 'condition_1' and 'condition_2'.

       if 'name_ft' is not set, it is equal to 'name_tf'

       'condition_n' are functions that takes as arguments a set of compound names
       (a tuple of strings) and a variant (integer)"""
    if not name_ft:
        name_ft = name_tf

    def kdual(names: Compset, k: Parameters, variant: int) -> float:  # ConstBuilder
        c1 = condition_1(names, variant)
        c2 = condition_2(names, variant)
        if c1:
            if c2:
                return k[name_tt]
            else:
                return k[name_tf]
        if c2:
            return k[name_ft]
        return k[name_ff]

    return kdual


# Reaction with no variants
novariant: VariantBuilder = lambda reactants: (invalidint,)


# Reactins with a single variant
def singlevariant(num: int) -> VariantBuilder:
    return lambda reactants: (num,)


def rangevariant(first_offset: int, last_offset: int, reacnum: int) -> VariantBuilder:
    return lambda reactants: range(first_offset, len(reactants[reacnum]) + last_offset)


def joiner(sep: str) -> ProdBuilder:
    """Generate a "joiner" ProdBuilder using a "sep" as a separator string.
       e.g. chainer=joiner("-") will give a ProdBuilder named chainer
       that will provide "A-B-C" from chainer(["A","B","C"])"""
    return lambda names, variant: (sep.join(names),)


def splitter(sep: str) -> ProdBuilder:
    """Generate a "splitter" ProdBuilder using a "sep" as a separator string.
       e.g. cutter=splitter("-") will give a ProdBuilder named cutter
       that will provide ["A","B","C"] from cutter("A-B-C")"""
    return lambda names, variant: tuple(names[0].split(sep))
