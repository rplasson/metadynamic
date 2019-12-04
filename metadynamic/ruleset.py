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

from typing import Callable, Any, Dict, KeysView, Tuple, Set, Iterable, List
from itertools import product
from importlib import import_module
from dataclasses import dataclass, field


# from functools import partial, cached_property from 3.8 only
# from functools import reduce


from metadynamic.ends import InitError


# Type alias (~~ data struct)
Compset = Tuple[str, ...]
Paramset = List[float]
Categorizer = Callable[[str], bool]
Propertizer = Callable[[str], Any]
# reactants, variant -> products
ProdBuilder = Callable[[Compset, int], Compset]
# reactants, parameters, variant -> constant
ConstBuilder = Callable[[Compset, Paramset, int], float]
# reactants -> variants
VariantBuilder = Callable[[Compset], Iterable[int]]
Builder = Tuple[ProdBuilder, ConstBuilder, VariantBuilder]
# rule, reactants, variant
ReacDescr = Tuple[str, Compset, int]
# products, constant
ReacProp = Tuple[Compset, float]
ChemDescr = str


class Descriptor:
    def __init__(self) -> None:
        self.cat_dict: Dict[str, Categorizer] = {}
        self.prop_dict: Dict[str, Propertizer] = {}

    @property
    def catlist(self) -> KeysView[str]:
        return self.cat_dict.keys()

    def prop(self, propname: str, name: str) -> Any:
        return self.prop_dict[propname](name)

    def categories(self, name: str) -> Set[str]:
        return {catname for catname, rule in self.cat_dict.items() if rule(name)}

    def __repr__(self) -> str:
        return f"Descriptor: {self.cat_dict.keys()}"

    def add_cat(self, catname: str, rule: Categorizer) -> None:
        self.cat_dict[catname] = rule

    def add_prop(self, propname: str, func: Propertizer) -> None:
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

    def build(self, reactants: Compset, variant: int) -> ReacProp:
        if not self.initialized:
            raise InitError("Rule {self} used before constant initialization")
        products: Compset = self.builder[0](reactants, variant)
        constant: float = self.builder[1](reactants, self.constants, variant)
        return products, constant

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

    def get_reaction(self, reacdescr: ReacDescr) -> ReacProp:
        rulename: str
        reactantnames: Compset
        variant: float
        rulename, reactantnames, variant = reacdescr
        return self.rules[rulename].build(reactantnames, variant)


class Model:
    def __init__(self):
        self.descriptor = Descriptor()
        self.rules: Dict[str, Rule]

    def add_cat(self, catname: str, rule: Categorizer) -> None:
        self.descriptor.add_cat(catname, rule)

    def add_prop(self, propname: str, func: Propertizer) -> None:
        self.descriptor.add_prop(propname, func)

    def add_rule(
        self, rulename: str, reactants: Compset, builder: Builder, descr: str,
    ) -> None:
        self.rules[rulename] = Rule(
            name=rulename, reactants=reactants, builder=builder, descr=descr,
        )


class Ruled:
    ruleset: Ruleset
    descriptor: Descriptor
    model: Model

    @classmethod
    def setrules(cls, modelpath: str, paramdict: Dict[str, Paramset]) -> None:
        """The set of rules must be set from a python file that can be reached
        from 'modelpath'.
        This file must contain a Model object named 'model'."""
        try:
            cls.model = import_module(modelpath).model
        except AttributeError:
            raise InitError(
                f"The given modelpath '{modelpath}' must contain a Model object named 'model'"
            )
        try:
            cls.descriptor = cls.model.descriptor
        except AttributeError:
            raise InitError(
                f"the object 'model' from {modelpath} must be of type Model"
            )
        cls.ruleset = Ruleset(cls.descriptor)
        for rulename in paramdict:
            try:
                cls.ruleset.add_rule(rulename, cls.model.rules[rulename])
            except KeyError:
                raise InitError(
                    f"The rule named '{rulename}' is not defined in '{modelpath}'"
                )
        cls.ruleset.initialize(paramdict)
