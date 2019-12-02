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

from typing import Callable, List, Any, Dict, KeysView, Tuple, Set, Iterable
from itertools import product
from dataclasses import dataclass, field


# from functools import partial, cached_property from 3.8 only
# from functools import reduce


from metadynamic.chemical import Collected
from metadynamic.ends import InitError
from metadynamic.inval import InvalidInt

invalidint = InvalidInt()


# Type alias
Categorizer = Callable[[str], bool]
Propertizer = Callable[[str], Any]
ProdBuilder = Callable[[List[str], int], List[str]]  # reactants, variant -> products
ConstBuilder = Callable[
    [List[str], List[float], int], float
]  # reactants, parameters, variant -> constant
VariantBuilder = Callable[[List[str]], Iterable[int]]  # reactants -> variants
Builder = Tuple[ProdBuilder, ConstBuilder, VariantBuilder]


class Descriptor:
    def __init__(self) -> None:
        self.cat_dict: Dict[str, Categorizer] = {}
        self.prop_dict: Dict[str, Propertizer] = {}

    @property
    def catlist(self) -> KeysView[str]:
        return self.cat_dict.keys()

    def prop(self, propname: str, name: str) -> Any:
        return self.prop_dict[propname](name)

    def categories(self, name: str) -> List[str]:
        return [catname for catname, rule in self.cat_dict.items() if rule(name)]

    def __repr__(self) -> str:
        return f"Descriptor: {self.cat_dict.keys()}"

    def add_cat(self, catname: str, rule: Categorizer) -> None:
        self.cat_dict[catname] = rule

    def add_prop(self, propname: str, func: Propertizer) -> None:
        self.prop_dict[propname] = func


class ReacDescr:
    def __init__(self, rule: "Rule", reactantnames: List[str]):
        self.rule: Rule = rule
        self.reactantnames: List[str] = reactantnames

    def builprod(self, variant: int) -> List[str]:
        return self.rule.builder[0](self.reactantnames, variant)

    # Which best place for storing k???
    def buildconst(self, k: List[float], variant: int) -> float:
        return self.rule.builder[1](self.reactantnames, k, variant)

    @property
    def name(self) -> str:
        return f"{self.rule.name}:{'+'.join(self.reactantnames)}"


@dataclass
class Rule:
    name: str
    reactants: List[str]
    builder: Builder
    descr: str
    constants: List[float] = field(default_factory=list)
    initialized: bool = False

    def set_constants(self, const_list: List[float]) -> None:
        self.constants = const_list
        self.initialized = True

    def __call__(self, reactants: List[str]) -> List[Tuple[List[str], float]]:
        if not self.initialized:
            raise InitError("Rule {self} used before constant initialization")
        res = []
        for variant in self.builder[2](reactants):
            res.append(
                (
                    self.builder[0](reactants, variant),
                    self.builder[1](reactants, self.constants, variant),
                )
            )
        return res

    def __str__(self) -> str:
        return self.descr


class Ruleset(Collected):
    def __init__(self, descriptor: Descriptor):
        self.descriptor: Descriptor = descriptor
        self.categories: Dict[str, Set[str]] = {}
        for catname in descriptor.catlist:
            self.categories[catname] = set()
        self.rules: Dict[str, Rule] = {}

    def add_rule(
        self, rulename: str, reactants: List[str], builder: Builder, descr: str,
    ) -> None:
        self.rules[rulename] = Rule(
            name=rulename, reactants=reactants, builder=builder, descr=descr,
        )
        for reac in reactants:
            try:
                self.categories[reac].add(rulename)
            except KeyError:
                raise ValueError(f"Unrecognize category {reac}")

    def get_related(self, comp_name: str) -> List[ReacDescr]:
        # Maybe memoize the list of rule for a given list of categories...
        # rule_related = reduce(
        #    lambda x, y: x | y,
        #    [self.categories[catname] for catname in self.descriptor.categories(comp_name)],
        # )
        # or simply use self.rules??? Far less overhead, maybe not much perf loss.
        categories = self.descriptor.categories(comp_name)
        # Will look fot the list of reactions for each rule
        result: List[ReacDescr] = []
        for rule in self.rules.values():
            res: Set[List[str]] = set()
            # get the list of reactant type for the rule
            for reaclist in rule.reactants:
                # Then scan all possible combinations, with fixing comp_name in each possible pos
                for pos, reacname in enumerate(reaclist):
                    if reacname in categories:
                        res |= set(
                            product(
                                [  # Check /!\ ... bad ...
                                    [comp_name] if pos2 == pos
                                    # expect comp_collect.categories to return str
                                    # collector will have to be adapted
                                    else self.comp_collect.categories[catname]
                                    for pos2, catname in enumerate(reaclist)
                                ]
                            )
                        )
            # Then create all possible reaction descriptions
            result += [ReacDescr(rule, reactantnames) for reactantnames in res]
        # Finally return all build reac decription for each rule
        return result


# Functions for descriptor


def ismono(name: str) -> bool:
    return len(name) == 1


def ispolym(name: str) -> bool:
    return name.isalpha()


def isact(name: str) -> bool:
    return name[-1] == "*" and name[:-1].isalpha()


def length(name: str) -> int:
    if ismono(name):
        return 1
    if ispolym(name):
        return len(name)
    if isact(name):
        return len(name) - 1
    return 0


# ProdBuilder


def joiner(names: List[str], sep: str) -> List[str]:
    return [sep.join(names)]


# joiner_direct = partial(joiner, sep="") from python 3.8 only
def joiner_direct(names: List[str], variant: int = -1) -> List[str]:
    return ["".join(names)]


def cut(names: List[str], pos: int) -> List[str]:
    tocut = names[0]
    # if other names, they are catalysts
    return [tocut[:pos], tocut[pos:]]


# ConstBuilder


def samecase(one: str, two: str) -> bool:
    return (one.islower() and two.islower()) or (one.isupper() and two.isupper())


def kpol(names: List[str], k: List[float], variant: int = invalidint) -> float:
    return k[0] if samecase(names[0][-1], names[1][0]) else k[1]


def khyd(names: List[str], k: List[float], pos: int) -> float:
    tocut: str = names[0]
    return k[0] if samecase(tocut[pos - 1], tocut[pos]) else k[1]


# VariantBuilder


def novariant(reactants: List[str]) -> Iterable[int]:
    return [invalidint]


def intervariant(reactants: List[str]) -> Iterable[int]:
    return range(len(reactants[0]) - 1)


def lenvariant(reactants: List[str]) -> Iterable[int]:
    return range(len(reactants[0]))


# Define a specific ruleset

chemdescriptor = Descriptor()
chemdescriptor.add_cat("mono", ismono)
chemdescriptor.add_cat("polym", ispolym)
chemdescriptor.add_cat("actpol", isact)
chemdescriptor.add_prop("length", length)

ruleset = Ruleset(chemdescriptor)
ruleset.add_rule(
    "P", ["polym", "polym"], (joiner_direct, kpol, novariant), "Polymerization"
)
ruleset.add_rule("H", ["polym"], (cut, khyd, intervariant), "Hydrolysis")
