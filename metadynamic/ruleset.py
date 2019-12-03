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
# reactants, variant -> products
ProdBuilder = Callable[[Tuple[str, ...], int], Tuple[str, ...]]
# reactants, parameters, variant -> constant
ConstBuilder = Callable[[Tuple[str, ...], Tuple[float, ...], int], float]
# reactants -> variants
VariantBuilder = Callable[[Tuple[str, ...]], Iterable[int]]
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
    def __init__(self, rule: "Rule", reactantnames: Tuple[str, ...]):
        self.rule: Rule = rule
        self.reactantnames: Tuple[str, ...] = reactantnames

    def builprod(self, variant: int) -> Tuple[str, ...]:
        return self.rule.builder[0](self.reactantnames, variant)

    # Which best place for storing k???
    def buildconst(self, k: Tuple[float, ...], variant: int) -> float:
        return self.rule.builder[1](self.reactantnames, k, variant)

    @property
    def name(self) -> str:
        return f"{self.rule.name}:{'+'.join(self.reactantnames)}"


@dataclass
class Rule:
    name: str
    reactants: Tuple[str, ...]
    builder: Builder
    descr: str
    constants: Tuple[float, ...] = field(default_factory=tuple)
    initialized: bool = False

    def set_constants(self, const_list: Tuple[float, ...]) -> None:
        self.constants = const_list
        self.initialized = True

    def __call__(
        self, reactants: Tuple[str, ...]
    ) -> List[Tuple[Tuple[str, ...], float]]:
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
        self, rulename: str, reactants: Tuple[str, ...], builder: Builder, descr: str,
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

        # get the categories to which belongs comp_name
        comp_categories = self.descriptor.categories(comp_name)
        # Will look fot the list of reactions for each rule
        result: List[ReacDescr] = []
        for rule in self.rules.values():
            res: Set[Tuple[str, ...]] = set()
            for pos, reactant_category in enumerate(rule.reactants):
                if reactant_category in comp_categories:
                    # Then scan all possible combinations, with fixing comp_name in each possible pos
                    for pos2, other_category in enumerate(rule.reactants):
                        for reactants in product(
                            *[
                                [comp_name] if pos2 == pos
                                #  need to adapt comp_collect ...
                                else self.comp_collect.categories[other_category]
                                for pos2, other_category in enumerate(rule.reactants)
                            ]
                        ):
                            res.add(tuple(reactants))
            # Then create all possible reaction descriptions (to be fixed)
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


def isactmono(name: str) -> bool:
    return isact(name) and len(name) == 2


# ProdBuilder


def joiner(sep: str) -> ProdBuilder:
    def joinersep(names: Tuple[str, ...], variant: int = invalidint) -> Tuple[str, ...]:
        return (sep.join(names),)

    return joinersep


def cut(names: Tuple[str, ...], variant: int) -> Tuple[str, ...]:
    tocut = names[0]
    # if other names, they are catalysts
    return tocut[:variant], tocut[variant:]


def act_polym(names: Tuple[str, ...], variant: int = invalidint) -> Tuple[str, ...]:
    return (names[0][:-1] + names[1],)


def activ(names: Tuple[str, ...], variant: int = invalidint) -> Tuple[str, ...]:
    return (names[0] + "*",)


def deactiv(names: Tuple[str, ...], variant: int = invalidint) -> Tuple[str, ...]:
    return (names[0][:-1],)


def epimer(names: Tuple[str, ...], variant: int) -> Tuple[str, ...]:
    return (
        names[0][:variant] + names[0][variant].swapcase() + names[0][variant + 1 :],
    )


# ConstBuilder


def samecase(one: str, two: str) -> bool:
    return (one.islower() and two.islower()) or (one.isupper() and two.isupper())


def kfastmono(
    names: Tuple[str, ...], k: Tuple[float, ...], variant: int = invalidint
) -> float:
    return k[0] * k[1] if length(names[0]) == 1 else k[0]


def khyd(names: Tuple[str, ...], k: Tuple[float, ...], pos: int) -> float:
    tocut: str = names[0]
    return k[0] if samecase(tocut[pos - 1], tocut[pos]) else k[1] * k[0]


def kactselect(
    names: Tuple[str, ...], k: Tuple[float, ...], variant: int = invalidint
) -> float:
    frag0 = names[0][-2]
    frag1 = names[1][0]
    return k[0] if samecase(frag0, frag1) else k[0] * k[1]


def kmidselect(names: Tuple[str, ...], k: Tuple[float, ...], variant: int) -> float:
    name = names[0]
    res = k[0]
    if variant < (length(name) - 1) and samecase(name[variant], name[variant + 1]):
        res *= k[1]
    if (variant > 0) and samecase(name[variant], name[variant - 1]):
        res *= k[1]
    return res


# VariantBuilder


def novariant(reactants: Tuple[str, ...]) -> Iterable[int]:
    return (invalidint,)


def intervariant(reactants: Tuple[str, ...]) -> Iterable[int]:
    return range(length(reactants[0]) - 1)


def lenvariant(reactants: Tuple[str, ...]) -> Iterable[int]:
    return range(length(reactants[0]))


def singlevariant(num: int) -> VariantBuilder:
    def singlevariantn(reactants: Tuple[str, ...]) -> Iterable[int]:
        return (num,)

    return singlevariantn


# Define a specific ruleset

chemdescriptor = Descriptor()
chemdescriptor.add_cat("mono", ismono)
chemdescriptor.add_cat("polym", ispolym)
chemdescriptor.add_cat("actpol", isact)
chemdescriptor.add_cat("actmono", isactmono)
chemdescriptor.add_prop("length", length)

ruleset = Ruleset(chemdescriptor)
ruleset.add_rule(
    rulename="P",
    reactants=("polym", "polym"),
    builder=(joiner(""), kfastmono, novariant),
    descr="Polymerization",
)
ruleset.add_rule(
    rulename="A",
    reactants=("actpol", "polym"),
    builder=(act_polym, kactselect, novariant),
    descr="Activated Polymerization",
)
ruleset.add_rule(
    rulename="M",
    reactants=("actmono", "polym"),
    builder=(act_polym, kactselect, novariant),
    descr="Activated Monomer Polymerization",
)
ruleset.add_rule(
    rulename="a",
    reactants=("polym",),
    builder=(activ, kfastmono, novariant),
    descr="Activation",
)
ruleset.add_rule(
    rulename="d",
    reactants=("actpol",),
    builder=(deactiv, kfastmono, novariant),
    descr="Deactivation",
)
ruleset.add_rule(
    rulename="H",
    reactants=("polym",),
    builder=(cut, khyd, intervariant),
    descr="Hydrolysis",
)
ruleset.add_rule(
    rulename="R",
    reactants=("polym",),
    builder=(epimer, kmidselect, lenvariant),
    descr="Epimerization",
)
ruleset.add_rule(
    rulename="E",
    reactants=("polym",),
    builder=(epimer, kmidselect, singlevariant(0)),
    descr="Epimerization at first end",
)
