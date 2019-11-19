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

from typing import Callable, List, Any
from itertools import product
from dataclasses import dataclass

# from functools import partial, cached_property from 3.8 only
# from functools import reduce


from metadynamic.collector import Collected


class Descriptor:
    def __init__(self):
        self.cat_dict = {}
        self.prop_dict = {}

    @property
    def catlist(self):
        return self.cat_dict.keys()

    @property
    def prop(self, propname):
        return self.prop_dict[propname]()

    def categories(self, name) -> List[str]:
        return [catname for catname, rule in self.cat_dict.items() if rule(name)]

    def __repr__(self) -> str:
        return f"Descriptor: {self.cat_dict.keys()}"

    def add_cat(self, catname: str, rule: Callable[[str], bool]) -> None:
        self.cat_dict[catname] = rule

    def add_prop(self, propname: str, func: Callable[[None], Any]) -> None:
        self.prop_dict[propname] = func


class ReacDescr:
    def __init__(self, rule, reactantnames):
        self.rule = rule
        self.reactantnames = reactantnames

    def builprod(self):
        return self.rule.prodbuilder(self.reactantnames)

    # Which best place for storing k???
    def buildconst(self, k):
        return self.rule.constbuilder(k)

    @property
    def name(self):
        return f"{self.rule.name}:{'+'.join(self.reactantnames)}"


@dataclass
class Rule:
    name: str
    reactants: List[str]
    prodbuilder: Callable[[List[str]], str]
    constbuilder: Callable[List[str], List[float]]
    descr: str


class Ruleset(Collected):
    def __init__(self, descriptor: Descriptor):
        self.descriptor = descriptor
        self.categories = {}
        for catname in descriptor.catlist:
            self.categories[catname] = set()
        self.rules = {}

    def add_rule(
        self,
        rulename: str,
        reactants: List[str],
        prodbuilder: Callable[[List[str]], str],
        constbuilder: Callable[List[str], List[float]],
        descr: str = "",
    ):
        self.rules[rulename] = Rule(
            name=rulename,
            reactants=reactants,
            prodbuilder=prodbuilder,
            constbuilder=constbuilder,
            descr=descr,
        )
        for reac in reactants:
            try:
                self.categories[reac].add(rulename)
            except KeyError:
                raise ValueError(f"Unrecognize category {reac}")

    def get_related(self, name: str):
        # Maybe memoize the list of rule for a given list of categories...
        # rule_related = reduce(
        #    lambda x, y: x | y,
        #    [self.categories[catname] for catname in self.descriptor.categories(name)],
        # )
        # or simply use self.rules??? Far less overhead, maybe not much perf loss.
        categories = self.descriptor.categories(name)
        # Will look fot the list of reactions for each rule
        result = []
        for rule in self.rules:
            res = set()
            # get the list of reactant type for the rule
            for reaclist in rule.reactants:
                # Then scan all possible combinations, with fixing name in each possible pos
                for pos, reacname in enumerate(reaclist):
                    if reacname in categories:
                        res |= set(
                            product(
                                *[
                                    [name] if pos2 == pos
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


# Functions for rulesets


def joiner(names, sep):
    return sep.join(names)


# joiner_direct = partial(joiner, sep="") from python 3.8 only
def joiner_direct(names):
    return "".join(names)


def samecase(one: str, two: str) -> bool:
    return (one.islower() and two.islower()) or (one.isupper() and two.isupper())


def kpol(names, k) -> float:
    return k[0] if samecase(names[0][-1], names[1][0]) else k[1]


def cut(names, pos):
    tocut = names[0]
    # if other names, they are catalysts
    return tocut[:pos], tocut[pos:]


def khyd(names, k, pos) -> float:
    tocut = names[0]
    return k[0] if samecase(tocut[pos - 1], tocut[pos]) else k[1]


# Define a specific ruleset

chemdescriptor = Descriptor()
chemdescriptor.add_cat("mono", ismono)
chemdescriptor.add_cat("polym", ispolym)
chemdescriptor.add_cat("actpol", isact)
chemdescriptor.add_prop("length", length)

ruleset = Ruleset(chemdescriptor)
ruleset.add_rule("P", ["polym", "polym"], joiner_direct, kpol, "Polymerization")
ruleset.add_rule("H", ["polym"], cut, khyd, "Hydrolysis")  # /!\ need parameter (pos)
