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

from metadynamic.inval import invalidint
from metadynamic.ruleset import (
    Categorizer,
    Propertizer,
    ProdBuilder,
    ConstBuilder,
    VariantBuilder,
    Compset,
    Paramset,
)

# Categorizer

ispolym: Categorizer = lambda name: name.isalpha()
ismono: Categorizer = lambda name: ispolym(name) and len(name) == 1
isact: Categorizer = lambda name: name[-1] == "*" and name[:-1].isalpha()
isactmono: Categorizer = lambda name: isact(name) and len(name) == 2
islongpol: Categorizer = lambda name: ispolym(name) and len(name) > 1

# Propertizer

length: Propertizer = lambda name: (
    1
    if ismono(name)
    else len(name)
    if ispolym(name)
    else len(name) - 1
    if isact(name)
    else 0
)


def asym(name: str) -> int:  # Propertizer
    res = 0
    for char in name:
        if char.isupper():
            res += 1
        elif char.islower():
            res -= 1
    return res


# asym: Propertizer = lambda name: (
#    sum(map(str.isupper, name)) - sum(map(str.islower, name))
# )  ### Slower.

isright: Categorizer = lambda name: asym(name) > 0
isleft: Categorizer = lambda name: asym(name) < 0

# ProdBuilder


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


merge: ProdBuilder = joiner("")
cut: ProdBuilder = lambda names, variant: (names[0][:variant], names[0][variant:])
act_polym: ProdBuilder = lambda names, variant: (names[0][:-1] + names[1],)
activ: ProdBuilder = lambda names, variant: (names[0] + "*",)
deactiv: ProdBuilder = lambda names, variant: (names[0][:-1],)
epimer: ProdBuilder = lambda names, variant: (
    names[0][:variant] + names[0][variant].swapcase() + names[0][variant + 1 :],
)


# ConstBuilder


def samecase(one: str, two: str) -> bool:
    return (one.islower() and two.islower()) or (one.isupper() and two.isupper())


# Invariant constant
kinvar: ConstBuilder = lambda names, k, variant: k[0]

kfastmono: ConstBuilder = lambda names, k, variant: (
    k[0] if length(names[0]) == 1 else k[1]
)
kslowmono: ConstBuilder = lambda names, k, variant: (
    k[0] if length(names[0]) != 1 else k[1]
)
khyd: ConstBuilder = lambda names, k, variant: (
    k[0] if samecase(names[0][variant - 1], names[0][variant]) else k[1]
)
kactselect: ConstBuilder = lambda names, k, variant: (
    k[0] if samecase(names[0][-2], names[1][0]) else k[1]
)


def kmidselect(names: Compset, k: Paramset, variant: int) -> float:  # ConstBuilder
    name = names[0]
    samebefore: bool = variant < (length(name) - 1) and samecase(
        name[variant], name[variant + 1]
    )
    sameafter: bool = (variant > 0) and samecase(name[variant], name[variant - 1])
    if samebefore:
        if sameafter:
            return k[2]
        else:
            return k[1]
    if sameafter:
        return k[1]
    return k[0]


# VariantBuilder

novariant: VariantBuilder = lambda reactants: (invalidint,)
intervariant: VariantBuilder = lambda reactants: range(1, int(length(reactants[0])))
lenvariant: VariantBuilder = lambda reactants: range(int(length(reactants[0])))


def singlevariant(num: int) -> VariantBuilder:
    return lambda reactants: (num,)


firstonly: VariantBuilder = singlevariant(0)
