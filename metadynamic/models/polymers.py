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

from metadynamic.ruleset import (
    Categorizer,
    Propertizer,
    ProdBuilder,
    ConstBuilder,
    VariantBuilder,
    Compset,
    Parameters,
    kalternate,
    kdualchoice,
    novariant,
    singlevariant,
    rangevariant,
    joiner,
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


isright: Categorizer = lambda name: asym(name) > 0
isleft: Categorizer = lambda name: asym(name) < 0

# ProdBuilder

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


def samebefore(names: Compset, variant: int) -> bool:
    name = names[0]
    return variant < (length(name) - 1) and samecase(name[variant], name[variant + 1])


def sameafter(names: Compset, variant: int) -> bool:
    name = names[0]
    return (variant > 0) and samecase(name[variant], name[variant - 1])


kpol: ConstBuilder = kalternate(
    condition=lambda names, variant: length(names[0]) == 1,
    name_t="kpol_mono",
    name_f="kpol_long",
)

kpola: ConstBuilder = kalternate(
    condition=lambda names, variant: samecase(names[0][-2], names[1][0]),
    name_t="kpola_same",
    name_f="kpola_diff",
)

kpola_mono: ConstBuilder = kalternate(
    condition=lambda names, variant: samecase(names[0][-2], names[1][0]),
    name_t="kpola_mono_same",
    name_f="kpola_mono_diff",
)

kact: ConstBuilder = kalternate(
    condition=lambda names, variant: length(names[0]) == 1,
    name_t="kact_mono",
    name_f="kact_pol",
)

kdeact: ConstBuilder = kalternate(
    condition=lambda names, variant: length(names[0]) == 1,
    name_t="kdeact_mono",
    name_f="kdeact_pol",
)

khyd: ConstBuilder = kalternate(
    condition=lambda names, variant: samecase(names[0][variant - 1], names[0][variant]),
    name_t="khyd_same",
    name_f="khyd_diff",
)

kepi: ConstBuilder = kdualchoice(
    condition_1=samebefore,
    condition_2=sameafter,
    name_tt="kepi_same",
    name_ff="kepi_diff",
    name_tf="kepi_mixed",
)

krac: ConstBuilder = kdualchoice(
    condition_1=samebefore,
    condition_2=sameafter,
    name_tt="krac_same",
    name_ff="krac_diff",
    name_tf="krac_mixed",
)

# VariantBuilder

intervariant: VariantBuilder = rangevariant(first_offset=1, last_offset=0, reacnum=0)
lenvariant: VariantBuilder = rangevariant(first_offset=0, last_offset=0, reacnum=0)
firstonly: VariantBuilder = singlevariant(num=0)
lastonly: VariantBuilder = singlevariant(num=-1)
