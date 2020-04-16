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

from typing import Dict, Any

from metadynamic.ruleset import (
    Categorizer,
    Propertizer,
    ProdBuilder,
    ConstBuilder,
    VariantBuilder,
    Compset,
    kalternate,
    kdualchoice,
    novariant_gen,
    singlevariant,
    rangevariant,
    joiner,
)

# Categorizer

polym: Categorizer = lambda name: name.isalpha()
mono: Categorizer = lambda name: polym(name) and len(name) == 1
actpol: Categorizer = lambda name: name[-1] == "*" and name[:-1].isalpha()
actmono: Categorizer = lambda name: actpol(name) and len(name) == 2
longpol: Categorizer = lambda name: polym(name) and len(name) > 1

# Propertizer

length: Propertizer = lambda name: (
    1 if mono(name) else len(name) if polym(name) else len(name) - 1 if actpol(name) else 0
)


def asym(name: str) -> int:  # Propertizer
    res = 0
    for char in name:
        if char.isupper():
            res += 1
        elif char.islower():
            res -= 1
    return res


right: Categorizer = lambda name: asym(name) > 0
left: Categorizer = lambda name: asym(name) < 0

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

novariant: VariantBuilder = novariant_gen()
intervariant: VariantBuilder = rangevariant(reacnum=0, first_offset=1)
lenvariant: VariantBuilder = rangevariant(reacnum=0)
firstonly: VariantBuilder = singlevariant(num=0)
lastonly: VariantBuilder = singlevariant(num=-1)


# Default Ruleset

default_ruleset: Dict[str, Any] = {
    "categories": ["mono", "polym", "longpol", "actpol", "actmono", "left", "right"],
    "properties": ["length", "asym"],
    "rules": {
        "P": {
            "reactants": ["polym", "polym"],
            "builder_func": "merge",
            "builder_const": "kpol",
            "builder_variant": "novariant",
            "descr": "Polymerization",
        },
        "A": {
            "reactants": ["actpol", "polym"],
            "builder_func": "act_polym",
            "builder_const": "kpola",
            "builder_variant": "novariant",
            "descr": "Activated Polymerization",
        },
        "M": {
            "reactants": ["actmono", "polym"],
            "builder_func": "act_polym",
            "builder_const": "kpola_mono",
            "builder_variant": "novariant",
            "descr": "Activated Monomer Polymerization",
        },
        "a": {
            "reactants": ["polym"],
            "builder_func": "activ",
            "builder_const": "kact",
            "builder_variant": "novariant",
            "descr": "Activation",
        },
        "d": {
            "reactants": ["actpol"],
            "builder_func": "deactiv",
            "builder_const": "kdeact",
            "builder_variant": "novariant",
            "descr": "Deactivation",
        },
        "H": {
            "reactants": ["polym"],
            "builder_func": "cut",
            "builder_const": "khyd",
            "builder_variant": "intervariant",
            "descr": "Hydrolysis",
        },
        "R": {
            "reactants": ["longpol"],
            "builder_func": "epimer",
            "builder_const": "krac",
            "builder_variant": "lenvariant",
            "descr": "Epimerization",
        },
        "E": {
            "reactants": ["longpol"],
            "builder_func": "epimer",
            "builder_const": "kepi",
            "builder_variant": "firstonly",
            "descr": "Epimerization at first end",
        },
    },
}
