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
    Paramrel,
    parmul,
    linproc,
    arrhenius,
    kinvar,
    kalternate,
    novariant,
    rangevariant,
    joiner,
)

# Categorizer #

# A polymer is a chain of characters, e.g. abc
polym: Categorizer = lambda name: name.isalpha()
# a source is a polymer postfixed by a '#', e.g. abc#
source: Categorizer = lambda name: name[-1] == "#" and name[:-1].isalpha()

# Propertizer #

length: Propertizer = lambda name: len(name)

# ProdBuilder #

# e.g. abc + def -> abcdef
merge: ProdBuilder = joiner("")
# e.g. abcdef -[3]-> abc + def
cut: ProdBuilder = lambda names, variant: (names[0][:variant], names[0][variant:])
# e.g. abc# -> abc# + abc
fromsource: ProdBuilder = lambda names, variant: (names[0], names[0][:-1])
# e.g. abc ->
destroy: ProdBuilder = lambda names, variant: ()


# Paramrel builder #

# input flux proportionale to process number
# (for parameter scan)
p_in: Paramrel = linproc("kin_min", "kin_max")

kpol0_2: Paramrel = parmul("kpol0", "alpha")
khyd0_2: Paramrel = parmul("khyd0", "alpha")
kpol_mono: Paramrel = arrhenius("kpol0", "Ea_pol")
khyd_head: Paramrel = arrhenius("khyd0", "Ea_hyd")
kpol_long: Paramrel = arrhenius("kpol0_2", "Ea_pol")
khyd_queue: Paramrel = arrhenius("khyd0_2", "Ea_hyd")


# ConstBuilder #

# Constants following arrhenius law

kin: ConstBuilder = kinvar("p_in")

kpol: ConstBuilder = kalternate(
    condition=lambda names, variant: length(names[0]) == 1,
    name_t="kpol_mono",
    name_f="kpol_long",
)

khyd: ConstBuilder = kalternate(
    condition=lambda names, variant: variant == 1,
    name_t="khyd_head",
    name_f="khyd_queue",
)

# constants proportional to reactant length
kout: ConstBuilder = lambda names, k, variant: (k["kout0"] * length(names[0]))


# VariantBuilder #

# (length-1) possible reactions from a given reaction
# (e.g. abc -[1]-> a+bc and abc -[2]->  ab+c)
intervariant: VariantBuilder = rangevariant(first_offset=1, last_offset=0, reacnum=0)


default_ruleset: Dict[str, Any] = {
    "categories": ["polym", "source"],
    "properties": ["length"],
    "relations": ["p_in", "kpol0_2", "khyd0_2", "kpol_mono", "khyd_head", "kpol_long", "khyd_queue"],
    "rules": {
        "P": {
            "reactants": ["polym", "polym"],
            "builder_func": "merge",
            "builder_const": "kpol",
            "builder_variant": "novariant",
            "descr": "Polymerization",
        },
        "H": {
            "reactants": ["polym"],
            "builder_func": "cut",
            "builder_const": "khyd",
            "builder_variant": "intervariant",
            "descr": "Hydrolysis",
        },
        "in": {
            "reactants": ["source"],
            "builder_func": "fromsource",
            "builder_const": "kin",
            "builder_variant": "novariant",
            "descr": "Constant input from source",
        },
        "out": {
            "reactants": ["polym"],
            "builder_func": "destroy",
            "builder_const": "kout",
            "builder_variant": "novariant",
            "descr": "Compound destruction (or output to sink)",
        },
    },
}
