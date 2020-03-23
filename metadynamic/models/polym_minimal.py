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
    kalternate,
    novariant_gen
)

from metadynamic.models import polymers


# Categorizer

# A polymer is a chain of characters, e.g. abc
polym: Categorizer = polymers.polym

# Propertizer

length: Propertizer = polymers.length

# ProdBuilder

# e.g. abc + def -> abcdef
merge: ProdBuilder = polymers.merge
# e.g. abcdef -[3]-> abc + def
cut: ProdBuilder = polymers.cut

# ConstBuilder

kpol: ConstBuilder = polymers.kpol

khyd: ConstBuilder = kalternate(
    condition=lambda names, variant: variant == 1,
    name_t="khyd_head",
    name_f="khyd_queue",
)

# VariantBuilder

novariant: VariantBuilder = novariant_gen()

# (length-1) possible reactions from a given reaction
# (e.g. abc -[1]-> a+bc and abc -[2]->  ab+c)
intervariant: VariantBuilder = polymers.intervariant

# Default Ruleset

default_ruleset: Dict[str, Any] = {
    "categories": ["polym"],
    "properties": ["length"],
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
    },
}
