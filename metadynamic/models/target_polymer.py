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

from fuzzywuzzy.fuzz import ratio

from metadynamic.ruleset import Categorizer, ProdBuilder, ConstBuilder, splitter, kinvar

from metadynamic.models.polymers import *


default_ruleset = default_ruleset.copy()

#  target compounds as '{abcd}' patterns
target: Categorizer = lambda name: name[0] == "{" and name[-1] == "}" and polym(
    name[1:-1]
)


def complex(name: str) -> bool:  # Categorizer
    """ complex compounds are the aggregation of a target and a polymer
    {abcd}:efgh
    """
    if "." in name:
        left, right = name.split(".")
        return polym(right) and target(left)
    return False


complexation: ProdBuilder = joiner(sep=".")
dissociation: ProdBuilder = splitter(sep=".")

k_complex: ConstBuilder = kinvar("k_c")
k_disso: ConstBuilder = lambda names, k, variant: k["k_d"] * k["k_aff"] ** (
    -float(ratio(*names[0].split("."))) / 100
)

default_ruleset["categories"].append("target")
default_ruleset["categories"].append("complex")
default_ruleset["rules"]["ct"] = {
    "reactants": ["target", "polym"],
    "builder_func": "complexation",
    "builder_const": "k_complex",
    "builder_variant": "novariant",
    "descr": "complexation between a target and a polymer",
}
default_ruleset["rules"]["dt"] = {
    "reactants": ["complex"],
    "builder_func": "dissociation",
    "builder_const": "k_disso",
    "builder_variant": "novariant",
    "descr": "dissociation of a complex into  a target and a polymer",
}
