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

from metadynamic.ruleset import Categorizer, ProdBuilder, ConstBuilder

from metadynamic.models.polymers import *

default_ruleset = default_ruleset.copy()


# Necessary/useful for avinding clash between polymers/catpolymers model?

# Example of catalized polymer by a dimer
# iff the dimer ends corresponds to polymers ends
# e.g.|   aaa + bbb + ab --> aaabbb + ab

dimer: Categorizer = lambda name: polym(name) and len(name) == 2

cat_polym: ProdBuilder = lambda names, variant: (names[0] + names[1], names[2])

k_cat_dimer_pol: ConstBuilder = lambda names, k, variant: (
    k["k_cat"] if names[0][-1] == names[2][0] and names[1][0] == names[2][-1] else 0.0
)

default_ruleset["categories"].append("dimer")

default_ruleset["rules"]["dP"] = {
    "reactants": ["polym", "polym", "dimer"],
    "builder_func": "cat_polym",
    "builder_const": "k_cat_dimer_pol",
    "builder_variant": "novariant",
    "descr": "Catalized polym by dimer if ends fits",
}
