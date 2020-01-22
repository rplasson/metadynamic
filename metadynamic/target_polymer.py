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

from copy import deepcopy

from fuzzywuzzy.fuzz import partial_ratio

# import * should be sufficient, but explicit everything used for better readability
# from metadynamic.polymers import *
from metadynamic.polymers import (
    model,
    Categorizer,
    ProdBuilder,
    ConstBuilder,
    novariant,
    ispolym,
    joiner,
    splitter,
    kinvar,
)

# Necessary/useful for avinding clash between polymers/catpolymers model?
# model = deepcopy(model)

#  target compounds as '{abcd}' patterns
istarget: Categorizer = lambda name: name[0] == "{" and name[-1] == "}" and ispolym(
    name[1:-1]
)


def iscomplex(name: str) -> bool:  # Categorizer
    """ complex compounds are the aggregation of a target and a polymer
    {abcd}:efgh
    """
    if ":" in name:
        left, right = name.split(":")
        return ispolym(right) and istarget(left)
    return False


complexation: ProdBuilder = joiner(sep=":")
dissociation: ProdBuilder = splitter(sep=":")

k_complex: ConstBuilder = kinvar
k_disso: ConstBuilder = lambda names, k, variant: k[0] * partial_ratio(
    *names[0].split(":")
)

model.add_cat("target", istarget)
model.add_cat("complex", iscomplex)

model.add_rule(
    rulename="ct",
    reactants=("target", "polym"),
    builder=(complexation, k_complex, novariant),
    descr="complexation between a target and a polymer",
)

model.add_rule(
    rulename="dt",
    reactants=("complex",),
    builder=(dissociation, k_disso, novariant),
    descr="dissociation of a complex into  a target and a polymer",
)
