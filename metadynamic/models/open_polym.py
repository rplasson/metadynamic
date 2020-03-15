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
    kinvar,
    novariant,
)

# Categorizer #

# A polymer is a chain of characters, e.g. abc
ispolym: Categorizer = lambda name: name.isalpha()
# a source is a polymer postfixed by a '#', e.g. abc#
issource: Categorizer = lambda name: name[-1] == "#" and name[:-1].isalpha()

# Propertizer #

length: Propertizer = lambda name: len(name)

# ProdBuilder #

# e.g. abc + def -> abcdef
merge: ProdBuilder = lambda names, variant: (names[0]+names[1],)
# e.g. abcdef -[3]-> abc + def
cut: ProdBuilder = lambda names, variant: (names[0][:variant], names[0][variant:])
# e.g. abc# -> abc# + abc
fromsource: ProdBuilder = lambda names, variant: (names[0], names[0][:-1])
# e.g. abc ->
destroy: ProdBuilder = lambda names, variant: ()


# ConstBuilder #

# different constant if the first compound is a monomer
kpol: ConstBuilder = lambda names, k, variant: (
    k[0] if length(names[0]) == 1 else k[1]
)

# different constant for hydrolysis at first position
khyd: ConstBuilder = lambda names, k, variant: (
    k[0] if variant == 1 else k[1]
)

# constants proportional to reactant length
kout: ConstBuilder = lambda names, k, variant: (
    k[0]*length(names[0])
)


# VariantBuilder #

# (length-1) possible reactions from a given reaction
# (e.g. abc -[1]-> a+bc and abc -[2]->  ab+c)
intervariant: VariantBuilder = lambda reactants: range(1, int(length(reactants[0])))
