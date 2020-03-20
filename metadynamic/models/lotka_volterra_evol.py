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

from random import randint
from typing import Dict, Any, List

from metadynamic.ruleset import (
    Categorizer,
    Propertizer,
    ProdBuilder,
    ConstBuilder,
    VariantBuilder,
    Paramrel,
    Parameters,
    Compset,
    parmul,
    kalternate,
    novariant,
    rangevariant,
    joiner,
)

# Categorizer #

rabbit: Categorizer = lambda name: name[0] == "R"
hungry_fox: Categorizer = lambda name: name[0] == "H"
fox: Categorizer = lambda name: name[0] == "F"


# ProdBuilder #


def merge_genes(gene1: str, gene2: str) -> str:
    cutat = randint(1, max(len(gene1), len(gene2)) - 1)
    return gene1[:cutat] + gene2[cutat:]


def rotate_gene(gene: str) -> str:
    startat = randint(0, len(gene) - 1)
    return gene[startat:] + gene[:startat]


def mutate_gene(gene: str) -> str:
    mutateat = randint(0, len(gene) - 2)
    return gene[:mutateat] + str(randint(0, 9)) + gene[mutateat+1:]


# 2 rabbits -> 3 rabbits, or 2 foxes -> 3 hungry foxes
def reproduce(names: Compset, variant: int) -> Compset:
    parent1 = names[0]
    parent2 = names[1]
    genes1 = parent1[1:].split("-")
    genes2 = parent1[1:].split("-")
    genes3 = [merge_genes(g1, g2) for g1, g2 in zip(genes1, genes2)]
    if variant == 1:
        genes3 = [rotate_gene(gene) for gene in genes3]
    if variant == 2:
        genes3 = [mutate_gene(gene) for gene in genes3]
    if parent1[0] == "R":
        child = "R" + "-".join(genes3)
    else:
        parent1 = "H" + parent1[1:]
        parent2 = "H" + parent2[1:]
        child = "H" + "-".join(genes3)
    return parent1, parent2, child


# hungry fox + rabbit -> fox
def eat(names: Compset, variant: int) -> Compset:
    return ("F" + names[0][1:],)


# animal ->
def die(names: Compset, variant: int) -> Compset:
    return ()


# Paramrel builder #

k_repro_r_rot: Paramrel = parmul("k_repro_r", "rot_ratio")
k_repro_f_rot: Paramrel = parmul("k_repro_f", "rot_ratio")
k_repro_r_mute: Paramrel = parmul("k_repro_r", "mute_ratio")
k_repro_f_mute: Paramrel = parmul("k_repro_f", "mute_ratio")
k_death_hf0: Paramrel = parmul("k_death_f0", "hunger_factor")

# ConstBuilder #


def translate_gene(name: str) -> List[float]:
    genes = [float("0." + g) for g in name[1:].split("-")]
    return genes


def k_repro_gen(name: str, name_rot: str, name_mute: str) -> Compset:
    def k_repro(names: Compset, k: Parameters, variant: int) -> float:  # ConstBuilder
        repro_power_a = translate_gene(names[0])[0]
        repro_power_b = translate_gene(names[1])[0]
        if variant == 1:
            goodname = name_rot
        elif variant == 2:
            goodname = name_mute
        else:
            goodname = name
        return k[goodname] * repro_power_a * repro_power_b

    return k_repro


k_repro_rabbit: ConstBuilder = k_repro_gen(
    "k_repro_r", "k_repro_r_rot", "k_repro_r_mute"
)
k_repro_fox: ConstBuilder = k_repro_gen("k_repro_f", "k_repro_f_rot", "k_repro_f_mute")


def k_eat(names: Compset, k: Parameters, variant: int) -> float:  # ConstBuilder
    hunt_power = translate_gene(names[0])[1]
    escape_power = translate_gene(names[1])[1]
    return k["k_eat_0"] * hunt_power * escape_power


def k_death_rabbit(
    names: Compset, k: Parameters, variant: int
) -> float:  # ConstBuilder
    powers = translate_gene(names[0])
    return k["k_death_r0"] * powers[0] * (1 - powers[1])


def k_death_fox_gen(name: str) -> Compset:
    def k_death(names: Compset, k: Parameters, variant: int) -> float:  # ConstBuilder
        powers = translate_gene(names[0])
        return k[name] * powers[0] * powers[1]

    return k_death


k_death_fox: ConstBuilder = k_death_fox_gen("k_death_f0")
k_death_hungry_fox: ConstBuilder = k_death_fox_gen("k_death_hf0")


threevariant: VariantBuilder = lambda reactants: range(3)

# Default Ruleset #

default_ruleset: Dict[str, Any] = {
    "categories": ["rabbit", "fox", "hungry_fox"],
    "relations": [
        "k_repro_r_rot",
        "k_repro_f_rot",
        "k_repro_r_mute",
        "k_repro_f_mute",
        "k_death_hf0",
    ],
    "rules": {
        "Rrepro": {
            "reactants": ["rabbit", "rabbit"],
            "builder_func": "reproduce",
            "builder_const": "k_repro_rabbit",
            "builder_variant": "threevariant",
        },
        "Frepro": {
            "reactants": ["fox", "fox"],
            "builder_func": "reproduce",
            "builder_const": "k_repro_fox",
            "builder_variant": "threevariant",
        },
        "Eat": {
            "reactants": ["hungry_fox", "rabbit"],
            "builder_func": "eat",
            "builder_const": "k_eat",
            "builder_variant": "novariant",
        },
        "Rdeath": {
            "reactants": ["rabbit"],
            "builder_func": "die",
            "builder_const": "k_death_rabbit",
            "builder_variant": "novariant",
        },
        "Fdeath": {
            "reactants": ["fox"],
            "builder_func": "die",
            "builder_const": "k_death_fox",
            "builder_variant": "novariant",
        },
        "Hdeath": {
            "reactants": ["hungry_fox"],
            "builder_func": "die",
            "builder_const": "k_death_hungry_fox",
            "builder_variant": "novariant",
        },
    },
}
