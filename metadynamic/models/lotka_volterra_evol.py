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
    kinvar,
    novariant,
    rangevariant,
    joiner,
)

# Categorizer #

grass: Categorizer = lambda name: name[0] == "G"
rabbit: Categorizer = lambda name: name[0] == "R" or name[0] == "r"
fox: Categorizer = lambda name: name[0] == "F" or name[0] == "f"
hungry: Categorizer = lambda name: name[0].islower()
hungry_fox: Categorizer = lambda name: name[0] == "f"
hungry_rabbit: Categorizer = lambda name: name[0] == "r"
satiated_fox: Categorizer = lambda name: name[0] == "F"
satiated_rabbit: Categorizer = lambda name: name[0] == "R"

energy: Propertizer = lambda name: 1 if grass(name) else 2 if rabbit(
    name
) else 4 if fox(name) else 0

# ProdBuilder #


def merge_genes(gene1: str, gene2: str) -> str:
    cutat = randint(1, max(len(gene1), len(gene2)) - 1)
    return gene1[:cutat] + gene2[cutat:]


def rotate_gene(gene: str) -> str:
    startat = randint(0, len(gene) - 1)
    return gene[startat:] + gene[:startat]


def mutate_gene(gene: str) -> str:
    mutateat = randint(0, len(gene) - 2)
    return gene[:mutateat] + str(randint(0, 9)) + gene[mutateat + 1 :]


# 2 satiated animals -> 3 hungry animals
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
    kind = parent1[0].lower()
    parent1 = kind + parent1[1:]
    parent2 = kind + parent2[1:]
    child = kind + "-".join(genes3)
    return parent1, parent2, child


# hungry animal + food -> satiated animal
def eat(names: Compset, variant: int) -> Compset:
    return (names[0][0].upper() + names[0][1:],)


# animal -> grass
def die_gen(size: int) -> ProdBuilder:
    def die(names: Compset, variant: int) -> Compset:
        return ("G",) * size

    return die


r_die: Compset = die_gen(2)
f_die: Compset = die_gen(4)


# Paramrel builder #

k_repro_f: Paramrel = parmul("k_repro_r", "repro_f_factor")
k_repro_r_rot: Paramrel = parmul("k_repro_r", "rot_ratio")
k_repro_f_rot: Paramrel = parmul("k_repro_f", "rot_ratio")
k_repro_r_mute: Paramrel = parmul("k_repro_r", "mute_ratio")
k_repro_f_mute: Paramrel = parmul("k_repro_f", "mute_ratio")

k_death_f0: Paramrel = parmul("k_death_r0", "death_f_factor")
k_death_hf0: Paramrel = parmul("k_death_f0", "hunger_factor")
k_death_hr0: Paramrel = parmul("k_death_r0", "hunger_factor")

k_eat_f0: Paramrel = parmul("k_eat_r0", "eat_f_factor")

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


def k_eat_f(names: Compset, k: Parameters, variant: int) -> float:  # ConstBuilder
    hunt_power = translate_gene(names[0])[1]
    escape_power = translate_gene(names[1])[1]
    return k["k_eat_f0"] * hunt_power * escape_power


k_eat_r: ConstBuilder = kinvar("k_eat_r0")


def k_death_rabbit(
    names: Compset, k: Parameters, variant: int
) -> float:  # ConstBuilder
    powers = translate_gene(names[0])
    const = k["k_death_r0"] if names[0][0].isupper() else k["k_death_hr0"]
    return const * powers[0] * (1 - powers[1])


def k_death_fox(names: Compset, k: Parameters, variant: int) -> float:  # ConstBuilder
    powers = translate_gene(names[0])
    const = k["k_death_f0"] if names[0][0].isupper() else k["k_death_hf0"]
    return const * powers[0] * powers[1]


threevariant: VariantBuilder = lambda reactants: range(3)

# Default Ruleset #

default_ruleset: Dict[str, Any] = {
    "categories": [
        "grass",
        "hungry",
        "rabbit",
        "hungry_rabbit",
        "satiated_rabbit",
        "fox",
        "hungry_fox",
        "satiated_fox",
    ],
    "properties": ["energy"],
    "relations": [
        "k_repro_f",
        "k_repro_r_rot",
        "k_repro_f_rot",
        "k_repro_r_mute",
        "k_repro_f_mute",
        "k_death_f0",
        "k_death_hf0",
        "k_death_hr0",
        "k_eat_f0",
    ],
    "rules": {
        "Rrepro": {
            "reactants": ["satiated_rabbit", "satiated_rabbit"],
            "builder_func": "reproduce",
            "builder_const": "k_repro_rabbit",
            "builder_variant": "threevariant",
        },
        "Frepro": {
            "reactants": ["satiated_fox", "satiated_fox"],
            "builder_func": "reproduce",
            "builder_const": "k_repro_fox",
            "builder_variant": "threevariant",
        },
        "Feat": {
            "reactants": ["hungry_fox", "rabbit"],
            "builder_func": "eat",
            "builder_const": "k_eat_f",
            "builder_variant": "novariant",
        },
        "Reat": {
            "reactants": ["hungry_rabbit", "grass"],
            "builder_func": "eat",
            "builder_const": "k_eat_r",
            "builder_variant": "novariant",
        },
        "Rdeath": {
            "reactants": ["rabbit"],
            "builder_func": "r_die",
            "builder_const": "k_death_rabbit",
            "builder_variant": "novariant",
        },
        "Fdeath": {
            "reactants": ["fox"],
            "builder_func": "f_die",
            "builder_const": "k_death_fox",
            "builder_variant": "novariant",
        },
    },
}
