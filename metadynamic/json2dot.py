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

from itertools import product
from numpy import array
from json import load
from typing import TextIO

from metadynamic.inputs import Json2dotParam


class Scaler:
    def __init__(self, data, minimal, maximal, cutoff=0, powerscale=1):
        self.minval, self.maxval = self.minmax(data, cutoff)
        self.minimal = minimal
        self.maximal = maximal
        self.powerscale = powerscale

    def __call__(self, value):
        n = self.powerscale
        return (
            self.maximal * (value ** n - self.minval ** n)
            + self.minimal * (self.maxval ** n - value ** n)
        ) / (self.maxval ** n - self.minval ** n)

    @staticmethod
    def minmax(data, cutoff=0):
        maxval = max(data)
        data = array(data)
        data = data[data > maxval * cutoff]
        minval = min(data)
        return minval, maxval


class Dotwriter:
    def __init__(self, out: TextIO) -> None:
        self.out = out

    def head(self) -> None:
        self.out.write("digraph {\n")

    def foot(self) -> str:
        self.out.write("}\n")

    def comment(self, text: str) -> None:
        self.out.write(f"    # {text}\n")

    def compound(self, name: str, width: float, fontsize: float, color: str) -> None:
        self.out.write(
            f'"{name}" [shape="circle", width={width}, fontsize={fontsize}, color={color}];\n'
        )

    def reaction(self, name: str, width: float, color: str) -> None:
        self.out.write(f'"{name}" [shape="point", width={width}, color={color}];\n')

    def edge(self, start: str, end: str, width: float, color: str) -> None:
        self.out.write(f'"{start}" -> "{end}" [penwidth={width}, color={color}];\n')


class Json2dot:
    def __init__(self, filename, parameterfile=""):
        with open(filename) as infile:
            data = load(infile)
        self.compounds = data["Compounds"]
        self.reactions = data["Reactions"]
        self.param = self.Json2dotParam.readfile(parameterfile) if parameterfile else Json2dotParam()

    def write(self, filename):
        with open(filename, "w") as out:
            io = Dotwriter(out)
            io.head()
            binode = self.param.binode
            compounds = set()
            reactions = set()
            io.comment("Flows")
            color = self.param.f_color
            scaler = Scaler(
                data=[rate for _, rate in self.reactions.values()],
                minimal=self.param.min_f_width,
                maximal=self.param.max_f_width,
                cutoff=self.param.cutoff,
                powerscale=self.param.f_powerscale,
            )
            for name, (_, rate) in self.reactions.items():
                if rate >= scaler.minval:
                    reactions.add(name)
                    width = scaler(rate)
                    reactants, products = name.split("->")
                    if not binode:
                        reaclist = []
                        prodlist = []
                    for reac in reactants.split("+"):
                        num, reacname = self.cutdown(reac)
                        compounds.add(reacname)
                        if binode:
                            for _ in range(num):
                                io.edge(start=reacname, end=name, width=width, color=color)
                        else:
                            reaclist.append(reacname)
                    for prod in products.split("+"):
                        num, prodname = self.cutdown(prod)
                        compounds.add(prodname)
                        if binode:
                            for _ in range(num):
                                io.edge(start=name, end=prodname, width=width, color=color)
                        else:
                            prodlist.append(prodname)
                    if not binode:
                        for reacname, prodname in product(reaclist, prodlist):
                            io.edge(start=reacname, end=prodname, width=width, color=color)
            io.comment("Compounds")
            color = self.param.c_color
            scaler = Scaler(
                data=list(self.compounds.values()),
                minimal=self.param.min_c_width,
                maximal=self.param.max_c_width,
                powerscale=self.param.c_powerscale,
            )
            f_scaler = Scaler(
                data=list(self.compounds.values()),
                minimal=self.param.min_fontsize,
                maximal=self.param.max_fontsize,
                powerscale=self.param.font_powerscale,
            )
            for name in compounds:
                try:
                    pop = self.compounds[name]
                    width = scaler(pop)
                    fontsize = f_scaler(pop)
                except KeyError:
                    width = 0
                    fontsize = 0
                io.compound(name=name, width=width, fontsize=fontsize, color=color)
            if binode:
                io.comment("Reactions")
                color = self.param.r_color
                scaler = Scaler(
                    data=[const for const, _ in self.reactions.values()],
                    minimal=self.param.min_r_width,
                    maximal=self.param.max_r_width,
                    powerscale=self.param.r_powerscale,
                )
                for name in reactions:
                    const, _ = self.reactions[name]
                    width = scaler(const)
                    io.reaction(name=name, width=width, color=color)
            io.foot()

    @staticmethod
    def minmax(data, cutoff=0):
        maxval = max(data)
        data = array(data)
        data = data[data > maxval * cutoff]
        minval = min(data)
        return minval, maxval

    @staticmethod
    def cutdown(name):
        num = ""
        comp = ""
        start = True
        for char in name:
            if start and char.isdigit():
                num += char
            else:
                comp += char
                start = False
        num = 1 if num == "" else int(num)
        return num, comp
