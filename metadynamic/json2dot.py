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

from numpy import array

from json import load
from metadynamic.inputs import Json2dotParam


class Scaler:
    def __init__(self, data, minimal, maximal, cutoff=0):
        self.minval, self.maxval = self.minmax(data, cutoff)
        self.minimal = minimal
        self.maximal = maximal

    def __call__(self, value):
        return (
            self.maximal * (value - self.minval) + self.minimal * (self.maxval - value)
        ) / (self.maxval - self.minval)

    @staticmethod
    def minmax(data, cutoff=0):
        maxval = max(data)
        data = array(data)
        data = data[data > maxval * cutoff]
        minval = min(data)
        return minval, maxval


class Json2dot:
    def __init__(self, filename, parameterfile=""):
        with open(filename) as infile:
            data = load(infile)
        self.compounds = data["Compounds"]
        self.reactions = data["Reactions"]
        self.param = Json2dotParam()
        if parameterfile:
            self.param.readfile(parameterfile)

    def write(self, filename):
        with open(filename, "w") as out:
            out.write("digraph {\n")
            compounds = set()
            reactions = set()
            out.write("# Flows\n")
            color = self.param.f_color
            scaler = Scaler(
                data=[rate for _, rate in self.reactions.values()],
                minimal=self.param.min_f_width,
                maximal=self.param.max_f_width,
                cutoff=self.param.cutoff,
            )
            for name, (_, rate) in self.reactions.items():
                if rate >= scaler.minval:
                    reactions.add(name)
                    width = scaler(rate)
                    reactants, products = name.split("->")
                    for reac in reactants.split("+"):
                        num, reacname = self.cutdown(reac)
                        compounds.add(reacname)
                        for _ in range(num):
                            out.write(
                                f'"{reacname}" -> "{name}" [penwidth={width}, color={color}];\n'
                            )
                    for prod in products.split("+"):
                        num, prodname = self.cutdown(prod)
                        compounds.add(reacname)
                        for _ in range(num):
                            out.write(
                                f'"{name}" -> "{prodname}" [penwidth={width}, color={color}];\n'
                            )
            out.write("# Compounds\n")
            color = self.param.c_color
            scaler = Scaler(
                data=list(self.compounds.values()),
                minimal=self.param.min_c_width,
                maximal=self.param.max_c_width,
            )
            f_scaler = Scaler(
                data=list(self.compounds.values()),
                minimal=self.param.min_fontsize,
                maximal=self.param.max_fontsize,
            )
            for name in compounds:
                pop = self.compounds[name]
                width = scaler(pop)
                fontsize = f_scaler(pop)
                out.write(
                    f'"{name}" [shape="circle", width={width}, fontsize={fontsize}, color={color}];\n'
                )
            out.write("# Reactions\n")
            color = self.param.r_color
            scaler = Scaler(
                data=[const for const, _ in self.reactions.values()],
                minimal=self.param.min_r_width,
                maximal=self.param.max_r_width,
            )
            for name in reactions:
                const, _ = self.reactions[name]
                width = scaler(const)
                out.write(f'"{name}" [shape="point", width={width}, color={color}];\n')
            out.write("}\n")

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
