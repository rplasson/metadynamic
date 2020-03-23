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

import numpy as np

from os import path
from itertools import product
from json import load
from graphviz import Digraph
from subprocess import CalledProcessError
from typing import Any, Tuple, Dict, List

from metadynamic.inputs import DotParam


class Scaler:
    def __init__(
        self,
        data: Any,
        minimal: float,
        maximal: float,
        cutoff: float = 0.0,
        powerscale: float = 1.0,
    ):
        self.minval, self.maxval = self.minmax(data, cutoff)
        self.minimal = minimal
        self.maximal = maximal
        self.powerscale = powerscale

    def __call__(self, value: float) -> float:
        n = self.powerscale
        return (
            self.maximal * (value ** n - self.minval ** n)
            + self.minimal * (self.maxval ** n - value ** n)
        ) / (self.maxval ** n - self.minval ** n)

    @staticmethod
    def minmax(data: Any, cutoff: float = 0.0) -> Tuple[float, float]:
        try:
            maxval = max(data)
            data = np.array(data)
            data = data[data > maxval * cutoff]
            minval = min(data)
            return minval, maxval
        except ValueError:
            return 0.0, 1.0  # generic values for empty data


class Graphwriter:
    def __init__(self, margin: float, concentrate: bool, maxsize: float) -> None:
        self.dot = Digraph(
            graph_attr={
                "margin": str(margin),
                "concentrate": str(concentrate).lower(),
                "size": str(maxsize),
            }
        )

    def compound(
        self,
        name: str,
        width: float,
        fontsize: float,
        color: str,
        margin: float,
        penwidth: float,
    ) -> None:
        self.dot.node(
            name,
            shape="square",
            margin=str(margin),
            style="rounded",
            width=str(width),
            fontsize=str(fontsize),
            color=color,
            penwidth=str(penwidth),
        )

    def reaction(self, name: str, width: float, color: str) -> None:
        self.dot.node(name, shape="point", width=str(width), color=color)

    def edge(self, start: str, end: str, width: float, color: str) -> None:
        self.dot.edge(start, end, penwidth=str(width), color=color)

    def render(self, filename: str, engine: str = "dot", view: bool = False) -> bool:
        try:
            _, export = path.splitext(filename)
            export = export[1:]
        except ValueError:
            export = "dot"
        if export == "dot":
            with open(filename, "w") as out:
                out.write(self.dot.source)
        else:
            self.dot.engine = engine
            self.dot.format = export
            try:
                self.dot.render(filename, view=view, cleanup=True)
            except CalledProcessError:
                return False
        return True


class Json2dot:
    def __init__(self, filename: str, parameterfile: str = ""):
        with open(filename) as infile:
            data = load(infile)
        compounds = data["Compounds"]
        reactions = data["Reactions"]
        self.converter = Data2dot(compounds, reactions, parameterfile)

    def write(self, outfilename: str) -> bool:
        return self.converter.write(outfilename)


class Data2dot:
    def __init__(
        self,
        compdict: Dict[str, int],
        reacdict: Dict[str, List[float]],
        parameterfile: str = "",
    ):
        self.compounds = compdict
        self.reactions = reacdict
        self.param = DotParam.readfile(parameterfile) if parameterfile else DotParam()
        self.crn = Graphwriter(
            margin=self.param.margin,
            concentrate=self.param.concentrate,
            maxsize=self.param.maxsize,
        )
        binode = self.param.binode
        compset = set()
        reacset = set()
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
                reacset.add(name)
                width = scaler(rate)
                reactants, products = name.split("->")
                if not binode:
                    reaclist = []
                    prodlist = []
                for reac in reactants.split("+"):
                    num, reacname = self.cutdown(reac)
                    compset.add(reacname)
                    for _ in range(num):
                        if binode:
                            self.crn.edge(
                                start=reacname, end=name, width=width, color=color
                            )
                        else:
                            reaclist.append(reacname)
                for prod in products.split("+"):
                    num, prodname = self.cutdown(prod)
                    compset.add(prodname)
                    for _ in range(num):
                        if binode:
                            self.crn.edge(
                                start=name, end=prodname, width=width, color=color
                            )
                        else:
                            prodlist.append(prodname)
                if not binode:
                    for reacname, prodname in product(reaclist, prodlist):
                        self.crn.edge(
                            start=reacname, end=prodname, width=width, color=color
                        )
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
        for name in compset:
            try:
                pop = self.compounds[name]
                width = scaler(pop)
                fontsize = f_scaler(pop)
            except KeyError:
                width = 0
                fontsize = 0
            self.crn.compound(
                name=name,
                width=width,
                fontsize=fontsize,
                color=color,
                margin=self.param.c_margin,
                penwidth=self.param.c_penwidth,
            )
        if binode:
            color = self.param.r_color
            scaler = Scaler(
                data=[const for const, _ in self.reactions.values()],
                minimal=self.param.min_r_width,
                maximal=self.param.max_r_width,
                powerscale=self.param.r_powerscale,
            )
            for name in reacset:
                const, _ = self.reactions[name]
                width = scaler(const)
                self.crn.reaction(name=name, width=width, color=color)

    def write(self, filename: str) -> bool:
        return self.crn.render(filename)

    def view(self, filename: str) -> None:
        self.crn.dot.view()

    @staticmethod
    def cutdown(name: str) -> Tuple[int, str]:
        num = ""
        comp = ""
        start = True
        for char in name:
            if start and char.isdigit():
                num += char
            else:
                comp += char
                start = False
        retnum = 1 if num == "" else int(num)
        return retnum, comp