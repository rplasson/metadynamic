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

"""
metadynamic.network
===================

converts a set of compounds and reaction as a graphvize network representation.


Provides:
---------

 - L{Data2dot}: generate a graphviz Digraph from dict description of compounds/reactions

"""

from os import path
from itertools import product
from json import load
from subprocess import CalledProcessError
from typing import Any, Tuple, Dict, List
from graphviz import Digraph

import numpy as np

from metadynamic.inputs import DotParam


class Scaler:
    """convert values from one scale to another"""

    def __init__(
        self,
        data: Any,
        minimal: float,
        maximal: float,
        cutoff: float = 0.0,
        powerscale: float = 1.0,
    ):
        """
        Create a scaler from data as the origin dataset
        to [minimal, maximal] scale with powerscale exponent

        the origin scale is set from the minimal and maximal values
        of data, once outlier values smaller than cutoff×max(data) are removed.

        @param data: origin data reference to scale
        @param minimal: minimal destination scale value
        @type minimal: float
        @param maximal: maximal destination scale value
        @type maximal: float
        @param cutoff: cutoff fraction for removing small outliers (Default value = 0.0)
        @type cutoff: float
        @param powerscale: scaling exponent (Default value = 1.0)
        @type powerscale: float
        """
        self.minimal: float = minimal
        """min value of the destination scale"""
        self.maximal: float = maximal
        """max value of the destination scale"""
        self.maxval: float
        """max value of the origin scale"""
        self.minval: float
        """min value of the origin scale"""
        self._lambda: float
        """scaling factor"""
        self._const: float
        """scaling constant"""
        self.n = powerscale
        """scaling exponent"""
        try:
            self.maxval = max(data)
            data = np.array(data)
            data = data[data > self.maxval * cutoff]
            self.minval = min(data)
        except ValueError:  # generic values for empty data
            self.minval = 0.0
            self.maxval = 1.0
        self._lambda = (self.maximal - self.minimal) / (
            self.maxval ** self.n - self.minval ** self.n
        )
        self._const = (
            +self.minimal * self.maxval ** self.n - self.maximal * self.minval ** self.n
        ) / (self.maxval ** self.n - self.minval ** self.n)

    def __call__(self, value: float) -> float:
        """
        scale as _factor × value^n + _const

        the scaled value will range from minimal to maximal

        original values outside of the [minval, maxval] range
        will be scaled to either minimal or maximal

        @param value: value to be scaled
        @param type: float
        @return: scaled value
        @rtype: float
        """
        if value <= self.minval:
            return self.minimal
        if value >= self.maxval:
            return self.maximal
        return self._lambda * value ** self.n + self._const


class Graphwriter:
    """Interface for building a graphviz Digraph from a CRN"""

    def __init__(self, margin: float, concentrate: bool, maxsize: float) -> None:
        """
        Create the network

        @param margin: margin size
        @type margin: float
        @param concentrate: if True, try to reduce graph size
        @type concentrate: bool
        @param maxsize: graph size
        @type maxsize: float
        """
        self.dot: Digraph = Digraph(
            graph_attr={
                "margin": str(margin),
                "concentrate": str(concentrate).lower(),
                "size": str(maxsize),
            }
        )
        """graphviz representation of the CRN"""

    def compound(
        self,
        name: str,
        width: float,
        fontsize: float,
        color: str,
        margin: float,
        penwidth: float,
    ) -> None:
        """
        Add a compound node to the network

        @param name: compound name
        @type name: str
        @param width: compound width
        @type width: float
        @param fontsize: compound font size
        @type fontsize: float
        @param color: compound color
        @type color: str
        @param margin: compound margin size
        @type margin: float
        @param penwidth: compound box pen width
        @type penwidth: float
        """
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
        """
        Add a reaction node to the network

        @param name: reaction name
        @type name: str
        @param width: reaction width
        @type width: float
        @param color: reaction color
        @type color: str
        """
        self.dot.node(name, shape="point", width=str(width), color=color)

    def edge(self, start: str, end: str, width: float, color: str) -> None:
        """
        Add a chemical flux edge to the network

        @param start: node start name
        @type start: str
        @param end: node end name
        @type end: str
        @param width: edge width
        @type width: float
        @param color: edge color
        @type color: str
        """
        self.dot.edge(start, end, penwidth=str(width), color=color)

    def render(self, filename: str, engine: str = "dot", view: bool = False) -> bool:
        """
        render the network in a file.

        The output format is guessed from filename extension.

        @param filename: name of the output file
        @type filename: str
        @param engine: output engine (Default value = "dot")
        @type engine: str
        @param view: if True, a viewer will be opened after file creation (Default value = False)
        @type view: bool
        @return: return True if rendering is OK, False if a problem occured
        @rtype: bool
        """
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


class Data2dot:
    """generate a graphviz Digraph from dict description of compounds/reaction"""

    @classmethod
    def fromjson(cls, filename: str, parameterfile: str = "") -> "Data2dot":
        """
        generate a graphviz Digraph from a json file describing compounds/reactions

        @param filename: name of json file containing compiunds/reactions
        @type filename: str
        @param parameterfile: Name of json Dotparam file (if empty, use default values) (Default value = "")
        @type parameterfile: str
        """
        with open(filename) as infile:
            data = load(infile)
        compounds = data["Compounds"]
        reactions = data["Reactions"]
        return cls(compounds, reactions, parameterfile)

    def __init__(
        self,
        compdict: Dict[str, int],
        reacdict: Dict[str, List[float]],
        parameterfile: str = "",
    ):
        """
        Generate the network

        @param compdict: Chemical compounds {name:population}
        @type compdict: Dict[str, int]
        @param reacdict: Chemical reaction {name:(constant,probability)}
        @type reactdict: Dict[str, List[float]]
        @param parameterfile: Name of json Dotparam file (if empty, use default values) (Default value = "")
        @type parameterfile: str
        """
        self.compounds: Dict[str, int] = compdict
        """Chemical compounds {name:population}"""
        self.reactions: Dict[str, List[float]] = reacdict
        """Chemical reaction {name:(constant,probability)}"""
        self.param: DotParam = DotParam.readfile(
            parameterfile
        ) if parameterfile else DotParam()
        """Rendering parameters"""
        self.crn: Graphwriter = Graphwriter(
            margin=self.param.margin,
            concentrate=self.param.concentrate,
            maxsize=self.param.maxsize,
        )
        """Chemical Reaction Network renderer"""
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
        """
        Write the network in a file
        (format guessed from file extension)

        @param filename: name of the output file
        @type filename: str
        @return: return True if rendering is OK, False if a problem occured
        @rtype: bool
        """
        return self.crn.render(filename)

    def view(self) -> None:
        """Open a viewer to visualize the network"""
        self.crn.dot.view()

    @staticmethod
    def cutdown(name: str) -> Tuple[int, str]:
        """
        Identify the stoechiometry+name from a str

           >>> Data2dot.cutdown("ab3")
           (1, 'ab3')
           >>> Data2dot.cutdown("23ab3")
           (23, 'ab3')

        @param name: string to confert
        @type name: str
        @return: stoechiometry, compound name
        @rtype: Tuple[int, str]
        """
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
