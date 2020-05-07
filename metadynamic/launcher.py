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

"""High level interface to direct launch a simulation run.

It provides L{launch}, a function for launching a simulation from a .json parameter file, or a .hdf5
result file from a previous run, with eventual additional parameters, and storing the results in a
new .hdf5 file.

"""

from os import path
from typing import Any

from metadynamic.system import System
from metadynamic.result import ResultReader


def launch(parameters: str, **kwd: Any) -> ResultReader:
    """Launch a metadynamic run.

    It saves the result in a .hdf5 file, and returns a ResulReader object.

    It can be launched either form a .json Param file
    or from a previous .hdf5 result file.

    @param parameters: name of the parameter file (.json or .hdf5)
    @type parameters: str
    @param kwd: additional parameters (override the one defined in parameters file)
    @return: Interface object to the generated .hdf5 file.
    @rtype: ResultReader

    """
    ext = path.splitext(parameters)[-1]
    if ext == ".json":
        syst = System.fromjson(parameters, **kwd)
    elif ext in (".hdf5", ".h5"):
        syst = System.fromhdf5(parameters, **kwd)
    syst.run()
    return ResultReader(syst.output.h5file)
