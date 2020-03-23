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

from metadynamic import System, ResultReader
from os import path
from typing import Any


def launch(parameters: str, **kwd: Any) -> ResultReader:
    """launch a metadynamic run, saves the result in a .hdf5 file,
       and returns a ResulReader object

    :param parameters: name of the parameter file (.json or .hdf5)
    :param **kwd: additional parameters (override the one defined in parameters file)
    :type parameters: str
    :return: Interface object to the generated .hdf5 file.
    :rtype: ResultReader

    """
    ext = path.splitext(parameters)[-1]
    if ext == ".json":
        syst = System.fromjson(parameters, **kwd)
    elif ext == ".hdf5" or ext == ".h5":
        syst = System.fromhdf5(parameters, **kwd)
    syst.run()
    return ResultReader(syst.output.h5file)
