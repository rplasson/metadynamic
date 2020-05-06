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

"""Simulation of chemical systems based on Gillespie's Stochastic algorithm.

This package is for modeling metadynamic systems, with on-the-fly creation/destruction of compounds
and reactions.

- A simulation can be simply performed using the function L{launch}, or totally controlled as a
  L{System} object.
- Results can be read from a L{ResultReader} object.
- Chemical reaction networks can be analyzed from L{Crn} objects.
- L{LOGGER} and L{MPI_STATUS} global objects are available for logging and MPI facilities.

"""


from metadynamic.version import __version__
from metadynamic.system import System
from metadynamic.chemical import Crn
from metadynamic.result import ResultReader
from metadynamic.launcher import launch
from metadynamic.logger import LOGGER
from metadynamic.mpi import MPI_STATUS


__all__ = [
    "System",
    "Crn",
    "ResultReader",
    "launch",
    "LOGGER",
    "MPI_STATUS",
    "__version__",
]
"""To be imported by 'from metadynamic import *'"""
