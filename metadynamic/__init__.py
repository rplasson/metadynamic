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
Metadynamic stochastic simulator
================================

Simulation of chemical systems based on Gillespie's Stochastic algorithm, with on-the-fly
creation/destruction of compounds and reactions.


Provides
--------

launch: function for launching a simulation from a .json parameter file, a .hdf5 result file from a
previous run, with eventual additional parameters, and storing the results in a new .hdf5 file.

ResultReader: class for reading and extracting data from a .hdf5 result file.

System: class for creating, running, and directly controlling a simmulation.

Crn: class describing a Chemical Reaction Network.

LOGGER: global object for logging messages

MPI_STATUS: global object for getting information about the MPI status


CLI interaction
---------------

A script named 'metarun' is provided for launching a simulation from the command line.

usage: metarun [-h] [--logdir [logdir]] [--comment [comment]] [--loglevel [loglevel]] [--compress [compress]] parameters

"""


from metadynamic.version import __version__
from metadynamic.system import System
from metadynamic.chemical import Crn
from metadynamic.result import ResultReader
from metadynamic.launcher import launch
from metadynamic.logger import LOGGER
from metadynamic.mpi import MPI_STATUS

# if somebody does "from somepackage import *", this is what they will
# be able to access:
__all__ = [
    "System",
    "Crn",
    "ResultReader",
    "launch",
    "LOGGER",
    "MPI_STATUS",
    "__version__",
]
