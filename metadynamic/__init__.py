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

    - L{launch}: function for launching a simulation from a .json parameter file, a .hdf5 result
      file from a previous run, with eventual additional parameters, and storing the results in a
      new .hdf5 file.

    - L{ResultReader}: class for reading and extracting data from a .hdf5 result file.

    - L{System}: class for creating, running, and directly controlling a simmulation.

    - L{Crn}: class describing a Chemical Reaction Network.

    - L{LOGGER}: global object for logging messages

    - L{MPI_STATUS}: global object for getting information about the MPI status


CLI interaction
---------------

A script named 'metarun' is provided for launching a simulation from the command line.

Usage::

    metarun [-h] [-c comment] [-d logdir] [-l log_level] [-x compress_level] param_file

    Launch a metadynamic run from a json or hdf5 file

    positional arguments:
      param_file            parameter json file, or previous hdf5 result file

    optional arguments:
      -h, --help            show this help message and exit
      -c comment, --comment comment
                            run comment string
      -d logdir, --logdir logdir
                            Directory for saving log text files
                            (if set, will override the value from param_file)
      -l log_level, --loglevel log_level
                            Logging level: DEBUG, INFO, WARNING, ERROR
                            (if set, will override the value from param_file)
      -x compress_level, --compress compress_level
                            compression level for the hdf5 output file

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
