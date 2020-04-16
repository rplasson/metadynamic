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
metadynamic.outputs
===================

deal with output files organization


Provides:
---------

 - L{Output}: access to output files and folders

"""

from os import path

from metadynamic.inputs import Param
from metadynamic.mpi import MPI_STATUS
from metadynamic.logger import LOGGER
from metadynamic.ends import NotAFolder


class Output:
    """Access to output files and folders"""

    def __init__(self, param: Param):
        """
        Generate files structures from parameters

        @param param: parameters
        @type param: Param
        """
        self.name: str = param.name
        """run name"""
        # self.timeformat = param.timeformat
        self.savedir: str = param.savedir if param.savedir else path.curdir
        """save folder"""
        if not path.isdir(self.savedir):
            raise NotAFolder("Bad folder name {self.savedir}")
        self.logdir: str = param.logdir
        """log folder"""
        if not path.isdir(self.logdir):
            LOGGER.debug(
                "Logs will be sent to standard output because {self.logdir} is not a folder."
            )

    @property
    def process(self) -> str:
        """
        return '-p<n>' where n is the process number in a MPI run,
        or an empty string if the run is not MPI_GATE

        @return: process string
        @rtype: str
        """
        if MPI_STATUS.ismpi:
            return f"-p{MPI_STATUS.rank}"
        return ""

    @property
    def basename(self) -> str:
        """
        file basename, as <run name>-<hostname>-<start time>

        @return: file basename
        @rtype: str
        """
        return f"{self.name}-{MPI_STATUS.hostname}-{MPI_STATUS.starttime}"

    @property
    def h5file(self) -> str:
        """
        hdf5 file name (with full path)

        @return: hdf5 file name
        @rtype: str
        """
        return path.join(self.savedir, f"{self.basename}.hdf5")

    @property
    def logfile(self) -> str:
        """
        log file name (with full path)
        empty string if log to standard output

        @return: log file name
        @rtype: str
        """
        if self.logdir and path.isdir(self.logdir):
            return path.join(self.logdir, f"{self.basename}{self.process}.log")
        return ""
