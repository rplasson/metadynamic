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


from os import path
from datetime import datetime
from socket import gethostname

from metadynamic.inputs import Param
from metadynamic.mpi import MPI_STATUS
from metadynamic.logger import LOGGER
from metadynamic.ends import NotAFolder


class Output:
    def __init__(self, param: Param):
        self.name = param.name
        self.timeformat = param.timeformat
        self.savedir = param.savedir if param.savedir else path.curdir
        if not path.isdir(self.savedir):
            raise NotAFolder("Bad folder name {self.savedir}")
        self.logdir = param.logdir
        self.loginfile = self.logdir and path.isdir(self.logdir)
        if not path.isdir(self.logdir):
            LOGGER.debug("Logs will be sent to standard output because {self.logdir} is not a folder.")

    @property
    def hostname(self) -> str:
        return gethostname()

    @property
    def datetime(self) -> str:
        if not hasattr(self, "_datetime"):
            self._datetime = datetime.now().strftime(self.timeformat)
        return self._datetime

    @property
    def process(self) -> str:
        if MPI_STATUS.ismpi:
            return f"-p{MPI_STATUS.rank}"
        return ""

    @property
    def basename(self) -> str:
        return f"{self.name}-{self.hostname}-{self.datetime}"

    @property
    def h5file(self) -> str:
        return path.join(self.savedir, f"{self.basename}.hdf5")

    @property
    def logfile(self) -> str:
        if self.loginfile:
            return path.join(self.logdir, f"{self.basename}{self.process}.log")
        return ""
