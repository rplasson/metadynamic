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

from metadynamic import System
from metadynamic.hdf5 import ResultReader


def launch(parameters: str, logfile: str, loglevel: str, comment: str = "") -> ResultReader:
    syst = System(parameters, logfile=logfile, loglevel=loglevel, comment=comment)
    syst.run()
    return ResultReader(syst.param.hdf5)