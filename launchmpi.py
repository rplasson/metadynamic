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

from mpi4py import MPI
import h5py
from argparse import ArgumentParser

from metadynamic import System

parser = ArgumentParser(description="Test multithread.")

parser.add_argument("output", type=str, help="output file (hdf5)")
parser.add_argument("--log", metavar="logfile", type=str, nargs="?", help="log file")
parser.add_argument(
    "--level", metavar="loglevel", type=str, nargs="?", help="log level", default="INFO"
)

args = parser.parse_args()

rank = MPI.COMM_WORLD.rank
size = MPI.COMM_WORLD.size

logfile = args.log if size == 1 else f"-{rank}_{size}.".join(args.log.split("."))

syst = System("aped-nosnap.json", logfile=logfile, loglevel=args.level)
syst.set_param(dropmode="drop")
res = syst.run()

data = res.table()
lines, cols = data.shape

out = h5py.File(args.output, "w", driver="mpio", comm=MPI.COMM_WORLD)

dataset = out.create_dataset("data", (size, lines, cols), maxshape=(size, lines, None))
dataset[rank] = res.table()

out.close()
print(res.end())
