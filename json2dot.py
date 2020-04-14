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

from metadynamic.network import Data2dot
from argparse import ArgumentParser

parser = ArgumentParser(description="Test multithread.")

parser.add_argument("file_in", type=str, help="input file (json)")
parser.add_argument("file_out", type=str, help="output file (dot)")
parser.add_argument("--param", metavar="param", type=str, help="output file (dot)", default="", nargs="?")

args = parser.parse_args()

Data2dot.fromjson(args.file_in, parameterfile=args.param).write(args.file_out)
