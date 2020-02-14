#!/usr/bin/env -S python3 -O

from argparse import ArgumentParser

from metadynamic import launch
from metadynamic.inputs import Param
from tempfile import NamedTemporaryFile

parser = ArgumentParser(description="Launch run from a json file")

parser.add_argument("parameters", type=str, help="parameter json file")
parser.add_argument("--log", metavar="logfile", type=str, nargs="?", help="log file")
parser.add_argument(
    "--level", metavar="loglevel", type=str, nargs="?", help="log level", default="INFO"
)

args = parser.parse_args()

param = Param.readfile(args.parameters)
param.set_param(dropmode="drop")
param.set_param(init={"a": 15000, "A": 30000}, maxsteps=100000)

paramfile = NamedTemporaryFile()

param.tojson(paramfile.name)

res = launch(paramfile.name, logfile=args.log, loglevel=args.level)

paramfile.close()

print(res.end[...])
