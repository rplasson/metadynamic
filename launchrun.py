#!/usr/bin/env -S python3 -O

from argparse import ArgumentParser
from subprocess import call
from os import rename

from metadynamic import launch
from metadynamic.inputs import Param
from metadynamic.mpi import Parallel
from tempfile import NamedTemporaryFile

parser = ArgumentParser(description="Launch run from a json file")

parser.add_argument("parameters", type=str, help="parameter json file")
parser.add_argument("--log", metavar="logfile", type=str, nargs="?", help="log file")
parser.add_argument(
    "--comment",
    metavar="comment",
    type=str,
    nargs="?",
    help="comments on the run",
    default="",
)
parser.add_argument(
    "--level", metavar="loglevel", type=str, nargs="?", help="log level", default="INFO"
)
parser.add_argument(
    "--compress", metavar="compress", type=str, nargs="?", help="log level", default="GZIP=9"
)

args = parser.parse_args()

param = Param.readfile(args.parameters)
# param.set_param(dropmode="drop")
# param.set_param(init={"a": 1500, "A": 3000}, maxsteps=100000)

paramfile = NamedTemporaryFile()

param.tojson(paramfile.name)

res = launch(
    paramfile.name, logfile=args.log, loglevel=args.level, comment=args.comment
)

if args.compress != "no" and Parallel.mpi.rank == 0:
    try:
        old = param.hdf5
        new = old+".repacked"
        call(["h5repack", "-f", args.compress, old, new])
        rename(new, old)
    except FileNotFoundError:
        print("Couldn't find 'h5repack' utility, outputfile is uncompressed")

paramfile.close()

if Parallel.mpi.rank == 0:
    print(res.printinfo)
