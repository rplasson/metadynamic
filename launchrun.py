#!/usr/bin/env -S python3 -O

from argparse import ArgumentParser
from subprocess import call
from os import rename

from metadynamic import launch
from metadynamic.inputs import Param
from metadynamic import MPI_STATUS
from tempfile import NamedTemporaryFile

parser = ArgumentParser(description="Launch run from a json file")

parser.add_argument("parameters", type=str, help="parameter json file")
parser.add_argument(
    "--logdir", metavar="logdir", type=str, nargs="?", help="log file", default=""
)
parser.add_argument(
    "--comment",
    metavar="comment",
    type=str,
    nargs="?",
    help="comments on the run",
    default="",
)
parser.add_argument(
    "--loglevel", metavar="loglevel", type=str, nargs="?", help="log level", default=""
)
parser.add_argument(
    "--compress",
    metavar="compress",
    type=str,
    nargs="?",
    help="log level",
    default="GZIP=9",
)

args = parser.parse_args()

if MPI_STATUS.root:
    print(f"Launched run '{args.comment}' on {MPI_STATUS.size} processes...")

param = Param.readfile(args.parameters)
# param.set_param(dropmode="drop")
# param.set_param(init={"a": 1500, "A": 3000}, maxsteps=100000)

if args.comment:
    param.set_param(comment=args.comment)
if args.logdir:
    param.set_param(logdir=args.logdir)
if args.loglevel:
    param.set_param(loglevel=args.loglevel)


paramfile = NamedTemporaryFile()

param.tojson(paramfile.name)

res = launch(paramfile.name)

if args.compress != "no" and MPI_STATUS.root:
    try:
        old = res.filename
        new = old + ".repacked"
        call(["h5repack", "-f", args.compress, old, new])
        rename(new, old)
    except FileNotFoundError:
        print("Couldn't find 'h5repack' utility, outputfile is uncompressed")

paramfile.close()

if MPI_STATUS.root:
    print(res.printinfo)
