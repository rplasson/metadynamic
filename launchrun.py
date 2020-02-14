#!/usr/bin/env -S python3 -O

from argparse import ArgumentParser

from metadynamic import launch


parser = ArgumentParser(description="Launch run from a json file")

parser.add_argument("parameters", type=str, help="parameter json file")
parser.add_argument("--log", metavar="logfile", type=str, nargs="?", help="log file")
parser.add_argument(
    "--level", metavar="loglevel", type=str, nargs="?", help="log level", default="INFO"
)

args = parser.parse_args()

# syst.set_param(dropmode="drop")
# syst.set_param(init={"a": 15000, "A": 30000}, maxsteps=100000)


res = launch(args.parameters, logfile=args.log, loglevel=args.level)

print(res.end[...])
