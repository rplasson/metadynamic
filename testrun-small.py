from metadynamic import System
from numpy import savetxt
from argparse import ArgumentParser

parser = ArgumentParser(description="Test multithread.")

parser.add_argument("output", type=str, help="output file (csv)")
parser.add_argument("--log", metavar="logfile", type=str, nargs="?", help="log file")
parser.add_argument(
    "--level", metavar="loglevel", type=str, nargs="?", help="log level", default="INFO"
)

args = parser.parse_args()

syst = System(
    {"a": 50, "A": 50, "b": 100, "c": 100},
    consts={"P": 10.0, "H": 1.0, "E": 1.0},
    altconsts={"P": 10.0},
    catconsts={"P": 10000.0, "H": 10000.0},
    dropreac=True,
    logfile=args.log,
    loglevel=args.level,
)

syst.set_param(conc=1.0, seed=1234)
syst.set_param(tend=2, tstep=0.1)
syst.set_param(save=["a", "b", "c", "A", "B", "C"])

table, lendist, pooldist, the_end = syst.run()

table.to_csv(args.output, sep=",")
print(the_end)
