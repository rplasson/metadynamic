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

syst = System("aped.json", logfile=args.log, loglevel=args.level)
table, lendist, pooldist, the_end = syst.multirun(4)

table.to_csv(args.output, sep=",")
print(the_end)
