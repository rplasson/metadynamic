from metadynamic import System
from numpy import savetxt
from argparse import ArgumentParser

parser = ArgumentParser(description="Test multithread.")

parser.add_argument("output", type=str, help="output file (csv)")
parser.add_argument("--log", metavar="logfile", type=str, nargs="?", help="log file")
parser.add_argument(
    "--level", metavar="loglevel", type=str, nargs="?", help="log level", default="DEBUG"
)

args = parser.parse_args()


syst = System("test.json", logfile=args.log, loglevel=args.level)
syst.runparam.set_param(dropmode="drop")
res = syst.run()

res.table().to_csv(args.output, sep=",")
print(res.end())