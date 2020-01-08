from metadynamic import System
from argparse import ArgumentParser

parser = ArgumentParser(description="Test multithread.")

parser.add_argument("output", type=str, help="output file (csv)")
parser.add_argument("--log", metavar="logfile", type=str, nargs="?", help="log file")
parser.add_argument(
    "--level", metavar="loglevel", type=str, nargs="?", help="log level", default="INFO"
)

args = parser.parse_args()

syst = System("aped.json", logfile=args.log, loglevel=args.level)
syst.runparam.set_param(dropmode="keep")
syst.set_param(init={"a": 15000, "A": 30000}, maxsteps=100000)
res = syst.run()

res.table().to_csv(args.output, sep=",")
print(res.end())