from metadynamic import System
from argparse import ArgumentParser

parser = ArgumentParser(description="Test multithread.")

parser.add_argument("--log", metavar="logfile", type=str, nargs="?", help="log file")
parser.add_argument(
    "--level", metavar="loglevel", type=str, nargs="?", help="log level", default="INFO"
)

args = parser.parse_args()

syst = System("aped.json", logfile=args.log, loglevel=args.level)
syst.param.set_param(dropmode="drop")
syst.param.set_param(nbthread=-1)
res = syst.run()

print(res)
