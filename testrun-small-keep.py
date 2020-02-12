from metadynamic import System
from argparse import ArgumentParser

parser = ArgumentParser(description="Test multithread.")

parser.add_argument("--log", metavar="logfile", type=str, nargs="?", help="log file")
parser.add_argument(
    "--level", metavar="loglevel", type=str, nargs="?", help="log level", default="DEBUG"
)

args = parser.parse_args()


syst = System("small.json", logfile=args.log, loglevel=args.level)
syst.set_param(dropmode="keep")
res = syst.run()

print(res)
