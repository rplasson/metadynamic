from polym import System
from numpy import savetxt
from argparse import ArgumentParser

parser = ArgumentParser(description="Test multithread.")

parser.add_argument("output", type=str, help="output file (csv)")
parser.add_argument("--log", metavar="logfile", type=str, nargs="?", help="log file")
parser.add_argument("--level", metavar="loglevel", type=str, nargs="?", help="log level", default="INFO")

args = parser.parse_args()

syst = System(
    {"a": 5, "A":5, "b": 5, "c": 5},
    consts={"P": 10.0, "H": 1.0, "E":1.0},
    altconsts={"P": 10.0},
    catconsts={"P":10000.0, "H": 10000.0},
    conc=1.0, seed = 1234, dropreac = False,
    logfile = args.log, loglevel = args.level,
)
syst.addkept("P.ba+ba.baba.")
syst.addkept("H.baba.baba.2")


syst.set_run(tend=1., tstep=0.05)
syst.set_run(save=["a", "b", "c", "A", "B", "C"])

table, lendist, pooldist, the_end = syst.run()

table.to_csv(args.output, sep=",")
print(the_end)
