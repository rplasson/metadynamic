from polym import System
from numpy import savetxt
from argparse import ArgumentParser

parser = ArgumentParser(description="Test multithread.")

parser.add_argument("output", metavar="O", type=str, help="output file (csv)")
parser.add_argument("--log", metavar="logfile", type=str, nargs="?", help="log file")
parser.add_argument("--level", metavar="loglevel", type=str, nargs="?", help="log level", default="INFO")

args = parser.parse_args()

syst = System(
    {"a": 500, "A":500, "b": 500, "B": 500},
    consts={"A": 100.0, "H": 1.0, "E":0.1, "a":1., "d":0.1},
    altconsts={"A": 10.0, "H":10.0, "E":10.0, "a": 0.},
    conc=1.0, seed = 1234, dropreac=False,
    logfile = args.log, loglevel = args.level,
)

syst.set_run(tend=10, tstep=1)
syst.set_run(save=["a", "b", "A", "B", "a*", "b*", "A*", "B*", "aa", "aA", "Aa", "AA", "aa*"])

table, lendist, pooldist, the_end = syst.run()

table.to_csv(args.output, sep=",")

print(syst.reac_collect.pool)
print(syst.reac_collect.active)
print(the_end)
for cat in ["A", "H", "E", "a", "d"]:
    print(cat,syst.reac_collect.cat_list(cat))

for cat in ["actpol", "polym", "mono"]:
    print(cat,syst.comp_collect.cat_list(cat))
