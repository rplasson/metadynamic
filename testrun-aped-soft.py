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


a = 1.0
p = 1.0
h = 1.0
e = 1.0
alpha = 0.1
beta = 0.0
gamma = 0.0
n0 = 4000
dn0 = 1000

syst = System(
    {"a": n0 / 2 + dn0, "A": n0 / 2 - dn0},
    consts={"A": p, "H": h, "E": e, "a": a},
    altconsts={"A": alpha, "H": beta, "E": gamma, "a": 0.0},
    dropmode="soft",
    logfile=args.log,
    loglevel=args.level,
)

max_time_min = 10
end = 50
step = 1
c0 = 30.0
rtlim = 60 * max_time_min

syst.set_param(conc=c0,  tend=end, tstep=step, rtlim=rtlim)
syst.set_param(save=["a", "aa", "aaa", "aaaa", "A", "Aa", "aA", "AA", "a*", "A*"])
syst.set_param(seed=1234)

table, lendist, pooldist, the_end = syst.run()

table.to_csv(args.output, sep=",")
print(the_end)
