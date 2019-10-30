from metadynamic import System


def test_system() -> None:
    syst = System(
        {"a": 3, "A": 6},
        consts={"P": 10.0, "A": 10.0, "H": 1.0, "E": 1.0, "a": 1.0, "d": 1.0, "R": 0.1},
        dropreac=True,
        logfile="testlog/test_logging.log",
        loglevel="DEBUG"
    )
    syst.set_param(conc=1.0, tend=1.0, tstep=0.05)
    syst.set_param(save=["a", "A", "aa", "aA", "Aa", "AA"])
    table, lendist, pooldist, the_end = syst.run()
    table.to_csv("testlog/test_logging-result.txt", sep=",")
    print(the_end)
