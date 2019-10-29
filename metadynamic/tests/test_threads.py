from metadynamic import System


def test_system() -> None:
    syst = System({"a": 1000, "A": 2000}, consts={"P": 10.0, "A": 10.0, "H": 1.0, "E": 1.0, "a": 1.0, "d": 1.0, "R": 0.1}, dropreac=False)
    syst.set_param(conc=1.0, tend=1.0, tstep=0.05)
    syst.set_param(save=["a", "A", "aa", "aA", "Aa", "AA"])
    res = syst.multirun(4)
    res.table().to_csv("testlog/test-thread.log", sep=",")
    print(res.the_end())
