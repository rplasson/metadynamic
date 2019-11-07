from metadynamic import System


def test_system() -> None:
    syst = System("metadynamic/tests/simplesyst-keepreac.json", logfile="testlog/test-keepreac.log")
    res = syst.run()
    res.table().to_csv("testlog/test-keepreac-result.txt", sep=",")
    print(res.end())
