from metadynamic import System


def test_system() -> None:
    syst = System("metadynamic/tests/simplesyst-keepreac.json", logfile="testlog/test-keepreac.log")
    table, lendist, pooldist, the_end = syst.run()
    table.to_csv("testlog/test-keepreac-result.txt", sep=",")
    print(the_end)
