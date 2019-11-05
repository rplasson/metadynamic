from metadynamic import System


def test_system() -> None:
    syst = System("metadynamic/tests/simplesyst.json", logfile="testlog/test.log")
    table, lendist, pooldist, the_end = syst.run()
    table.to_csv("testlog/test-result.txt", sep=",")
    print(the_end)
