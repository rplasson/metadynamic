from metadynamic import System


def test_system() -> None:
    syst = System("metadynamic/tests/verysimple.json", logfile="testlog/test-logging.log", loglevel="DEBUG")
    table, lendist, pooldist, the_end = syst.run()
    table.to_csv("testlog/test-logging-result.txt", sep=",")
    print(the_end)
