from metadynamic import System


def test_system() -> None:
    syst = System("metadynamic/tests/verysimple.json", logfile="testlog/test-logging.log", loglevel="DEBUG")
    res = syst.run()
    res.table().to_csv("testlog/test-logging-result.txt", sep=",")
    syst.log.info(f"Finished: {res.end()}")
    syst.log.disconnect()
