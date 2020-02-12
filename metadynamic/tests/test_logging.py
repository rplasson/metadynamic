from metadynamic import System


def test_system() -> None:
    syst = System(
        "metadynamic/tests/verysimple.json",
        logfile="testlog/test-logging.log",
        loglevel="DEBUG",
    )
    res = syst.run()
    syst.log.info(f"Finished: {res}")
    syst.log.disconnect()
