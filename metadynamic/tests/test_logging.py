from metadynamic import System, LOGGER


def test_system() -> None:
    syst = System(
        "metadynamic/tests/verysimple.json",
        logfile="testlog/test-logging.log",
        loglevel="DEBUG",
    )
    res = syst.run()
    LOGGER.info(f"Finished: {res}")
    LOGGER.disconnect()
