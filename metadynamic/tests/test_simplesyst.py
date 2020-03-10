from metadynamic import System, LOGGER


def test_system() -> None:
    syst = System(
        "metadynamic/tests/simplesyst.json", logfile="testlog/test-keepreac.log"
    )
    syst.set_param(dropmode="keep")
    res = syst.run()
    LOGGER.info(f"Finished: {res}")
    LOGGER.disconnect()
