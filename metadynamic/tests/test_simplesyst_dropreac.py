from metadynamic import System, LOGGER


def test_system() -> None:
    syst = System(
        "metadynamic/tests/simplesyst.json", logfile="testlog/test-dropreac.log"
    )
    syst.set_param(dropmode="drop")
    syst.set_param(name="testlog/simplesyst-dropreac")
    res = syst.run()
    LOGGER.info(f"Finished: {res}")
    LOGGER.disconnect()
