from metadynamic import System, LOGGER


def test_system() -> None:
    syst = System(
        "metadynamic/tests/simple.json",
    )
    res = syst.run()
    LOGGER.info(f"Finished: {res}")
    LOGGER.disconnect()
