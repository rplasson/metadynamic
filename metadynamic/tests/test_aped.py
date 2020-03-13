from metadynamic import System, LOGGER


def test_system() -> None:
    syst = System(
        "metadynamic/tests/aped.json",
    )
    res = syst.run()
    LOGGER.info(f"Finished: {res}")
    LOGGER.disconnect()
