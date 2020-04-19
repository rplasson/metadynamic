from metadynamic import System, LOGGER


def test_system() -> None:
    syst = System.fromjson("metadynamic/tests/aped.json")
    res = syst.run()
    LOGGER.info(f"Finished: {res}")
    LOGGER.disconnect()
