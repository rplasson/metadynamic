from metadynamic import System, LOGGER


def test_system() -> None:
    syst = System.fromjson("metadynamic/tests/polym_minimal.json")
    res = syst.run()
    LOGGER.info(f"Finished: {res}")
    LOGGER.disconnect()
