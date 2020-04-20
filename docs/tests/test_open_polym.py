from metadynamic import System, LOGGER


def test_system() -> None:
    syst = System.fromjson("docs/tests/open_polym.json")
    res = syst.run()
    LOGGER.info(f"Finished: {res}")
    LOGGER.disconnect()
