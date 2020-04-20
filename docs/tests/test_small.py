from metadynamic import System, LOGGER


def test_system() -> None:
    syst = System.fromjson("docs/tests/small.json")
    res = syst.run()
    LOGGER.info(f"Finished: {res}")
    LOGGER.disconnect()
