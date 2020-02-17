from metadynamic import System


def test_system() -> None:
    syst = System(
        "metadynamic/tests/simplesyst.json", logfile="testlog/test-keepreac.log"
    )
    syst.set_param(dropmode="keep")
    res = syst.run()
    syst.log.info(f"Finished: {res}")
    syst.log.disconnect()
