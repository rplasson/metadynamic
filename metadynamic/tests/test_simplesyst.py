from metadynamic import System


def test_system() -> None:
    syst = System(
        "metadynamic/tests/simplesyst.json", logfile="testlog/test-keepreac.log"
    )
    syst.set_param(dropmode="keep")
    res = syst.run()
    res.table().to_csv("testlog/test-keepreac-result.txt", sep=",")
    syst.log.info(f"Finished: {res.end()}")
    syst.log.disconnect()
