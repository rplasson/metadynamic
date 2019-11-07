from metadynamic import System


def test_system() -> None:
    syst = System(
        "metadynamic/tests/simplesyst.json", logfile="testlog/test-dropreac.log"
    )
    syst.runparam.set_param(dropmode="drop")
    res = syst.run()
    res.table().to_csv("testlog/test-dropreac-result.txt", sep=",")
    syst.log.info(f"Finished: {res.end()}")
    syst.log.disconnect()
