from metadynamic import System


def test_system() -> None:
    syst = System("metadynamic/tests/simplesyst.json", logfile="testlog/test-thread.log")
    syst.runparam.set_param(nbthread=-1)
    res = syst.run()
    res.table().to_csv("testlog/test-thread-result.txt", sep=",")
    syst.log.info(f"Finished: {res.end()}")
    syst.log.disconnect()
