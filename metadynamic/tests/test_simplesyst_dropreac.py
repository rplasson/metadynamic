from metadynamic import System


def test_system() -> None:
    syst = System(
        "metadynamic/tests/simplesyst.json", logfile="testlog/test-dropreac.log"
    )
    res = syst.run()
    res.table().to_csv("testlog/test-dropreac-result.txt", sep=",")
    print(res.end())
