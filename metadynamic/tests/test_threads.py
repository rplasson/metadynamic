from metadynamic import System


def test_system() -> None:
    syst = System("metadynamic/tests/simplesyst.json", logfile="testlog/test-thread.log")
    res = syst.multirun(4)
    res.table().to_csv("testlog/test-thread-result.txt", sep=",")
    print(res.end())
