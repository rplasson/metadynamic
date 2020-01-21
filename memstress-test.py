from metadynamic import System

print("Loading...")
syst = System("memstress.json", logfile="testlog/memstress-test.log")
print("Set parameters...")
syst.log.info("Set parameters...")
syst.set_param(tend=4.0)
syst.set_param(dropmode="keep")
print("Running...")
syst.log.info("Running...")
res = syst.run()
syst.log.info("Finished")
print("Finished")
print(res.end())
