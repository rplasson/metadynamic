from metadynamic import System


syst = System("memstress.json", logfile="testlog/memstress-test.log")
syst.set_param(tend=1.)
syst.runparam.set_param(dropmode="keep")
res = syst.run()
print(res.end())

