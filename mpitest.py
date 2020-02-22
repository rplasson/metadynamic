from metadynamic.mpi import MpiStatus
from random import random, choice
from time import sleep


def logmsg(msg: str):
    print(f"#{mpi.rank} -> {msg}")


def op1():
    logmsg("Operation 1")


def op2():
    logmsg("Operation 2")


mpi = MpiStatus()
mpi.gate.register_function("op1", op1)
mpi.gate.register_function("op2", op2)

logmsg("Entering")
for i in range(30):
    sleep(1 + random())
    if random() > 0.9:
        msg = choice(["op1", "op2"])
        logmsg(f"Send message {msg}")
        mpi.gate.close(msg)
    sleep(random())
    mpi.gate.checkpoint()
logmsg("Exiting...")
mpi.gate.exit()
logmsg("... I am free!")
