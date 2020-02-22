from metadynamic.mpi import MpiGate
from random import random, choice
from time import sleep


def logmsg(msg: str):
    print(f"#{gate.rank} -> {msg}")


def op1():
    logmsg("Operation 1")


def op2():
    logmsg("Operation 2")


def dostuff():
    sleep(1 + random())
    if random() > 0.9:
        msg = choice(["op1", "op2"])
        logmsg(f"Send message {msg}")
        gate.close(msg)
    sleep(random())


with MpiGate(taginit=100, operations={"op1": op1, "op2": op2}) as gate:
    logmsg("Entering")
    for i in range(30):
        dostuff()
        gate.checkpoint()
logmsg("... I am free!")
