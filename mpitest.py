from metadynamic.mpi import MpiGate
from random import random, choice
from time import sleep


def logmsg(msg: str):
    print(f"#{gate.rank} -> {msg}")


def op1():
    logmsg("Operation 1")


def op2():
    logmsg("Operation 2")


gate = MpiGate(taginit=100, operations={"op1": op1, "op2": op2})

with gate.context():
    logmsg("Entering")
    while gate.cont.ok:
        # Do stuff!
        sleep(1 + random())
        if random() > 0.9:
            msg = choice(["op1", "op2"])
            logmsg(f"Send message {msg}")
            gate.close(msg)
        if random() > 0.98:
            logmsg(f"Made a mistake")
            raise ValueError("Oups")
        if random() > 0.95:
            logmsg(f"Kill half !")
            gate.close("oom")
        sleep(random())
        # Stuff done
        gate.checkpoint()
logmsg("... I am free!")
