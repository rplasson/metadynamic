import gc
from multiprocessing import get_context
from itertools import repeat
from os import getpid
from typing import List, Optional, Dict, Tuple
from psutil import Process
from dataclasses import dataclass, replace, field
from json import load

from pandas import DataFrame

from metadynamic.ends import (
    Finished,
    TimesUp,
    RuntimeLim,
    NotFound,
    NoMore,
    HappyEnding,
    BadEnding,
)
from metadynamic.logger import Logged
from metadynamic.proba import Probalistic
from metadynamic.processing import Result
from metadynamic.chemical import (
    Collected,
    CollectofCompound,
    CollectofReaction,
    Compound,
    Reaction,
)


class Readerclass:
    _default_section: str = "Parameter"

    @classmethod
    def readfile(cls, filename: str, section: str = "") -> "Readerclass":
        """Return a SysParam object, updated by the data from filename"""
        if section == "":
            section = cls._default_section
        with open(filename) as json_data:
            parameters = load(json_data)[section]
        return cls(**parameters)


@dataclass
class SysParam(Readerclass):
    _default_section = "System"
    conc: float = 0.1
    tend: float = 1.0
    tstep: float = 0.01
    rtlim: float = 900.0  # Limit to 15min runtime
    maxsteps: int = 10000
    ptot: int = field(init=False)
    vol: float = field(init=False)
    seed: Optional[int] = None
    save: List[str] = field(default_factory=list)
    init: Dict[str, int] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.ptot = sum([pop * len(comp) for comp, pop in self.init.items()])
        self.vol = self.ptot / self.conc


@dataclass
class RunParam(Readerclass):
    _default_section = "Run"
    dropreac: bool = True
    autoclean: bool = True
    minprob: float = 1e-10
    dropmode: str = ""
    gcperio: bool = True
    context: str = "fork"
    consts: Dict[str, float] = field(default_factory=dict)
    altconsts: Dict[str, float] = field(default_factory=dict)
    catconsts: Dict[str, float] = field(default_factory=dict)


class System(Logged, Probalistic, Collected):
    def __init__(
        self, filename: str, logfile: Optional[str] = None, loglevel: str = "INFO"
    ):
        Logged.setlogger(logfile, loglevel)
        self.param = SysParam.readfile(filename)
        self.log.info(f"Parameters loaded: {self.param}")
        self.runparam = RunParam.readfile(filename)
        self.log.info(f"Run Parameters loaded: {self.runparam}")
        if self.runparam.gcperio:
            gc.disable()
        self.time = 0.0
        # If seed set, always restart from the same seed. For timing/debugging purpose
        self.step = 0
        Probalistic.setprobalist(vol=self.param.vol, minprob=self.runparam.minprob)
        CollectofCompound()
        CollectofReaction(
            self.runparam.dropreac, categorize=False, dropmode=self.runparam.dropmode
        )
        self.reac_collect.init_ruleset(
            self.runparam.consts, self.runparam.altconsts, self.runparam.catconsts
        )
        self.initialized = False
        self.log.info("System created")

    @property
    def memuse(self) -> float:
        return float(Process(getpid()).memory_info().rss) / 1024 / 1024

    def conc_of(self, compound: str) -> float:
        try:
            return self.comp_collect.active[compound].pop / self.param.vol
        except KeyError:  # compounds doesn't exist => conc=0
            return 0.0

    @property
    def poplist(self) -> Dict[str, int]:
        return {
            comp.description.name: comp.pop
            for comp in self.comp_collect.active.values()
        }

    @property
    def lendist(self) -> Dict[int, int]:
        return self.comp_collect.dist(lenweight=True, full=False)

    @property
    def pooldist(self) -> Dict[int, int]:
        return self.comp_collect.dist(lenweight=False, full=True)

    def statlist(self) -> Tuple[Dict[str, int], Dict[str, Dict[int, int]]]:
        stat = {}
        dist = {}
        dist["lendist"] = self.lendist
        dist["pooldist"] = self.pooldist
        stat["nbcomp"] = len(self.comp_collect.active)
        stat["nbreac"] = len(self.reac_collect.active)
        stat["poolsize"] = len(self.comp_collect.pool)
        stat["poolreac"] = len(self.reac_collect.pool)
        stat["maxlength"] = max(dist["lendist"])
        return stat, dist

    def initialize(self) -> None:
        if not self.initialized:
            for compound, pop in self.param.init.items():
                self.comp_collect[compound].init_pop(pop)
            Compound.trigger_update()
            Reaction.trigger_update()

    def _process(self, tstop: float) -> bool:
        # Check if a cleanup should be done
        if self.runparam.autoclean:
            self.probalist.clean()
        # Check if end of time is nigh
        if self.time >= self.param.tend:
            raise TimesUp(f"t={self.time}")
        if self.log.runtime() >= self.param.rtlim:
            raise RuntimeLim(f"t={self.time}")
        # Then process self.maxsteps times
        for _ in repeat(None, self.param.maxsteps):
            # ... but stop is step ends
            if self.time >= tstop:
                return True
            # choose a random event
            chosen, dt = self.probalist.choose()
            if chosen.description.name not in self.reac_collect.active:
                self.log.error(
                    f"Weird! {chosen} is not active... Is it in the pool? {chosen.description.name in self.reac_collect.pool}"
                )
            # check if there even was an event to choose
            if chosen is None:
                raise NotFound(f"t={self.time}")
            # perform the (chosen one) event
            chosen.process()
            if self.probalist.probtot == 0:
                raise NoMore(f"t={self.time}")
            # update time for next step
            self.time += dt
            self.step += 1
        self.log.warning(f"maxsteps per process (={self.param.maxsteps}) too low")
        return False

    def run(self, num: int = 0) -> Tuple[DataFrame, int, int, str]:
        lines = (
            ["thread", "ptime", "memuse", "step", "time"]
            + self.param.save
            + ["maxlength", "nbcomp", "poolsize", "nbreac", "poolreac"]
        )
        self.log.connect(f"Reconnected from thread {num+1}", num + 1)
        self.log.info("Will initialize")
        self.probalist.seed(self.param.seed)
        self.initialize()
        table = DataFrame(index=lines)
        lendist = DataFrame()
        pooldist = DataFrame()
        self.time = 0.0
        tnext = 0.0
        step = 0
        # Process(getpid()).cpu_affinity([num % cpu_count()])
        while True:
            try:
                self.log.info(f"#{step}: {self.time} -> {tnext}")
                finished = self._process(tnext)
                stat, dist = self.statlist()
                res = (
                    [num, self.log.runtime(), self.memuse, self.step, self.time]
                    + [self.conc_of(comp) for comp in self.param.save]
                    + [
                        stat["maxlength"],
                        stat["nbcomp"],
                        stat["poolsize"],
                        stat["nbreac"],
                        stat["poolreac"],
                    ]
                )
                table[step] = res
                self.log.debug(str(res))
                lendist = lendist.join(
                    DataFrame.from_dict({step: dist["lendist"]}), how="outer"
                ).fillna(0)
                pooldist = pooldist.join(
                    DataFrame.from_dict({step: dist["pooldist"]}), how="outer"
                ).fillna(0)
                if self.runparam.gcperio:
                    gc.collect()
                if finished:
                    tnext += self.param.tstep
                step += 1
            except Finished as the_end:
                end = f"{the_end} ({self.log.runtime()} s)"
                if isinstance(the_end, HappyEnding):
                    self.log.info(str(the_end))
                elif isinstance(the_end, BadEnding):
                    self.log.error(str(the_end))
                else:
                    self.log.warning(str(the_end))
                break
            except KeyboardInterrupt:
                end = f"Interrupted! ({self.time} s)"
                self.log.error(end)
                break
        self.log.disconnect(f"Disconnected from thread {num}")
        return (table, lendist.astype(int), pooldist.astype(int), end)

    def multirun(self, nbthread: Optional[int] = None) -> Result:
        ctx = get_context(self.context)
        if nbthread is None:
            nbthread = ctx.cpu_count()
        self.log.disconnect(reason="Entering Multithreading envoronment")
        with ctx.Pool(nbthread) as pool:
            res = pool.map(self.run, range(nbthread))
        self.log.connect(reason="Leaving Multithreading environment")
        self.log.info("Multithread run finished")
        return Result(res)

    def set_param(self, **kw) -> None:
        self.param = replace(self.param, **kw)
        self.probalist.vol = self.param.vol  # Need to update volume if changed!
        self.log.info(f"Parameters changed: {self.param}")

    def addkept(self, reac: str) -> None:
        """reactions that are kept"""
        self.reac_collect[reac].addkept()
