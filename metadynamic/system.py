from multiprocessing import get_context, cpu_count
from itertools import repeat, chain
from os import getpid
from typing import List, Optional, Dict, Any, Tuple, Type
from psutil import Process
from dataclasses import dataclass, replace, field

from numpy import log, random
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
from metadynamic.logger import Log
from metadynamic.proba import Probalist
from metadynamic.collector import Collect
from metadynamic.processing import Result
from metadynamic.chemical import Compound, Reaction
from metadynamic.description import ReacDescr, CompDescr

# D = TypeVar("D", "Descr", "ReacDescr", "CompDescr")
# C = TypeVar("C", "Chemical", "Reaction", "Compound")
# T = TypeVar("T")
# A = TypeVar("A")
# R = TypeVar("R")


class Ruleset:
    def __init__(
        self,
        consts: Dict[str, float],
        altconsts: Dict[str, float],
        catconsts: Dict[str, float],
    ):
        # Test if bad keys were entered
        self.testkinds(consts)
        self.testkinds(altconsts)
        self.testkinds(catconsts)
        # set class variables
        # dictionary copy may be useless... but only called once at startup thus is a no cost safety
        self.consts = consts.copy()
        kinds = set(self.consts.keys())
        self.descrs: Dict[str, Type[ReacDescr]] = {
            cls.kind: cls for cls in ReacDescr.subclasses() if cls.kind in kinds
        }
        self.altconsts = altconsts.copy()
        self.catconsts = catconsts.copy()
        # Missing keys in alt/cat are set to 1
        for kind in kinds - altconsts.keys():
            self.altconsts[kind] = 1
        for kind in kinds - catconsts.keys():
            self.catconsts[kind] = 1

    def reac_from_name(self, description: str) -> ReacDescr:
        try:
            kind, reactants, catal, pos = description.split(".")
        except ValueError:
            raise ValueError(f"Badly formatted description {description}")
        return self._newreac(
            description,
            kind,
            reactants.split("+"),
            catal,
            -1 if pos == "" else int(pos),
        )

    def reac_from_descr(
        self, kind: str, reactants: List[str], catal: str, pos: int
    ) -> ReacDescr:
        return self._newreac("", kind, reactants, catal, pos)

    def full_reac(self, reactant: "Compound") -> List[ReacDescr]:
        return list(
            chain.from_iterable(
                [
                    descriptor.fullfrom(
                        reactant,
                        self.consts[kind],
                        self.altconsts[kind],
                        self.catconsts[kind],
                    )
                    for kind, descriptor in self.descrs.items()
                ]
            )
        )

    def _newreac(
        self, name: str, kind: str, reactants: List[str], catal: str, pos: int
    ) -> ReacDescr:
        return self.descrs[kind](
            name,
            [CompDescr(x) for x in reactants],
            CompDescr(catal),
            pos,
            const=self.consts[kind],
            altconst=self.altconsts[kind],
            catconst=self.catconsts[kind],
        )

    @classmethod
    def testkinds(cls, dic: Dict[str, Any]) -> None:
        badkeys = list(dic.keys() - ReacDescr.kinds())
        if badkeys:
            raise ValueError(f"'{badkeys}' are not recognized reaction types")


class CollectofCompound(Collect[Compound]):
    def _create(self, name: str) -> Compound:
        return Compound(CompDescr(name), self.system)

    def _categorize(self, obj: "Compound") -> List[str]:
        return obj.description.categories

    def dist(self, lenweight: bool = False, full: bool = False) -> Dict[int, int]:
        res: Dict[int, int] = {}
        search = self.pool if False else self.active
        for comp in search.values():
            length = comp.length
            inc = comp.pop if lenweight else 1
            if length not in res:
                res[length] = 0
            res[length] += inc
        return res


class CollectofReaction(Collect[Reaction]):
    def _create(self, name: str) -> Reaction:
        assert isinstance(name, str), f"{name} is not a string..."
        return Reaction(self.ruleset.reac_from_name(name), self.system)

    def _categorize(self, obj: Reaction) -> List[str]:
        return obj.description.categories

    def init_ruleset(
        self,
        consts: Dict[str, float],
        altconsts: Dict[str, float],
        catconsts: Dict[str, float],
    ) -> None:
        self.ruleset = Ruleset(consts, altconsts, catconsts)

    def from_descr(self, description: ReacDescr) -> Reaction:
        try:
            return self.pool[description.name]
        except KeyError:
            newobj = Reaction(description, self.system)
            self.pool[description.name] = newobj
            return newobj

    def describe(
        self, kind: str, reactants: List[str], catal: str, pos: int
    ) -> Reaction:
        return self.from_descr(
            self.ruleset.reac_from_descr(kind, reactants, catal, pos)
        )

    def get_related(self, reactant: "Compound") -> List[Reaction]:
        return [
            self.from_descr(description)
            for description in self.ruleset.full_reac(reactant)
        ]


@dataclass
class SysParam:
    conc: float = 0.1
    tend: float = 1.0
    tstep: float = 0.01
    rtlim: float = 900.0  # Limit to 15min runtime
    maxsteps: int = 10000
    ptot: int = 100
    vol: float = field(init=False)
    seed: Optional[int] = None
    save: List[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.vol = self.ptot / self.conc


class System:
    def __init__(
        self,
        init: Dict[str, int],
        consts: Dict[str, float],
        altconsts: Optional[Dict[str, float]] = None,
        catconsts: Optional[Dict[str, float]] = None,
        logfile: Optional[str] = None,
        loglevel: str = "INFO",
        dropreac: bool = True,
        autoclean: bool = True,
        minprob: float = 1e-10,
    ):
        self.log = Log(logfile, loglevel)
        if altconsts is None:
            altconsts = {}
        if catconsts is None:
            catconsts = {}
        self.time = 0.0
        self.dropreac = dropreac
        self.autoclean = autoclean
        # If seed set, always restart from the same seed. For timing/debugging purpose
        self.step = 0
        self.param = SysParam(ptot=sum([pop * len(comp) for comp, pop in init.items()]))
        self.comp_collect = CollectofCompound(self)
        self.reac_collect = CollectofReaction(self, dropreac)  # set of active reactions
        self.reac_collect.init_ruleset(consts, altconsts, catconsts)
        self.probalist = Probalist(minprob=minprob)
        for compound, pop in init.items():
            self.comp_collect[compound].init_pop(pop)
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

    def _process(self, tstop: float) -> bool:
        # Check if a cleanup should be done
        if self.autoclean:
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
            chosen = self.probalist.choose()
            # check if there even was an event to choose
            if chosen is None:
                raise NotFound(f"t={self.time}")
            # perform the (chosen one) event
            chosen.process()
            if self.probalist.probtot == 0:
                raise NoMore(f"t={self.time}")
            # update time for next step
            self.time += log(1 / random.rand()) / self.probalist.probtot
            self.step += 1
        self.log.warning(f"maxsteps per process (={self.param.maxsteps}) too low")
        return False

    def run(self, num: int = 0) -> Tuple[DataFrame, int, int, str]:
        lines = (
            ["thread", "ptime", "memuse", "step", "time"]
            + self.param.save
            + ["maxlength", "nbcomp", "poolsize", "nbreac", "poolreac"]
        )
        table = DataFrame(index=lines)
        lendist = DataFrame()
        pooldist = DataFrame()
        random.seed(self.param.seed)
        self.log.connect(f"Reconnected from thread {num+1}", num + 1)
        self.time = 0.0
        tnext = 0.0
        step = 0
        #Process(getpid()).cpu_affinity([num % cpu_count()])
        while True:
            try:
                self.log.info(f"#{step}: {self.time} -> {tnext}")
                finished = self._process(tnext)
                stat, dist = self.statlist()
                table[step] = (
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
                lendist = lendist.join(
                    DataFrame.from_dict({step: dist["lendist"]}), how="outer"
                ).fillna(0)
                pooldist = pooldist.join(
                    DataFrame.from_dict({step: dist["pooldist"]}), how="outer"
                ).fillna(0)
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
        self.log.info("Start multithread run")
        self.log.disconnect(reason="Entering Multithreading environment")
        ctx = get_context("fork")
        if nbthread is None:
            nbthread = ctx.cpu_count()
        with ctx.Pool(nbthread) as pool:
            res = pool.map(self.run, range(nbthread))
        self.log.connect(reason="Leaving Multithreading environment")
        self.log.info("Multithread run finished")
        return Result(res)

    def set_param(self, **kw):
        self.param = replace(self.param, **kw)

    def addkept(self, reac: str) -> None:
        """reactions that are kept"""
        self.reac_collect[reac].addkept()
