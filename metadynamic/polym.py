from __future__ import annotations
from multiprocessing import get_context, cpu_count, current_process
from itertools import repeat, chain
from os import getpid
from typing import *
from psutil import Process
from collections import deque
from dataclasses import dataclass, replace, field

import numpy as N
import pandas as P

from metadynamic.utils import *
from metadynamic.proba import Probaobj, Probalist
from metadynamic.collector import Collect

D = TypeVar("D", "Descr", "ReacDescr", "CompDescr")
C = TypeVar("C", "Chemical", "Reaction", "Compound")
T = TypeVar("T")
A = TypeVar("A")
R = TypeVar("R")


class Result:
    def __init__(self, data):
        res = N.array(data)
        self.data = {
            "table": res[:, 0],
            "lendist": res[:, 1],
            "pooldist": res[:, 2],
            "end": res[:, 3],
        }

    def _format(self, name, field, num, mean=None):
        if num is None or num == "m":
            res = self.data[name].mean()
        elif num == "s":
            res = self.data[name].std()
        elif num == "+":
            res = self.data[name].mean() + self.data[name].std()
        elif num == "-":
            res = self.data[name].mean() - self.data[name].std()
        else:
            res = self.data[name][num]
        if field is not None:
            res = res.loc[field]
        if mean:
            res = self.running_mean(res, mean)
        return res

    def table(self, field=None, num=None, mean=None):
        return self._format("table", field, num, mean)

    def lendist(self, field=None, num=None, mean=None):
        return self._format("lendist", field, num, mean)

    def pooldist(self, field=None, num=None, mean=None):
        return self._format("pooldist", field, num, mean)

    def end(self, num=None):
        if num is None:
            return self.data["end"]
        return self.data["end"][num]

    @staticmethod
    def running_mean(data, length):
        """Running mean from https://stackoverflow.com/questions/13728392/moving-average-or-running-mean"""
        return N.convolve(data, N.ones((length,)) / length, mode="valid")


class Descr:
    _descrtype = "Chemical Description"

    def __init__(self, name: str):
        self.name: str = name

    @property
    def categories(self) -> List[str]:
        return []

    def __repr__(self) -> str:
        return f"{self._descrtype}: self.name"

    def __str__(self) -> str:
        return self.name


class CompDescr(Descr):
    _descrtype = "Compound Description"

    # Check memoize effects
    @property
    # @memoize_property
    def length(self) -> int:
        return len(self.name)

    @property
    # @memoize_property
    def ismono(self) -> bool:
        return self.length == 1

    @property
    # @memoize_property
    def ispolym(self) -> bool:
        return self.name.isalpha()

    @property
    # @memoize_property
    def isact(self) -> bool:
        return self.name[-1] == "*" and self.name[:-1].isalpha()

    @property
    # @memoize_property
    def categories(self) -> List[str]:
        res = []
        if self.ispolym:
            res.append("polym")
            if self.ismono:
                res.append("mono")
        elif self.isact:
            res.append("actpol")
        return res

    @property
    # @memoize_property
    def activate(self) -> str:
        assert not self.isact, f"{self} is already activated"
        return self.name + "*"

    @property
    # @memoize_property
    def unactivate(self) -> str:
        assert self.isact, f"{self} is not activated"
        return self.name[:-1]

    # @memoize_oneparam
    def cut(self, pos: int) -> List[str]:
        name = self.name
        assert (
            1 <= pos < len(name)
        ), f"cut position should be at least 1, at most the chain length minus one"
        return [name[:pos], name[pos:]]

    # @memoize_oneparam
    def epimer(self, pos: int) -> str:
        name = self.name
        assert (
            0 <= pos < len(name)
        ), f"epimer position should be at least 0, at most the chain length minus one"
        return name[:pos] + name[pos].swapcase() + name[pos + 1:]

    # @memoize_oneparam
    def extract(self, pos: int) -> str:
        return self.name[pos]


class ReacDescr(Descr):
    _descrtype = "Reaction Description"
    kind = "X"
    nbreac = 0

    def __init__(
        self,
        name: str,
        reactants: List[CompDescr],
        catal: Optional[CompDescr],
        pos: Optional[int],
        const: float,
        altconst: float,
        catconst: float,
    ):
        assert len(reactants) == self.nbreac
        if not name:
            reactantnames = (
                ""
                if self.nbreac == 0
                else reactants[0].name
                if self.nbreac == 1
                else reactants[0].name + "+" + reactants[1].name
            )
            name = f"{self.kind}.{reactantnames}.{catal if catal else ''}.{pos if pos else ''}"
        super().__init__(name)
        assert self.name[0] == self.kind, f"{self.name} created as {self.kind}"
        self.reactants: List[CompDescr] = reactants
        self.catal: Optional[CompDescr] = catal
        self._const: float = const
        self._altconst: float = altconst
        self._catconst: float = catconst
        self.pos: Optional[int] = pos
        self.checkdimer()

    @property
    def categories(self) -> List[str]:
        return [self.kind]

    def checkdimer(self) -> None:
        self.dimer = False

    def build_products(self) -> List[str]:
        return []

    def build_const(self) -> float:
        const = (self._altconst ** self.nbalt()) * self._const
        if self.catal:
            const *= self._catconst
        return const

    def nbalt(self) -> int:
        return 0

    @classmethod
    def fullfrom(
        cls, reactant: Compound, const: float, altconst: float, catconst: float
    ) -> List[ReacDescr]:
        return []

    @classmethod
    def subclasses(cls):
        for subclass in cls.__subclasses__():
            yield from subclass.subclasses()
            yield subclass

    @classmethod
    def kinds(cls):
        return (subcls.kind for subcls in cls.subclasses())


class PolymDescr(ReacDescr):
    kind = "P"
    nbreac = 2

    def checkdimer(self) -> None:
        self.dimer = self.reactants[0] == self.reactants[1]

    def build_products(self) -> List[str]:
        return [self.reactants[0].name + self.reactants[1].name]

    def nbalt(self) -> int:
        return 1 if self.reactants[0].ismono else 0

    @classmethod
    def fullfrom(
        cls, reactant: Compound, const: float, altconst: float, catconst: float
    ) -> List[ReacDescr]:
        first: Callable[[Compound], ReacDescr] = lambda other: cls(
            "",
            [reactant.description, other.description],
            None,
            None,
            const,
            altconst,
            catconst,
        )
        last: Callable[[Compound], ReacDescr] = lambda other: cls(
            "",
            [other.description, reactant.description],
            None,
            None,
            const,
            altconst,
            catconst,
        )
        return [
            create(other)
            for other in reactant.system.comp_collect.cat_list("polym")
            for create in (first, last)
        ]


class ActpolymDescr(ReacDescr):
    kind = "A"
    nbreac = 2

    def build_products(self) -> List[str]:
        return [self.reactants[0].unactivate + self.reactants[1].name]

    def nbalt(self) -> int:
        name0 = self.reactants[0].extract(-2)
        name1 = self.reactants[1].extract(0)
        return 0 if samecase(name0, name1) else 1

    @classmethod
    def fullfrom(
        cls, reactant: Compound, const: float, altconst: float, catconst: float
    ) -> List[ReacDescr]:
        if reactant.description.ispolym:
            return [
                cls(
                    "",
                    [other.description, reactant.description],
                    None,
                    None,
                    const,
                    altconst,
                    catconst,
                )
                for other in reactant.system.comp_collect.cat_list("actpol")
            ]
        if reactant.description.isact:
            return [
                cls(
                    "",
                    [reactant.description, other.description],
                    None,
                    None,
                    const,
                    altconst,
                    catconst,
                )
                for other in reactant.system.comp_collect.cat_list("polym")
            ]
        return []


class ActivationDescr(ReacDescr):
    kind = "a"
    nbreac = 1

    def build_products(self) -> List[str]:
        return [self.reactants[0].activate]

    def nbalt(self) -> int:
        return 0 if self.reactants[0].ismono else 1

    @classmethod
    def fullfrom(
        cls, reactant: Compound, const: float, altconst: float, catconst: float
    ) -> List[ReacDescr]:
        if reactant.description.ispolym:
            return [
                cls("", [reactant.description], None, None, const, altconst, catconst)
            ]
        return []


class DectivationDescr(ReacDescr):
    kind = "d"
    nbreac = 1

    def build_products(self) -> List[str]:
        return [self.reactants[0].unactivate]

    def nbalt(self) -> int:
        return 0 if self.reactants[0].length == 2 else 1

    @classmethod
    def fullfrom(
        cls, reactant: Compound, const: float, altconst: float, catconst: float
    ) -> List[ReacDescr]:
        if reactant.description.isact:
            return [
                cls("", [reactant.description], None, None, const, altconst, catconst)
            ]
        return []


class HydroDescr(ReacDescr):
    kind = "H"
    nbreac = 1

    def build_products(self) -> List[str]:
        assert self.pos is not None
        return self.reactants[0].cut(self.pos)

    def nbalt(self) -> int:
        assert self.pos is not None
        left = self.reactants[0].extract(self.pos - 1)
        right = self.reactants[0].extract(self.pos)
        return 0 if samecase(right, left) else 1

    @classmethod
    def fullfrom(
        cls, reactant: Compound, const: float, altconst: float, catconst: float
    ) -> List[ReacDescr]:
        return [
            cls("", [reactant.description], None, i, const, altconst, catconst)
            for i in range(1, reactant.length)
        ]


class RacemDescr(ReacDescr):
    kind = "R"
    nbreac = 1

    def build_products(self) -> List[str]:
        assert self.pos is not None
        return [self.reactants[0].epimer(self.pos)]

    def nbalt(self) -> int:
        assert self.pos is not None, f"{self.name} as no position..."
        name = self.reactants[0]
        pos = self.pos
        res = 0
        if pos < (name.length - 1) and samecase(
            name.extract(pos), name.extract(pos + 1)
        ):
            res += 1
        if (pos > 0) and samecase(name.extract(pos), name.extract(pos - 1)):
            res += 1
        return res

    @classmethod
    def fullfrom(
        cls, reactant: Compound, const: float, altconst: float, catconst: float
    ) -> List[ReacDescr]:
        return [
            cls("", [reactant.description], None, i, const, altconst, catconst)
            for i in range(reactant.length)
        ]


class EpimDescr(RacemDescr):
    kind = "E"

    @classmethod
    def fullfrom(
        cls, reactant: Compound, const: float, altconst: float, catconst: float
    ) -> List[ReacDescr]:
        if reactant.description.ismono:
            return []
        return [cls("", [reactant.description], None, 0, const, altconst, catconst)]


class Chemical(Generic[D]):
    _descrtype = "Chemical"

    def __init__(self, description: D, system: System):
        self.system = system
        self.description: D = description
        # self.system.log.debug(f"Creating {description} Chemical")

    def __repr__(self) -> str:
        return f"{self._descrtype}: {self}"

    def __str__(self) -> str:
        return self.description.name


class Reaction(Chemical[ReacDescr]):
    _descrtype = "Reaction"

    def __init__(self, description: ReacDescr, system: System):
        super().__init__(description, system)
        self._vol = system.param.vol
        self._probaobj = Probaobj(self, system.probalist)
        self._reactants: List[Compound] = [
            system.comp_collect[i.name] for i in description.reactants
        ]
        self._productnames = description.build_products()
        self.const = description.build_const()
        self._started = False
        catal = description.catal
        if catal:
            self._catalized = True
            self._catal = system.comp_collect[catal.name]
        else:
            self._catalized = False
        self._reaccalc: Optional[Callable[[], float]] = None
        self._uncatcalc: Optional[Callable[[], float]] = None
        self._set_reaccalc()
        self.register()
        # self.system.log.debug(f"{self} created from {description}, transforms {self._reactants} into {self._productnames}")

    def destroy(self) -> None:
        if self._probaobj.registered:
            # self.system.log.debug(f"Destroying {self}")
            self._probaobj.unregister()
            for comp in self._reactants:
                comp.unregister_reaction(self)
            if self._catalized:
                self._catal.unregister_reaction(self)
            self.system.reac_collect.remove(self.description.name)

    def register(self) -> None:
        # self.system.log.debug(f"Registering {self}")
        for comp in self._reactants:
            comp.register_reaction(self)
        if self._catalized:
            self._catal.register_reaction(self)

    def process(self) -> None:
        for reac in self._reactants:
            try:
                reac.dec()
            except Finished as f:
                # Add details to the message
                detail = f" from {self} (p={self.proba}, "
                for comp in self._reactants:
                    detail += f"[{comp.description}]={comp.pop} ,"
                if self._catalized:
                    detail += f"catal[{self._catal.description}]={self._catal.pop})"
                else:
                    detail += ")"
                raise DecrZero(f.detail + detail)
        if not self._started:
            self._products = [
                self.system.comp_collect[name] for name in self._productnames
            ]
            self._started = True
        for prod in self._products:
            prod.inc()

    def _order0(self) -> float:
        # Reaction ->
        return self.const

    def _order1(self) -> float:
        # Reaction A ->
        return self.const * self._reactants[0].pop

    def _order2dim(self) -> float:
        pop0 = self._reactants[0].pop
        # Reaction A + A ->
        return self.const * pop0 * (pop0 - 1)

    def _order2(self) -> float:
        # Reaction A + B ->
        return self.const * self._reactants[0].pop * self._reactants[1].pop

    def _autocatdim(self) -> float:
        # Reaction A + A -[A]->
        assert self._uncatcalc is not None
        return self._uncatcalc() * (self._catal.pop - 2)

    def _autocat(self) -> float:
        # Reaction A (+other...)  -[A}->
        assert self._uncatcalc is not None
        return self._uncatcalc() * (self._catal.pop - 1)

    def _cat(self) -> float:
        # Reaction (other...) -[B]->
        assert self._uncatcalc is not None
        return self._uncatcalc() * self._catal.pop

    def _set_reaccalc(self) -> None:
        order = self.description.nbreac
        dimer = self.description.dimer
        if order == 0:
            self._reaccalc = self._order0
        elif order == 1:
            self._reaccalc = self._order1
        elif order == 2:
            if dimer:
                self._reaccalc = self._order2dim
                self.const = self.const / self._vol / 2.0
            else:
                self._reaccalc = self._order2
                self.const = self.const / self._vol
        else:
            raise ValueError(
                "Kinetic orders different from 0,1 or 2 are not supported (yet)"
            )
        if self._catalized:
            self._uncatcalc = self._reaccalc
            self.const = self.const / self._vol
            if self._catal in self._reactants:
                if dimer:
                    self._reaccalc = self._autocatdim
                else:
                    self._reaccalc = self._autocat
            else:
                self._reaccalc = self._cat

    def calcproba(self) -> float:
        assert self._reaccalc is not None
        return 0.0 if self._reactants[0].pop == 0 else self._reaccalc()

    @property
    def proba(self) -> float:
        return self._probaobj.proba

    def update(self) -> None:
        self._probaobj.update(self.calcproba())

    def addkept(self) -> None:
        if self._catal:
            self._catal.addkept(self)
        for comp in self._reactants:
            comp.addkept(self)


class Compound(Chemical[CompDescr]):
    _descrtype = "Compound"

    def __init__(self, description: CompDescr, system: System):
        super().__init__(description, system)
        self._reactions: Set[Reaction] = set()
        self.system = system
        self._kept_descr: List[Reaction] = []
        self.pop = 0
        self.length = description.length

    def inc(self) -> None:
        self.pop += 1
        if self.pop == 1:
            self._activate()
        self._upd_reac()

    def dec(self) -> None:
        if self.pop > 0:
            self.pop -= 1
        else:
            raise DecrZero(self.description.name)
        if self.pop == 0:
            self._unactivate()
        else:
            self._upd_reac()

    def init_pop(self, start: int) -> None:
        if start != 0:
            if self.pop == 0:
                self.pop = start
                self._activate()
            else:
                self.pop = start
            self._upd_reac()
        else:
            self.pop = 0
            self._unactivate()

    def _activate(self) -> None:
        self.system.comp_collect.activate(self.description.name)
        self.reac_set()
        # for reac in self.reac_set():
        #    self.system.reac_collect[reac.name]

    def _unactivate(self) -> None:
        for reac in self.reac_set():
            reac.destroy()
        self.system.comp_collect.unactivate(self.description.name)

    def _upd_reac(self) -> None:
        # copy as list because self._reactions to be changed in loop. Elegant?...
        for reac in list(self._reactions):
            reac.update()

    def reac_set(self) -> List[Reaction]:
        return self._kept_descr + self.system.reac_collect.get_related(self)

    def addkept(self, reac: Reaction) -> None:
        self._kept_descr.append(reac)

    def register_reaction(self, reaction: Reaction) -> None:
        self._reactions.add(reaction)

    def unregister_reaction(self, reaction: Reaction) -> None:
        try:
            self._reactions.remove(reaction)
        except KeyError:
            pass


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
        self.descrs = {
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

    def full_reac(self, reactant: Compound) -> List[ReacDescr]:
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

    def _categorize(self, obj: Compound):
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

    def _categorize(self, obj: Reaction):
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

    def get_related(self, reactant: Compound) -> List[Reaction]:
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

    def __post_init__(self):
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
        self.timer = Timer()
        self.log = Log(self, logfile, loglevel)
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
        return Process(getpid()).memory_info().rss / 1024 / 1024

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
        if self.timer.time >= self.param.rtlim:
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
            self.time += N.log(1 / N.random.rand()) / self.probalist.probtot
            self.step += 1
        self.log.warning(f"maxsteps per process (={self.param.maxsteps}) too low")
        return False

    def run(self, num: int = 0) -> Tuple[P.DataFrame, int, int, str]:
        lines = (
            ["thread", "ptime", "memuse", "step", "time"]
            + self.param.save
            + ["maxlength", "nbcomp", "poolsize", "nbreac", "poolreac"]
        )
        table = P.DataFrame(index=lines)
        lendist = P.DataFrame()
        pooldist = P.DataFrame()
        N.random.seed(
            self.param.seed
        )  # necessary for multiprocessing from different seeds
        self.log.connect(f"Reconnected from thread {num+1}", num + 1)
        self.time = 0.0
        tnext = 0.0
        step = 0
        self.timer.reset()
        Process(getpid()).cpu_affinity([num % cpu_count()])
        while True:
            try:
                self.log.info(f"#{step}: {self.time} -> {tnext}")
                finished = self._process(tnext)
                stat, dist = self.statlist()
                table[step] = (
                    [num, self.timer.time, self.memuse, self.step, self.time]
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
                    P.DataFrame.from_dict({step: dist["lendist"]}), how="outer"
                ).fillna(0)
                pooldist = pooldist.join(
                    P.DataFrame.from_dict({step: dist["pooldist"]}), how="outer"
                ).fillna(0)
                if finished:
                    tnext += self.param.tstep
                step += 1
            except Finished as the_end:
                end = f"{the_end} ({self.timer.time} s)"
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
