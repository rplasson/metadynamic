from __future__ import annotations
from multiprocessing import get_context, cpu_count
from itertools import repeat, chain
from time import process_time
from os import getpid
from typing import List, Set, Dict, Tuple
from psutil import Process

import numpy as N
import pandas as P


def samecase(one, two):
    return (one.islower() and two.islower()) or (one.isupper() and two.isupper())


def get_all_subclasses(cls):
    all_subclasses = []
    for subclass in cls.__subclasses__():
        all_subclasses.append(subclass)
        all_subclasses.extend(get_all_subclasses(subclass))
    return all_subclasses


class Finished(Exception):
    num = -1
    error_message = "Stopped for unknown reason"

    def __init__(self, detail=""):
        self.detail = detail

    @property
    def message(self) -> str:
        msg = self.error_message
        if self.detail != "":
            msg = msg + " -> " + self.detail
        return msg

    def __str__(self):
        return f"End ({self.num}): {self.message}"


class TimesUp(Finished):
    num = 0
    error_message = "Time is up"


class NoMore(Finished):
    num = 1
    error_message = "No more reactions can be processed"


class NotFound(Finished):
    num = 2
    error_message = "No reaction could be find"


class RoundError(Finished):
    num = 3
    error_message = "Rounding problem/negative probability detected"


class DecrZero(Finished):
    num = 4
    error_message = "Tried to decrement unpopulated species"


class RuntimeLim(Finished):
    num = 5
    error_message = "Runtime limit exceeded"


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
        return self._format("lendist", field, num)

    def pooldist(self, field=None, num=None, mean=None):
        return self._format("pooldist", field, num, mean)

    def end(self, num=None, mean=None):
        if num is None:
            return self.data["end"]
        return self.data["end"][num]

    @staticmethod
    def running_mean(data, length):
        """Running mean from https://stackoverflow.com/questions/13728392/moving-average-or-running-mean"""
        return N.convolve(data, N.ones((length,)) / length, mode="valid")


class Probalist:
    def __init__(self, minprob: float = 1e-10):
        self._minprob = minprob
        self._map = [N.array([])]
        self._mapobj = [N.array([])]
        self._nblist = 1
        self._actlist = 0
        self._problist = N.array([0.0])
        self.probtot = 0.0
        self._queue: List[Tuple[int, int]] = []

    def _updateprob(self, nlist: int, delta: float):
        newproba = self._problist[nlist] + delta
        # reset to 0 probabilities that are too low (rounding errors correction)
        if newproba < self._minprob:
            delta += newproba
            newproba = 0.0
        self._problist[nlist] = newproba
        self.probtot += delta

    def unregister(self, obj: Reaction):
        nlist = obj.nlist
        npos = obj.npos
        oldproba = self._map[nlist][npos]
        self._map[nlist][npos] = 0
        self._mapobj[nlist][npos] = None
        self._queue.append((nlist, npos))
        self._updateprob(nlist, -oldproba)
        obj.registered = False

    def update(self, obj: Reaction, proba: float):
        if obj.registered:
            delta = proba - self._map[obj.nlist][obj.npos]
            self._map[obj.nlist][obj.npos] = proba
            self._updateprob(obj.nlist, delta)
        else:
            if self._queue:
                obj.nlist, obj.npos = self._addfromqueue(obj, proba)
            else:
                obj.nlist, obj.npos = self._addfrommap(obj, proba)
            obj.registered = True
            self._updateprob(obj.nlist, proba)

    def getproba(self, obj: Reaction) -> float:
        return self._map[obj.nlist][obj.npos] if obj.registered else 0.0

    def _addfrommap(self, newobj: Reaction, proba: float) -> Tuple[int, int]:
        rpos = len(self._map[self._actlist])
        if rpos <= self._nblist:
            rlist = self._actlist
            self._map[rlist] = N.append(self._map[rlist], proba)
            self._mapobj[rlist] = N.append(self._mapobj[rlist], newobj)
            return rlist, rpos
        if len(self._map[0]) <= self._nblist:
            self._actlist = 0
        else:
            self._actlist += 1
            if self._actlist >= self._nblist:
                self._map.append(N.array([]))
                self._mapobj.append(N.array([]))
                self._problist = N.append(self._problist, 0.0)
                self._nblist += 1
        return self._addfrommap(newobj, proba)

    def _addfromqueue(self, newobj: Reaction, proba: float) -> Tuple[int, int]:
        rlist, rpos = self._queue.pop(0)
        self._map[rlist][rpos] = proba
        self._mapobj[rlist][rpos] = newobj
        return rlist, rpos

    def choose(self) -> Reaction:
        try:
            # First choose a random line in the probability map
            nlist = N.random.choice(self._nblist, p=self._problist / self.probtot)
            # Then choose a random column in the chosen probability line
            return N.random.choice(
                self._mapobj[nlist], p=self._map[nlist] / self._problist[nlist]
            )
        except ValueError:
            raise RoundError

    def clean(self):
        """(Re-)Compute probality for sums. 
           It is far slower than just updating individual probabilities,
           but the fast methods used in update function leads to accumulate small rounding
           errors.
           This functions is intended to be called regularly for cleaning
           these rounding errors."""
        self._problist = N.array([data.sum() for data in self._map])
        self.probtot = self._problist.sum()


class ReacDescr:
    kind = "X"
    nbreac = 0

    def __init__(
        self,
        name: str,
        reactants: List[str],
        catal: str,
        pos: int,
        const: float,
        altconst: float,
        catconst: float,
    ):
        self.reactants = reactants  # .copy()  Necessary? has parformance impact...
        self.name = name
        self.catal = catal
        self._const = const
        self._altconst = altconst
        self._catconst = catconst
        self.pos = pos
        self.checkdimer()

    def checkdimer(self):
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

    def __str__(self):
        return self.name

    @staticmethod
    def fullfrom(reactant: Chemical, system: System) -> List[ReacDescr]:
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

    def checkdimer(self):
        self.dimer = self.reactants[0] == self.reactants[1]

    def build_products(self) -> List[str]:
        return [self.reactants[0] + self.reactants[1]]

    def nbalt(self) -> int:
        return 1 if len(self.reactants[0]) == 1 else 0

    @staticmethod
    def fullfrom(reactant: Chemical, system: System) -> List[ReacDescr]:
        first = lambda other: system.ruleset.reac_from_descr(
            "P", [reactant.name, other], "", -1
        )
        last = lambda other: system.ruleset.reac_from_descr(
            "P", [other, reactant.name], "", -1
        )
        return [
            create(other) for other in system.comp_collect for create in (first, last)
        ]


class ActpolymDescr(PolymDescr):
    kind = "A"

    def nbalt(self) -> int:
        name0 = self.reactants[0][0]
        name1 = self.reactants[1][0]
        return 0 if samecase(name0, name1) else 1

    @staticmethod
    def fullfrom(reactant: Chemical, system: System) -> List[ReacDescr]:
        res = [
            system.ruleset.reac_from_descr("A", [other, reactant.name], "", -1)
            for other in system.mono_collect
        ]
        if reactant.length == 1:
            res += [
                system.ruleset.reac_from_descr("A", [reactant.name, other], "", -1)
                for other in system.comp_collect
            ]
        return res


class HydroDescr(ReacDescr):
    kind = "H"
    nbreac = 1

    def build_products(self) -> List[str]:
        name = self.reactants[0]
        if not (1 <= self.pos < len(name)):
            raise ValueError(
                "Hydrolysis cut position should be at least 1, at most the chain length minus one"
            )
        return [name[: self.pos], name[self.pos :]]

    def nbalt(self) -> int:
        left = self.reactants[0][self.pos - 1]
        right = self.reactants[0][self.pos]
        return 0 if samecase(right, left) else 1

    @staticmethod
    def fullfrom(reactant: Chemical, system: System) -> List[ReacDescr]:
        return [
            system.ruleset.reac_from_descr("H", [reactant.name], "", i)
            for i in range(1, reactant.length)
        ]


class RacemDescr(ReacDescr):
    kind = "R"
    nbreac = 1

    def build_products(self) -> List[str]:
        name = self.reactants[0]
        pos = self.pos
        if not (0 <= pos < len(name)):
            raise ValueError(
                "Racemization position should be at least 0, at most the chain length minus one"
            )
        return [name[:pos] + name[pos].swapcase() + name[pos + 1 :]]

    def nbalt(self) -> int:
        name = self.reactants[0]
        pos = self.pos
        res = 0
        if (pos < (len(name) - 1)) and samecase(name[pos], name[pos + 1]):
            res += 1
        if (pos > 0) and samecase(name[pos], name[pos - 1]):
            res += 1
        return res

    @staticmethod
    def fullfrom(reactant: Chemical, system: System) -> List[ReacDescr]:
        return [
            system.ruleset.reac_from_descr("R", [reactant.name], "", i)
            for i in range(reactant.length)
        ]


class EpimDescr(RacemDescr):
    kind = "E"

    @staticmethod
    def fullfrom(reactant: Chemical, system: System) -> List[ReacDescr]:
        if reactant.length == 1:
            return []
        return [system.ruleset.reac_from_descr("E", [reactant.name], "", 0)]


class Ruleset:
    def __init__(self, consts: dict, altconsts: dict, catconsts: dict):
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
        except (ValueError, AttributeError):
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
        name = "{kind}.{reactants}.{catal}.{pos}".format(
            kind=kind,
            reactants="+".join(reactants),
            catal=catal,
            pos=pos if pos >= 0 else "",
        )  # keep .format() to avoid too much code in f"..."
        return self._newreac(name, kind, reactants, catal, pos)

    def full_reac(self, reactant: Chemical, system: System) -> List[ReacDescr]:
        return list(
            chain.from_iterable(
                [
                    descriptor.fullfrom(reactant, system)
                    for descriptor in self.descrs.values()
                ]
            )
        )

    def _newreac(
        self, name: str, kind: str, reactants: list, catal: str, pos: int
    ) -> ReacDescr:
        return self.descrs[kind](
            name,
            reactants,
            catal,
            pos,
            const=self.consts[kind],
            altconst=self.altconsts[kind],
            catconst=self.catconsts[kind],
        )

    @classmethod
    def testkinds(cls, dic):
        badkeys = list(dic.keys() - ReacDescr.kinds())
        if badkeys:
            raise ValueError(f"'{badkeys}' are not recognized reaction types")


class Reaction:
    def __init__(self, description: ReacDescr, system: System):
        self._system = system
        self._description = description
        self._reactants = [self._system.get_compound(i) for i in description.reactants]
        self._products = description.build_products()
        self.const = description.build_const()
        self._started = False
        catal = description.catal
        if catal:
            self._catalized = True
            self._catal = system.get_compound(catal)
        else:
            self._catalized = False
        # self.nlist, self.npos = -1, -1 # Should be useless now (tested with self.registered)
        self.registered = False
        self._set_reaccalc()
        self.register()

    def destroy(self):
        self._system.probalist.unregister(self)
        # self.nlist, self.npos = -1, -1 # Should be useless now (tested with self.registered)
        for comp in self._reactants:
            comp.unregister_reaction(self)
        if self._catalized:
            self._catal.unregister_reaction(self)
        if self._system.dropreac:
            del self._system.reac_collect[str(self)]

    def register(self):
        for comp in self._reactants:
            comp.register_reaction(self)
        if self._catalized:
            self._catal.register_reaction(self)

    def process(self):
        for reac in self._reactants:
            try:
                reac.dec()
            except Finished as f:
                # Add details to the message
                detail = f" from {self} (p={self.proba}, "
                for comp in self._reactants:
                    detail += f"[{comp.name}]={comp.pop} ,"
                if self._catal:
                    detail += f"catal[{self._catal.name}]={self._catal.pop})"
                else:
                    detail += ")"
                raise DecrZero(f.detail + detail)
        if not self._started:
            self._products = [
                self._system.get_compound(name) for name in self._products
            ]
            self._started = True
        for prod in self._products:
            prod.inc()

    def __repr__(self):
        return self._description.name

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
        return self._uncatcalc() * (self._catal.pop - 2)

    def _autocat(self) -> float:
        # Reaction A (+other...)  -[A}->
        return self._uncatcalc() * (self._catal.pop - 1)

    def _cat(self) -> float:
        # Reaction (other...) -[B]->
        return self._uncatcalc() * self._catal.pop

    def _set_reaccalc(self):
        order = self._description.nbreac
        dimer = self._description.dimer
        if order == 0:
            self._reaccalc = self._order0
        elif order == 1:
            self._reaccalc = self._order1
        elif order == 2:
            if dimer:
                self._reaccalc = self._order2dim
                self.const = self.const / self._system.vol / 2.0
            else:
                self._reaccalc = self._order2
                self.const = self.const / self._system.vol
        else:
            raise ValueError(
                "Kinetic orders different from 0,1 or 2 are not supported (yet)"
            )
        if self._catalized:
            self._uncatcalc = self._reaccalc
            self.const = self.const / self._system.vol
            if self._catal in self._reactants:
                if dimer:
                    self._reaccalc = self._autocatdim
                else:
                    self._reaccalc = self._autocat
            else:
                self._reaccalc = self._cat

    def calcproba(self) -> float:
        return 0.0 if self._reactants[0].pop == 0 else self._reaccalc()

    @property
    def proba(self) -> float:
        return self._system.probalist.getproba(self)

    def update(self):
        self._system.probalist.update(self, self.calcproba())

    def addkept(self):
        if self._catal:
            self._catal.addkept(self._description)
        for comp in self._reactants:
            comp.addkept(self._description)


class Chemical:
    def __init__(self, name: str, system: System):
        if name == "":
            self._real = False
        else:
            self._real = True
            self._reactions: Set[Reaction] = set()
            self._system = system
            self._kept_descr: List[ReacDescr] = []
            self.pop = 0
            self.name = name
            self.length = len(name)

    def __repr__(self):
        return self.name

    def __bool__(self):
        return self._real

    def inc(self):
        self.pop += 1
        if self.pop == 1:
            self._generate()
        self._upd_reac()

    def dec(self):
        if self.pop > 0:
            self.pop -= 1
        else:
            raise DecrZero(self.name)
        if self.pop == 0:
            self._destroy()
        else:
            self._upd_reac()

    def init_pop(self, start: int):
        if start != 0:
            if self.pop == 0:
                self.pop = start
                self._generate()
            else:
                self.pop = start
            self._upd_reac()
        else:
            self.pop = 0
            self._destroy()

    def _generate(self):
        if self.name not in self._system.comp_collect:
            self._system.comp_collect[self.name] = self
            if self.length == 1:
                self._system.mono_collect[self.name] = self
        for reac in self._reac_descr():
            self._system.get_reaction(reac)

    def _upd_reac(self):
        # copy as list because self._reactions to be changed in loop
        for reac in list(self._reactions):
            reac.update()

    def _destroy(self):
        # copy as list because self._reactions to be changed in loop
        for reac in list(self._reactions):
            reac.destroy()
        del self._system.comp_collect[self.name]
        if self.length == 1:
            del self._system.mono_collect[self.name]

    def _reac_descr(self) -> List[ReacDescr]:
        return self._kept_descr + self._system.ruleset.full_reac(self, self._system)

    def addkept(self, description: ReacDescr):
        self._kept_descr.append(description)

    def register_reaction(self, reaction: Reaction):
        self._reactions.add(reaction)

    def unregister_reaction(self, reaction: Reaction):
        try:
            self._reactions.remove(reaction)
        except KeyError:
            pass


class System:
    def __init__(
        self,
        init: dict,
        consts: dict,
        altconsts: dict = None,
        catconsts: dict = None,
        conc: float = 0.1,
        tend: float = 1.0,
        tstep: float = 0.01,
        rtlim: float = 900.0,  # Limit to 15min runtime
        save: list = None,
        maxsteps: int = 10000,
        minprob: float = 1e-10,
        dropreac: bool = True,
        autoclean: bool = True,
        seed: int = None,
    ):
        if altconsts is None:
            altconsts = {}
        if catconsts is None:
            catconsts = {}
        self.ruleset = Ruleset(consts, altconsts, catconsts)
        self.runtime_init()
        # self.proc_init()
        self.time = 0.0
        self.conc = conc
        self.dropreac = dropreac
        self.autoclean = autoclean
        # If seed set, always restart from the same seed. For timing/debugging purpose
        self._seed = seed
        self.step = 0
        self._ptot = sum([pop * len(comp) for comp, pop in init.items()])
        self.vol = self._ptot / self.conc
        self.pool_collect: Dict[str, Chemical] = {}  # full set of created compound
        self.comp_collect: Dict[str, Chemical] = {}  # set of active compounds
        self.mono_collect: Dict[str, Chemical] = {}  # set of monomer compounds
        self.reac_collect: Dict[str, Reaction] = {}  # set of active reactions
        self.probalist = Probalist(minprob=minprob)
        for compound, pop in init.items():
            self.get_compound(compound).init_pop(pop)
        if save is None:
            save = []
        self.set_run(tend, tstep, save, maxsteps, rtlim)

    def runtime_init(self):
        self._ptime0 = process_time()

    @property
    def runtime(self) -> float:
        return process_time() - self._ptime0

    # def proc_init(self):
    # pass
    # self._proc = Process(getpid()) ###Broken, cannot pickle anymore !!

    @property
    def memuse(self) -> float:
        return Process(getpid()).memory_info().rss / 1024 / 1024
        ###Broken, cannot Process  pickle anymore !!
        # return self._proc.memory_info().rss / 1024 / 1024

    def conc_of(self, compound: str) -> float:
        try:
            return self.comp_collect[compound].pop / self.vol
        except KeyError:  # compounds doesn't exist => conc=0
            return 0.0

    @property
    def poplist(self) -> dict:
        return {comp: comp.pop for comp in self.comp_collect.values()}

    @property
    def lendist(self) -> dict:
        res: Dict[int, int] = {}
        for comp in self.comp_collect.values():
            length = comp.length
            pop = comp.pop
            if length not in res:
                res[length] = 0
            res[length] += pop
        return res

    @property
    def pooldist(self) -> dict:
        res: Dict[int, int] = {}
        for comp in self.pool_collect.values():
            length = comp.length
            if length not in res:
                res[length] = 0
            res[length] += 1
        return res

    def statlist(self) -> tuple:
        stat = {}
        dist = {}
        dist["lendist"] = self.lendist
        dist["pooldist"] = self.pooldist
        stat["nbcomp"] = len(self.comp_collect)
        stat["nbreac"] = len(self.reac_collect)
        stat["poolsize"] = len(self.pool_collect)
        stat["maxlength"] = max(dist["lendist"])
        return stat, dist

    def process(self, tstop: float):
        # Check if a cleanup should be done
        if self.autoclean:
            self.probalist.clean()
        # Check if end of time is nigh
        if self.time >= self.tend:
            raise TimesUp(f"t={self.time}")
        if self.runtime >= self.rtlim:
            raise RuntimeLim(f"t={self.time}")
        # Then process self.maxsteps times
        for _ in repeat(None, self.maxsteps):
            # ... but stop is step ends
            if self.time >= tstop:
                break
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

    def run(self, num: int = 0, output: bool = False) -> tuple:
        lines = (
            ["thread", "ptime", "memuse", "step", "time"]
            + self.save
            + ["maxlength", "nbcomp", "poolsize", "nbreac"]
        )
        table = P.DataFrame(index=lines)
        lendist = P.DataFrame()
        pooldist = P.DataFrame()
        # if self._seed:
        N.random.seed(self._seed)  # necessary for multiprocessing from different seeds
        self.time = 0.0
        tnext = 0.0
        step = 0
        self.runtime_init()
        # self.proc_init()
        Process(getpid()).cpu_affinity(
            [num % cpu_count()]
        )  ###Broken, cannot Process  pickle anymore !!
        # self._proc.cpu_affinity([num % cpu_count()])
        while True:
            try:
                if output:
                    print(f"#{step} -> {tnext}")
                self.process(tnext)
                stat, dist = self.statlist()
                table[step] = (
                    [num, self.runtime, self.memuse, self.step, self.time]
                    + [self.conc_of(comp) for comp in self.save]
                    + [
                        stat["maxlength"],
                        stat["nbcomp"],
                        stat["poolsize"],
                        stat["nbreac"],
                    ]
                )
                lendist = lendist.join(
                    P.DataFrame.from_dict({step: dist["lendist"]}), how="outer"
                ).fillna(0)
                pooldist = pooldist.join(
                    P.DataFrame.from_dict({step: dist["pooldist"]}), how="outer"
                ).fillna(0)
                tnext += self.tstep
                step += 1
            except Finished as the_end:
                end = f"{the_end} ({self.runtime} s)"
                break
            except KeyboardInterrupt:
                end = f"Interrupted! ({self.time} s)"
                break
        return (table, lendist.astype(int), pooldist.astype(int), end)

    def multirun(self, nbthread: int = None) -> Result:
        ctx = get_context("fork")
        if nbthread is None:
            nbthread = ctx.cpu_count()
        with ctx.Pool(nbthread) as pool:
            ### Warning! self.run cannot be pickeled anymore...
            res = pool.map(self.run, range(nbthread))
        return Result(res)

    def set_run(
        self,
        tend: float = None,
        tstep: float = None,
        save: list = None,
        maxsteps: int = None,
        rtlim: float = None,
    ):
        if tend is not None:
            self.tend = tend
        if tstep is not None:
            self.tstep = tstep
        if save is not None:
            self.save = save
        if maxsteps is not None:
            self.maxsteps = maxsteps
        if rtlim is not None:
            self.rtlim = rtlim

    def addkept(self, reac: str):
        """reactions that are kept"""
        self.get_reaction(self.ruleset.reac_from_name(reac)).addkept()

    def get_compound(self, name: str) -> Chemical:
        if name == "":
            return Chemical("", self)
        try:
            return self.pool_collect[name]
        except KeyError:
            new_comp = Chemical(name, self)
            self.pool_collect[name] = new_comp
            return new_comp

    def get_reaction(self, description: ReacDescr) -> Reaction:
        try:
            return self.reac_collect[description.name]
        except KeyError:
            new_reac = Reaction(description, self)
            self.reac_collect[description.name] = new_reac
            return new_reac


# def multirun(compute, nbthread: int = None):  ### Check if still used somewhere...
#    ctx = get_context("fork")
#    if nbthread is None:
#        nbthread = ctx.cpu_count()
#    with ctx.Pool(nbthread) as pool:
#        res = pool.map(compute, range(nbthread))
#    return Result(res)
