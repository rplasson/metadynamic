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
        super().__init__()

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
        #print(f"Update {nlist} by {delta}")
        newproba = self._problist[nlist] + delta
        # reset to 0 probabilities that are too low (rounding errors correction)
        if newproba < self._minprob:
            delta += newproba
            newproba = 0.0
            #print(f"...but reset to 0")
        self._problist[nlist] = newproba
        self.probtot += delta

    def unregister(self, obj: Reaction):
        #print(f"kill {obj}")
        nlist, npos = obj.proba_pos
        oldproba = self._map[nlist][npos]
        self._map[nlist][npos] = 0
        self._mapobj[nlist][npos] = None
        self._queue.append((nlist, npos))
        self._updateprob(nlist, -oldproba)
        obj.unset_proba_pos()

    def update(self, obj: Reaction, proba: float):
        #print(f"update {obj}")
        if obj.registered:
            nlist, npos = obj.proba_pos
            delta = proba - self._map[nlist][npos]
            self._map[nlist][npos] = proba
            self._updateprob(nlist, delta)
        else:
            if self._queue:
                obj.set_proba_pos(self._addfromqueue(obj, proba))
            else:
                obj.set_proba_pos(self._addfrommap(obj, proba))
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
        # First choose a random line in the probability map
        try:
            nlist = N.random.choice(self._nblist, p=self._problist / self.probtot)
        except ValueError as v:
            raise RoundError(
                f"(reason: {v}; probtot={self.probtot}=?={self._problist.sum()}; problist={self._problist})"
            )
        # Then choose a random column in the chosen probability line
        try:
            return N.random.choice(
                self._mapobj[nlist], p=self._map[nlist] / self._problist[nlist]
            )
        except ValueError as v:
            raise RoundError(
                f"(reason: {v}; probtot={self._problist[nlist]}=?={self._map[nlist].sum()}; problist={self._map[nlist]})"
            )

    def clean(self):
        """(Re-)Compute probality for sums.
           It is far slower than just updating individual probabilities,
           but the fast methods used in update function leads to accumulate small rounding
           errors.
           This functions is intended to be called regularly for cleaning
           these rounding errors."""
        #old_problist = self._problist
        #old_probtot = self.probtot
        self._problist = N.array([data.sum() for data in self._map])
        self.probtot = self._problist.sum()
        #print(self.probtot - old_probtot, self._problist - old_problist)


class Descr(str):
    def __new__(cls, name, *vals, **kvals):
        return str.__new__(cls, name)

    def __init__(self, name: str):
        self.name = name
        super().__init__()

    def __add__(self, other):
        return self.__class__(super().__add__(other))

    def __radd__(self, other):
        return self.__add__(other)

    def __getitem__(self, key):
        return self.__class__(super().__getitem__(key))

    @property
    def categories(self) -> List[str]:
        return []


class CompDescr(Descr):
    # Memoize properties? May be called often
    @property
    def ismono(self) -> bool:
        return len(self) == 1

    @property
    def ispolym(self) -> bool:
        return self.isalpha()

    @property
    def isact(self) -> bool:
        return self[-1] == "*" and self[:-1].isalpha()

    @property
    def categories(self) -> List[str]:
        res = []
        if self.ispolym:
            res.append("polym")
            if self.ismono:
                res.append("mono")
        elif self.isact:
            res.append("actpol")
        return res


class ReacDescr(Descr):
    kind = "X"
    nbreac = 0

    def __init__(
        self,
        name: str,
        reactants: List[CompDescr],
        catal: CompDescr,
        pos: int,
        const: float,
        altconst: float,
        catconst: float,
    ):
        super().__init__(name)
        self.reactants: List[CompDescr] = reactants
        self.catal: CompDescr = catal
        self._const: float = const
        self._altconst: float = altconst
        self._catconst: float = catconst
        self.pos: int = pos
        self.checkdimer()

    @property
    def categories(self) -> List[str]:
        return [self.kind]

    def checkdimer(self):
        self.dimer = False

    def build_products(self) -> List[CompDescr]:
        return []

    def build_const(self) -> float:
        const = (self._altconst ** self.nbalt()) * self._const
        if self.catal:
            const *= self._catconst
        return const

    def nbalt(self) -> int:
        return 0

    @staticmethod
    def fullfrom(reactant: Compound, system: System) -> List[ReacDescr]:
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

    def build_products(self) -> List[CompDescr]:
        return [self.reactants[0] + self.reactants[1]]

    def nbalt(self) -> int:
        return 1 if self.reactants[0].ismono else 0

    @staticmethod
    def fullfrom(reactant: Compound, system: System) -> List[ReacDescr]:
        first = lambda other: system.ruleset.reac_from_descr(
            "P", [reactant.description, other], "", -1
        )
        last = lambda other: system.ruleset.reac_from_descr(
            "P", [other, reactant.description], "", -1
        )
        return [
            create(other)
            for other in system.comp_collect.cat_list("polym")
            for create in (first, last)
        ]


class ActpolymDescr(ReacDescr):
    kind = "A"

    def build_products(self) -> List[CompDescr]:
        return [self.reactants[0][:-1] + self.reactants[1]]

    def nbalt(self) -> int:
        name0 = self.reactants[0][-2]
        name1 = self.reactants[1][0]
        return 0 if samecase(name0, name1) else 1

    @staticmethod
    def fullfrom(reactant: Compound, system: System) -> List[ReacDescr]:
        if reactant in system.comp_collect.cat_list("polym"):
            return [
                system.ruleset.reac_from_descr(
                    "A", [other, reactant.description], "", -1
                )
                for other in system.comp_collect.cat_list("actpol")
            ]
        if reactant in system.comp_collect.cat_list("actpol"):
            return [
                system.ruleset.reac_from_descr(
                    "A", [reactant.description, other], "", -1
                )
                for other in system.comp_collect.cat_list("polym")
            ]
        return []


class ActivationDescr(ReacDescr):
    kind = "a"
    nbreac = 1

    def build_products(self) -> List[CompDescr]:
        return [self.reactants[0] + "*"]

    def nbalt(self) -> int:
        return 0 if len(self.reactants[0]) == 1 else 1

    @staticmethod
    def fullfrom(reactant: Compound, system: System) -> List[ReacDescr]:
        if reactant in system.comp_collect.cat_list("polym"):
            return [system.ruleset.reac_from_descr("a", [reactant.description], "", "")]
        return []


class DectivationDescr(ReacDescr):
    kind = "d"
    nbreac = 1

    def build_products(self) -> List[CompDescr]:
        return [self.reactants[0][:-1]]

    def nbalt(self) -> int:
        return 0 if len(self.reactants[0]) == 2 else 1

    @staticmethod
    def fullfrom(reactant: Compound, system: System) -> List[ReacDescr]:
        if reactant in system.comp_collect.cat_list("actpol"):
            return [system.ruleset.reac_from_descr("d", [reactant.description], "", "")]
        return []


class HydroDescr(ReacDescr):
    kind = "H"
    nbreac = 1

    def build_products(self) -> List[CompDescr]:
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
    def fullfrom(reactant: Compound, system: System) -> List[ReacDescr]:
        return [
            system.ruleset.reac_from_descr("H", [reactant.description], "", i)
            for i in range(1, reactant.length)
        ]


class RacemDescr(ReacDescr):
    kind = "R"
    nbreac = 1

    def build_products(self) -> List[CompDescr]:
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
    def fullfrom(reactant: Compound, system: System) -> List[ReacDescr]:
        return [
            system.ruleset.reac_from_descr("R", [reactant.description], "", i)
            for i in range(reactant.length)
        ]


class EpimDescr(RacemDescr):
    kind = "E"

    @staticmethod
    def fullfrom(reactant: Compound, system: System) -> List[ReacDescr]:
        if reactant.description.ismono:
            return []
        return [system.ruleset.reac_from_descr("E", [reactant.description], "", 0)]


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

    def full_reac(self, reactant: Compound, system: System) -> List[ReacDescr]:
        return list(
            chain.from_iterable(
                [
                    descriptor.fullfrom(reactant, system)
                    for descriptor in self.descrs.values()
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
    def testkinds(cls, dic):
        badkeys = list(dic.keys() - ReacDescr.kinds())
        if badkeys:
            raise ValueError(f"'{badkeys}' are not recognized reaction types")


class Chemical:
    def __init__(self, description: Descr, system: System):
        self._system = system
        self.description = description

    def __repr__(self):
        return self.description.name


class Reaction(Chemical):
    def __init__(self, description: ReacDescr, system: System):
        super().__init__(description, system)
        self._reactants = [self._system.comp_collect[i] for i in description.reactants]
        self._products = description.build_products()
        self.const = description.build_const()
        self._started = False
        catal = description.catal
        if catal:
            self._catalized = True
            self._catal = system.comp_collect[Descr(catal)]
        else:
            self._catalized = False
        self._set_reaccalc()
        self.register()
        self.unset_proba_pos()

    def set_proba_pos(self, proba_pos: Tuple[int, int], register: bool = True):
        self.nlist, self.npos = proba_pos
        self.registered = register

    def unset_proba_pos(self):
        self.set_proba_pos((-1, -1), False)

    @property
    def proba_pos(self) -> Tuple[int, int]:
        return self.nlist, self.npos

    def destroy(self):
        if self.registered:
            self._system.probalist.unregister(self)
            for comp in self._reactants:
                comp.unregister_reaction(self)
            if self._catalized:
                self._catal.unregister_reaction(self)
            self._system.reac_collect.remove(self.description)

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
                    detail += f"[{comp.description}]={comp.pop} ,"
                if self._catalized:
                    detail += f"catal[{self._catal.description}]={self._catal.pop})"
                else:
                    detail += ")"
                raise DecrZero(f.detail + detail)
        if not self._started:
            self._products = [
                self._system.comp_collect[name] for name in self._products
            ]
            self._started = True
        for prod in self._products:
            prod.inc()

    def __repr__(self):
        return self.description.name

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
        order = self.description.nbreac
        dimer = self.description.dimer
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
            self._catal.addkept(self.description)
        for comp in self._reactants:
            comp.addkept(self.description)


class Compound(Chemical):
    def __init__(self, description: CompDescr, system: System):
        super().__init__(description, system)
        if description == "":
            self._real = False
        else:
            self._real = True
            self._reactions: Set[Reaction] = set()
            self._system = system
            self._kept_descr: List[ReacDescr] = []
            self.pop = 0
            self.length = len(description)

    def __bool__(self):
        return self._real

    def inc(self):
        self.pop += 1
        if self.pop == 1:
            self._activate()
        self._upd_reac()

    def dec(self):
        if self.pop > 0:
            self.pop -= 1
        else:
            raise DecrZero(self.description)
        if self.pop == 0:
            self._unactivate()
        else:
            self._upd_reac()

    def init_pop(self, start: int):
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

    def _activate(self): # /!\ Clean activation/unactivation/register/unregister mess /!\ Check with v4
        self._system.comp_collect.activate(self.description)
        for reac in self.reac_descr():
            self._system.reac_collect[reac]

    def _unactivate(self):
        for reac in self.reac_descr():
            self._system.reac_collect[reac].destroy()
        self._system.comp_collect.unactivate(self.description)

    def _upd_reac(self):
        # copy as list because self._reactions to be changed in loop. Elegant?...
        for reac in list(self._reactions):
            reac.update()

    def reac_descr(self) -> List[ReacDescr]:
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


class Collection:
    _type = Chemical
    _descrtype = Descr

    def __init__(self, system: System, drop=False):
        self._system = system
        self.pool: Dict[Descr, Chemical] = {}
        self.active: Dict[Descr, Chemical] = {}
        self.categories: Dict[str, Set[Descr]] = {}
        self.drop = drop

    def __getitem__(self, description: str) -> Chemical:
        if not isinstance(description, self._descrtype):
            description = self._descrtype(description)
        try:
            return self.pool[description]
        except KeyError:
            newobj = self._type(description, self._system)
            self.pool[description] = newobj
            return newobj

    def cat_list(self, category: str) -> Set[Descr]:
        try:
            return self.categories[category]
        except KeyError:
            return set()

    def remove(self, description: Descr):
        if self.drop:
            del self.pool[description]

    def activate(self, description: Descr):
        if description not in self.active:
            self.active[description] = self[description]
        for catname in description.categories:
            try:
                self.categories[catname].add(description)
            except KeyError:
                self.categories[catname] = {description}

    def unactivate(self, description: Descr):
        try:
            del self.active[description]
        except KeyError:
            pass
        for cat in self.categories.values():
            try:
                cat.remove(description)
            except KeyError:
                pass


class CollectofReaction(Collection):
    _type = Reaction
    _descrtype = ReacDescr


class CollectofCompound(Collection):
    _type = Compound
    _descrtype = CompDescr


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
        self.time = 0.0
        self.conc = conc
        self.dropreac = dropreac
        self.autoclean = autoclean
        # If seed set, always restart from the same seed. For timing/debugging purpose
        self._seed = seed
        self.step = 0
        self._ptot = sum([pop * len(comp) for comp, pop in init.items()])
        self.vol = self._ptot / self.conc
        self.comp_collect = CollectofCompound(self)
        self.reac_collect = CollectofReaction(self, dropreac)  # set of active reactions
        self.probalist = Probalist(minprob=minprob)
        for compound, pop in init.items():
            self.comp_collect[compound].init_pop(pop)
        if save is None:
            save = []
        self.set_run(tend=tend, tstep=tstep, save=save, maxsteps=maxsteps, rtlim=rtlim)

    def runtime_init(self):
        self._ptime0 = process_time()

    @property
    def runtime(self) -> float:
        return process_time() - self._ptime0

    @property
    def memuse(self) -> float:
        return Process(getpid()).memory_info().rss / 1024 / 1024

    def conc_of(self, compound: str) -> float:
        try:
            return self.comp_collect.active[CompDescr(compound)].pop / self.vol
        except KeyError:  # compounds doesn't exist => conc=0
            return 0.0

    @property
    def poplist(self) -> dict:
        return {comp: comp.pop for comp in self.comp_collect.active.values()}

    @property
    def lendist(self) -> dict:
        res: Dict[int, int] = {}
        for comp in self.comp_collect.active.values():
            length = comp.length
            pop = comp.pop
            if length not in res:
                res[length] = 0
            res[length] += pop
        return res

    @property
    def pooldist(self) -> dict:
        res: Dict[int, int] = {}
        for comp in self.comp_collect.pool.values():
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
        stat["nbcomp"] = len(self.comp_collect.active)
        stat["nbreac"] = len(self.reac_collect.pool)
        stat["poolsize"] = len(self.comp_collect.pool)
        stat["maxlength"] = max(dist["lendist"])
        return stat, dist

    def _process(self, tstop: float):
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
        N.random.seed(self._seed)  # necessary for multiprocessing from different seeds
        self.time = 0.0
        tnext = 0.0
        step = 0
        self.runtime_init()
        Process(getpid()).cpu_affinity([num % cpu_count()])
        while True:
            try:
                if output:
                    print(f"#{step} -> {tnext}")
                self._process(tnext)
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
            res = pool.map(self.run, range(nbthread))
        return Result(res)

    def set_run(
        self,
        *,  # force to name arguments
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
        self.reac_collect[self.ruleset.reac_from_name(reac)].addkept()
