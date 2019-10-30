from itertools import chain
from typing import Generic, List, Optional, Callable, Set, TypeVar, Dict, Type, Any

from metadynamic.collector import Collect
from metadynamic.proba import Probaobj, Probalist
from metadynamic.ends import DecrZero
from metadynamic.description import ReacDescr, CompDescr
from metadynamic.logger import Logged

D = TypeVar("D", "ReacDescr", "CompDescr")
C = TypeVar("C", "CollectofReaction", "CollectofCompound")


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


class CollectofCompound(Collect["Compound"], Logged):
    _colltype = "Compound"
    
    def _create(self, name: str) -> "Compound":
        newcomp = Compound(CompDescr(name), self)
        return newcomp

    def _initialize(self, comp: "Compound") -> None:
        comp.initialize(self.system.reac_collect)

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


class CollectofReaction(Collect["Reaction"], Logged):
    _colltype = "Reaction"
    
    def _create(self, name: str) -> "Reaction":
        assert isinstance(name, str), f"{name} is not a string..."
        newreac = Reaction(self.ruleset.reac_from_name(name), self)
        return newreac

    def _initialize(self, reac: "Reaction") -> None:
        reac.initialize(
            self.system.param.vol, self.system.probalist, self.system.comp_collect
        )

    def _categorize(self, obj: "Reaction") -> List[str]:
        return obj.description.categories

    def init_ruleset(
        self,
        consts: Dict[str, float],
        altconsts: Dict[str, float],
        catconsts: Dict[str, float],
    ) -> None:
        self.ruleset = Ruleset(consts, altconsts, catconsts)

    def from_descr(self, description: ReacDescr) -> "Reaction":
        # try:
        return self[description.name]
        # except KeyError:
        #    newobj = self._create(description.name)
        #    self.pool[description.name] = newobj
        #    return newobj

    def describe(
        self, kind: str, reactants: List[str], catal: str, pos: int
    ) -> "Reaction":
        return self.from_descr(
            self.ruleset.reac_from_descr(kind, reactants, catal, pos)
        )

    def get_related(self, reactant: "Compound") -> List["Reaction"]:
        return [
            self.from_descr(description)
            for description in self.ruleset.full_reac(reactant)
        ]


class Chemical(Generic[D, C], Logged):
    _descrtype = "Chemical"

    def __init__(self, description: D, collect: C):
        self.description: D = description
        self.collect: C = collect
        self.activated: bool = False
        self.log.debug(f"Creating {self}")

    def __repr__(self) -> str:
        return f"{self._descrtype}: {self}"

    def __str__(self) -> str:
        return self.description.name

    def activate(self) -> None:
        self.log.debug(f"Trying to activate {self}...")
        if not self.activated:
            if self._start_activation():
                self.collect.activate(self.description.name)
                self.activated = True
                self.log.debug(f"{self} activation...done")
                self._finish_activation()
                self.log.debug(f"{self} activation...finished")

    def _start_activation(self) -> bool:
        """To be implemented in derived class (if needed)
        Will be processed *before* common activation procedure is performed
        Return True if activation can be processed"""
        return True

    def _finish_activation(self):
        """To be implemented in derived class (if needed)
        Will be processed *after* common activation procedure is performed"""
        pass

    def unactivate(self) -> None:
        self.log.debug(f"Trying to unactivate {self}...")
        if self.activated:
            if self._start_unactivation():
                self.collect.unactivate(self.description.name)
                self.activated = False
                self.log.debug(f"{self} unactivation...done")

    def _start_unactivation(self) -> bool:
        """To be implemented in derived class (if needed)
        Will be processed *before* common unactivation procedure is performed"""
        return True


class Reaction(Chemical[ReacDescr, CollectofReaction], Logged):
    _descrtype = "Reaction"

    def initialize(
        self, vol: float, probalist: Probalist, comp_collect: CollectofCompound
    ) -> None:
        self._vol: float = vol
        self._probaobj: Probaobj = probalist.get_probaobj(self)
        self.comp_collect = comp_collect
        self._reactants: List[Compound] = [
            self.comp_collect[i.name] for i in self.description.reactants
        ]
        self._productnames = self.description.build_products()
        self._products: Optional[List[Compound]] = None
        self.const = self.description.build_const()
        catal = self.description.catal
        if catal:
            self._catalized = True
            self._catal = self.comp_collect[catal.name]
        else:
            self._catalized = False
        self._reaccalc: Optional[Callable[[], float]] = None
        self._uncatcalc: Optional[Callable[[], float]] = None
        self._set_reaccalc()
        for comp in self._reactants:
            comp.register_reaction(self)
        if self._catalized:
            self._catal.register_reaction(self)
        #self.activate()

    def _start_unactivation(self) -> bool:
        assert self.calcproba() == 0.0
        self._probaobj.unregister()
        return True

    def _start_activation(self) -> bool:
        # Only activate if probability >0
        return self.calcproba() > 0.0

    @property
    def proba(self) -> float:
        return self._probaobj.proba

    def update(self) -> None:
        newproba = self.calcproba()
        if newproba > 0.0:
            self._probaobj.update(newproba)
        else:
            self.unactivate()
        assert self.proba == self.calcproba()

    def process(self) -> None:
        # If first time processed, activate
        # if not self.activated:
        #    self.activate()
        # Increment products
        self.log.debug(f"Processing {self}, p={self.proba}={self.calcproba}, concs={[c.pop for c in self._reactants]}")
        if self._products is None:
            self._products = [self.comp_collect[name] for name in self._productnames]
        for prod in self._products:
            prod.inc()
        # Decrement reactants
        for reac in self._reactants:
            try:
                reac.dec()
            except DecrZero as end:
                # Decremented from 0...
                # Thus exit with max of information
                detail = f" from {self}, that is activated? ({self.activated}) (p={self.proba}={self.calcproba()}, "
                for comp in self._reactants:
                    detail += f"[{comp.description}]={comp.pop} ,"
                if self._catalized:
                    detail += f"catal[{self._catal.description}]={self._catal.pop})"
                else:
                    detail += ")"
                raise DecrZero(end.detail + detail)

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

    def addkept(self) -> None:
        if self._catal:
            self._catal.addkept(self)
        for comp in self._reactants:
            comp.addkept(self)


class Compound(Chemical[CompDescr, CollectofCompound], Logged):
    _descrtype = "Compound"

    def initialize(self, reac_collect: CollectofReaction) -> None:
        self.reac_collect = reac_collect
        self._reactions: Set[Reaction] = set()
        self._kept_descr: Set[Reaction] = set()
        self.pop = 0
        self.length = self.description.length

    def _finish_activation(self) -> bool:
        self._reactions = self._kept_descr | set(self.reac_collect.get_related(self))
        for reac in self._reactions:
            reac.activate()
        return True

    #Only unactivate reaction from update
    #def _start_unactivation(self) -> bool:
    #    for reac in self._reactions.copy():
    #        reac.unactivate()
    #    return True

    def _upd_reac(self) -> None:
        # copy as list because self._reactions to be changed in loop. Elegant?...
        for reac in list(self._reactions):
            reac.update()

    def addkept(self, reaction: Reaction) -> None:
        self._kept_descr.add(reaction)

    def register_reaction(self, reaction: Reaction) -> None:
        self._reactions.add(reaction)

    def inc(self) -> None:
        self.pop += 1
        if self.pop == 1:
            self.activate()
        self._upd_reac()

    def dec(self) -> None:
        if self.pop > 0:
            self.pop -= 1
        else:
            raise DecrZero(self.description.name)
        if self.pop == 0:
            self.unactivate()
        else:
            self._upd_reac()

    def init_pop(self, start: int) -> None:
        if start != 0:
            if self.pop == 0:
                self.pop = start
                self.activate()
            else:
                self.pop = start
            self._upd_reac()
        else:
            self.pop = 0
            self.unactivate()
