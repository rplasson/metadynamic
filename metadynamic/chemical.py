from itertools import chain
from typing import Generic, List, Optional, Callable, TypeVar, Dict, Type, Any, Set

from metadynamic.collector import Collect
from metadynamic.proba import Probaobj, Probalistic
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

    def __init__(self, drop: bool = False, categorize: bool = True, dropmode: str = ""):
        super().__init__(drop, categorize)
        Collected.setcompcoll(self)

    def _create(self, name: str) -> "Compound":
        newcomp = Compound(CompDescr(name))
        return newcomp

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


class CollectofReaction(Collect["Reaction"], Logged, Probalistic):
    _colltype = "Reaction"

    def __init__(self, drop: bool = False, categorize: bool = True, dropmode: str = ""):
        super().__init__(drop, categorize, dropmode)
        Collected.setreaccoll(self)

    def _create(self, name: str) -> "Reaction":
        assert isinstance(name, str), f"{name} is not a string..."
        newreac = Reaction(self.ruleset.reac_from_name(name))
        return newreac

    def _categorize(self, obj: "Reaction") -> List[str]:
        return obj.description.categories

    def init_ruleset(
        self,
        consts: Dict[str, float],
        altconsts: Dict[str, float],
        catconsts: Dict[str, float],
    ) -> None:
        self.ruleset = Ruleset(consts, altconsts, catconsts)

    def describe(
        self, kind: str, reactants: List[str], catal: str, pos: int
    ) -> "Reaction":
        return self[self.ruleset.reac_from_descr(kind, reactants, catal, pos).name]

    def get_related(self, reactant: "Compound") -> List["Reaction"]:
        return [
            self[description.name] for description in self.ruleset.full_reac(reactant)
        ]


class Collected:
    comp_collect: CollectofCompound
    reac_collect: CollectofReaction

    @classmethod
    def setcompcoll(cls, comp_collect: CollectofCompound) -> None:
        cls.comp_collect = comp_collect

    @classmethod
    def setreaccoll(cls, reac_collect: CollectofReaction) -> None:
        cls.reac_collect = reac_collect


class Chemical(Generic[D, C], Logged, Collected):
    _descrtype = "Chemical"
    _updatelist: Dict["Chemical[D, C]", int]

    @property
    def collect(self) -> C:
        """Should return the collection the object belongs to. ti be implemented in subclasses"""
        raise NotImplementedError

    def __init__(self, description: D):
        self.description: D = description
        self.activated: bool = False
        # self.log.debug(f"Creating {self}")

    def __repr__(self) -> str:
        return f"{self._descrtype}: {self}"

    def __str__(self) -> str:
        return self.description.name

    def activate(self) -> None:
        if not self.activated:
            self.collect.activate(self.description.name)
            self.activated = True

    def unactivate(self) -> None:
        self.collect.unactivate(self.description.name)
        self.activated = False

    @classmethod
    def toupdate(cls, obj: "Chemical[D, C]", change: int = 0) -> None:
        """Add object to update"""
        try:
            cls._updatelist[obj] += change
        except KeyError:
            cls._updatelist[obj] = change
        # cls.log.debug(f"Should update {cls._updatelist} for {cls}")

    @classmethod
    def trigger_update(cls) -> None:
        """Trigger update events"""
        for obj, change in cls._updatelist.items():
            obj.update(change)
        cls._updatelist = {}

    def update(self, change: int = 0) -> None:
        """to be implemented in subclasses"""
        raise NotImplementedError


class Reaction(Chemical[ReacDescr, CollectofReaction], Logged, Collected, Probalistic):
    _descrtype = "Reaction"
    _updatelist: Dict[Chemical[ReacDescr, CollectofReaction], int] = {}

    @property
    def collect(self) -> CollectofReaction:
        return Chemical.reac_collect

    def __init__(self, description: ReacDescr):
        super().__init__(description)
        self._vol: float = self.probalist.vol
        self._probaobj: Probaobj = self.probalist.get_probaobj(self)
        self.proba: float = 0.0
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

    def update(self, change: int = 0) -> None:
        oldproba = self.proba
        self.proba = self.calcproba()
        # only perform update if something changed
        if oldproba != self.proba:
            if self.proba != 0.0:
                if oldproba == 0.0:
                    # was unactivated, thus activate
                    self.activate()
                    self._probaobj.register()
                    for comp in self._reactants:
                        comp.register_reaction(self)
                    if self._catalized:
                        self._catal.register_reaction(self)
                self._probaobj.update(self.proba)
                # assert self.proba == self.calcproba()
            else:
                if oldproba != 0.0:
                    # was activated, thus deactivate
                    self._probaobj.unregister()
                    self.unactivate()
                    for comp in self._reactants:
                        comp.unregister_reaction(self)
                    if self._catalized:
                        self._catal.unregister_reaction(self)

    def process(self) -> None:
        if self._products is None:
            self._products = [self.comp_collect[name] for name in self._productnames]
        for prod in self._products:
            prod.inc()
        # Decrement reactants
        for reac in self._reactants:
            reac.dec()
        trigger_changes(self)

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
    _updatelist: Dict[Chemical[CompDescr, CollectofCompound], int] = {}

    @property
    def collect(self) -> CollectofCompound:
        return Chemical.comp_collect

    def __init__(self, description: CompDescr):
        super().__init__(description)
        self.reactions: Set[Reaction] = set()
        self.pop = 0
        self.length = self.description.length

    def addkept(self, reaction: Reaction) -> None:
        pass  # to be removed. Will create an "explicit" reac description instead

    def register_reaction(self, reaction: Reaction) -> None:
        self.reactions.add(reaction)

    def unregister_reaction(self, reaction: Reaction) -> None:
        try:
            self.reactions.remove(reaction)
        except KeyError:
            if not reaction.description.dimer:
                # Only sensible reason for a double unregistration would be for dimer
                self.log.debug(
                    f"Tried to unregister twice {reaction} from {self} (p={self.pop})"
                )

    def scan_reaction(self) -> None:  # Check if is useful!
        self.reactions = set(self.reac_collect.get_related(self))

    def update(self, change: int = 0) -> None:
        if change != 0:
            # self.log.debug(f"Really updating {self}")
            pop0 = self.pop
            self.pop = pop0 + change
            if self.pop < 0:
                raise DecrZero(self.description.name)
            if self.pop == 0:
                self.unactivate()
            elif pop0 == 0:
                self.activate()
                self.scan_reaction()
            for reac in self.reactions:  # impactedreac:
                Reaction.toupdate(reac)

    def inc(self) -> None:
        Compound.toupdate(self, +1)

    def dec(self) -> None:
        Compound.toupdate(self, -1)

    def init_pop(self, start: int) -> None:
        Compound.toupdate(self, start)


def trigger_changes(fromreac: Optional[Reaction] = None) -> None:
    try:
        Compound.trigger_update()
    except DecrZero as end:
        if fromreac is None:
            detail = end.detail
        else:
            # Decremented from 0...
            # Thus exit with max of information
            detail = (
                end.detail
                + f" from {fromreac}, that is activated? ({fromreac.activated}) (p={fromreac.proba}={fromreac.calcproba()}, "
            )
            for comp in fromreac._reactants:
                detail += f"[{comp.description}]={comp.pop} ,"
            if fromreac._catalized:
                detail += f"catal[{fromreac._catal.description}]={fromreac._catal.pop})"
            else:
                detail += ")"
        raise DecrZero(detail)
    Reaction.trigger_update()
