from typing import Generic, List, Optional, Callable, Set, TypeVar

from metadynamic.proba import Probaobj
from metadynamic.ends import Finished, DecrZero
from metadynamic.description import ReacDescr, CompDescr

D = TypeVar("D", "ReacDescr", "CompDescr")


class Chemical(Generic[D]):
    _descrtype = "Chemical"

    def __init__(self, description: D, system: "System"):
        self.system = system
        self.description: D = description

    def __repr__(self) -> str:
        return f"{self._descrtype}: {self}"

    def __str__(self) -> str:
        return self.description.name


class Reaction(Chemical[ReacDescr]):
    _descrtype = "Reaction"

    def __init__(self, description: ReacDescr, system: "System"):
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

    def destroy(self) -> None:
        if self._probaobj.registered:
            self._probaobj.unregister()
            for comp in self._reactants:
                comp.unregister_reaction(self)
            if self._catalized:
                self._catal.unregister_reaction(self)
            self.system.reac_collect.remove(self.description.name)

    def register(self) -> None:
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

    def __init__(self, description: CompDescr, system: "System"):
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
