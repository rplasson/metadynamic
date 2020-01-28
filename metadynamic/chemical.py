#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2019 by Raphaël Plasson
#
# This file is part of metadynamic
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.

# from itertools import chain
from typing import (
    Generic,
    List,
    Callable,
    TypeVar,
    Dict,
    Any,
    Set,
    Hashable,
    Tuple,
)
from math import factorial
from itertools import repeat

from metadynamic.collector import Collect, Collectable
from metadynamic.proba import Probalistic
from metadynamic.ends import DecrZero
from metadynamic.ruleset import Ruled, ReacDescr
from metadynamic.inval import isvalid, Invalid, invalidint


class Memcalc:
    def __init__(self, func: Callable[[int], int]):
        self.results: Dict[int, int] = {}
        self.func = func

    def __call__(self, val: int) -> int:
        try:
            return self.results[val]
        except KeyError:
            self.results[val] = self.func(val)
            return self.results[val]


# only use the value of order! (but often), with order rarely above 3...
fact = Memcalc(factorial)


# For naming a reaction... to be moved in chemical (instead of reacdescr.name ?)
#    @property
#    def name(self) -> str:
#        return f"{self.rule.name}:{'+'.join(self.reactantnames)}"


K = TypeVar("K", bound=Hashable)
C = TypeVar("C", "CollectofCompound", "CollectofReaction")


class CollectofCompound(Collect[str, "Compound"]):
    _colltype = "Compound"

    def _create(self, name: str) -> "Compound":
        newcomp = Compound(name)
        return newcomp

    def _categorize(self, obj: "Compound") -> Set[str]:
        return self.descriptor.categories(obj.description)

    def dist(
        self, prop: str, lenweight: bool = False, full: bool = False
    ) -> Dict[Any, int]:
        res: Dict[Any, int] = {}
        search = self.pool if full else self.active
        for comp in search.values():
            prop_value = self.descriptor.prop(prop, comp.description)
            inc = comp.pop if lenweight else 1
            if prop_value not in res:
                res[prop_value] = inc
            else:
                res[prop_value] += inc
        return res


class CollectofReaction(Collect[ReacDescr, "Reaction"]):
    _colltype = "Reaction"

    def _create(self, description: ReacDescr) -> "Reaction":
        newreac = Reaction(description)
        return newreac

    def _categorize(self, obj: "Reaction") -> Set[str]:
        return {obj.description[0]}


class Collected(Ruled):
    comp_collect: CollectofCompound
    reac_collect: CollectofReaction

    @classmethod
    def setcollections(
        cls,
        categorize_comp: bool = True,
        dropmode_comp: str = "",
        categorize_reac: bool = False,
        dropmode_reac: str = "",
    ) -> None:
        cls.comp_collect = CollectofCompound(categorize_comp, dropmode_comp)
        cls.reac_collect = CollectofReaction(categorize_reac, dropmode_reac)


class Chemical(Generic[K], Collected, Collectable):
    _descrtype = "Chemical"
    _updatelist: Dict["Chemical[K]", int]

    def __init__(self, description: K):
        self.description: K = description
        self.activated: bool = False
        # self.log.debug(f"Creating {self}")

    def __repr__(self) -> str:
        return f"{self._descrtype}: {self}"

    def __str__(self) -> str:
        return str(self.description)

    def activate(self) -> None:
        if not self.activated:
            self._activate()
            self.activated = True

    def _activate(self) -> None:
        pass

    def unactivate(self) -> None:
        if self.activated:
            self._unactivate()
            self.activated = False

    def _unactivate(self) -> None:
        pass

    @classmethod
    def toupdate(cls, obj: "Chemical[K]", change: int = 0) -> None:
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

    def delete(self) -> None:
        """to be implemented in subclasses"""
        raise NotImplementedError


class Reaction(Chemical[ReacDescr], Probalistic):
    _descrtype = "Reaction"
    _updatelist: Dict[Chemical[ReacDescr], int] = {}

    def __init__(self, description: ReacDescr):
        super().__init__(description)
        self.name: str = ""
        # If name is empty => invalid reaction, no process to be done
        if description[0] != "":
            self.proba: float = 0.0
            stoechreac, self._stoechproduct, const, = self.ruleset.buildreac(
                self.description
            )
            self.stoechio: List[Tuple[Compound, int]] = []
            order: int = 0
            # stochastic rate between n reactions must be divided by V^(n-1)
            # (k <-> concentration, stoch rate <-> number of particles)
            self.const = const
            for reacname, stoechnum in stoechreac:
                self.stoechio.append((self.comp_collect[reacname], stoechnum))
                order += stoechnum
                # stochastic rate implying 2 distinct compounds is to be diveded by two
                # as it is not proportional to X.X, nor X.(X-1), but X.(X-1)/2,
                # for order N, it is X.(X-1)...(X-n+1)/n!
                # the n1 part is only computed here.
                if stoechnum > 1:
                    self.const /= fact(stoechnum)
            self.const /= self.probalist.vol ** (order - 1)
            self.tobeinitialized = True
            self._unset_proba_pos()

    def _unset_proba_pos(self) -> None:
        self.proba_pos: int = invalidint
        self.registered: bool = False

    def _activate(self) -> None:
        self.reac_collect.activate(self.description)
        for comp, _ in self.stoechio:
            comp.register_reaction(self)

    def _unactivate(self) -> None:
        self.reac_collect.unactivate(self.description)
        for comp, _ in self.stoechio:
            comp.unregister_reaction(self)

    def update(self, change: int = 0) -> None:
        newproba, changed = self.updateproba()
        if changed:
            # only perform update if something changed
            if newproba != 0.0:
                if not self.registered:
                    # was unactivated, thus activate
                    self.activate()
                    self.proba_pos = self.probalist.register(self)
                    self.registered = True
                self.probalist.update(self.proba_pos, newproba)
            elif self.registered:
                # was activated, thus deactivate
                self.delete()

    def process(self) -> None:
        if self.tobeinitialized:
            self.products = [
                (self.comp_collect[name], order) for name, order in self._stoechproduct
            ]
            self.tobeinitialized = False
        for prod, order in self.products:
            prod.change_pop(order)
        # Decrement reactants
        for reac, order in self.stoechio:
            reac.change_pop(-order)
        trigger_changes(self)

    def _ordern(self, pop: int, order: int) -> int:
        return pop if order == 1 else pop * self._ordern(pop - 1, order - 1)

    def updateproba(self) -> Tuple[float, bool]:
        oldproba = self.proba
        self.proba = self.const
        for reactant, stoechnum in self.stoechio:
            # self.proba *= self._ordern(reactant.pop, stoechnum)
            pop = reactant.pop
            self.proba *= pop
            for _ in repeat(None, stoechnum - 1):
                pop -= 1
                self.proba *= pop
        return self.proba, self.proba != oldproba

    def delete(self) -> None:
        if self.registered:
            self.probalist.unregister(self.proba_pos)
            self._unset_proba_pos()
        self.unactivate()

    def serialize(self) -> Any:
        return self.const, self.proba

    def __str__(self) -> str:
        if not self.name:
            self.name = "->".join(
                [
                    "+".join(
                        [
                            f"{num}{name}" if num > 1 else str(name)
                            for name, num in stoechio
                        ]
                    )
                    for stoechio in [self.stoechio, self._stoechproduct]
                ]
            )
        return self.name


class Compound(Chemical[str]):
    _descrtype = "Compound"
    _updatelist: Dict[Chemical[str], int] = {}

    def __str__(self) -> str:
        #  Already a string, conversion useless, thus overload
        return self.description

    def __init__(self, description: str):
        super().__init__(description)
        if self.description == "":
            self.log.error("Created empty compound!!!")
        self.reactions: Set[Reaction] = set()
        self.pop = 0
        # self.length = self.descriptor("length", description)

    def _activate(self) -> None:
        self.comp_collect.activate(self.description)
        self.scan_reaction()

    def _unactivate(self) -> None:
        self.comp_collect.unactivate(self.description)

    def register_reaction(self, reaction: Reaction) -> None:
        self.reactions.add(reaction)

    def unregister_reaction(self, reaction: Reaction) -> None:
        try:
            self.reactions.remove(reaction)
        except KeyError:
            self.log.debug(
                f"Tried to unregister twice {reaction} from {self} (p={self.pop})"
            )

    def scan_reaction(self) -> None:
        self.reactions = {
            self.reac_collect[descr]
            for descr in self.ruleset.get_related(
                self.description, self.comp_collect.categories
            )
        }

    def update(self, change: int = 0) -> None:
        if change != 0:
            # self.log.debug(f"Really updating {self}")
            pop0 = self.pop
            self.pop = pop0 + change
            if self.pop < 0:
                #  population cannot be negative
                raise DecrZero(self.description)
            if self.pop == 0:
                #  unactivate compounds whose population reaches 0
                self.unactivate()
            elif pop0 == 0:
                #  activate compounds whose population switched from 0 to nonzero
                self.activate()
            for reac in self.reactions:  # impactedreac:
                Reaction.toupdate(reac)

    def change_pop(self, start: int) -> None:
        Compound.toupdate(self, start)

    def delete(self) -> None:
        self._unactivate()
        self.change_pop(-self.pop)

    def serialize(self) -> Any:
        return self.pop


class InvalidReaction(Invalid, Reaction):
    _invalrepr = "Invalid Reaction"


invalidreaction = InvalidReaction(("", (), 0))


def trigger_changes(fromreac: Reaction = invalidreaction) -> None:
    try:
        Compound.trigger_update()
    except DecrZero as end:
        if not isvalid(fromreac):
            detail = end.detail
        else:
            # Decremented from 0...
            # Thus exit with max of information
            detail = (
                end.detail
                + f" from {fromreac}, that is activated? ({fromreac.activated})"
                + f" (p={fromreac.proba}, "
            )
            for comp, _ in fromreac.stoechio:
                detail += f"[{comp.description}]={comp.pop} ,"
            else:
                detail += ")"
        raise DecrZero(detail)
    Reaction.trigger_update()
