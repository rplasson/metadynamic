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

"""
metadynamic.chemical
==================

Objects specific to chemical data.


Provides
--------

CollectofCompound: Collect class for storing a pool of Compound

CollectofReaction: Collect class for storing a pool of Reaction

Chemical: Collectable class, for objects with chemical properties

Reaction: Chemical class, for the description of specific reactions

Compound: Chemical class, for the description of specific compounds

Crn: Chemical Reaction Network class: high level interface containing
     a CollectofCompound, a CollectofReaction, and a Model related to a
     given run.

"""


import numpy as np

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
    Iterable,
    Union,
)
from math import factorial
from itertools import repeat

from metadynamic.collector import Collect, Collectable
from metadynamic.proba import Probalist
from metadynamic.ends import DecrZero, NoMore, NotFound
from metadynamic.ruleset import Model, ReacDescr
from metadynamic.inval import invalidint
from metadynamic.logger import LOGGER
from metadynamic.inputs import Param


class Memcalc:
    """Function memoizer"""

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


def entro(x: Union[float, int]) -> float:
    """Return x×ln(x), 0 if x=0"""
    return 0.0 if x == 0 else float(x * np.log(x))


K = TypeVar("K", bound=Hashable)
C = TypeVar("C", "CollectofCompound", "CollectofReaction")


class CollectofCompound(Collect[str, "Compound"]):
    """Collect class for storing a pool of Compound"""

    _colltype = "Compound"

    def _create(self, name: str) -> "Compound":
        """
        Create a Compound from its key

        :param name: name of the compound to be created
        :type name: str
        :return: newly created Compound object
        :rtype: Compound
        """
        newcomp = Compound(name, self.crn)
        return newcomp

    def _categorize(self, obj: "Compound") -> Set[str]:
        """
        List the categories of the Compound

        :param obj: object to catagorize
        :type obj: Compound
        :return: set of categories
        :rtype: Set[str]
        """
        return self.model.descriptor.categories(obj.description)

    def set_crn(self, crn: "Crn") -> None:
        """
        set the parent Crn

        :param crn: parent Crn
        :type crn: Crn
        """
        self.crn: Crn = crn

    def _getprop(self, prop: str, obj: "Compound") -> float:
        """
        Return a given property of a Compound object.
        Shouldn't be directly used (used by proplist method only)

        The property may be:
          - "count": returns 1.0
          - "pop": returns the compound population
          - "entropy": returns pop×log(pop)
          - other values: a property defined in model.descriptor

        :param prop: property name
        :type prop: str
        :param obj: object to be probed
        :type obj: "Compound"
        """
        return float(
            1.0
            if prop == "count"
            else obj.pop
            if prop == "pop"
            else entro(obj.pop)
            if prop == "entropy"
            else self.model.descriptor.prop(prop, obj.description)
        )


class CollectofReaction(Collect[ReacDescr, "Reaction"]):
    """Collect class for storing a pool of Reaction"""

    _colltype = "Reaction"

    def _create(self, description: ReacDescr) -> "Reaction":
        """
        Create a Reaction from its key

        :param description: description of the compound to be created
        :type name: ReacDescr
        :return: newly created Reaction object
        :rtype: Reaction
        """
        newreac = Reaction(description, self.crn)
        return newreac

    def _categorize(self, obj: "Reaction") -> Set[str]:
        """
        List the categories of the Reaction.
        A reaction belongs to a sole category, its reaction type
        (i.e. the name of the rule that created it)

        :param obj: object to categorize
        :type obj: Reaction
        :return: set of categories
        :rtype: Set[str]
        """
        return {obj.description[0]}

    def set_crn(self, crn: "Crn") -> None:
        """
        set the parent Crn

        :param crn: parent Crn
        :type crn: Crn
        """
        self.crn: Crn = crn

    def _getprop(self, prop: str, obj: "Reaction") -> float:
        """
        Return a given property of a Reaction object.
        Shouldn't be directly used (used by proplist method only)

        The property may be:
          - "count": returns 1.0
          - "pop": returns the compound probability
          - "entropy": returns proba×log(proba)
          - other values: a property defined in model.descriptor

        :param prop: property name
        :type prop: str
        :param obj: object to be probed
        :type obj: "Compound"
        """
        return float(
            1.0
            if prop == "count"
            else obj.proba
            if prop == "rate"
            else entro(obj.proba)
            if prop == "entropy"
            else self.model.descriptor.prop(
                prop, str(obj.description)
            )  # Check here... properties of reactions?
        )


class Chemical(Generic[K], Collectable):
    """ """

    _descrtype = "Chemical"

    def __init__(self, description: K, crn: "Crn"):
        self.crn: Crn = crn
        self.description: K = description
        self.activated: bool = False

    def __repr__(self) -> str:
        return f"{self._descrtype}: {self}"

    def __str__(self) -> str:
        return str(self.description)

    def activate(self) -> None:
        """ """
        if not self.activated:
            self._activate()
            self.activated = True

    def _activate(self) -> None:
        """ """
        pass

    def unactivate(self) -> None:
        """ """
        if self.activated:
            self._unactivate()
            self.activated = False

    def _unactivate(self) -> None:
        """ """
        pass

    def update(self, change: int = 0) -> None:
        """to be implemented in subclasses

        :param change: (Default value = 0)
        :type change: int :
        """
        raise NotImplementedError

    def delete(self) -> None:
        """to be implemented in subclasses"""
        raise NotImplementedError


class Reaction(Chemical[ReacDescr]):
    """ """

    _descrtype = "Reaction"
    _updatelist: Dict[Chemical[ReacDescr], int] = {}

    def __init__(self, description: ReacDescr, crn: "Crn"):
        super().__init__(description, crn)
        self.name: str = ""
        # If name is empty => invalid reaction, no process to be done
        if description[0] != "":
            self.proba: float = 0.0
            (
                stoechreac,
                self._stoechproduct,
                const,
                self.robust,
            ) = self.crn.model.ruleset.buildreac(self.description)
            self.stoechio: List[Tuple[Compound, int]] = []
            order: int = 0
            # stochastic rate between n reactions must be divided by V^(n-1)
            # (k <-> concentration, stoch rate <-> number of particles)
            self.const = const
            for reacname, stoechnum in stoechreac:
                self.stoechio.append((self.crn.comp_collect[reacname], stoechnum))
                order += stoechnum
                # stochastic rate implying 2 distinct compounds is to be diveded by two
                # as it is not proportional to X.X, nor X.(X-1), but X.(X-1)/2,
                # for order N, it is X.(X-1)...(X-n+1)/n!
                # the n1 part is only computed here.
                #
                # /!\  Check if this cannot be computed earlier using ruleset.Parameters facilities
                #
                if stoechnum > 1:
                    self.const /= fact(stoechnum)
            self.const /= self.crn.probalist.vol ** (order - 1)
            self.tobeinitialized = True
            self._unset_proba_pos()

    def _unset_proba_pos(self) -> None:
        """ """
        self.proba_pos: int = invalidint
        self.registered: bool = False

    def _activate(self) -> None:
        """ """
        self.crn.reac_collect.activate(self.description)
        for comp, _ in self.stoechio:
            comp.register_reaction(self)

    def _unactivate(self) -> None:
        """ """
        self.crn.reac_collect.unactivate(self.description)
        for comp, _ in self.stoechio:
            comp.unregister_reaction(self)

    def update(self, change: int = 0) -> None:
        """

        :param change: (Default value = 0)
        :type change: int :
        :param change: int:  (Default value = 0)

        
        """
        newproba, changed = self.updateproba()
        if changed:
            # only perform update if something changed
            if newproba != 0.0:
                if not self.registered:
                    # was unactivated, thus activate
                    self.activate()
                    self.proba_pos = self.crn.probalist.register(self)
                    self.registered = True
                self.crn.probalist.update(self.proba_pos, newproba)
            elif self.registered:
                # was activated, thus deactivate
                self.delete()

    def process(self) -> None:
        """ """
        if self.tobeinitialized:
            self.products = [
                (self.crn.comp_collect[name], order)
                for name, order in self._stoechproduct
            ]
            if self.robust:
                self.tobeinitialized = False
            else:
                self._stoechproduct = self.crn.model.ruleset.rebuild_prod(
                    self.description
                )
        for prod, order in self.products:
            prod.change_pop(order)
        # Decrement reactants
        for reac, order in self.stoechio:
            reac.change_pop(-order)
        self.crn.update()
        if self.crn.probalist.probtot == 0:
            raise NoMore(f"after processing {self}")

    def _ordern(self, pop: int, order: int) -> int:
        """

        :param pop: 
        :type pop: int :
        :param order: 
        :type order: int :
        :param pop: int: 
        :param order: int: 

        
        """
        return pop if order == 1 else pop * self._ordern(pop - 1, order - 1)

    def updateproba(self) -> Tuple[float, bool]:
        """ """
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
        """ """
        if self.registered:
            self.crn.probalist.unregister(self.proba_pos)
            self._unset_proba_pos()
        self.unactivate()

    def serialize(self) -> Any:
        """ """
        return self.const, self.proba

    @staticmethod
    def _join_compounds(stoechio: Iterable[Tuple[Any, int]]) -> str:
        """

        :param stoechio: 
        :type stoechio: Iterable[Tuple[Any :
        :param int]]: 
        :param stoechio: Iterable[Tuple[Any: 

        
        """
        return "+".join(
            [f"{num}{name}" if num > 1 else str(name) for name, num in stoechio]
        )

    def __str__(self) -> str:
        if not self.name:
            self.name = "->".join(
                [
                    self._join_compounds(self.stoechio),
                    self._join_compounds(self._stoechproduct),
                ]
            )
        return self.name


class Compound(Chemical[str]):
    """ """

    _descrtype = "Compound"
    _updatelist: Dict[Chemical[str], int] = {}

    def __str__(self) -> str:
        #  Already a string, conversion useless, thus overload
        return self.description

    def __init__(self, description: str, crn: "Crn"):
        super().__init__(description, crn)
        if self.description == "":
            LOGGER.error("Created empty compound!!!")
        self.reactions: Set[Reaction] = set()
        self.pop: int = 0
        # self.length = self.descriptor("length", description)

    def _activate(self) -> None:
        """ """
        self.crn.comp_collect.activate(self.description)
        self.scan_reaction()

    def _unactivate(self) -> None:
        """ """
        self.crn.comp_collect.unactivate(self.description)

    def register_reaction(self, reaction: Reaction) -> None:
        """

        :param reaction: 
        :type reaction: Reaction :
        :param reaction: Reaction: 

        
        """
        self.reactions.add(reaction)

    def unregister_reaction(self, reaction: Reaction) -> None:
        """

        :param reaction: 
        :type reaction: Reaction :
        :param reaction: Reaction: 

        
        """
        try:
            self.reactions.remove(reaction)
        except KeyError:
            LOGGER.debug(
                f"Tried to unregister twice {reaction} from {self} (p={self.pop})"
            )

    def scan_reaction(self) -> None:
        """ """
        self.reactions = {
            self.crn.reac_collect[descr]
            for descr in self.crn.model.ruleset.get_related(
                self.description, self.crn.comp_collect.categories
            )
        }

    def update(self, change: int = 0) -> None:
        """

        :param change: (Default value = 0)
        :type change: int :
        :param change: int:  (Default value = 0)

        
        """
        if change != 0:
            # LOGGER.debug(f"Really updating {self}")
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
            # for reac in self.reactions:  # impactedreac:
            #     Reaction.toupdate(reac)

    def change_pop(self, start: int) -> None:
        """

        :param start: 
        :type start: int :
        :param start: int: 

        
        """
        self.crn.comp_toupdate(self, start)

    def delete(self) -> None:
        """ """
        self._unactivate()
        self.change_pop(-self.pop)

    def serialize(self) -> Any:
        """ """
        return self.pop


class Crn:
    """ """

    def __init__(self, param: Param):
        # update trackers
        self._reac_update: Set[Reaction] = set()
        self._comp_update: Dict[Compound, int] = {}
        # Create Crn objects
        self.model = Model(param.rulemodel, param.reactions, param.parameters)
        self.vol = param.vol
        self.dropmode = param.dropmode
        self.init_collect()
        for compound, pop in param.init.items():
            self.comp_collect[compound].change_pop(pop)
        self.update()
        LOGGER.debug(f"Initialized with {param}")

    def init_collect(self) -> None:
        """ """
        self.probalist = Probalist(self.vol)
        self.comp_collect = CollectofCompound(self.model, dropmode="keep")
        self.comp_collect.set_crn(self)
        self.reac_collect = CollectofReaction(
            self.model, dropmode=self.dropmode
        )  # set categorize to False or True?
        self.reac_collect.set_crn(self)

    def close(self) -> None:
        """ """
        #  (to be tested): create new data for cleaning memory
        self.init_collect()

    def clean(self) -> None:
        """ """
        self.probalist.clean()

    def stepping(self) -> float:
        """ """
        # choose a random event
        chosen, dt = self.probalist.choose()
        # check if there even was an event to choose
        if chosen is None:
            raise NotFound(f"dt={dt}")
        # perform the (chosen one) event
        chosen.process()
        return dt

    def reac_toupdate(self, reac: Reaction) -> None:
        """Add Reaction 'reac' to the updatelist

        :param reac: 
        :type reac: Reaction :
        :param reac: Reaction: 

        
        """
        self._reac_update.add(reac)

    def comp_toupdate(self, comp: Compound, change: int) -> None:
        """Add Compound 'comp' to the updatelist,
            por a population variation of 'change'

        :param comp: 
        :type comp: Compound :
        :param change: 
        :type change: int :
        :param comp: Compound: 
        :param change: int: 

        
        """
        try:
            self._comp_update[comp] += change
        except KeyError:
            self._comp_update[comp] = change

    def update(self) -> None:
        """Perform full Crn update"""
        for comp, change in self._comp_update.items():
            # Update compounds
            comp.update(change)
            # List impacted reactions
            for reac in comp.reactions:
                self.reac_toupdate(reac)
        self._comp_update = {}
        for reac in self._reac_update:
            reac.update()
        self._reac_update = set()

    def collstat(
        self, collection: str, prop: str, weight: str, method: str, full: bool
    ) -> float:
        """Get statistics

        :param collection: 
        :type collection: str :
        :param prop: 
        :type prop: str :
        :param weight: 
        :type weight: str :
        :param method: 
        :type method: str :
        :param full: 
        :type full: bool :
        
        """
        return (
            self.comp_collect.stat(prop, weight, method, full)
            if collection == "compounds"
            else self.reac_collect.stat(prop, weight, method, full)
        )

    def collmap(
        self,
        collection: str,
        prop: str,
        weight: str,
        sort: str,
        method: str,
        full: bool,
    ) -> Dict[float, float]:
        """Get a statistic map

        :param collection: 
        :type collection: str :
        :param prop: 
        :type prop: str :
        :param weight: 
        :type weight: str :
        :param sort: 
        :type sort: str :
        :param method: 
        :type method: str :
        :param full: 
        """
        return (
            self.comp_collect.map(prop, weight, sort, method, full)
            if collection == "compounds"
            else self.reac_collect.map(prop, weight, sort, method, full)
        )
