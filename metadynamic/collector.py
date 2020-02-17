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

from typing import Generic, TypeVar, Dict, Set, Union, Hashable, List, Any
from weakref import WeakValueDictionary
from collections import defaultdict
from numpy import array, sum, average, ndarray

# from metadynamic.logger import Logged
from metadynamic.ruleset import Ruled
from metadynamic.ends import BadFile


class Collectable:
    def delete(self) -> None:
        """
        Clean memory of the object.
        It is only intended for a final purge when the thread exits.
        Manual call may break many things!

        To be implemented in subclasses"""
        raise NotImplementedError

    def serialize(self) -> Any:
        """
        Output the value corresponding to the collected object
        In order to be dumped in a json file

        To be implemented in subclasses"""
        raise NotImplementedError


K = TypeVar("K", bound=Hashable)
T = TypeVar("T", bound=Collectable)


class WeakDict(Generic[K, T], WeakValueDictionary):
    pass


WDict = Union[Dict[K, T], WeakDict[K, T]]


class Collect(Generic[K, T], Ruled):
    _colltype = "Generic"

    def __init__(self, categorize: bool = True, dropmode: str = "drop"):
        self.dropmode = dropmode
        self.pool: WDict[K, T]
        if self.dropmode == "soft":
            self.pool = WeakDict()
        else:
            self.pool = {}
        self.categories: Dict[str, Set[K]] = defaultdict(set)
        self.active: WDict[K, T] = self.pool if self.dropmode == "drop" else {}
        self.categorize = categorize
        self.log.info(
            f"Created {self} as drop={self.dropmode}, with pool of type {type(self.pool)}"
        )

    def __repr__(self) -> str:
        return f"<Collect of {len(self.pool)} {self._colltype}>"

    def __getitem__(self, key: K) -> T:
        """ Return the object as described by its key
            If it is the first call of the object, create it
            Else, return the already created one"""
        try:
            return self.pool[key]
        except KeyError:
            newobj = self._create(key)
            if self.dropmode == "soft":
                self.active[key] = newobj
            self.pool[key] = newobj
            return newobj

    def activate(self, key: K) -> None:
        """Put the object 'key' in the active section, then categorize it"""
        obj = self[key]
        # will fail if activate an duplicated object (i.e. 2 different objects/same key exists)
        assert obj is self[key]
        if key not in self.active:
            self.active[key] = obj
        if self.categorize:
            for catkey in self._categorize(obj):
                try:
                    self.categories[catkey].add(key)
                except KeyError:
                    self.categories[catkey] = {key}

    def unactivate(self, key: K) -> None:
        """Remove the object 'key' from the active section, then
           uncategorize it"""
        try:
            # del self.active[key]
            self.active.pop(key)
        except KeyError:
            self.log.debug(f"Tried to unactivate twice {key}")
        if self.categorize:
            for cat in self.categories.values():
                try:
                    cat.remove(key)
                except KeyError:
                    pass

    def cat_list(self, category: str) -> Set[K]:
        """Return all (active) objects from the specified 'categories'"""
        try:
            return self.categories[category]
        except KeyError:
            return set()

    def _create(self, key: K) -> T:
        """Create the object <T> from its key.
        Must be implemented in subclasses"""
        raise NotImplementedError

    def _categorize(self, obj: T) -> Set[str]:
        """List the categories of the object.
        Must be implemented in subclasses"""
        raise NotImplementedError

    def purge(self) -> None:
        keys: List[K] = list(self.pool.keys())
        for key in keys:
            self.pool[key].delete()
            try:
                del self.pool[key]
            except KeyError:
                # Already purged by delete process
                pass

    def save(self, full: bool = False) -> Dict[str, T]:
        dataset = self.pool if full else self.active
        return {str(val): val for val in dataset.values()}

    def getprop(self, prop: str, obj: T) -> float:
        if prop != "count":
            raise BadFile(f"Unknown property {prop}")
        return 1.0

    def proplist(self, prop: str, full: bool = False) -> ndarray:
        search = self.pool if full else self.active
        return array([self.getprop(prop, obj) for obj in search.values()], dtype=float)

    def stat(self, prop: str, weight: str, method: str, full: bool = False) -> float:
        values = self.proplist(prop, full)
        weights = self.proplist(weight, full)
        if method == "+":
            return float(sum(values * weights))
        if method == "m":
            return float(average(values, weights=weights))
        if method == "max":
            return float(max(values * weights))
        if method == "min":
            return float(min(values * weights))
        raise BadFile(f"the method {method} is not recognized")

    def map(self, prop: str, weight: str, sort: str, method: str, full: bool = False) -> Dict[float, float]:
        res: Dict[float, float] = {}
        tot: Dict[float, float] = {}
        values = self.proplist(prop, full)
        weights = self.proplist(weight, full)
        sorts = self.proplist(sort, full)
        for v, w, s in zip(values, weights, sorts):
            try:
                res[s] += v*w
                tot[s] += w
            except KeyError:
                res[s] = v*w
                tot[s] = w
        if method == "+":
            return res
        for s in res:
            res[s] /= tot[s]
        return res
