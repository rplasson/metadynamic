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
"""Classes for collecting objects into a collection.

It provides L{Collectable} for classes of objects that can be collected, and L{Collect}: for the
collection Collectable objects.

"""

from typing import Generic, TypeVar, Dict, Set, Hashable, Any
from collections import defaultdict

import numpy as np

from metadynamic.ruleset import Model
from metadynamic.ends import BadFile
from metadynamic.logger import LOGGER


class Collectable:
    """Generic class, to be derived in classes that will be collected in Collect."""

    def delete(self) -> None:
        """Clean memory of the object.

        It is only intended for a final purge when the thread exits.  Manual call may break many
        things!

        To be implemented in subclasses

        """
        raise NotImplementedError

    def serialize(self) -> Any:
        """Output the value corresponding to the collected object for being dumped in a json file.

        To be implemented in subclasses

        """
        raise NotImplementedError


K = TypeVar("K", bound=Hashable)
"""generic Hashable type"""
T = TypeVar("T", bound=Collectable)
"""generic Collectable type"""


class Collect(Generic[K, T]):
    """Object collector.

    Generic class to be implemented for Collectable objects of type T,
    indexed with Hashable key of type K.
    Keeps two dictionaries {K:T} (a pool with all objects, an activ for only active ones

    """

    _colltype = "Generic"
    """string description of the collected object"""

    def __init__(self, model: Model, categorize: bool = True, dropmode: str = "drop"):
        """Create a collection from a model (for categorizing and building objects).

        @param model: Model containing all rules
        @type model: Model
        @param categorize: shall the objects be automatically categorized (Default value = True)
        @type categorize: bool
        @param dropmode: how the unactivated object are dealt
            ("drop"=> removed, "keep"=> kept in pool)
            (Default value = "drop")
        @type dropmode: str

        """
        self.model: Model = model
        """rule model"""
        self.dropmode: str = dropmode
        """drop mode 'keep' or 'drop'"""
        self.pool: Dict[K, T] = {}
        """pool of all objects"""
        self.categories: Dict[str, Set[K]] = defaultdict(set)
        """categories of objects"""
        self.active: Dict[K, T] = self.pool if self.dropmode == "drop" else {}
        """pool of active objects"""
        self.categorize: bool = categorize
        """shall the objects be categorized"""
        LOGGER.debug(
            f"Created {self} as drop={self.dropmode}, with pool of type {type(self.pool)}"
        )

    def __repr__(self) -> str:
        """Format the collection representation as <Collect of <n> T>.

        @return: formatted dtring
        @rtype: str

        """
        return f"<Collect of {len(self.pool)} {self._colltype}>"

    def __getitem__(self, key: K) -> T:
        """Return the object as described by its key.

        If it is the first call of the object, create it Else, return the already created one.

        """
        try:
            return self.pool[key]
        except KeyError:
            newobj = self._create(key)
            if self.dropmode == "soft":
                self.active[key] = newobj
            self.pool[key] = newobj
            return newobj

    def activate(self, key: K) -> None:
        """Put the object 'key' in the active section, then categorize it.

        @param key: object key
        @type key: K

        """
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
        """Remove the object 'key' from the active section, then uncategorize it.

        @param key: object key
        @type key: K

        """
        try:
            # del self.active[key]
            self.active.pop(key)
        except KeyError:
            pass  # Ok to try to unactivate already unactivated object
        if self.categorize:
            for cat in self.categories.values():
                try:
                    cat.remove(key)
                except KeyError:
                    pass

    def cat_list(self, category: str) -> Set[K]:
        """Return all (active) objects from the specified 'categories'.

        @param category: category name
        @type category: str
        @return: set of compounds in the corresponding category
        @rtype: Set[K]

        """
        try:
            return self.categories[category]
        except KeyError:
            return set()

    def _create(self, key: K) -> T:
        """Create the object <T> from its key.

        Must be implemented in subclasses

        @param key: object key
        @type key: K

        """
        raise NotImplementedError

    def _categorize(self, obj: T) -> Set[str]:
        """List the categories the object belongs to.

        Must be implemented in subclasses

        @param obj: object to categorize
        @type obj: T
        @return: list of the categories
        @rtype: Set[str]

        """
        raise NotImplementedError

    def save(self, full: bool = False) -> Dict[str, T]:
        """Return the collected data as a dict without value change.

        @param full: collect all pool data if True, else all active data
        @type full: bool
        @return: data dictionary
        @rtype: Dict[str, T]

        """
        dataset = self.pool if full else self.active
        return {str(val): val for val in dataset.values()}

    def asdict(self, full: bool = False) -> Dict[str, Any]:
        """Return the collected data as a dict, serializing the values.

        @param full: collect all pool data if True, else all active data
        @type full: bool
        @return: data dictionary
        @rtype: Dict[str, Any]

        """
        dataset = self.pool if full else self.active
        return {str(val): val.serialize() for val in dataset.values()}

    def _getprop(self, prop: str, obj: T) -> float:
        """Return a given property of a <T> object.

        Shouldn't be directly used (used by proplist method only)

        If not implemented in derived classes, only returns
        1.0 to "count" property

        @param prop: property name
        @type prop: str
        @param obj: object to be probed
        @type obj: T

        """
        if prop != "count":
            raise BadFile(f"Unknown property {prop}")
        return 1.0

    def proplist(self, prop: str, full: bool = False) -> np.ndarray:
        """Return a list of all property values of the collection.

        The data order is not important, used for statistics.

        @param prop: name of the property to be collected
        @type prop: str
        @param full: if False, only collect properties of active objects
            (Default value = False)
        @type full: bool
        @return: array of properties
        @rtype: nv.ndarray

        """
        search = self.pool if full else self.active
        return np.array(
            [self._getprop(prop, obj) for obj in search.values()], dtype=float
        )

    def stat(self, prop: str, weight: str, method: str, full: bool = False) -> float:
        """Return statistics on all property values of the collection.

        @param prop: name of the property to be collected
        @type prop: str
        @param weight: name of the property to be used as weight
        @type weight: str
        @param method: '+' returns the weighted sum, 'm' the weighted average,
            'max' the weighted max value, 'min' the weighted min value
        @type method: str
        @param full: if False, only collect properties of active objects
            (Default value = False)
        @type full: bool
        @return: stat value
        @rtype: float

        """
        values = self.proplist(prop, full)
        # FIX check use of "single"... still dubious...
        weights = 1 / values if weight == "single" else self.proplist(weight, full)
        if method == "+":
            return float(np.nansum(values * weights))
        if method == "m":
            try:
                return float(np.average(values, weights=weights))
            except ZeroDivisionError:
                return float(np.nan)  # float conversion for mypy...
        if method == "max":
            return float(max(values * weights))
        if method == "min":
            return float(min(values * weights))
        raise BadFile(f"the method {method} is not recognized")

    def map(
        self, prop: str, weight: str, sort: str, method: str, full: bool = False
    ) -> Dict[float, float]:
        """Return a statistics maps on all property values of the collection.

        @param prop: name of the property to be collected
        @type prop: str
        @param weight: name of the property to be used as weight
        @type weight: str
        @param sort: name of the property to be used as a sorting key
        @type sort: str
        @param method: '+' returns the weighted sum, 'm' the weighted average
        @type method: str
        @param full: if False, only collect properties of active objects
            (Default value = False)
        @type full: bool
        @return: stat values
        @rtype: Dict[float, float]

        """
        res: Dict[float, float] = {}
        tot: Dict[float, float] = {}
        values = self.proplist(prop, full)
        weights = 1 / values if weight == "single" else self.proplist(weight, full)
        sorts = self.proplist(sort, full)
        for val, wgh, srt in zip(values, weights, sorts):
            try:
                res[srt] += val * wgh
                tot[srt] += wgh
            except KeyError:
                res[srt] = val * wgh
                tot[srt] = wgh
        if method == "+":
            return res
        for srt in res:
            res[srt] /= tot[srt]
        return res
