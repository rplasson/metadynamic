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

from typing import Generic, TypeVar, Dict, Set, List, Union
from weakref import WeakValueDictionary

from metadynamic.logger import Logged

K = TypeVar("K")
T = TypeVar("T")


class WeakDict(Generic[K, T], WeakValueDictionary):
    pass


WDict = Union[Dict[str, T], WeakDict[str, T]]


class Collect(Generic[T], Logged):
    _colltype = "Generic"

    def __init__(self, categorize: bool = True, dropmode: str = "drop"):
        self.dropmode = dropmode
        self.pool: WDict[T]
        if self.dropmode == "soft":
            self.pool = WeakDict()
        else:
            self.pool = {}
        self.categories: Dict[str, Set[str]] = {}
        self.active: WDict[T] = self.pool if self.dropmode == "drop" else {}
        self.categorize = categorize
        self.log.info(
            f"Created {self} as drop={self.dropmode}, with pool of type {type(self.pool)}"
        )

    def __repr__(self) -> str:
        return f"<Collect of {len(self.pool)} {self._colltype}>"

    def __getitem__(self, name: str) -> T:
        """ Return the object as described by its name
            If it is the first call of the object, create it
            Else, return the already created one"""
        try:
            return self.pool[name]
        except KeyError:
            newobj = self._create(name)
            if self.dropmode == "soft":
                self.active[name] = newobj
            self.pool[name] = newobj
            return newobj

    def activate(self, name: str) -> None:
        """Put the object 'name' in the active section, then categorize it"""
        obj = self[name]
        # will fail if activate an duplicated object (i.e. 2 different objects/same name exists)
        assert obj is self[name]
        if name not in self.active:
            self.active[name] = obj
        if self.categorize:
            for catname in self._categorize(obj):
                try:
                    self.categories[catname].add(name)
                except KeyError:
                    self.categories[catname] = {name}

    def unactivate(self, name: str) -> None:
        """Remove the object 'name' from the active section, then
           uncategorize it"""
        try:
            # del self.active[name]
            self.active.pop(name)
        except KeyError:
            self.log.debug(f"Tried to unactivate twice {name}")
        if self.categorize:
            for cat in self.categories.values():
                try:
                    cat.remove(name)
                except KeyError:
                    pass

    def cat_list(self, category: str) -> Set[str]:
        """Return all (active) objects from the specified 'categories'"""
        try:
            return self.categories[category]
        except KeyError:
            return set()

    def _create(self, name: str) -> T:
        """Create the object <T> from its name.
        Must be implemented in subclasses"""
        raise NotImplementedError

    def _categorize(self, obj: T) -> List[str]:
        """List the categories of the object.
        Must be implemented in subclasses"""
        raise NotImplementedError
