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

from typing import Callable, TypeVar, Any

T = TypeVar("T")
A = TypeVar("A")


def memoize_property(func: Callable[[T], Any]) -> Callable[[T], Any]:
    def memoized(self: T) -> Any:
        try:
            return self.__getattribute__(f"_{func.__name__}_memoized")
        except AttributeError:
            self.__setattr__(f"_{func.__name__}_memoized", func(self))
            return self.__getattribute__(f"_{func.__name__}_memoized")

    return memoized


def memoize_oneparam(func: Callable[[T, A], Any]) -> Callable[[T, A], Any]:
    def memoized(self: T, pos: A) -> Any:
        try:
            return self.__getattribute__(f"_{func.__name__}_{pos}_memoized")
        except AttributeError:
            self.__setattr__(f"_{func.__name__}_{pos}_memoized", func(self, pos))
            return self.__getattribute__(f"_{func.__name__}_{pos}_memoized")

    return memoized
