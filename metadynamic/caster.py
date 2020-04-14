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
metadynamic.caster
==================

generic type caster

Provides:
---------

 - L{Caster}: generic type-caster generator

"""

from typing import Type, Any, List


class Caster:
    """generic type-caster generator"""
    def __init__(self, target: Type[Any]):
        """
        creates an callable object that will try to cast any value to
        the given target.

        >>> caster = Caster(List[int])
        >>> caster([1, "2", 4.5])
        [1,2,4]

        @param target: the target type. may be complexe types based on typing module
        @type target: Type[Any]
        """
        self.dest: Type[Any]
        """destination conversion type"""
        self.args: List[Caster]
        """destination type arguments (e.g. for list or dict)"""
        if hasattr(target, "__origin__"):
            self.dest = target.__origin__
            args = target.__args__
            self.args = [Caster(args[0])] if Ellipsis in args else [Caster(arg) for arg in args]
        else:
            self.dest = target
            self.args = []

    def __call__(self, value: Any) -> Any:
        """
        cast the value to the pre-defined target type

        @param value: the value to be casted
        @return: the casted value
        """
        if isinstance(value, bytes):
            value = value.decode()
        if self.dest is dict:
            return {
                self.args[0](key): self.args[1](val) for key, val in dict(value).items()
            }
        if self.dest is list:
            return [self.args[0](val) for val in value]
        if self.dest is tuple:
            if len(self.args) == 1:
                return tuple([self.args[0](val) for val in value])
            return tuple([conv(val) for conv, val in zip(self.args, value)])
        return self.dest(value)

    def __repr__(self) -> str:
        """represent a caster as 'TO<destination type>'"""
        args = "" if len(self.args) == 0 else " , ".join([str(arg) for arg in self.args])
        return f"TO{self.dest}[{args}]"
