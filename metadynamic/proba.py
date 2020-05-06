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

"""Probalistic calculations (Gillespie's algorithm to be found here !).

It provides L{Probalist}, a class of object storing and performing all probalistic events

"""

from typing import Tuple, Deque, Iterable, Any
from collections import deque
from secrets import SystemRandom
from numba import jit, float64, int32

import numpy as np

from metadynamic.ends import RoundError, NoMore
from metadynamic.logger import LOGGER
from metadynamic.inval import isvalid


@jit(int32(float64[:], float64), nopython=True, cache=True)  # type: ignore
def choice(data: Iterable[float], proba: float) -> int:
    """Return the index i where proba < sum(data) from 0 to i.

    compiled with numba for performance gain

    Assumes that the total sum of probabilities is 1 !

    @param data: list of probalities
    @type data: Iterable[float]
    @param proba: probability (randomly chosed somewhere else)
    @type proba: float
    @return: index of chosen event
        (-1 if failed to find one... only occurs if sum of probabilities is < 1)
    @rtype: int

    """
    res: float = 0.0
    for index, val in enumerate(data):
        res += val
        if proba < res:
            return index
    return -1


class Probalist:
    """Class of object storing all probalistic events."""

    def __init__(self, vol: float = 1, maxlength: int = 100):
        """Create a Probalist object.

        @param vol: system volume for constant conversion
        @type vol: float
        @param maxlength: maximum probalistic events to be stored
            (when run out of space, size is increased by maxlength increments)
        @type maxlength: int

        """
        # FIX check if vol can be removed and placed somewhere else !
        self.vol: float = vol
        """system volume"""
        # Check replacement with list instead of numpy array !
        self._mapobj: np.ndarray = np.zeros(maxlength, dtype=np.dtype("O"))
        """List of probalistic objects"""
        self._problist: np.ndarray = np.zeros(maxlength, dtype=np.float64)
        """List of probas"""
        self._actlist: int = 0
        """Next available position"""
        self._maxlength: int = maxlength
        """maximum probalistic event increment to be stored"""
        self.probtot: float = 0.0
        """total probability"""
        self._queue: Deque[int] = deque()
        """Queue of freed positions"""
        self.sysrand: SystemRandom = SystemRandom()
        """random number generator (uses system generator for optimal entropy)"""

    def register(self, obj: Any) -> int:
        """Register a new probabilistic event.

        @param obj: object to be registered
        @return: storage index Probalist object
        @rtype: int

        """
        # Still room left from previously removed object
        if self._queue:
            nlist = self._queue.popleft()
            self._mapobj[nlist] = obj
            return nlist
        # No more free room, add at the end
        try:
            self._mapobj[self._actlist] = obj
        except IndexError:
            # arrays too small, extend them
            self._mapobj = np.append(self._mapobj, np.full(self._maxlength, None))
            self._problist = np.append(self._problist, np.zeros(self._maxlength))
            self._mapobj[self._actlist] = obj
        self._actlist += 1
        return self._actlist - 1

    def unregister(self, proba_pos: int) -> None:
        """Unregister object located at the given storage index.

        @param proba_pos: storage index of the object
        @param proba_pos: int

        """
        if isvalid(proba_pos):
            self.probtot -= self._problist[proba_pos]
            self._problist[proba_pos] = 0.0
            self._mapobj[proba_pos] = None
            self._queue.append(proba_pos)

    def update(self, proba_pos: int, proba: float) -> None:
        """Update the probability of an object.

        @param proba_pos: storage index of the object
        @param proba_pos: int
        @param proba: new probability
        @type proba: float

        """
        # assertion shall greatly reduce perf for non-optimized python code!
        assert isvalid(proba_pos)
        #  get proba change from the event
        delta = proba - self._problist[proba_pos]
        #  Set the new probability of the event
        self._problist[proba_pos] = proba
        #  Update the probability of the proba sum
        self.probtot += delta

    def choose(self) -> Tuple[Any, float]:
        """Choose a probabilistic event according to Gillespie's algorithm.

        @return: the chosen object, and the delta_t after which the event occured
        @rtype: Tuple[Any, float]

        @raise RoundError: if something bad happened (most likely being a rounding error somewhere)

        """
        # First choose a random line in the probability map
        try:
            # chosen = random.choice(self._mapobj, p=self._problist / self.probtot)
            chosen_num = choice(self._problist, self.probtot * self.sysrand.random())
            chosen = self._mapobj[chosen_num]
            if chosen_num == -1:
                if self.probtot == 0.0:
                    raise NoMore("choice() had nothing to choose!")
                raise RoundError(
                    "choice() couldn't find a suitable object."
                    f"probtot={self.probtot}=?={self._problist.sum()}; "
                    f"problist={self._problist})"
                )
            delta_t = np.log(1 / self.sysrand.random()) / self.probtot
            if chosen is None:
                LOGGER.error(
                    f"Badly destroyed reaction from {self._mapobj} with proba {self._problist}"
                )
            return (chosen, delta_t)
        except ValueError as err:
            raise RoundError(
                f"(reason: {err}; "
                f"probtot={self.probtot}=?={self._problist.sum()}; "
                f"problist={self._problist})"
            )

    def clean(self) -> None:
        """Recompute total probality.

        It is far slower than just updating individual probabilities, but the fast methods used in
        update function leads to accumulate small rounding errors.  This functions is intended to be
        called regularly for cleaning these rounding errors.

        """
        self.probtot = self._problist.sum()
