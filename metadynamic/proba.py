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
    """Return the index i where proba < sum(data) from 0 to i"""
    res: float = 0.0
    for index, val in enumerate(data):
        res += val
        if proba < res:
            return index
    return -1


class Probalist:
    def __init__(self, vol: float = 1, maxlength: int = 100):
        self.vol = vol
        # List of objects ### Check replacement with list instead of numpy array !
        self._mapobj = np.zeros(maxlength, dtype=np.dtype("O"))
        # List of probas
        self._problist = np.zeros(maxlength, dtype=np.float64)
        # Next available position
        self._actlist = 0
        self._maxlength = maxlength
        self.probtot = 0.0
        self._queue: Deque[int] = deque()
        self.sysrand = SystemRandom()

    def register(self, obj: Any) -> int:
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
        if isvalid(proba_pos):
            self.probtot -= self._problist[proba_pos]
            self._problist[proba_pos] = 0.0
            self._mapobj[proba_pos] = None
            self._queue.append(proba_pos)

    def update(self, proba_pos: int, proba: float) -> None:
        # assertion shall greatly reduce perf for non-optimized python code!
        assert isvalid(proba_pos)
        #  get proba change from the event
        delta = proba - self._problist[proba_pos]
        #  Set the new probability of the event
        self._problist[proba_pos] = proba
        #  Update the probability of the proba sum
        self.probtot += delta

    def choose(self) -> Tuple[Any, float]:
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
            dt = np.log(1 / self.sysrand.random()) / self.probtot
            if chosen is None:
                LOGGER.error(
                    f"Badly destroyed reaction from {self._mapobj} with proba {self._problist}"
                )
            return (chosen, dt)
        except ValueError as v:
            raise RoundError(
                f"(reason: {v}; "
                f"probtot={self.probtot}=?={self._problist.sum()}; "
                f"problist={self._problist})"
            )

    def clean(self) -> None:
        """(Re-)Compute probality for sums.
           It is far slower than just updating individual probabilities,
           but the fast methods used in update function leads to accumulate small rounding
           errors.
           This functions is intended to be called regularly for cleaning
           these rounding errors."""
        self.probtot = self._problist.sum()
