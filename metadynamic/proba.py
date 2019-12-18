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

from typing import List, Tuple, Deque
from collections import deque
from numpy import array, append, log, random, zeros, dtype

from metadynamic.ends import RoundError
from metadynamic.logger import Logged


class Activable:
    def __init__(self) -> None:
        self.activated: bool = False

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

    def process(self) -> None:
        pass

    def choose(self) -> "Activable":
        return self


class Probalistic(Logged):
    probalist: "Probalist"

    @classmethod
    def setprobalist(cls, vol: float = 1, minprob: float = 1e-10) -> None:
        cls.probalist = Probalist(vol, minprob)


class Probaobj(Probalistic):
    """Probalistic object, that can be used in a Probalist object
    for Gillespie-like computation.

    Each object to be stored in the probalist must contain one Probaobj object
    """

    def __init__(self, obj: Activable):
        self.proba_pos: List[int]
        self.obj: Activable = obj
        self.unset_proba_pos()

    def unset_proba_pos(self) -> None:
        self.proba_pos = []
        self.registered = False

    def update(self, newproba: float) -> None:
        # only perform update if something changed
        if newproba != 0.0:
            if not self.registered:
                # was unactivated, thus activate
                self.obj.activate()
                self.proba_pos = self.probalist.register(self.obj)
                self.registered = True
            self.probalist.update(self.proba_pos, newproba)
        elif self.registered:
            # was activated, thus deactivate
            self.probalist.unregister(self.proba_pos)
            self.unset_proba_pos()
            self.obj.unactivate()
        else:
            self.log.warning("Updated a reaction from 0 to 0...")


class Probalist(Logged):
    def __init__(self, vol: float = 1, minprob: float = 1e-10, maxlength=100):
        self.vol = vol
        self._minprob = minprob
        # List of objects
        self._mapobj = zeros([], dtype("O"))
        # List of probas
        self._problist = zeros([])
        # Next available position
        self._actlist = 0
        self._maxlength = maxlength
        self.probtot = 0.0
        self._queue: Deque[int] = deque()

    def register(self, obj: Activable) -> List[int]:
        # Still room left from previously removed object
        if self._queue:
            nlist = self._queue.popleft()
            # self.log.debug(f"Registering {obj} in position {nlist} from queue")
            self._mapobj[nlist] = obj
            return [nlist]
        # No more free room, add at the end
        try:
            self._mapobj[self._actlist] = obj
            # self.log.debug(f"Registering {obj} in position {self._actlist} at end")
        except IndexError:
            # arrays too small, extend them
            self._mapobj = append(self._mapobj, zeros(self._maxlength, dtype("O")))
            self._problist = append(self._problist, zeros(self._maxlength))
            self._mapobj[self._actlist] = obj
            # self.log.debug(f"Registering {obj} in position {self._actlist} at end after extension")
        self._actlist += 1
        return [self._actlist - 1]

    def unregister(self, proba_pos: List[int]) -> None:
        assert len(proba_pos) == 1
        nlist = proba_pos[0]
        # self.log.debug(f"Unregistering {self._mapobj[nlist]} in position {nlist}")
        self.probtot -= self._problist[nlist]
        self._problist[nlist] = 0.
        self._mapobj[nlist] = None
        self._queue.append(nlist)

    def update(self, proba_pos: List[int], proba: float) -> None:
        # assertion shall greatly reduce perf for non-optimized python code!
        assert len(proba_pos) == 1
        nlist = proba_pos[0]
        #  get proba change from the event
        delta = proba - self._problist[nlist]
        #  Set the new probability of the event
        self._problist[nlist] = proba
        #  Update the probability of the proba sum
        self.probtot += delta

    def choose(self) -> Tuple[Activable, float]:
        # First choose a random line in the probability map
        try:
            chosen = random.choice(self._mapobj, p=self._problist / self.probtot)
            if chosen is None or chosen == 0:
                self.log.error(f"Badly destroyed reaction from {self._mapobj} with proba {self._problist}")
            return (chosen, log(1 / random.rand()) / self.probtot)
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

    def get_probaobj(self, obj: Activable) -> Probaobj:
        return Probaobj(obj)

    @staticmethod
    def seed(seed: int) -> None:
        if seed == 0:
            random.seed()
        else:
            random.seed(seed)
