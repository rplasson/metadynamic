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

from typing import Any, Tuple, Deque
from collections import deque
from numpy import array, append, log, random

from metadynamic.ends import RoundError
from metadynamic.logger import Logged
from metadynamic.inval import invalidint, isvalid


class Activable:
    def __init__(self) -> None:
        self.activated: bool = False

    def activate(self) -> None:
        if not self.activated:
            self._activate()
            self.activated = True

    def _activate(self) -> None:
        raise NotImplementedError

    def unactivate(self) -> None:
        if self.activated:
            self._unactivate()
            self.activated = False

    def _unactivate(self) -> None:
        raise NotImplementedError


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

    def __init__(self, obj: Any):
        self.nlist: int
        self.npos: int
        self.obj: Activable = obj
        self.unset_proba_pos()

    def unset_proba_pos(self) -> None:
        self.nlist = invalidint
        self.npos = invalidint
        self.registered = False

    @property
    def proba_pos(self) -> Tuple[int, int]:
        #  assert for perf ??? Assume always OK?
        if self.registered:
            return self.nlist, self.npos
        else:
            raise ValueError("Unregistered")

    def update(self, oldproba: float, newproba: float) -> None:
        # only perform update if something changed
        if oldproba != newproba:
            if newproba != 0.0:
                if oldproba == 0.0:
                    # was unactivated, thus activate
                    self.obj.activate()
                    self.nlist, self.npos = self.probalist.register(self)
                    self.registered = True
                self.probalist.update(self, newproba)
            else:
                if oldproba != 0.0:
                    # was activated, thus deactivate
                    self.probalist.unregister(self)
                    self.obj.unactivate()

    @property
    def proba(self) -> float:
        return self.probalist.getproba(self)


class Probalist(Logged):
    def __init__(self, vol: float = 1, minprob: float = 1e-10):
        self.vol = vol
        self._minprob = minprob
        #  /!\  Start from pre-determined array sized that do not grow instead of append???
        #  May gain performance...
        # ... but need to track correctly the effective size
        self._map = [array([])]
        self._mapobj = [array([])]
        self.npblist = 1
        self._actlist = 0
        self._problist = array([0.0])
        self.probtot = 0.0
        self._queue: Deque[Tuple[int, int]] = deque()

    def _updateprob(self, nlist: int, delta: float) -> None:
        # updating by delta avoid to recompute proba sums
        newproba = self._problist[nlist] + delta
        # reset to 0 probabilities that are too low
        # for rounding errors correction
        if newproba < self._minprob:
            delta += newproba
            newproba = 0.0
        #  update column proba
        self._problist[nlist] = newproba
        #  update total proba
        self.probtot += delta

    def register(self, probaobj: Probaobj) -> Tuple[int, int]:
        # free place from queue list
        if self._queue:
            return self._addfromqueue(probaobj)
        # No free place, create a new one
        return self._addfrommap(probaobj)

    def unregister(self, probaobj: Probaobj) -> None:
        nlist, npos = probaobj.proba_pos
        oldproba = self._map[nlist][npos]
        self._map[nlist][npos] = 0
        self._mapobj[nlist][npos] = None
        self._queue.append((nlist, npos))
        self._updateprob(nlist, -oldproba)
        probaobj.unset_proba_pos()

    def update(self, probaobj: Probaobj, proba: float) -> None:
        assert probaobj.registered
        nlist, npos = probaobj.proba_pos
        #  get proba change from the event
        delta = proba - self._map[nlist][npos]
        #  Set the new probability of the event
        self._map[nlist][npos] = proba
        #  Update the probability of the proba sums (column and tot)
        self._updateprob(nlist, delta)

    def getproba(self, probaobj: Probaobj) -> float:
        if probaobj.registered:
            nlist, npos = probaobj.proba_pos
            return self._map[nlist][npos]
        return 0.0

    def _addfrommap(self, newobj: Probaobj) -> Tuple[int, int]:
        rpos = len(self._map[self._actlist])
        if rpos <= self.npblist:
            rlist = self._actlist
            self._map[rlist] = append(self._map[rlist], 0)
            self._mapobj[rlist] = append(self._mapobj[rlist], newobj.obj)
            return rlist, rpos
        if len(self._map[0]) <= self.npblist:
            self._actlist = 0
        else:
            self._actlist += 1
            if self._actlist >= self.npblist:
                self._map.append(array([]))
                self._mapobj.append(array([]))
                self._problist = append(self._problist, 0.0)
                self.npblist += 1
        return self._addfrommap(newobj)

    def _addfromqueue(self, newobj: Probaobj) -> Tuple[int, int]:
        rlist, rpos = self._queue.popleft()
        self._map[rlist][rpos] = 0
        self._mapobj[rlist][rpos] = newobj.obj
        return rlist, rpos

    def choose(self) -> Tuple[Any, float]:
        # First choose a random line in the probability map
        try:
            nlist = random.choice(self.npblist, p=self._problist / self.probtot)
        except ValueError as v:
            raise RoundError(
                f"(reason: {v}; "
                f"probtot={self.probtot}=?={self._problist.sum()}; "
                f"problist={self._problist})"
            )
        # Then choose a random column in the chosen probability line
        try:
            chosen = random.choice(
                self._mapobj[nlist], p=self._map[nlist] / self._problist[nlist]
            )
            if chosen is None:
                explanation = [
                    (a, b)
                    for a, b in zip(self._mapobj[nlist], self._map[nlist])
                    if a is not None and a() is None
                ]
                self.log.error(f"Badly destroyed reaction in {explanation}")
            return (chosen, log(1 / random.rand()) / self.probtot)
        except ValueError as v:
            raise RoundError(
                f"(reason: {v}; "
                f"probtot={self._problist[nlist]}=?={self._map[nlist].sum()}; "
                f"problist={self._map[nlist]})"
            )

    def clean(self) -> None:
        """(Re-)Compute probality for sums.
           It is far slower than just updating individual probabilities,
           but the fast methods used in update function leads to accumulate small rounding
           errors.
           This functions is intended to be called regularly for cleaning
           these rounding errors."""
        # old_problist = self._problist
        # old_probtot = self.probtot
        self._problist = array([data.sum() for data in self._map])
        self.probtot = self._problist.sum()

    def get_probaobj(self, obj: Any) -> Probaobj:
        return Probaobj(obj)

    @staticmethod
    def seed(seed: int) -> None:
        if seed == 0:
            random.seed()
        else:
            random.seed(seed)
