from typing import Any, Tuple, Optional, Deque
from collections import deque
from numpy import array, append, log, random

from metadynamic.ends import RoundError
from metadynamic.logger import Logged


class Probalistic:
    probalist: "Probalist"

    @classmethod
    def setprobalist(cls, vol: float = 1, minprob: float = 1e-10) -> None:
        cls.probalist = Probalist(vol, minprob)


class Probaobj(Logged, Probalistic):
    """Probalistic object, that can be used in a Probalist object
    for Gillespie-like computation.

    Each object to be stored in the probalist must contain one Probaobj object
    """

    def __init__(self, obj: Any):
        self.nlist: Optional[int]
        self.npos: Optional[int]
        self.obj = obj
        self.unset_proba_pos()

    def set_proba_pos(self, nlist: int, npos: int) -> None:
        self.nlist = nlist
        self.npos = npos
        self.registered = True

    def unset_proba_pos(self) -> None:
        self.nlist = None
        self.npos = None
        self.registered = False

    @property
    def proba_pos(self) -> Tuple[int, int]:
        if self.registered:
            assert self.nlist is not None and self.npos is not None
            return self.nlist, self.npos
        else:
            raise ValueError("Unregistered")

    def register(self, proba: float = 0) -> None:
        self.probalist.register(self, proba)

    def update(self, proba: float) -> None:
        self.probalist.update(self, proba)

    def unregister(self) -> None:
        self.probalist.unregister(self)

    @property
    def proba(self) -> float:
        return self.probalist.getproba(self)


class Probalist:
    def __init__(self, vol: float = 1, minprob: float = 1e-10):
        self.vol = vol
        self._minprob = minprob
        self._map = [array([])]
        self._mapobj = [array([])]
        self.npblist = 1
        self._actlist = 0
        self._problist = array([0.0])
        self.probtot = 0.0
        self._queue: Deque[Tuple[int, int]] = deque()

    def _updateprob(self, nlist: int, delta: float) -> None:
        newproba = self._problist[nlist] + delta
        # reset to 0 probabilities that are too low
        # for rounding errors correction
        if newproba < self._minprob:
            delta += newproba
            newproba = 0.0
        self._problist[nlist] = newproba
        self.probtot += delta

    def register(self, probaobj: Probaobj, proba: float = 0.0) -> None:
        if self._queue:
            probaobj.set_proba_pos(*self._addfromqueue(probaobj, proba))
        else:
            probaobj.set_proba_pos(*self._addfrommap(probaobj, proba))
        assert probaobj.nlist is not None
        self._updateprob(probaobj.nlist, proba)

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
        delta = proba - self._map[nlist][npos]
        self._map[nlist][npos] = proba
        self._updateprob(nlist, delta)

    def getproba(self, probaobj: Probaobj) -> float:
        if probaobj.registered:
            nlist, npos = probaobj.proba_pos
            return self._map[nlist][npos]
        return 0.0

    def _addfrommap(self, newobj: Probaobj, proba: float) -> Tuple[int, int]:
        rpos = len(self._map[self._actlist])
        if rpos <= self.npblist:
            rlist = self._actlist
            self._map[rlist] = append(self._map[rlist], proba)
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
        return self._addfrommap(newobj, proba)

    def _addfromqueue(self, newobj: Probaobj, proba: float) -> Tuple[int, int]:
        rlist, rpos = self._queue.popleft()
        self._map[rlist][rpos] = proba
        self._mapobj[rlist][rpos] = newobj.obj
        return rlist, rpos

    def choose(self) -> Tuple[Any, float]:
        # First choose a random line in the probability map
        try:
            nlist = random.choice(self.npblist, p=self._problist / self.probtot)
        except ValueError as v:
            raise RoundError(
                f"(reason: {v}; probtot={self.probtot}=?={self._problist.sum()}; problist={self._problist})"
            )
        # Then choose a random column in the chosen probability line
        try:
            return (
                random.choice(
                    self._mapobj[nlist], p=self._map[nlist] / self._problist[nlist]
                ),
                log(1 / random.rand()) / self.probtot,
            )
        except ValueError as v:
            raise RoundError(
                f"(reason: {v}; probtot={self._problist[nlist]}=?={self._map[nlist].sum()}; problist={self._map[nlist]})"
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

    def get_probaobj(self, obj: Any):
        return Probaobj(obj)

    @staticmethod
    def seed(seed: Optional[int]) -> None:
        random.seed(seed)
