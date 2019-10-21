from multiprocessing import get_context, cpu_count
from itertools import repeat
from time import process_time, time
from psutil import Process
from os import getpid

import numpy as N
import pandas as P


class Finished(Exception):
    def __init__(self, msgnum, detail=""):
        self.num = msgnum
        self.detail = detail

    @property
    def message(self):
        if self.num == 0:
            return "Time is up!"
        elif self.num == 1:
            return "No more reactions can be processed"
        elif self.num == 2:
            return "No reaction could be find"
        elif self.num == 3:
            return "Rounding problem, negative proba (likely)!"
        elif self.num == 4:
            return "Wanted to decrement unpopulated species " + self.detail
        else:
            return "Unknown error!"

    def __str__(self):
        return "End ({num}): {mess}".format(num=self.num, mess=self.message)


# May be simply subtituted by a set, with correct hash for objects?
class ObjCollect(dict):
    """Dictionnary tracking objects, removing duplicates.
       Object are stored as str(obj):obj
       Objects must thus be differentiated by their conversion to string as a hash
       Objects should only be added by .add method (no direct dictionnary assignment)
       Objects cannot be overriden. if an object with the same str hash is added, the old
       one will be returned (Use 'add' method for this functionality, no direct assignment).
       Modifications implies to first delete the entry
       Note: direct dict assignment is not forbiden,
             as it is needed by multiprocessing (need to copy object)
    """

    def add(self, new):
        name = str(new)
        if name not in self:  # is new object
            super().__setitem__(name, new)
            return new
        return self[name]

    def remove(self, old):
        if old in self:
            super().__delitem__(old)


class Result:
    def __init__(self, data):
        res = N.array(data)
        self.data = {
            "table": res[:, 0],
            "lendist": res[:, 1],
            "pooldist": res[:, 2],
            "end": res[:, 3],
        }

    def _format(self, name, field, num):
        if num is None or num == "m":
            res = self.data[name].mean()
        elif num == "s":
            res = self.data[name].std()
        elif num == "+":
            res = self.data[name].mean() + self.data[name].std()
        elif num == "-":
            res = self.data[name].mean() - self.data[name].std()
        else:
            res = self.data[name][num]
        if field is not None:
            res = res.loc[field]
        return res

    def table(self, field=None, num=None):
        return self._format("table", field, num)

    def lendist(self, field=None, num=None):
        return self._format("lendist", field, num)

    def pooldist(self, field=None, num=None):
        return self._format("pooldist", field, num)

    def end(self, num=None):
        if num is None:
            return self.data["end"]
        else:
            return self.data["end"][num]


class Probalist:
    def __init__(self, minprob=1e-10):
        self._minprob = minprob
        self._map = [N.array([])]
        self._mapobj = [N.array([])]
        self._nblist = 1
        self._actlist = 0
        self._problist = N.array([0.0])
        self.probtot = 0.0
        self._queue = []

    def register(self, newobj, proba):
        if self._queue:
            rlist, rpos = self._addfromqueue(newobj, proba)
        else:
            rlist, rpos = self._addfrommap(newobj, proba)
        self._updateprob(rlist, proba)
        return rlist, rpos

    def _updateprob(self, nlist, delta):
        self._problist[nlist] += delta
        if abs(self._problist[nlist]) < self._minprob:
            self.probtot -= self._problist[nlist]
            self._problist[nlist] = 0
        self.probtot += delta

    def unregister(self, nlist, npos):
        if nlist is not None:
            oldproba = self._map[nlist][npos]
            self._map[nlist][npos] = 0
            self._mapobj[nlist][npos] = None
            self._queue.append((nlist, npos))
            self._updateprob(nlist, -oldproba)

    def update(self, nlist, npos, proba):
        delta = proba - self._map[nlist][npos]
        self._map[nlist][npos] = proba
        self._updateprob(nlist, delta)

    def getproba(self, nlist, npos):
        return self._map[nlist][npos]

    def _addfrommap(self, newobj, proba):
        rpos = len(self._map[self._actlist])
        if rpos <= self._nblist:
            rlist = self._actlist
            self._map[rlist] = N.append(self._map[rlist], proba)
            self._mapobj[rlist] = N.append(self._mapobj[rlist], newobj)
            return rlist, rpos
        if len(self._map[0]) <= self._nblist:
            self._actlist = 0
        else:
            self._actlist += 1
            if self._actlist >= self._nblist:
                self._map.append(N.array([]))
                self._mapobj.append(N.array([]))
                self._problist = N.append(self._problist, 0.0)
                self._nblist += 1
        return self._addfrommap(newobj, proba)

    def _addfromqueue(self, newobj, proba):
        rlist, rpos = self._queue.pop(0)
        self._map[rlist][rpos] = proba
        self._mapobj[rlist][rpos] = newobj
        return rlist, rpos

    def choose(self):
        try:
            nlist = N.random.choice(self._nblist, p=self._problist / self.probtot)
            return N.random.choice(
                self._mapobj[nlist], p=self._map[nlist] / self._problist[nlist]
            )
        except ValueError:
            raise Finished(3)

    def clean(self, report=False):
        """Compute probality for sums. It is slower than just updating individual probabilities,
           but the fast methods used in update function leads to accumulate small rounding
           errors.
           This functions is intended to be called regularly for cleaning 
           these rounding errors."""
        if report:
            old_problist = self._problist
            old_probtot = self.probtot
        self._problist = N.array([data.sum() for data in self._map])
        self.probtot = self._problist.sum()
        if report:
            return self.probtot - old_probtot, self._problist - old_problist


class System:
    def __init__(
        self,
        init: dict,
        kpol: float = 1.0,
        khyd: float = 1.0,
        act: float = 1.0,
        kcat: float = 1000.0,
        conc: float = 0.1,
        tend: float = 1.0,
        tstep: float = 0.01,
        save: list = None,
        minprob: float = 1e-10,
        dropreac: bool = True,
        autoclean: bool = True,
    ):
        self.runtime_init()
        #self.proc_init()
        self.time = 0.0
        self.kpol = kpol
        self.khyd = khyd
        self.kcat = kcat
        self.act = act
        self.conc = conc
        self.tend = tend
        self.tstep = tstep
        self.save = save
        self.dropreac = dropreac
        self.autoclean = autoclean
        self.step = 0
        self._ptot = sum([pop * len(comp) for comp, pop in init.items()])
        self.vol = self._ptot / self.conc
        self.pool_collect = ObjCollect()
        self.comp_collect = ObjCollect()
        self.reac_collect = ObjCollect()
        self.probalist = Probalist(minprob=minprob)
        for compound, pop in init.items():
            self.get_compound(compound).init_pop(pop)
        self.set_run(tend=1.0, tstep=0.1, save=[], maxsteps=10000)

    def runtime_init(self):
        self._ptime0 = process_time()

    @property
    def runtime(self):
        return process_time() - self._ptime0

    #def proc_init(self):
        #pass
        #self._proc = Process(getpid()) ###Broken, cannot pickle anymore !!

    @property
    def memuse(self):
        return Process(getpid()).memory_info().rss / 1024 / 1024
        #return self._proc.memory_info().rss / 1024 / 1024   ###Broken, cannot Process  pickle anymore !!

    def concentration(self, compound: str):
        try:
            return self.comp_collect[compound].pop / self.vol
        except KeyError:
            return 0.0

    @property
    def poplist(self):
        return {comp: comp.pop for comp in self.comp_collect.values()}

    @property
    def lendist(self):
        res = {}
        for comp in self.comp_collect.values():
            length = comp.length
            pop = comp.pop
            if length not in res:
                res[length] = 0
            res[length] += pop
        return res

    @property
    def pooldist(self):
        res = {}
        for comp in self.pool_collect.values():
            length = comp.length
            if length not in res:
                res[length] = 0
            res[length] += 1
        return res

    def statlist(self):
        stat = {}
        dist = {}
        dist["lendist"] = self.lendist
        dist["pooldist"] = self.pooldist
        stat["nbcomp"] = len(self.comp_collect)
        stat["nbreac"] = len(self.reac_collect)
        stat["poolsize"] = len(self.pool_collect)
        stat["maxlength"] = max(dist["lendist"])
        return stat, dist

    def process(self, tstop):
        # Check if a cleanup should be done
        if self.autoclean:
            self.probalist.clean()
        # Check if end of time is nigh
        if self.time >= self.tend:
            raise Finished(0)
        # Then process self.maxsteps times
        for _ in repeat(None, self.maxsteps):
            # ... but stop is step ends
            if self.time >= tstop:
                break
            # choose a random event
            chosen = self.probalist.choose()
            # check if there even was an event to choose
            if chosen is None:
                raise Finished(2)
            # perform the (chosen one) event
            chosen.process()
            if self.probalist.probtot == 0:
                raise Finished(1)
            # update time for next step
            self.time += N.log(1 / N.random.rand()) / self.probalist.probtot
            self.step += 1

    def run(self, num: int = 0):
        lines = (
            ["thread", "ptime", "memuse", "step", "time"]
            + self.save
            + ["maxlength", "nbcomp", "poolsize", "nbreac"]
        )
        table = P.DataFrame(index=lines)
        lendist = P.DataFrame()
        pooldist = P.DataFrame()
        N.random.seed()  # necessary for multiprocessing from different seeds
        self.time = 0.0
        tnext = 0.0
        step = 0
        self.runtime_init()
        #self.proc_init()
        Process(getpid()).cpu_affinity([num % cpu_count()])   ###Broken, cannot Process  pickle anymore !!
        #self._proc.cpu_affinity([num % cpu_count()])   ###Broken, cannot Process  pickle anymore !!
        while True:
            try:
                self.process(tnext)
                stat, dist = self.statlist()
                table[step] = (
                    [num, self.runtime, self.memuse, self.step, self.time]
                    + [self.concentration(comp) for comp in self.save]
                    + [
                        stat["maxlength"],
                        stat["nbcomp"],
                        stat["poolsize"],
                        stat["nbreac"],
                    ]
                )
                lendist = lendist.join(
                    P.DataFrame.from_dict({step: dist["lendist"]}), how="outer"
                ).fillna(0)
                pooldist = pooldist.join(
                    P.DataFrame.from_dict({step: dist["pooldist"]}), how="outer"
                ).fillna(0)
                tnext += self.tstep
                step += 1
            except Finished as the_end:
                end = "{mess} ({time} s)".format(mess=the_end, time=self.runtime)
                break
        return (table, lendist.astype(int), pooldist.astype(int), end)

    def multirun(self, nbthread: int = None):
        ctx = get_context("fork")
        if nbthread is None:
            nbthread = ctx.cpu_count()
        with ctx.Pool(nbthread) as pool:
            res = pool.map(self.run, range(nbthread)) ### Warning! self.run cannot be pickeled anymore...
        return Result(res)

    def set_run(
        self,
        tend: float = None,
        tstep: float = None,
        save: list = None,
        maxsteps: int = None,
    ):
        if tend is not None:
            self.tend = tend
        if tstep is not None:
            self.tstep = tstep
        if save is not None:
            self.save = save
        if maxsteps is not None:
            self.maxsteps = maxsteps

    def addkept(self, reac: str):
        """reactions that are kept"""
        self.get_reaction(reac).addkept()

    def get_compound(self, name: str):
        try:
            return self.pool_collect[name]
        except KeyError:
            return self.pool_collect.add(Polymer(name, self))

    def get_reaction(self, name: str):
        try:
            reac = self.reac_collect[name]
        except KeyError:
            return self.reac_collect.add(newreaction(name, self))
        reac.register()
        return reac



class Reaction:
    def __init__(self, description: str, reactants: str, catal:str, pos:str, system: System):
        self._system = system
        self._reactantnames = reactants.split("+")
        self.name = description
        self._pos = int(pos) if pos else ""
        self._create()
        self._started = False
        self._catinreac = False
        if catal:
            catal = system.get_compound(catal)
            if catal in self._reactants:
                self._catinreac = True
            self.const *= self._system.kcat
        self._catal = catal
        self._nlist, self._npos = None, None
        self.register()

    def _create(self): # Must be implemented in children
        pass

    def destroy(self):
        self._system.probalist.unregister(self._nlist, self._npos)
        self._nlist, self._npos = None, None
        for comp in self._reactants:
            comp.unregister_reaction(self)
        if self._catal:
            self._catal.unregister_reaction(self)
        if self._system.dropreac:
            self._system.reac_collect.remove(self.name)

    def register(self):
        for comp in self._reactants:
            comp.register_reaction(self)
        if self._catal:
            self._catal.register_reaction(self)

    def process(self):
        for reac in self._reactants:
            try:
                reac.dec()
            except Finished as f:
                detail = " from {} (p={}, ".format(self.name, self.proba)
                for comp in self._reactants:
                    detail += "[{}]={} ,".format(comp.name, comp.pop)
                if self._catal:
                    detail += "cat[{}]={})".format(self._catal.name, self._catal.pop)
                else:
                    detail += ")"
                raise Finished(4, f.detail + detail)
        if not self._started:
            self._products = [
                self._system.get_compound(name) for name in self._products
            ]
            self._started = True
        for prod in self._products:
            prod.inc()

    def __repr__(self):
        return self.name

    def calcproba(self):
        pop0 = self._reactants[0].pop
        if pop0 == 0:
            return 0
        proba = self._uncatal(pop0)
        if self._catal:
            if self._catal in self._reactants:
                if self._dimer:
                    proba *= (self._catal.pop - 2) / self._system.vol
                else:
                    proba *= (self._catal.pop - 1) / self._system.vol
            else:
                proba *= self._catal.pop / self._system.vol
        return proba

    def _uncatal(self):
        pass # must be implemented in childrens
    
    @property
    def proba(self):
        if self._nlist is None:
            return 0.0
        return self._system.probalist.getproba(self._nlist, self._npos)

    def update(self):
        proba = self.calcproba()
        if self._nlist is None:
            self._nlist, self._npos = self._system.probalist.register(self, proba)
        else:
            self._system.probalist.update(self._nlist, self._npos, proba)

    def addkept(self):
        if self._catal:
            self._catal.addkept(self)
        for comp in self._reactants:
            comp.addkept(self)


class Polymerization(Reaction):
    def _create(self):
        try:
            self._products = [self._reactantnames[0] + self._reactantnames[1]]
            self._reactants = [
                self._system.get_compound(self._reactantnames[0]),
                self._system.get_compound(self._reactantnames[1]),
            ]
            self._dimer = self._reactants[0] is self._reactants[1]
            self.const = self._system.kpol
            if len(self._reactantnames[0]) == 1:
                self.const *= self._system.act
        except IndexError:
            raise ValueError("Bad name for polymerization reaction")

    def _uncatal(self, pop0:int):
        if self._dimer:
            return self.const * pop0 * (pop0 - 1) / self._system.vol / 2.0
        else:
            return self.const * pop0 * self._reactants[1].pop / self._system.vol

        
class Hydrolysis(Reaction):
    def _create(self):
        try:
            self._dimer = False
            name = self._reactantnames[0]
            self._reactants = [self._system.get_compound(name)]
            if self._pos < 1 or self._pos >= len(name):
                raise ValueError(
                    "Hydrolysis cut position should be at least 1, at most the chain length minus one"
                )
            self._products = [name[: self._pos], name[self._pos :]]
            self.const = self._system.khyd
        except IndexError:
            raise ValueError("Bad name for hydrolysis reaction")

    def _uncatal(self, pop0:int):
        return self.const * pop0

    
def newreaction(description: str, system: System):
        kind, reactants, catal, pos = description.split(".")
        if kind == "P":  # Polymerisation
            actual = Polymerization
        elif kind == "H":  # Hydrolysis
            actual = Hydrolysis
        else:
            raise ValueError("Unknown reaction kind")
        return actual(description, reactants, catal, pos, system)


class Polymer:
    def __init__(self, name: str, system: System):
        self._reactions = set()
        self._system = system
        self._keptname = set()
        self.pop = 0
        self.name = name
        self.length = len(name)

    def __repr__(self):
        return self.name

    def inc(self):
        self.pop += 1
        if self.pop == 1:
            self._generate()
        self._upd_reac()

    def dec(self):
        if self.pop > 0:
            self.pop -= 1
        else:
            raise Finished(4, self.name)
        if self.pop == 0:
            self._destroy()
        else:
            self._upd_reac()

    def init_pop(self, start: int):
        if start != 0:
            if self.pop == 0:
                self.pop = start
                self._generate()
            else:
                self.pop = start
            self._upd_reac()
        else:
            self.pop = 0
            self._destroy()

    def _generate(self):
        self._system.comp_collect.add(self)
        for reac in self._reacnames():
            self._system.get_reaction(reac)

    def _upd_reac(self):
        for reac in list(self._reactions):
            reac.update()

    def _destroy(self):
        for reac in list(self._reactions):
            reac.destroy()
        self._system.comp_collect.remove(self.name)

    def _build_hydro_name(self):
        return {"H." + self.name + ".." + str(i) for i in range(1, self.length)}

    def _build_pol_name(self):
        res = set()
        for other in self._system.comp_collect:
            res.add("P." + self.name + "+" + other + "..")
            res.add("P." + other + "+" + self.name + "..")
        return res

    def _reacnames(self):
        return self._build_pol_name() | self._build_hydro_name() | self._keptname

    def addkept(self, reaction: Reaction):
        self._keptname.add(reaction.name)

    def register_reaction(self, reaction: Reaction):
        self._reactions.add(reaction)

    def unregister_reaction(self, reaction: Reaction):
        try:
            self._reactions.remove(reaction)
        except KeyError:
            pass


def multirun(compute, nbthread: int = None):
    ctx = get_context("fork")
    if nbthread is None:
        nbthread = ctx.cpu_count()
    with ctx.Pool(nbthread) as pool:
        res = pool.map(compute, range(nbthread))
    return Result(res)
