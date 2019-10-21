from multiprocessing import get_context
from itertools import repeat
import numpy as N


class Finished(Exception):
    pass


class ObjCollect(dict):
    """Dictionnary tracking objects, removing duplicates.
       Object are stored as str(obj):obj
       Objects must thus be differentiated by their conversion to string as a hash
       Objects should only be added by .add method (no direct dictionnary assignment)
       Objects cannot be overriden. if an object with the same str hash is added, the old
       one will be returned (Use 'add' method for this functionality, no direct assignment).
       Modifications implies to first delete the entry
       Note: direct dict assignment is not forbiden, as it is needed by multiprocessing (need to copy object)
    """

    def add(self, new):
        name = str(new)
        if name not in self:  # is new object
            super().__setitem__(name, new)
            return new
        return self[name]

    def remove(self, old):
        if old not in self:
            raise ValueError("{key} not in collection".format(key=old))
        super().__delitem__(old)


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
        if len(self._queue) > 0:
            rlist, rpos = self._addfromqueue(newobj, proba)
        else:
            rlist, rpos = self._addfrommap(newobj, proba)
        self._updateprob(rlist, proba)
        return rlist, rpos

    def _updateprob(self, nlist, delta):
        self._problist[nlist] += delta
        # if self._problist[nlist] <0:
        #    self.probtot -= self._problist[nlist]
        #    self._problist[nlist] = 0
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

    def update(self, nlist, npos, delta):
        self._map[nlist][npos] += delta
        self._updateprob(nlist, delta)

    def getproba(self, nlist, npos):
        return self._map[nlist][npos]

    def _addfrommap(self, newobj, proba):
        # Possible to optimize by avoiding N.append???
        rpos = len(self._map[self._actlist])
        if rpos <= self._nblist:
            rlist = self._actlist
            self._map[rlist] = N.append(self._map[rlist], proba)
            self._mapobj[rlist] = N.append(self._mapobj[rlist], newobj)
            return rlist, rpos
        else:
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
            raise Finished("Oups, rounding problem, negative proba!")


class System:
    def __init__(
        self,
        init,
        kpol=1.0,
        khyd=1.0,
        act=1.0,
        kcat=1000.0,
        conc=0.1,
        tend=1.0,
        tstep=0.01,
        save=None,
        minprob=1e-10,
    ):
        self.time = 0.0
        self.kpol = kpol
        self.khyd = khyd
        self.kcat = kcat
        self.act = act
        self.conc = conc
        self.tend = tend
        self.tstep = tstep
        self.save = save
        self._ptot = sum([pop * len(comp) for comp, pop in init.items()])
        self.vol = self._ptot / self.conc
        self.pool = ObjCollect()
        self.compounds = ObjCollect()
        self.reactions = ObjCollect()
        self.created = ObjCollect()
        self.probalist = Probalist(minprob=minprob)
        for compound, pop in init.items():
            self.get_compound(compound).init_pop(pop)

    def concentration(self, compound):
        try:
            return self.compounds[compound].pop / self.vol
        except KeyError:
            return 0.0

    @property
    def poplist(self):
        return {comp: comp.pop for comp in self.compounds.values()}

    @property
    def lendist(self):
        res = {}
        for comp in self.compounds.values():
            length = comp.length
            pop = comp.pop
            if length not in res:
                res[length] = 0
            res[length] += pop
        return res

    @property
    def pooldist(self):
        res = {}
        for comp in self.pool.values():
            length = comp.length
            if length not in res:
                res[length] = 0.0
            res[length] += 1.0
        return res

    @property
    def createddist(self):
        res = {}
        for comp in self.created.values():
            length = comp.length
            if length not in res:
                res[length] = 0.0
            res[length] += 1.0
        return res

    @property
    def statlist(self):
        res = {}
        res["lendist"] = self.lendist
        res["nbcomp"] = len(self.compounds)
        res["nbreac"] = len(self.reactions)
        res["poolsize"] = len(self.pool)
        res["pooldist"] = self.pooldist
        res["createddist"] = self.createddist
        res["createdsize"] = len(self.created)
        res["maxlength"] = max(res["lendist"])
        return res

    def process(self, tend, maxsteps=1000):
        for _ in repeat(None, maxsteps):
            try:
                chosen = self.probalist.choose()
                if chosen is None:
                    raise Finished
                chosen.process()
                if self.probalist.probtot == 0:
                    raise Finished
                self.time += N.log(1 / N.random.rand()) / self.probalist.probtot
                if self.time >= tend:
                    break
            except Finished:
                raise Finished("No more reactions can be processed")

    def run(self, num=0):
        N.random.seed()  # necessary for multiprocessing from different seeds
        self.time = 0.0
        tend = self.tend
        tstep = self.tstep
        save = self.save
        nbstep = int(tend / tstep) + 1
        res = N.empty((nbstep, len(save) + 2))
        res.fill(N.nan)
        res[:, 0] = num
        tnext = 0.0
        for step in range(nbstep):
            try:
                self.process(tnext)
                res[step, 1] = self.time
                for col, comp in enumerate(save):
                    res[step, col + 2] = self.concentration(comp)
                tnext += tstep
            except Finished:
                break
        return res

    def multirun(self, tend, tstep, save, nbthread=None):
        self.tend = tend
        self.tstep = tstep
        self.save = save
        ctx = get_context("fork")
        if nbthread is None:
            nbthread = ctx.cpu_count()
        with ctx.Pool(nbthread) as p:
            res = p.map(self.run, range(nbthread))
        res = N.array(res)
        return N.nanmean(res, axis=0), N.nanstd(res, axis=0), res

    def addkept(self, reac):
        """reactions that are kept"""
        reaction = self.get_reaction(reac)
        if reaction._catal != "":
            reaction._catal.addkept(reaction)
        for comp in reaction._reactants:
            comp.addkept(reaction)

    def get_compound(self, name):
        try:
            return self.pool[name]
        except KeyError:
            return self.pool.add(Polymer(name, self))

    def get_reaction(self, name):
        try:
            return self.reactions[name]
        except KeyError:
            return self.reactions.add(Reaction(name, self))


class Polymer:
    def __init__(self, reactants, system):
        self._system = system
        self._kept = set()
        self._keptname = set()
        self.pop = 0
        self.name = reactants
        self.length = len(reactants)

    def __repr__(self):
        return self.name

    def inc(self):
        self.pop += 1
        if self.pop == 1:
            self._generate()
        self._upd_reac(+1)

    def dec(self):
        if self.pop > 0:
            self.pop -= 1
        else:
            raise ValueError("Cannot decrement unpopulated species !")
        if self.pop == 0:
            self._destroy()
        else:
            self._upd_reac(-1)

    def init_pop(self, start):
        if start != 0:
            if self.pop == 0:
                self._generate()
            if start > self.pop:
                for _ in repeat(
                    None, start - self.pop
                ):  # Inefficient, but rarely called...
                    self.inc()
            else:
                for _ in repeat(
                    None, self.pop - start
                ):  # Inefficient, but rarely called...
                    self.dec()
        else:
            self._destroy()

    def _generate(self):
        self._system.compounds.add(self)
        self._system.created.add(self)
        for reac in self._reacnames():
            self._system.get_reaction(reac)

    def _upd_reac(self, change):
        for reac in self._reacnames():
            self._system.reactions[reac].update(self, change)

    def _destroy(self):
        for reac in self._reacnames():
            self._system.reactions[reac].destroy()
        self._system.compounds.remove(self.name)

    def _build_hydro_name(self):
        return {"H." + self.name + ".." + str(i) for i in range(1, self.length)}

    def _build_pol_name(self):
        res = set()
        for other in self._system.compounds:
            res.add("P." + self.name + "+" + other + "..")
            res.add("P." + other + "+" + self.name + "..")
        return res

    def _reacnames(self):
        return self._build_pol_name() | self._build_hydro_name() | self._keptname

    def addkept(self, reaction):
        self._kept.add(reaction)
        self._keptname.add(reaction.name)


class Reaction:
    def __init__(self, description, system):
        self._system = system
        kind, reactants, catal, pos = description.split(".")
        self._reactants = [system.get_compound(comp) for comp in reactants.split("+")]
        catal = system.get_compound(catal) if catal else ""
        self.name = description
        self._pos = int(pos) if pos else ""
        if kind[0] == "P":  # Polymerisation
            self._create_polym()
        elif kind[0] == "H":  # Hydrolysis
            self._create_hydrol()
        else:
            raise ValueError("Unknown reaction kind")
        self._catinreac = False
        if catal != "":
            if catal in self._reactants:
                self._catinreac = True
            self.const *= self._system.kcat
        self._catal = catal
        self._nlist, self._npos = None, None

    def _create_polym(self):
        if len(self._reactants) != 2:
            raise ValueError("Two reactants needed for polymerization")
        self._kind = "P"
        self._products = [
            self._system.get_compound(self._reactants[0].name + self._reactants[1].name)
        ]
        if self._reactants[0] is self._reactants[1]:
            self._dimer = True
        else:
            self._dimer = False
        self.const = self._system.kpol
        if self._reactants[0].length == 1:
            self.const *= self._system.act

    def _create_hydrol(self):
        if len(self._reactants) != 1:
            raise ValueError("One reactant for hydrolysis")
        self._kind = "H"
        self._dimer = False
        name = self._reactants[0].name
        if self._pos < 1 or self._pos >= len(name):
            raise ValueError(
                "Hydrolysis cut position should be at least 1, at most the chain length minus one"
            )
        self._products = [
            self._system.get_compound(name[: self._pos]),
            self._system.get_compound(name[self._pos :]),
        ]
        self.const = self._system.khyd

    def destroy(self):
        self._system.probalist.unregister(self._nlist, self._npos)
        self._nlist, self._npos = None, None

    def process(self):
        for reac in self._reactants:
            reac.dec()
        for prod in self._products:
            prod.inc()

    def __repr__(self):
        return self.name

    def calcproba(self):
        c0 = self._reactants[0].pop
        if c0 == 0:
            return False
        if self._kind == "H":
            proba = self.const * c0
        elif self._dimer:
            proba = self.const * c0 * (c0 - 1) / self._system.vol / 2.0
        else:
            proba = self.const * c0 * self._reactants[1].pop / self._system.vol
        if self._catal != "":
            if self._catal in self._reactants:
                if self._dimer:
                    proba *= (self._catal.pop - 2) / self._system.vol
                else:
                    proba *= (self._catal.pop - 1) / self._system.vol
            else:
                proba *= self._catal.pop / self._system.vol
        return proba

    @property
    def proba(self):
        if self._nlist is None:
            return 0.0
        else:
            return self._system.probalist.getproba(self._nlist, self._npos)

    def update(self, compound, change):
        proba = self.calcproba()
        if proba is False:
            self.destroy()
        else:
            if self._nlist is None:
                self._nlist, self._npos = self._system.probalist.register(self, proba)
            else:
                self._system.probalist.update(
                    self._nlist, self._npos, proba - self.proba
                )
