from typing import Callable, List
from itertools import product

# from functools import partial, cached_property from 3.8 only
from functools import reduce


from metadynamic.collector import Collected


def ismono(name: str) -> bool:
    return len(name) == 1


def ispolym(name: str) -> bool:
    return name.isalpha()


def isact(name: str) -> bool:
    return name[-1] == "*" and name[:-1].isalpha()


class Descriptor:
    def __init__(self):
        self.cat_dict = {}

    @property
    def catlist(self):
        return self.cat_dict.keys()

    def categories(self, name) -> List[str]:
        return [catname for catname, rule in self.cat_dict.items() if rule(name)]

    def __repr__(self) -> str:
        return f"Descriptor: {self.cat_dict.keys()}"

    def add_cat(self, catname: str, rule: Callable[[str], bool]) -> None:
        self.cat_dict[catname] = rule


chemdescriptor = Descriptor()
chemdescriptor.add_cat("mono", ismono)
chemdescriptor.add_cat("polym", ispolym)
chemdescriptor.add_cat("actpol", isact)


class ReacDescr:
    pass


class Ruleset(Collected):
    def __init__(self, descriptor: Descriptor):
        self.descriptor = descriptor
        self.categories = {}
        for catname in descriptor.catlist:
            self.categories[catname] = set()
        self.reactants = {}
        self.prodbuilder = {}
        self.constbuilder = {}
        self.rules = set()

    def add_rule(
        self,
        rulename: str,
        reactants: List[str],
        prodbuilder: Callable[[List[str]], str],
        constbuilder: Callable[List[str], List[float]],
    ):
        self.rules.add(rulename)
        self.reactants[rulename] = reactants
        self.prodbuilder[rulename] = prodbuilder
        self.constbuilder[rulename] = constbuilder
        for reac in reactants:
            try:
                self.categories[reac].add(rulename)
            except KeyError:
                raise ValueError(f"Unrecognize category {reac}")

    def get_related(self, name: str):
        # Maybe memoize the list of rule for a given list of categories...
        # rule_related = reduce(
        #    lambda x, y: x | y,
        #    [self.categories[catname] for catname in self.descriptor.categories(name)],
        # )
        # or simply use self.rules??? Far less overhead, maybe not much perf loss.
        categories = self.descriptor.categories(name)
        res = set()
        for reaclist in self.reactants:
            for pos, reacname in enumerate(reaclist):
                if reacname in categories:
                    reactantnames = product(
                        *[
                            [name]
                            if pos2 == pos
                            else self._comp_collect.categories[catname]
                            for pos2, catname in enumerate(reaclist)
                        ]
                    )
                    res.add(
                        ReacDescr(
                            reacname,
                            reactantnames,
                            self.prodbuilder[reacname],
                            self.constbuilder[reacname],
                        )
                    )


def joiner(names, sep):
    return sep.join(names)


# joiner_direct = partial(joiner, sep="") from python 3.8 only
def joiner_direct(names):
    return "".join(names)


def samecase(one: str, two: str) -> bool:
    return (one.islower() and two.islower()) or (one.isupper() and two.isupper())


def kpol(names, k) -> float:
    return k[0] if samecase(names[0][-1], names[1][0]) else k[1]


ruleset = Ruleset()
ruleset.add_rule("Polym", ["polym", "polym"], joiner_direct, kpol)
