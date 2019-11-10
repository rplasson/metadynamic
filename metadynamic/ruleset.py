from typing import Callable, List, Any
from itertools import product

# from functools import partial, cached_property from 3.8 only
# from functools import reduce


from metadynamic.collector import Collected


class Descriptor:
    def __init__(self):
        self.cat_dict = {}
        self.prop_dict = {}

    @property
    def catlist(self):
        return self.cat_dict.keys()

    @property
    def prop(self, propname):
        return self.prop_dict[propname]()

    def categories(self, name) -> List[str]:
        return [catname for catname, rule in self.cat_dict.items() if rule(name)]

    def __repr__(self) -> str:
        return f"Descriptor: {self.cat_dict.keys()}"

    def add_cat(self, catname: str, rule: Callable[[str], bool]) -> None:
        self.cat_dict[catname] = rule

    def add_prop(self, propname: str, func: Callable[[None], Any]) -> None:
        self.prop_dict[propname] = func


class ReacDescr:
    def __init__(self, rule, reactantnames, prodbuilder, constbuilder):
        self.rule = rule
        self.reactantnames = reactantnames
        self._prodbuilder = prodbuilder
        self._constbuilder = constbuilder

    def builprod(self):
        return self._prodbuilder(self.reactantnames)

    # Which best place for storing k???
    def buildconst(self, k):
        return self._constbuilder(k)

    @property
    def name(self):
        return f"{self.rule}:{'+'.join(self.reactantnames)}"


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
        # Will look fot the list of reactions for each rule
        result = []
        for rule in self.rules:
            res = set()
            # get the list of reactant type for the rule
            for reaclist in self.reactants[rule]:
                # Then scan all possible combinations, with fixing name in each possible pos
                for pos, reacname in enumerate(reaclist):
                    if reacname in categories:
                        res |= set(
                            product(
                                *[
                                    [name]
                                    if pos2 == pos
                                    # expect comp_collect.categories to return str
                                    # collector will have to be adapted
                                    else self.comp_collect.categories[catname]
                                    for pos2, catname in enumerate(reaclist)
                                ]
                            )
                        )
            # Then create all possible reaction descriptions
            result += [
                ReacDescr(
                    rule, reactantnames, self.prodbuilder[rule], self.constbuilder[rule]
                )
                for reactantnames in res
            ]
        # Finally return all build reac decription for each rule
        return result


# Functions for descriptor


def ismono(name: str) -> bool:
    return len(name) == 1


def ispolym(name: str) -> bool:
    return name.isalpha()


def isact(name: str) -> bool:
    return name[-1] == "*" and name[:-1].isalpha()


def length(name: str) -> int:
    if ismono(name):
        return 1
    if ispolym(name):
        return len(name)
    if isact(name):
        return len(name) - 1


# Functions for rulesets


def joiner(names, sep):
    return sep.join(names)


# joiner_direct = partial(joiner, sep="") from python 3.8 only
def joiner_direct(names):
    return "".join(names)


def samecase(one: str, two: str) -> bool:
    return (one.islower() and two.islower()) or (one.isupper() and two.isupper())


def kpol(names, k) -> float:
    return k[0] if samecase(names[0][-1], names[1][0]) else k[1]


# Define a specific ruleset

chemdescriptor = Descriptor()
chemdescriptor.add_cat("mono", ismono)
chemdescriptor.add_cat("polym", ispolym)
chemdescriptor.add_cat("actpol", isact)
chemdescriptor.add_prop("length", length)

ruleset = Ruleset()
ruleset.add_rule("Polym", ["polym", "polym"], joiner_direct, kpol)
