from typing import List, Optional, Iterator, Callable, Type

from metadynamic.utils import samecase


class Descr:
    _descrtype = "Chemical Description"

    def __init__(self, name: str):
        self.name: str = name

    @property
    def categories(self) -> List[str]:
        return []

    def __repr__(self) -> str:
        return f"{self._descrtype}: {self.name}"

    def __str__(self) -> str:
        return self.name


class CompDescr(Descr):
    _descrtype = "Compound Description"

    def __bool__(self):
        return bool(self.name)

    # Check memoize effects
    @property
    # @memoize_property
    def length(self) -> int:
        return len(self.name)

    @property
    # @memoize_property
    def ismono(self) -> bool:
        return self.length == 1

    @property
    # @memoize_property
    def ispolym(self) -> bool:
        return self.name.isalpha()

    @property
    # @memoize_property
    def isact(self) -> bool:
        return self.name[-1] == "*" and self.name[:-1].isalpha()

    @property
    # @memoize_property
    def categories(self) -> List[str]:
        res = []
        if self.ispolym:
            res.append("polym")
            if self.ismono:
                res.append("mono")
        elif self.isact:
            res.append("actpol")
        return res

    @property
    # @memoize_property
    def actreac(self) -> str:
        assert not self.isact, f"{self} is already activated"
        return self.name + "*"

    @property
    # @memoize_property
    def unactreac(self) -> str:
        assert self.isact, f"{self} is not activated"
        return self.name[:-1]

    # @memoize_oneparam
    def cut(self, pos: int) -> List[str]:
        name = self.name
        assert (
            1 <= pos < len(name)
        ), f"cut position should be at least 1, at most the chain length minus one"
        return [name[:pos], name[pos:]]

    # @memoize_oneparam
    def epimer(self, pos: int) -> str:
        name = self.name
        assert (
            0 <= pos < len(name)
        ), f"epimer position should be at least 0, at most the chain length minus one"
        return name[:pos] + name[pos].swapcase() + name[pos + 1 :]

    # @memoize_oneparam
    def extract(self, pos: int) -> str:
        return self.name[pos]


class ReacDescr(Descr):
    _descrtype = "Reaction Description"
    kind = "X"
    nbreac = 0

    def __init__(
        self,
        name: str,
        reactants: List[CompDescr],
        catal: Optional[CompDescr],
        pos: Optional[int],
        const: float,
        altconst: float,
        catconst: float,
    ):
        assert len(reactants) == self.nbreac
        if not name:
            reactantnames = (
                ""
                if self.nbreac == 0
                else reactants[0].name
                if self.nbreac == 1
                else reactants[0].name + "+" + reactants[1].name
            )
            name = f"{self.kind}.{reactantnames}.{catal if catal else ''}.{'' if pos is None else pos}"
        super().__init__(name)
        assert self.name[0] == self.kind, f"{self.name} created as {self.kind}"
        self.reactants: List[CompDescr] = reactants
        self.catal: Optional[CompDescr] = catal
        self._const: float = const
        self._altconst: float = altconst
        self._catconst: float = catconst
        self.pos: Optional[int] = pos
        self.checkdimer()

    @property
    def categories(self) -> List[str]:
        return [self.kind]

    def checkdimer(self) -> None:
        self.dimer = False

    def build_products(self) -> List[str]:
        return []

    def build_const(self) -> float:
        const = (self._altconst ** self.nbalt()) * self._const
        if self.catal:
            const *= self._catconst
        return const

    def nbalt(self) -> int:
        return 0

    @classmethod
    def fullfrom(
        cls, reactant: "Compound", const: float, altconst: float, catconst: float
    ) -> List["ReacDescr"]:
        return []

    @classmethod
    def subclasses(cls) -> Iterator[Type["ReacDescr"]]:
        for subclass in cls.__subclasses__():
            yield from subclass.subclasses()
            yield subclass

    @classmethod
    def kinds(cls) -> Iterator[str]:
        return (subcls.kind for subcls in cls.subclasses())


class PolymDescr(ReacDescr):
    kind = "P"
    nbreac = 2

    def checkdimer(self) -> None:
        self.dimer = self.reactants[0].name == self.reactants[1].name

    def build_products(self) -> List[str]:
        return [self.reactants[0].name + self.reactants[1].name]

    def nbalt(self) -> int:
        return 1 if self.reactants[0].ismono else 0

    @classmethod
    def fullfrom(
        cls, reactant: "Compound", const: float, altconst: float, catconst: float
    ) -> List[ReacDescr]:
        first: Callable[["Compound"], ReacDescr] = lambda other: cls(
            "",
            [reactant.description, other.description],
            None,
            None,
            const,
            altconst,
            catconst,
        )
        last: Callable[["Compound"], ReacDescr] = lambda other: cls(
            "",
            [other.description, reactant.description],
            None,
            None,
            const,
            altconst,
            catconst,
        )
        if reactant.description.ispolym:
            return [
                create(other)
                for other in reactant.collect.cat_list("polym")
                for create in (first, last)
            ]
        return []


class ActpolymDescr(ReacDescr):
    kind = "A"
    nbreac = 2

    def build_products(self) -> List[str]:
        return [self.reactants[0].unactreac + self.reactants[1].name]

    def nbalt(self) -> int:
        name0 = self.reactants[0].extract(-2)
        name1 = self.reactants[1].extract(0)
        return 0 if samecase(name0, name1) else 1

    @classmethod
    def fullfrom(
        cls, reactant: "Compound", const: float, altconst: float, catconst: float
    ) -> List[ReacDescr]:
        if reactant.description.ispolym:
            return [
                cls(
                    "",
                    [other.description, reactant.description],
                    None,
                    None,
                    const,
                    altconst,
                    catconst,
                )
                for other in reactant.collect.cat_list("actpol")
            ]
        if reactant.description.isact:
            return [
                cls(
                    "",
                    [reactant.description, other.description],
                    None,
                    None,
                    const,
                    altconst,
                    catconst,
                )
                for other in reactant.collect.cat_list("polym")
            ]
        return []


class ActivationDescr(ReacDescr):
    kind = "a"
    nbreac = 1

    def build_products(self) -> List[str]:
        return [self.reactants[0].actreac]

    def nbalt(self) -> int:
        return 0 if self.reactants[0].ismono else 1

    @classmethod
    def fullfrom(
        cls, reactant: "Compound", const: float, altconst: float, catconst: float
    ) -> List[ReacDescr]:
        if reactant.description.ispolym:
            return [
                cls("", [reactant.description], None, None, const, altconst, catconst)
            ]
        return []


class DectivationDescr(ReacDescr):
    kind = "d"
    nbreac = 1

    def build_products(self) -> List[str]:
        return [self.reactants[0].unactreac]

    def nbalt(self) -> int:
        return 0 if self.reactants[0].length == 2 else 1

    @classmethod
    def fullfrom(
        cls, reactant: "Compound", const: float, altconst: float, catconst: float
    ) -> List[ReacDescr]:
        if reactant.description.isact:
            return [
                cls("", [reactant.description], None, None, const, altconst, catconst)
            ]
        return []


class HydroDescr(ReacDescr):
    kind = "H"
    nbreac = 1

    def build_products(self) -> List[str]:
        assert self.pos is not None
        return self.reactants[0].cut(self.pos)

    def nbalt(self) -> int:
        assert self.pos is not None
        left = self.reactants[0].extract(self.pos - 1)
        right = self.reactants[0].extract(self.pos)
        return 0 if samecase(right, left) else 1

    @classmethod
    def fullfrom(
        cls, reactant: "Compound", const: float, altconst: float, catconst: float
    ) -> List[ReacDescr]:
        if reactant.description.ispolym:
            return [
                cls("", [reactant.description], None, i, const, altconst, catconst)
                for i in range(1, reactant.length)
            ]
        return []


class RacemDescr(ReacDescr):
    kind = "R"
    nbreac = 1

    def build_products(self) -> List[str]:
        assert self.pos is not None
        return [self.reactants[0].epimer(self.pos)]

    def nbalt(self) -> int:
        assert self.pos is not None, f"{self.name} as no position..."
        name = self.reactants[0]
        pos = self.pos
        res = 0
        if pos < (name.length - 1) and samecase(
            name.extract(pos), name.extract(pos + 1)
        ):
            res += 1
        if (pos > 0) and samecase(name.extract(pos), name.extract(pos - 1)):
            res += 1
        return res

    @classmethod
    def fullfrom(
        cls, reactant: "Compound", const: float, altconst: float, catconst: float
    ) -> List[ReacDescr]:
        if reactant.description.ispolym:
            return [
                cls("", [reactant.description], None, i, const, altconst, catconst)
                for i in range(reactant.length)
            ]
        return []


class EpimDescr(RacemDescr):
    kind = "E"

    @classmethod
    def fullfrom(
        cls, reactant: "Compound", const: float, altconst: float, catconst: float
    ) -> List[ReacDescr]:
        if reactant.description.ispolym and not reactant.description.ismono:
            return [cls("", [reactant.description], None, 0, const, altconst, catconst)]
        return []
