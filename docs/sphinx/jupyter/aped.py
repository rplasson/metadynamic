"""APED model definition module

This system is defined by:
 - Activation reaction; e.g.  A -> A*
 - Activated polymerisation reaction; e.g.  A* + BcD -> ABcD
 - Epimerization reaction; e.g.  ABcD -> aBcD
 - Depolymerization reaction; e.g.  aBcD -> aBc + D
"""

from typing import Dict, Any

from metadynamic.ruleset import (
    Categorizer,
    Propertizer,
    ProdBuilder,
    ConstBuilder,
    VariantBuilder,
    Compset,
    kalternate,
    kdualchoice,
    novariant_gen,
    singlevariant,
    rangevariant,
)

# Definition of Categorizer and Propertizer functions

polym: Categorizer = lambda name: name.isalpha()
"""Definition of a polymer

Return True if the name is an alphabetic string:
>>> polym("AaDcs")
True
>>> polym("A2")
False
"""

mono: Categorizer = lambda name: polym(name) and len(name) == 1
"""Definition of a monomer

Return True if the name is polymer of length one:
>>> mono("Abc")
False
>>> mono("b")
True
"""

actpol: Categorizer = lambda name: name[-1] == "*" and name[:-1].isalpha()
"""Definition of an activated polymer

Return True if the name is an alphabetic string + a '*' char at its end:
>>> actpol("AbcD")
False
>>> actpol("Abcd*")
True
"""


actmono: Categorizer = lambda name: actpol(name) and len(name) == 2
"""Definition of an activated monomer

Return True is the name is an activated polymer of length 2:
>>> actmono("Abcd*")
False
>>> actmono("A*")
True
"""

longpol: Categorizer = lambda name: polym(name) and len(name) > 1
"""Definition of a long polymer

Return True if the name is a polymer of length 2 or more:
>>> longpol("A")
False
>>> longpol("Abcd")
True
"""


length: Propertizer = lambda name: (
    1
    if mono(name)
    else len(name)
    if polym(name)
    else len(name) - 1
    if actpol(name)
    else 0
)
"""Definition of length property

Return the length of a monomer, polymer or activated polymer:
>>> length('A')
1
>>> length('AbCD')
4
>>> length('Ab*')
3
>>> length('A2')
0
"""


def asym(name: str) -> int:  # Propertizer
    """Definition of asymmetry property

    Each upper char in the name is counted as +1, each lower char as -1:
    >>> asym('AAAa')
    2
    >>> asym('AAAA')
    4
    >>>> asym('aaaa')
    -4
    >>>> asym('aAAa')
    0
    """
    res = 0
    for char in name:
        if char.isupper():
            res += 1
        elif char.islower():
            res -= 1
    return res


right: Categorizer = lambda name: asym(name) > 0
"""Definition of right category

Return True if the name has a positive asymmetry property.
"""

left: Categorizer = lambda name: asym(name) < 0
"""Definition of left category

Return True if the name has a negative asymmetry property.
"""


# Definition of ProdBuilder functions

cut: ProdBuilder = lambda names, variant: (names[0][:variant], names[0][variant:])
act_polym: ProdBuilder = lambda names, variant: (names[0][:-1] + names[1],)
activ: ProdBuilder = lambda names, variant: (names[0] + "*",)
epimer: ProdBuilder = lambda names, variant: (
    names[0][:variant] + names[0][variant].swapcase() + names[0][variant + 1 :],
)


# ConstBuilder


def samecase(one: str, two: str) -> bool:
    return (one.islower() and two.islower()) or (one.isupper() and two.isupper())


def samebefore(names: Compset, variant: int) -> bool:
    name = names[0]
    return variant < (length(name) - 1) and samecase(name[variant], name[variant + 1])


def sameafter(names: Compset, variant: int) -> bool:
    name = names[0]
    return (variant > 0) and samecase(name[variant], name[variant - 1])


kpola: ConstBuilder = kalternate(
    condition=lambda names, variant: samecase(names[0][-2], names[1][0]),
    name_t="kpola_same",
    name_f="kpola_diff",
)

kact: ConstBuilder = kalternate(
    condition=lambda names, variant: length(names[0]) == 1,
    name_t="kact_mono",
    name_f="kact_pol",
)

khyd: ConstBuilder = kalternate(
    condition=lambda names, variant: samecase(names[0][variant - 1], names[0][variant]),
    name_t="khyd_same",
    name_f="khyd_diff",
)

kepi: ConstBuilder = kdualchoice(
    condition_1=samebefore,
    condition_2=sameafter,
    name_tt="kepi_same",
    name_ff="kepi_diff",
    name_tf="kepi_mixed",
)

# VariantBuilder

novariant: VariantBuilder = novariant_gen()
intervariant: VariantBuilder = rangevariant(reacnum=0, first_offset=1)
firstonly: VariantBuilder = singlevariant(num=0)


# Default Ruleset

default_ruleset: Dict[str, Any] = {
    "categories": ["mono", "polym", "longpol", "actpol", "actmono", "left", "right"],
    "properties": ["length", "asym"],
    "rules": {
        "A": {
            "reactants": ["actpol", "polym"],
            "builder_func": "act_polym",
            "builder_const": "kpola",
            "builder_variant": "novariant",
            "descr": "Activated Polymerization",
        },
        "a": {
            "reactants": ["polym"],
            "builder_func": "activ",
            "builder_const": "kact",
            "builder_variant": "novariant",
            "descr": "Activation",
        },
        "H": {
            "reactants": ["polym"],
            "builder_func": "cut",
            "builder_const": "khyd",
            "builder_variant": "intervariant",
            "descr": "Hydrolysis",
        },
        "E": {
            "reactants": ["longpol"],
            "builder_func": "epimer",
            "builder_const": "kepi",
            "builder_variant": "firstonly",
            "descr": "Epimerization at first end",
        },
    },
}
