from metadynamic.inval import invalidint
from metadynamic.ruleset import (
    Categorizer,
    Propertizer,
    ProdBuilder,
    ConstBuilder,
    VariantBuilder,
    Compset,
    Paramset,
    Model,
)

# Categorizer

ismono: Categorizer = lambda name: len(name) == 1
ispolym: Categorizer = lambda name: name.isalpha()
isact: Categorizer = lambda name: name[-1] == "*" and name[:-1].isalpha()
isactmono: Categorizer = lambda name: isact(name) and len(name) == 2

# Propertizer

length: Propertizer = lambda name: (
    1
    if ismono(name)
    else len(name)
    if ispolym(name)
    else len(name) - 1
    if isact(name)
    else 0
)


# ProdBuilder


def joiner(sep: str) -> ProdBuilder:
    """Generate a "joiner" ProdBuilder using a "sep" as a separator string.
       e.g. chainer=joiner("-") will give a ProdBuilder named chainer
       that will provide "A-B-C" from chainer(["A","B","C"])"""
    return lambda names, variant: (sep.join(names),)


cut: ProdBuilder = lambda names, variant: (names[0][:variant], names[0][variant:])
act_polym: ProdBuilder = lambda names, variant: (names[0][:-1] + names[1],)
activ: ProdBuilder = lambda names, variant: (names[0] + "*",)
deactiv: ProdBuilder = lambda names, variant: (names[0][:-1],)
epimer: ProdBuilder = lambda names, variant: (
    names[0][:variant] + names[0][variant].swapcase() + names[0][variant + 1 :],
)


# ConstBuilder


def samecase(one: str, two: str) -> bool:
    return (one.islower() and two.islower()) or (one.isupper() and two.isupper())


kfastmono: ConstBuilder = lambda names, k, variant: (
    k[0] * k[1] if length(names[0]) == 1 else k[0]
)
khyd: ConstBuilder = lambda names, k, variant: (
    k[0] if samecase(names[0][variant - 1], names[0][variant]) else k[1] * k[0]
)
kactselect: ConstBuilder = lambda names, k, variant: (
    k[0] if samecase(names[0][-2], names[1][0]) else k[0] * k[1]
)


def kmidselect(names: Compset, k: Paramset, variant: int) -> float:  # ConstBuilder
    name = names[0]
    res = k[0]
    if variant < (length(name) - 1) and samecase(name[variant], name[variant + 1]):
        res *= k[1]
    if (variant > 0) and samecase(name[variant], name[variant - 1]):
        res *= k[1]
    return res


# VariantBuilder

novariant: VariantBuilder = lambda reactants: (invalidint,)
intervariant: VariantBuilder = lambda reactants: range(length(reactants[0]) - 1)
lenvariant: VariantBuilder = lambda reactants: range(length(reactants[0]))


def singlevariant(num: int) -> VariantBuilder:
    return lambda reactants: (num,)


# Define a specific ruleset model


model = Model()

model.add_cat("mono", ismono)
model.add_cat("polym", ispolym)
model.add_cat("actpol", isact)
model.add_cat("actmono", isactmono)

model.add_prop("length", length)


model.add_rule(
    rulename="P",
    reactants=("polym", "polym"),
    builder=(joiner(""), kfastmono, novariant),
    descr="Polymerization",
)
model.add_rule(
    rulename="A",
    reactants=("actpol", "polym"),
    builder=(act_polym, kactselect, novariant),
    descr="Activated Polymerization",
)
model.add_rule(
    rulename="M",
    reactants=("actmono", "polym"),
    builder=(act_polym, kactselect, novariant),
    descr="Activated Monomer Polymerization",
)
model.add_rule(
    rulename="a",
    reactants=("polym",),
    builder=(activ, kfastmono, novariant),
    descr="Activation",
)
model.add_rule(
    rulename="d",
    reactants=("actpol",),
    builder=(deactiv, kfastmono, novariant),
    descr="Deactivation",
)
model.add_rule(
    rulename="H",
    reactants=("polym",),
    builder=(cut, khyd, intervariant),
    descr="Hydrolysis",
)
model.add_rule(
    rulename="R",
    reactants=("polym",),
    builder=(epimer, kmidselect, lenvariant),
    descr="Epimerization",
)
model.add_rule(
    rulename="E",
    reactants=("polym",),
    builder=(epimer, kmidselect, singlevariant(0)),
    descr="Epimerization at first end",
)
