from copy import deepcopy

# import * should be sufficient, but explicit everything used for better readability
# from metadynamic.polymers import *
from metadynamic.polymers import (
    model,
    Categorizer,
    ProdBuilder,
    ConstBuilder,
    novariant,
    ispolym,
)

# Necessary/useful for avinding clash between polymers/catpolymers model???
model = deepcopy(model)

# Example of catalized polymer by a dimer
# iff the dimer ends corresponds to polymers ends
# e.g.|   aaa + bbb + ab --> aaabbb + ab

isdimer: Categorizer = lambda name: ispolym(name) and len(name) == 2
cat_polym: ProdBuilder = lambda names, variant: (names[0] + names[1], names[2])
k_cat_dimer_pol: ConstBuilder = lambda names, k, variant: (
    k[0] if names[0][-1] == names[2][0] and names[1][0] == names[2][-1] else 0.0
)

model.add_cat("dimer", isdimer)

model.add_rule(
    rulename="dP",
    reactants=("polym", "polym", "dimer"),
    builder=(cat_polym, k_cat_dimer_pol, novariant),
    descr="Catalized polym by dimer if ends fits",
)
