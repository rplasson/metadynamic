from typing import Dict, Any

from metadynamic.ruleset import (
    Categorizer,
    Propertizer,
    ProdBuilder,
    ConstBuilder,
    VariantBuilder,
    kalternate,
    rangevariant,
)


# Categorizer

# A polymer is a chain of characters, e.g. abc
polym: Categorizer = lambda name: name.isalpha()

# Propertizer
length: Propertizer = lambda name: len(name)

# ProdBuilder

# e.g. abcdef -[3]-> abc + def
cut: ProdBuilder = lambda names, variant: (names[0][:variant], names[0][variant:])

# ConstBuilder

khyd: ConstBuilder = kalternate(
    condition=lambda names, variant: variant == 1,
    name_t="khyd_head",
    name_f="khyd_queue",
)

# (length-1) possible reactions from a given reaction
# (e.g. abc -[1]-> a+bc and abc -[2]->  ab+c)
intervariant: VariantBuilder = rangevariant(reacnum=0, first_offset=1)

# Default Ruleset

default_ruleset: Dict[str, Any] = {
    "categories": ["polym"],
    "properties": ["length"],
    "rules": {
        "H": {
            "reactants": ["polym"],
            "builder_func": "cut",
            "builder_const": "khyd",
            "builder_variant": "intervariant",
            "descr": "Hydrolysis",
        },
    },
}
