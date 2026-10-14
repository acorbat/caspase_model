"""Observables of the SimBio models.

Counterparts of :func:`caspase_model.shared.observe_caspases`
and :func:`caspase_model.shared.observe_biosensors`,
which sum each protein in a given state over all the complexes it is part of.

Use them as the output of a simulation:

>>> from simbio import Simulator
>>> from caspase_model.simbio_model import ARM
>>> sim = Simulator(ARM).with_transform(observables(ARM))
"""

from functools import reduce
from operator import add

from simbio import System

# Paths to the variables containing each protein state.
# Paths absent in a model are skipped.
CASPASES = {
    "Cas3_inactive": (
        "caspases.C3_pro",
        "caspases.C8_A_C3_pro",
        "apoptosome.Apop_C3_pro",
        "apoptosome.Apaf_A_C3_pro",
    ),
    "Cas3_active": (
        "caspases.C3_A",
        "caspases.XIAP_C3_A",
        "caspases.C3_A_PARP_U",
        "caspases.C3_A_C6_pro",
        "apoptosome.C3_A_Apaf_A",
        "sCas3.cleavage.enzyme_dimer",
    ),
    "Cas8_inactive": (
        "rec_to_bid.C8_pro",
        "rec_to_bid.DISC_C8_pro",
        "caspases.C6_A_C8_pro",
    ),
    "Cas8_active": (
        "rec_to_bid.C8_A",
        "rec_to_bid.C8_A_Bid_U",
        "rec_to_bid.BAR_C8_A",
        "caspases.C8_A_C3_pro",
        "sCas8.cleavage.enzyme_dimer",
    ),
    "Cas6_inactive": ("caspases.C6_pro", "caspases.C3_A_C6_pro"),
    "Cas6_active": ("caspases.C6_A", "caspases.C6_A_C8_pro"),
    "Cas9_inactive": ("apoptosome.C9",),
    "Apaf_inactive": ("apoptosome.Apaf_I", "apoptosome.CytoC_A_Apaf_I"),
    "Apaf_active": (
        "apoptosome.Apaf_A",
        "apoptosome.Apaf_A_C3_pro",
        "apoptosome.C3_A_Apaf_A",
        "apoptosome.Apaf_A_XIAP",
        "sCas9_apaf.enzyme_dimer",
    ),
    "Apop_active": (
        "apoptosome.Apop",
        "apoptosome.Apop_C3_pro",
        "apoptosome.Apop_XIAP",
        "sCas9.cleavage.enzyme_dimer",
    ),
}

BIOSENSORS = ("sCas3", "sCas8", "sCas9")


def _get(model, path: str):
    obj = model
    for name in path.split("."):
        obj = getattr(obj, name)
    return obj


def _total(model, paths: tuple[str, ...]):
    terms = []
    for path in paths:
        try:
            terms.append(_get(model, path))
        except AttributeError:
            continue
    if len(terms) == 0:
        return None
    return reduce(add, terms)


def caspase_observables(model: System | type[System]) -> dict:
    """Caspases (and Apaf) in their inactive and active states."""
    observables = {}
    for name, paths in CASPASES.items():
        total = _total(model, paths)
        if total is not None:
            observables[name] = total
    return observables


def biosensor_observables(model: System | type[System]) -> dict:
    """Biosensors in their monomeric and dimeric states."""
    observables = {}
    for name in BIOSENSORS:
        try:
            sensor = getattr(model, name)
        except AttributeError:
            continue
        observables[f"{name}_monomer"] = sensor.monomer
        observables[f"{name}_dimer"] = sensor.dimer
    return observables


def observables(model: System | type[System]) -> dict:
    """All observables present in the model."""
    return caspase_observables(model) | biosensor_observables(model)
