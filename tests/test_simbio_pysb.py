"""Compare the SimBio models against their PySB counterparts."""

from collections.abc import Callable

import numpy as np
import pandas as pd
import pytest
from earm import albeck_modules
from poincare.solvers import LSODA
from pysb import Model
from pysb.simulator import ScipyOdeSimulator
from simbio import Simulator, System

from caspase_model import models as pysb_models
from caspase_model.simbio_model import models as simbio_models

from .name_mapping import name_mapping

TIME = np.linspace(0, 20_000, 1_001)
ATOL = 1e-6
RTOL = 1e-8


def earm_model(momp: str) -> Callable[[], Model]:
    """EARM 1.0 model with the given MOMP variant (and pore transport)."""

    def load() -> Model:
        model = Model()
        albeck_modules.all_monomers()
        albeck_modules.rec_to_bid()
        albeck_modules.pore_to_parp()
        getattr(albeck_modules, f"albeck_{momp}")(do_pore_transport=True)
        return model

    return load


def stimulated(model: Callable[..., Model], stimuli: str) -> Callable[[], Model]:
    return lambda: model(stimuli=stimuli)


MODELS = {
    "Albeck11b": (simbio_models.Albeck11b, earm_model("11b")),
    "Albeck11c": (simbio_models.Albeck11c, earm_model("11c")),
    "Albeck11d": (simbio_models.Albeck11d, earm_model("11d")),
    "Albeck11e": (simbio_models.Albeck11e, earm_model("11e")),
    "Albeck11f": (simbio_models.Albeck11f, earm_model("11f")),
    "AlbeckAsMatlab": (simbio_models.AlbeckAsMatlab, pysb_models.albeck_as_matlab),
    "Corbat2018-extrinsic": (
        simbio_models.Corbat2018,
        stimulated(pysb_models.corbat_2018, "extrinsic"),
    ),
    "Corbat2018-intrinsic": (
        simbio_models.intrinsic(simbio_models.Corbat2018),
        stimulated(pysb_models.corbat_2018, "intrinsic"),
    ),
    "ARM-extrinsic": (
        simbio_models.ARM,
        stimulated(pysb_models.arm, "extrinsic"),
    ),
    "ARM-intrinsic": (
        simbio_models.intrinsic(simbio_models.ARM),
        stimulated(pysb_models.arm, "intrinsic"),
    ),
}


def solve_simbio(model: System | type[System]) -> pd.DataFrame:
    sim = Simulator(model).with_solver(LSODA(atol=ATOL, rtol=RTOL))
    return sim.solve(save_at=TIME).to_dataframe()


def solve_pysb(model: Model, names: set[str]) -> pd.DataFrame:
    """Solve the PySB model and rename its species to the SimBio names."""
    sim = ScipyOdeSimulator(
        model,
        integrator="lsoda",
        integrator_options={"atol": ATOL, "rtol": RTOL},
        compiler="python",
    )
    result = sim.run(TIME)
    columns = []
    for species in model.species:
        (name,) = names.intersection(name_mapping[str(species)])
        columns.append(name)
    return pd.DataFrame(result.species, index=TIME, columns=columns)


@pytest.mark.parametrize("simbio_model, pysb_model", MODELS.values(), ids=MODELS)
def test_model(simbio_model, pysb_model):
    """Run models with SimBio and PySB, and compare results at each timepoint."""
    simbio = solve_simbio(simbio_model)
    pysb = solve_pysb(pysb_model(), set(simbio.columns))

    # Each PySB species is a different SimBio species,
    assert pysb.columns.is_unique
    # and SimBio species absent in PySB remain at zero.
    extra = simbio.columns.difference(pysb.columns)
    assert (simbio[extra] == 0).all(axis=None), extra

    simbio = simbio[pysb.columns]
    # Compare relative to the range of each species,
    # but at least to 1 molecule, as some species remain at ~0.
    scale = pysb.abs().max().clip(lower=1)
    np.testing.assert_allclose(simbio / scale, pysb / scale, rtol=0, atol=1e-4)
