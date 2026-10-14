"""Compare the SimBio observables against their PySB counterparts."""

import numpy as np
import pandas as pd
import pytest
from poincare.solvers import LSODA
from pysb.simulator import ScipyOdeSimulator
from simbio import Simulator

from caspase_model import models as pysb_models
from caspase_model import shared
from caspase_model.simbio_model import models as simbio_models
from caspase_model.simbio_model.observables import observables

from .test_simbio_pysb import ATOL, RTOL, TIME, stimulated


def pysb_model(model, biosensors: bool):
    def load():
        m = model()
        shared.observe_caspases()
        if biosensors:
            shared.observe_biosensors()
        return m

    return load


MODELS = {
    "AlbeckAsMatlab": (
        simbio_models.AlbeckAsMatlab,
        pysb_model(pysb_models.albeck_as_matlab, biosensors=False),
    ),
    "Corbat2018-intrinsic": (
        simbio_models.intrinsic(simbio_models.Corbat2018),
        pysb_model(stimulated(pysb_models.corbat_2018, "intrinsic"), biosensors=True),
    ),
    "ARM-extrinsic": (
        simbio_models.ARM,
        pysb_model(stimulated(pysb_models.arm, "extrinsic"), biosensors=True),
    ),
    "ARM-intrinsic": (
        simbio_models.intrinsic(simbio_models.ARM),
        pysb_model(stimulated(pysb_models.arm, "intrinsic"), biosensors=True),
    ),
}


@pytest.mark.parametrize("simbio_model, pysb_model", MODELS.values(), ids=MODELS)
def test_observables(simbio_model, pysb_model):
    sim = Simulator(simbio_model).with_solver(LSODA(atol=ATOL, rtol=RTOL))
    simbio = sim.with_transform(observables(simbio_model)).solve(save_at=TIME)
    simbio = simbio.to_dataframe()

    model = pysb_model()
    result = ScipyOdeSimulator(
        model,
        integrator="lsoda",
        integrator_options={"atol": ATOL, "rtol": RTOL},
        compiler="python",
    ).run(TIME)
    pysb = pd.DataFrame(
        {obs.name: result.observables[obs.name] for obs in model.observables},
        index=TIME,
    )

    # Observables absent in SimBio have no species in PySB.
    missing = pysb.columns.difference(simbio.columns)
    assert (pysb[missing] == 0).all(axis=None), missing
    assert simbio.columns.difference(pysb.columns).empty

    pysb = pysb[simbio.columns]
    scale = pysb.abs().max().clip(lower=1)
    np.testing.assert_allclose(simbio / scale, pysb / scale, rtol=0, atol=1e-4)
