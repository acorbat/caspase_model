import numpy as np
import pytest
from earm import albeck_modules
from pysb import *  # noqa: F403
from pysb.simulator import ScipyOdeSimulator
from simbio import Simulator
from simbio.simulator.solvers.scipy import ODEint

from caspase_model.simbio_model import albeck
from caspase_model.tests.name_mapping import name_mapping


def load_pysb_model(func, pore: bool):
    """Load an Albeck model.

    func: Callable
        Specific Albeck model function, such as albeck_11b.
    pore: bool
        If True, do pore transport in the model.
    """

    model = Model()  # noqa: F405
    # Declare monomers
    albeck_modules.all_monomers()
    # Generate the upstream and downstream sections
    albeck_modules.rec_to_bid()
    albeck_modules.pore_to_parp()
    # The specific MOMP model to use
    func(do_pore_transport=pore)
    return model


def pysb_dataframe(result, model):
    """Rename dataframe from PySB"""
    names = dict(zip(result.all.dtype.names, map(str, model.species)))
    return result.dataframe.rename(columns=names).rename(columns=name_mapping)


# Mapping of models between SimBio and PySB
MODELS = [
    (albeck.Albeck11b, load_pysb_model(albeck_modules.albeck_11b, pore=False)),
    (
        albeck.Albeck11bPoreTransport,
        load_pysb_model(albeck_modules.albeck_11b, pore=True),
    ),
    (albeck.Albeck11c, load_pysb_model(albeck_modules.albeck_11c, pore=False)),
    (
        albeck.Albeck11cPoreTransport,
        load_pysb_model(albeck_modules.albeck_11c, pore=True),
    ),
    (albeck.Albeck11d, load_pysb_model(albeck_modules.albeck_11d, pore=False)),
    (
        albeck.Albeck11dPoreTransport,
        load_pysb_model(albeck_modules.albeck_11d, pore=True),
    ),
    (albeck.Albeck11e, load_pysb_model(albeck_modules.albeck_11e, pore=False)),
    (
        albeck.Albeck11ePoreTransport,
        load_pysb_model(albeck_modules.albeck_11e, pore=True),
    ),
    (albeck.Albeck11f, load_pysb_model(albeck_modules.albeck_11f, pore=False)),
    (
        albeck.Albeck11fPoreTransport,
        load_pysb_model(albeck_modules.albeck_11f, pore=True),
    ),
]


@pytest.mark.parametrize("simbio_model, pysb_model", MODELS)
def test_model(simbio_model, pysb_model):
    """Run models with SimBio and PySB, and compare results at each timepoint."""

    t = np.linspace(0, 20_000, 1_000)
    solver_options = {"atol": 1e-6, "rtol": 1e-6}

    # SimBio
    sim = Simulator(
        simbio_model, builder="numpy", solver=ODEint, solver_kwargs=solver_options
    )
    _, df_simbio = sim.run(t)

    # PySB
    sim_pysb = ScipyOdeSimulator(
        pysb_model, integrator="lsoda", integrator_options=solver_options
    )
    df_pysb = pysb_dataframe(sim_pysb.run(t), pysb_model)

    assert np.allclose(df_pysb, df_simbio[df_pysb.columns], rtol=1e-2, atol=1e-2)
