from pysb import *
from pysb.macros import *
from pysb.util import alias_model_components
from .shared import add_biosensors

def albeck_as_matlab():
    """Loads model as stated in Albeck's 2008 paper."""
    from earm.albeck_11f import model

    # Add citoplasmic Bcl2 as it was in Albeck's model because it's absent in
    # EARM implementation.
    Monomer('Bcl2c', ['b'])
    Parameter('Bcl2c_0', 2e4)

    alias_model_components()
    Initial(Bcl2c(b=None), Bcl2c_0)

    bind(Bid(bf=None, state='T'), 'bf', Bcl2c(b=None), 'b', [1e-6, 0.001])

    # There is another discrepancy in C8A and Bid kf
    model.parameters['bind_C8A_BidU_to_C8ABidU_kf'].value = 1e-7


def corbat_2018():
    """EARM model as modified by Corbat et al. (2018)."""
    albeck_as_matlab()
    add_biosensors()

    alias_model_components()

    XIAP_0.value = 1e2
    L_0.value = 1e3
    R_0.value = 1e3

    
