from pysb import *
from pysb.macros import *
from pysb.util import alias_model_components
from .modules import pore_to_parp_double_apop
from .shared import add_biosensors, choose_stimuli, add_apaf_biosensor_cleavage


def albeck_as_matlab():
    """Loads model as stated in Albeck's 2008 paper."""
    from earm import albeck_modules
    model = Model()

    # Declare monomers
    albeck_modules.all_monomers()

    # Generate the upstream and downstream sections
    albeck_modules.rec_to_bid()
    albeck_modules.pore_to_parp()

    # The specific MOMP model to use
    albeck_modules.albeck_11e()

    # Add citoplasmic Bcl2 as it was in Albeck's model because it's absent in
    # EARM implementation.
    Monomer('Bcl2c', ['b'])
    Parameter('Bcl2c_0', 2e4)

    alias_model_components()
    Initial(Bcl2c(b=None), Bcl2c_0)

    bind(Bid(bf=None, state='T'), 'bf', Bcl2c(b=None), 'b', [1e-6, 0.001])

    # There is another discrepancy in C8A and Bid kf
    model.parameters['bind_C8A_BidU_to_C8ABidU_kf'].value = 1e-7

    return model


def corbat_2018(stimuli='extrinsic', add_CASPAM=True):
    """EARM model as modified by Corbat et al. (2018)."""
    model = albeck_as_matlab()
    if add_CASPAM:
        add_biosensors()

    alias_model_components()

    XIAP_0.value = 1e2
    L_0.value = 1e3
    R_0.value = 1e3

    model = choose_stimuli(model, stimuli)

    return model


def arm(stimuli='extrinsic', add_CASPAM=True):
    """Returns the new apoptotic reaction model. stimuli can be extrinsic
    (default) or intrinsic.

    Parameters
    ----------
    stimuli: string (default: extrinsic)
        Either intrinsic or extrinsic to choose stimuli
    add_CASPAM: bool (default: True)
        if True, CASPAM biosensors are added to the model
    """
    from earm import albeck_modules
    model = Model()

    # Declare monomers
    albeck_modules.all_monomers()

    # Generate the upstream and downstream sections
    albeck_modules.rec_to_bid()
    pore_to_parp_double_apop()

    # The specific MOMP model to use
    albeck_modules.albeck_11e()

    # Add citoplasmic Bcl2 as it was in Albeck's model because it's absent in
    # EARM implementation.
    Monomer('Bcl2c', ['b'])
    Parameter('Bcl2c_0', 2e4)

    alias_model_components()
    Initial(Bcl2c(b=None), Bcl2c_0)

    bind(Bid(bf=None, state='T'), 'bf', Bcl2c(b=None), 'b', [1e-6, 0.001])

    # There is another discrepancy in C8A and Bid kf
    model.parameters['bind_C8A_BidU_to_C8ABidU_kf'].value = 1e-7

    if add_CASPAM:
        # Add biosensors and their perturbation
        add_biosensors()

        model.parameters['dsCas3_0'].value = 7.5e5
        model.parameters['dsCas8_0'].value = 7.5e5
        model.parameters['dsCas9_0'].value = 7.5e5

        model.parameters['bind_sCas3sCas3_C3A_kf'].value = 2.8e-07
        model.parameters['bind_sCas8sCas8_C8A_kf'].value = 5.4e-08
        model.parameters['bind_sCas9sCas9_Apop_kf'].value = 2.8e-07

        # Add interaction between single activation
        add_apaf_biosensor_cleavage()

    # Other parameters that need to be corrected
    model.parameters['L_0'].value = 1e3
    model.parameters['Smac_0'].value = 1e5
    model.parameters['CytoC_0'].value = 1e5
    
    return choose_stimuli(model, stimuli)


def albeck_apoptosome_corrected(stimuli='extrinsic'):
    """to_deprecate: Loads model as stated in Albeck's 2008 paper and modifies
    apoptosome dyanmics."""
    from earm import albeck_modules
    model = Model()

    # Declare monomers
    albeck_modules.all_monomers()

    # Generate the upstream and downstream sections
    albeck_modules.rec_to_bid()
    pore_to_parp_double_apop()

    # The specific MOMP model to use
    albeck_modules.albeck_11e()

    # Add citoplasmic Bcl2 as it was in Albeck's model because it's absent in
    # EARM implementation.
    Monomer('Bcl2c', ['b'])
    Parameter('Bcl2c_0', 2e4)

    alias_model_components()
    Initial(Bcl2c(b=None), Bcl2c_0)

    bind(Bid(bf=None, state='T'), 'bf', Bcl2c(b=None), 'b', [1e-6, 0.001])

    # There is another discrepancy in C8A and Bid kf
    model.parameters['bind_C8A_BidU_to_C8ABidU_kf'].value = 1e-7

    return choose_stimuli(model, stimuli)
