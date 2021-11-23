from pysb import *
from pysb.macros import *
from pysb.util import alias_model_components
from .macros import cleave_dimer


def add_biosensors(cc=1e7):
    """Add biosensors monomers to model. All start with an initial
    concentration of 1e7 and the cleaving reactions proposed in the paper by
    Corbat et al. (2018)."""
    sensors = ['sCas3', 'sCas8', 'sCas9']
    sensors_monomer = {}
    for sensor in sensors:
        sensors_monomer[sensor] = Monomer(sensor, ['sl', 'bf']) # sl for sensitive linker and es for enzyme site

    sensor_cc = {sensor: cc for sensor in sensors}
    sensor_ini = {}

    for sensor in sensors:
        sensor_ini[sensor] = Parameter('d' + sensor + '_0', sensor_cc[sensor])

    alias_model_components()

    for sensor in sensors:
        Initial(
            sensors_monomer[sensor](sl=1, bf=None) % sensors_monomer[sensor](
                sl=1, bf=None),
            sensor_ini[sensor])

    # Forward rates have been halved because there are two binding sites for
    # enzymes and this causes the rate to be doubled.
    sensor_cleavers = {'sCas3': (C3(state='A'), (1e-6 / 2, 1e-2, 1)),
                       'sCas9': (Apop(bf=None), (5e-9 / 2, 1e-3, 1)),
                       'sCas8': (C8(state='A'), (1e-7 / 2, 1e-3, 1))}

    for sensor in sensors:
        dimer = sensors_monomer[sensor](sl=1, bf=None) % sensors_monomer[
            sensor](sl=1, bf=None)
        cleave_dimer(sensor_cleavers[sensor][0], 'bf', dimer, 'bf', 'sl',
                     sensor_cleavers[sensor][1])


def add_apaf_biosensor_cleavage():
    """In double apoptosome activation we need the first activation of apaf to
    cleave biosensor for cas9 at a slower rate."""
    # We need to halve the forward rate because it has two binding sites
    klist = [2e-10, 1e-3, 1]
    dimer = sCas9(sl=1, bf=None) % sCas9(sl=1, bf=None)
    cleave_dimer(Apaf(state='A'), 'bf', dimer, 'bf', 'sl', klist)


def add_effector_bid_feedback():
    """Add a feedback from effector caspase to bid as reported by Slee et al. (2000, doi: 10.1038/sj.cdd.4400689)."""
    klist = [1e-6, 1e-3, 1]  # Taken from caspase 3
    catalyze(C3(state='A'), 'bf', Bid(state='U'), 'bf', Bid(state='T'), klist)


def observe_biosensors():
    """Add biosensors in monomeric and dimeric state as observables."""
    alias_model_components()
    sensor_dict = {'sCas3': sCas3,
                   'sCas8': sCas8,
                   'sCas9': sCas9}

    for sensor_name, sensor in sensor_dict.items():
        Observable(sensor_name + '_monomer', sensor(sl=None, bf=None))
        Observable(sensor_name + '_dimer',
                   sensor(sl=1, bf=None) % sensor(sl=1, bf=None),
                   match='species')


def observe_caspases():
    """Add caspases in inactive and active state as observables."""
    alias_model_components()
    caspase_dict = {'Cas3': C3,
                   'Cas8': C8,
                   'Cas6': C6}

    for caspase_name, caspase in caspase_dict.items():
        Observable(caspase_name + '_inactive', caspase(state='pro'))
        Observable(caspase_name + '_active', caspase(state='A'))

    Observable('Cas9_inactive', C9)
    Observable('Apaf_inactive', Apaf(state='I'))
    Observable('Apaf_active', Apaf(state='A'))
    Observable('Apop_active', Apop)


def remove_extrinsic_stimuli(model):
    """Sets ligand to 0 to remove extrinsic stimuli."""
    model.parameters['L_0'].value = 0


def intrinsic_stimuli(model=None):
    """Add instrinsic stimuli activation through Bid truncation. If model is
    given, then extrinsic activation is removed from it."""

    Monomer('IntrinsicStimuli', ['bf'])
    Parameter('IntrinsicStimuli_0', 1e2)

    alias_model_components()

    Initial(IntrinsicStimuli(bf=None), IntrinsicStimuli_0)

    # =====================
    # tBid Intrinsic Activation Rules
    # ---------------------
    #        Bid + IntrinsicStimuli <--> Bid:IS --> tBid + IS
    # ---------------------
    catalyze(IntrinsicStimuli(bf=None), 'bf',  Bid(state='U'), 'bf',
             Bid(state='T'), [1e-6, 1e-3, 1])

    if model:
        # Remove extrinsic stimuli
        remove_extrinsic_stimuli(model)


def choose_stimuli(model, stimuli):
    """Adapts the model to the corresponding stimuli, checking that intrinsic
    (modify model) and extrinsic are the only possibilities."""
    if stimuli == 'intrinsic':
        intrinsic_stimuli(model)
        return model
    elif stimuli == 'extrinsic':
        return model
    else:
        raise ValueError('Stimuli can be either extrinsic or intrinsic.')
