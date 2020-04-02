from pysb import *
from pysb.util import alias_model_components
from .macros import cleave_dimer


def add_biosensors():
    """Add biosensors monomers to model. All start with an initial
    concentration of 1e7 and the cleaving reactions proposed in the paper by
    Corbat et al. (2018)."""
    sensors = ['BFP', 'Cit', 'mKate']
    sensors_monomer = {}
    for sensor in sensors:
        sensors_monomer[sensor] = Monomer(sensor, ['sl', 'bf']) # sl for sensitive linker and es for enzyme site

    sensor_cc = {sensor: 1e7 for sensor in sensors}
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
    sensor_cleavers = {'BFP': (C3(state='A'), (1e-6 / 2, 1e-2, 1)),
                       'Cit': (Apop(bf=None), (5e-9 / 2, 1e-3, 1)),
                       'mKate': (C8(state='A'), (1e-7 / 2, 1e-3, 1))}

    for sensor in sensors:
        dimer = sensors_monomer[sensor](sl=1, bf=None) % sensors_monomer[
            sensor](sl=1, bf=None)
        cleave_dimer(sensor_cleavers[sensor][0], 'bf', dimer, 'bf', 'sl',
                     sensor_cleavers[sensor][1])


def observe_biosensors():
    """Add biosensors in monomeric and dimeric state as observables."""
    alias_model_components()
    sensor_dict = {'BFP': BFP,
                   'Cit': Cit,
                   'mKate': mKate}

    for sensor_name, sensor in sensor_dict.items():
        Observable(sensor_name + '_monomer', sensor(sl=None, bf=None))
        Observable(sensor_name + '_dimer',
                   sensor(sl=1, bf=None) % sensor(sl=1, bf=None),
                   match='species')
