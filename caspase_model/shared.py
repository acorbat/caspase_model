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

    sensor_cleavers = {'BFP': (C3(state='A'), (1e-6, 1e-2, 1)),
                       'Cit': (Apop(bf=None), (5e-9, 1e-3, 1)),
                       'mKate': (C8(state='A'), (1e-7, 1e-3, 1))}

    for sensor in sensors:
        dimer = sensors_monomer[sensor](sl=1, bf=None) % sensors_monomer[
            sensor](sl=1, bf=None)
        cleave_dimer(sensor_cleavers[sensor][0], 'bf', dimer, 'bf', 'sl',
                     sensor_cleavers[sensor][1])
