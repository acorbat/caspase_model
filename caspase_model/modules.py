from pysb import *
from pysb.util import alias_model_components
from earm.shared import *
from pysb.macros import equilibrate


# Default forward, reverse, and catalytic rates:

KF = 1e-6
KR = 1e-3
KC = 1


def pore_to_parp_double_apop():
    """Defines what happens after the pore is activated and Cytochrome C and
    Smac are released.
    Uses CytoC, Smac, Apaf, Apop, C3, C6, C8, PARP, XIAP monomers and their
    associated parameters to generate the rules that describe apoptosome
    formation, caspase-3 further activation of apoptosome, XIAP inhibition,
    activation of caspases (including caspase-6-mediated feedback), and
    cleavage of effector caspase substrates as specified in EARM 1.0.
    Declares initial conditions for CytoC, Smac, Apaf-1, caspases
    3 and 6, XIAP, and PARP.
    """

    # Declare initial conditions:

    Parameter('Apaf_0'  , 1.0e3) # Apaf-1
    Parameter('C3_0'    , 1.0e4) # procaspase-3 (pro-C3)
    Parameter('C6_0'    , 1.0e4) # procaspase-6 (pro-C6)
    Parameter('XIAP_0'  , 1.0e4) # X-linked inhibitor of apoptosis protein
    Parameter('PARP_0'  , 1.0e6) # C3* substrate

    alias_model_components()

    Initial(Apaf(bf=None, state='I'), Apaf_0)
    Initial(C3(bf=None, state='pro'), C3_0)
    Initial(C6(bf=None, state='pro'), C6_0)
    Initial(PARP(bf=None, state='U'), PARP_0)
    Initial(XIAP(bf=None), XIAP_0)

    # CytoC and Smac activation after release
    # --------------------------------------

    equilibrate(Smac(bf=None, state='C'), Smac(bf=None, state='A'),
                          transloc_rates)

    equilibrate(CytoC(bf=None, state='C'), CytoC(bf=None, state='A'),
                          transloc_rates)

    # Apoptosome formation
    # --------------------
    #   Apaf + cCytoC <-->  Apaf:cCytoC --> aApaf + cCytoC
    #   aApaf + pC3 <-->  aApaf:pC3 --> aApaf + C3
    #   C3 + aApaf <-->  C3:aApaf --> Apop
    #   Apop + pC3 <-->  Apop:pC3 --> Apop + C3

    catalyze(CytoC(state='A'), Apaf(state='I'), Apaf(state='A'), [5e-7, KR, KC])
    catalyze(Apaf(state='A'), C3(state='pro'), C3(bf=None, state='A'), [5e-09, KR, KC])
    catalyze(C3(bf=None, state='A'), Apaf(state='A'), Apop(bf=None), [1.3e-06, KR, KC])
    catalyze(Apop(), C3(state='pro'), C3(bf=None, state='A'), [5e-9, KR, KC])

    # Apoptosome-related inhibitors
    # -----------------------------
    #   Apaf + XIAP <-->  Apaf:XIAP
    #   cSmac + XIAP <-->  cSmac:XIAP

    bind(Apaf(state='A'), XIAP(), [2e-6, KR])
    bind(Smac(state='A'), XIAP(), [7e-6, KR])

    # Caspase reactions
    # -----------------
    # Includes effectors, inhibitors, and feedback initiators:
    #
    #   pC3 + C8 <--> pC3:C8 --> C3 + C8 CSPS
    #   pC6 + C3 <--> pC6:C3 --> C6 + C3 CSPS
    #   XIAP + C3 <--> XIAP:C3 --> XIAP + C3_U CSPS
    #   PARP + C3 <--> PARP:C3 --> CPARP + C3 CSPS
    #   pC8 + C6 <--> pC8:C6 --> C8 + C6 CSPS
    catalyze(C8(state='A'), C3(state='pro'), C3(state='A'), [1e-7, KR, KC])
    catalyze(XIAP(), C3(state='A'), C3(state = 'ub'), [2e-6, KR, 1e-1])
    catalyze(C3(state='A'), PARP(state='U'), PARP(state='C'), [KF, 1e-2, KC])
    catalyze(C3(state='A'), C6(state='pro'), C6(state='A'), [KF, KR, KC])
    catalyze(C6(state='A'), C8(state='pro'), C8(state='A'), [3e-8, KR, KC])
