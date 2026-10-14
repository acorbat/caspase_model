"""Subsystem modules of the apoptotic network.

Each module is a :class:`simbio.System` describing a part of the network.
A module owns the species it introduces,
with their default initial conditions,
and exposes them, and its rate constants, as keyword arguments.
Species shared with other modules are linked on instantiation:

>>> class Model(System):
...     momp = MOMP11e()
...     rec_to_bid = ReceptorToBid(Bid_U=momp.Bid_U, Bid_T=momp.Bid_T)

Initial conditions and rate constants are changed in the same way:

>>> class Model(System):
...     momp = MOMP11e(Smac_M=1e5)
...     rec_to_bid = ReceptorToBid(
...         Bid_U=momp.Bid_U, Bid_T=momp.Bid_T, bid_cleavage_kf=1e-7
...     )

The modules follow the PySB implementation in :mod:`earm.albeck_modules`
and :mod:`caspase_model.modules`,
using the same default rate constants.
"""

from simbio import Parameter, System, Variable, assign, initial
from simbio.reactions import (
    CatalyzeConvert,
    Conversion,
    Equilibration,
    MichaelisMenten,
    ReversibleSynthesis,
)

# Default forward, reverse, and catalytic rates
KF = 1e-6
KR = 1e-3
KC = 1

# Rate scaling for reactions in the mitochondrial membrane,
# the inverse of its fractional volume.
MITO_RATE_SCALING = 1 / 0.07

# Note on homo-oligomerization rate constants:
# simbio's mass action for 2 A -> A2 with rate k gives d[A]/dt = -2 k [A]^2.
# BioNetGen applies symmetry factors to these reactions,
# which EARM compensates for in its PySB parameter values
# (see earm.albeck_modules.Bax_tetramerizes).
# The rate constants here are the effective ones,
# matching the original Albeck et al. (2008) ODEs.


class ReceptorToBid(System):
    """Receptor ligation to Bid truncation.

    - L + R <-> L:R --> DISC
    - DISC + C8_pro <-> DISC:C8_pro --> DISC + C8_A
    - C8_A + Bid_U <-> C8_A:Bid_U --> C8_A + Bid_T
    - DISC + flip <-> DISC:flip
    - BAR + C8_A <-> BAR:C8_A
    """

    L: Variable = initial(default=3000)  # 50 ng/ml SK-TRAIL
    R: Variable = initial(default=200)
    DISC: Variable = initial(default=0)
    flip: Variable = initial(default=1e2)
    C8_pro: Variable = initial(default=2e4)
    C8_A: Variable = initial(default=0)
    BAR: Variable = initial(default=1e3)
    Bid_U: Variable = initial(default=0)
    Bid_T: Variable = initial(default=0)

    L_R = initial(default=0)
    DISC_C8_pro = initial(default=0)
    C8_A_Bid_U = initial(default=0)
    DISC_flip = initial(default=0)
    BAR_C8_A = initial(default=0)

    disc_formation_kf: Parameter = assign(default=4e-7)
    disc_formation_kr: Parameter = assign(default=KR)
    disc_formation_kc: Parameter = assign(default=1e-5)
    c8_activation_kf: Parameter = assign(default=KF)
    c8_activation_kr: Parameter = assign(default=KR)
    c8_activation_kc: Parameter = assign(default=KC)
    bid_cleavage_kf: Parameter = assign(default=KF)
    bid_cleavage_kr: Parameter = assign(default=KR)
    bid_cleavage_kc: Parameter = assign(default=KC)
    flip_binding_kf: Parameter = assign(default=KF)
    flip_binding_kr: Parameter = assign(default=KR)
    bar_binding_kf: Parameter = assign(default=KF)
    bar_binding_kr: Parameter = assign(default=KR)

    disc_formation = CatalyzeConvert(
        A=L,
        B=R,
        AB=L_R,
        P=DISC,
        forward_rate=disc_formation_kf,
        reverse_rate=disc_formation_kr,
        conversion_rate=disc_formation_kc,
    )
    c8_activation = MichaelisMenten(
        E=DISC,
        S=C8_pro,
        ES=DISC_C8_pro,
        P=C8_A,
        forward_rate=c8_activation_kf,
        reverse_rate=c8_activation_kr,
        catalytic_rate=c8_activation_kc,
    )
    bid_cleavage = MichaelisMenten(
        E=C8_A,
        S=Bid_U,
        ES=C8_A_Bid_U,
        P=Bid_T,
        forward_rate=bid_cleavage_kf,
        reverse_rate=bid_cleavage_kr,
        catalytic_rate=bid_cleavage_kc,
    )
    flip_binding = ReversibleSynthesis(
        A=DISC,
        B=flip,
        AB=DISC_flip,
        forward_rate=flip_binding_kf,
        reverse_rate=flip_binding_kr,
    )
    bar_binding = ReversibleSynthesis(
        A=BAR,
        B=C8_A,
        AB=BAR_C8_A,
        forward_rate=bar_binding_kf,
        reverse_rate=bar_binding_kr,
    )


# Mitochondrial outer membrane permeabilization (MOMP)
# ====================================================
#
# Building blocks shared by the MOMP variants of Albeck et al. (2008), Figure 11.


class BaxOligomerization(System):
    """Bax dimerization and tetramerization.

    - 2 monomer <-> dimer
    - 2 dimer <-> tetramer
    """

    monomer: Variable = initial(default=0)
    dimer: Variable = initial(default=0)
    tetramer: Variable = initial(default=0)

    dimerization_kf: Parameter = assign(default=KF)
    dimerization_kr: Parameter = assign(default=KR)
    tetramerization_kf: Parameter = assign(default=KF)
    tetramerization_kr: Parameter = assign(default=KR)

    dimerization = Equilibration(
        A=2 * monomer,
        B=dimer,
        forward_rate=dimerization_kf,
        reverse_rate=dimerization_kr,
    )
    tetramerization = Equilibration(
        A=2 * dimer,
        B=tetramer,
        forward_rate=tetramerization_kf,
        reverse_rate=tetramerization_kr,
    )


class Bcl2Inhibition(System):
    """Bcl2 binds and inhibits Bax monomers, dimers and tetramers.

    - Bax + Bcl2 <-> Bax:Bcl2, for Bax in {monomer, dimer, tetramer}
    """

    Bcl2: Variable = initial(default=2e4)
    monomer: Variable = initial(default=0)
    dimer: Variable = initial(default=0)
    tetramer: Variable = initial(default=0)

    monomer_Bcl2 = initial(default=0)
    dimer_Bcl2 = initial(default=0)
    tetramer_Bcl2 = initial(default=0)

    kf: Parameter = assign(default=KF)
    kr: Parameter = assign(default=KR)

    monomer_binding = ReversibleSynthesis(
        A=monomer, B=Bcl2, AB=monomer_Bcl2, forward_rate=kf, reverse_rate=kr
    )
    dimer_binding = ReversibleSynthesis(
        A=dimer, B=Bcl2, AB=dimer_Bcl2, forward_rate=kf, reverse_rate=kr
    )
    tetramer_binding = ReversibleSynthesis(
        A=tetramer, B=Bcl2, AB=tetramer_Bcl2, forward_rate=kf, reverse_rate=kr
    )


class PoreTransport(System):
    """A pore transports Smac and CytoC from the mitochondria to the cytosol.

    - pore + Smac_M <-> pore:Smac_M --> pore + Smac_C
    - pore + CytoC_M <-> pore:CytoC_M --> pore + CytoC_C
    """

    pore: Variable = initial(default=0)
    Smac_M: Variable = initial(default=1e6)
    Smac_C: Variable = initial(default=0)
    CytoC_M: Variable = initial(default=1e6)
    CytoC_C: Variable = initial(default=0)

    pore_Smac_M = initial(default=0)
    pore_CytoC_M = initial(default=0)

    smac_kf: Parameter = assign(default=KF)
    smac_kr: Parameter = assign(default=KR)
    smac_kc: Parameter = assign(default=10)
    cytoc_kf: Parameter = assign(default=KF)
    cytoc_kr: Parameter = assign(default=KR)
    cytoc_kc: Parameter = assign(default=10)

    smac_transport = MichaelisMenten(
        E=pore,
        S=Smac_M,
        ES=pore_Smac_M,
        P=Smac_C,
        forward_rate=smac_kf,
        reverse_rate=smac_kr,
        catalytic_rate=smac_kc,
    )
    cytoc_transport = MichaelisMenten(
        E=pore,
        S=CytoC_M,
        ES=pore_CytoC_M,
        P=CytoC_C,
        forward_rate=cytoc_kf,
        reverse_rate=cytoc_kr,
        catalytic_rate=cytoc_kc,
    )


class MOMP11b(System):
    """Minimal MOMP model (Albeck et al. 2008, Figure 11b).

    - tBid activates Bax
    - Active Bax is inhibited by Bcl2
    - Active Bax transports Smac and CytoC to the cytosol
    """

    Bid_U: Variable = initial(default=1e5)
    Bid_T: Variable = initial(default=0)
    Bax_C: Variable = initial(default=1e5)
    Bax_A: Variable = initial(default=0)
    Bcl2: Variable = initial(default=2e4)
    Smac_M: Variable = initial(default=1e6)
    Smac_C: Variable = initial(default=0)
    CytoC_M: Variable = initial(default=1e6)
    CytoC_C: Variable = initial(default=0)

    Bid_T_Bax_C = initial(default=0)
    Bax_A_Bcl2 = initial(default=0)

    bax_activation_kf: Parameter = assign(default=1e-7)
    bax_activation_kr: Parameter = assign(default=KR)
    bax_activation_kc: Parameter = assign(default=KC)
    bcl2_binding_kf: Parameter = assign(default=KF)
    bcl2_binding_kr: Parameter = assign(default=KR)

    bax_activation = MichaelisMenten(
        E=Bid_T,
        S=Bax_C,
        ES=Bid_T_Bax_C,
        P=Bax_A,
        forward_rate=bax_activation_kf,
        reverse_rate=bax_activation_kr,
        catalytic_rate=bax_activation_kc,
    )
    bcl2_binding = ReversibleSynthesis(
        A=Bax_A,
        B=Bcl2,
        AB=Bax_A_Bcl2,
        forward_rate=bcl2_binding_kf,
        reverse_rate=bcl2_binding_kr,
    )
    transport = PoreTransport(
        pore=Bax_A, Smac_M=Smac_M, Smac_C=Smac_C, CytoC_M=CytoC_M, CytoC_C=CytoC_C
    )


class MOMP11c(System):
    """MOMP model with Bax oligomerization (Albeck et al. 2008, Figure 11c).

    - tBid activates Bax
    - Active Bax dimerizes, and dimers dimerize into tetramers
    - Bcl2 inhibits Bax monomers, dimers and tetramers
    - Bax tetramers transport Smac and CytoC to the cytosol
    """

    Bid_U: Variable = initial(default=4e4)
    Bid_T: Variable = initial(default=0)
    Bax_C: Variable = initial(default=1e5)
    Bax_A: Variable = initial(default=0)
    Bax_A2: Variable = initial(default=0)
    Bax_A4: Variable = initial(default=0)
    Bcl2: Variable = initial(default=2e4)
    Smac_M: Variable = initial(default=1e6)
    Smac_C: Variable = initial(default=0)
    CytoC_M: Variable = initial(default=1e6)
    CytoC_C: Variable = initial(default=0)

    Bid_T_Bax_C = initial(default=0)

    bax_activation_kf: Parameter = assign(default=1e-7)
    bax_activation_kr: Parameter = assign(default=KR)
    bax_activation_kc: Parameter = assign(default=KC)

    bax_activation = MichaelisMenten(
        E=Bid_T,
        S=Bax_C,
        ES=Bid_T_Bax_C,
        P=Bax_A,
        forward_rate=bax_activation_kf,
        reverse_rate=bax_activation_kr,
        catalytic_rate=bax_activation_kc,
    )
    oligomerization = BaxOligomerization(monomer=Bax_A, dimer=Bax_A2, tetramer=Bax_A4)
    bcl2 = Bcl2Inhibition(Bcl2=Bcl2, monomer=Bax_A, dimer=Bax_A2, tetramer=Bax_A4)
    transport = PoreTransport(
        pore=Bax_A4,
        Smac_M=Smac_M,
        Smac_C=Smac_C,
        CytoC_M=CytoC_M,
        CytoC_C=CytoC_C,
        smac_kf=2 * KF,
    )


class MOMP11d(System):
    """MOMP model with mitochondrial transport (Albeck et al. 2008, Figure 11d).

    - tBid activates Bax
    - Active Bax translocates to the mitochondria
    - Reactions in the mitochondrial membrane have increased association rates
    - Mitochondrial Bax dimerizes, and dimers dimerize into tetramers
    - Bcl2 inhibits Bax monomers, dimers and tetramers
    - Bax tetramers transport Smac and CytoC to the cytosol
    """

    Bid_U: Variable = initial(default=4e4)
    Bid_T: Variable = initial(default=0)
    Bax_C: Variable = initial(default=1e5)
    Bax_A: Variable = initial(default=0)
    Bax_M: Variable = initial(default=0)
    Bax_M2: Variable = initial(default=0)
    Bax_M4: Variable = initial(default=0)
    Bcl2: Variable = initial(default=2e4)
    Smac_M: Variable = initial(default=1e6)
    Smac_C: Variable = initial(default=0)
    CytoC_M: Variable = initial(default=1e6)
    CytoC_C: Variable = initial(default=0)

    Bid_T_Bax_C = initial(default=0)

    bax_activation_kf: Parameter = assign(default=1e-7)
    bax_activation_kr: Parameter = assign(default=KR)
    bax_activation_kc: Parameter = assign(default=KC)
    bax_translocation_kf: Parameter = assign(default=1e-2)
    bax_translocation_kr: Parameter = assign(default=1e-2)

    bax_activation = MichaelisMenten(
        E=Bid_T,
        S=Bax_C,
        ES=Bid_T_Bax_C,
        P=Bax_A,
        forward_rate=bax_activation_kf,
        reverse_rate=bax_activation_kr,
        catalytic_rate=bax_activation_kc,
    )
    bax_translocation = Equilibration(
        A=Bax_A,
        B=Bax_M,
        forward_rate=bax_translocation_kf,
        reverse_rate=bax_translocation_kr,
    )
    oligomerization = BaxOligomerization(
        monomer=Bax_M,
        dimer=Bax_M2,
        tetramer=Bax_M4,
        dimerization_kf=KF * MITO_RATE_SCALING,
        tetramerization_kf=KF * MITO_RATE_SCALING,
    )
    bcl2 = Bcl2Inhibition(
        Bcl2=Bcl2,
        monomer=Bax_M,
        dimer=Bax_M2,
        tetramer=Bax_M4,
        kf=KF * MITO_RATE_SCALING,
    )
    transport = PoreTransport(
        pore=Bax_M4,
        Smac_M=Smac_M,
        Smac_C=Smac_C,
        CytoC_M=CytoC_M,
        CytoC_C=CytoC_C,
        smac_kf=2 * KF * MITO_RATE_SCALING,
    )


class MOMP11e(System):
    """MOMP model with mitochondrial pore insertion (Albeck et al. 2008, Figure 11e).

    As :class:`MOMP11d`, but Bax tetramers bind to mitochondrial sites (Mito),
    which become active pores that transport Smac and CytoC to the cytosol.

    The model with cooperative pore assembly (Albeck et al. 2008, Figure 11f)
    is obtained by changing the rate constants,
    see :data:`MOMP11F_PARAMETERS`.
    """

    Bid_U: Variable = initial(default=4e4)
    Bid_T: Variable = initial(default=0)
    Bax_C: Variable = initial(default=1e5)
    Bax_A: Variable = initial(default=0)
    Bax_M: Variable = initial(default=0)
    Bax_M2: Variable = initial(default=0)
    Bax_M4: Variable = initial(default=0)
    Bcl2: Variable = initial(default=2e4)
    Mito_I: Variable = initial(default=5e5)
    Mito_A: Variable = initial(default=0)
    Smac_M: Variable = initial(default=1e6)
    Smac_C: Variable = initial(default=0)
    CytoC_M: Variable = initial(default=1e6)
    CytoC_C: Variable = initial(default=0)

    Bid_T_Bax_C = initial(default=0)
    Bax_M4_Mito_I = initial(default=0)

    bax_activation_kf: Parameter = assign(default=1e-7)
    bax_activation_kr: Parameter = assign(default=KR)
    bax_activation_kc: Parameter = assign(default=KC)
    bax_translocation_kf: Parameter = assign(default=1e-2)
    bax_translocation_kr: Parameter = assign(default=1e-2)
    bax_dimerization_kf: Parameter = assign(default=KF * MITO_RATE_SCALING)
    bax_tetramerization_kf: Parameter = assign(default=KF * MITO_RATE_SCALING)
    mito_binding_kf: Parameter = assign(default=KF * MITO_RATE_SCALING)
    mito_binding_kr: Parameter = assign(default=KR)
    mito_activation_kc: Parameter = assign(default=KC)

    bax_activation = MichaelisMenten(
        E=Bid_T,
        S=Bax_C,
        ES=Bid_T_Bax_C,
        P=Bax_A,
        forward_rate=bax_activation_kf,
        reverse_rate=bax_activation_kr,
        catalytic_rate=bax_activation_kc,
    )
    bax_translocation = Equilibration(
        A=Bax_A,
        B=Bax_M,
        forward_rate=bax_translocation_kf,
        reverse_rate=bax_translocation_kr,
    )
    oligomerization = BaxOligomerization(
        monomer=Bax_M,
        dimer=Bax_M2,
        tetramer=Bax_M4,
        dimerization_kf=bax_dimerization_kf,
        tetramerization_kf=bax_tetramerization_kf,
    )
    bcl2 = Bcl2Inhibition(
        Bcl2=Bcl2,
        monomer=Bax_M,
        dimer=Bax_M2,
        tetramer=Bax_M4,
        kf=KF * MITO_RATE_SCALING,
    )
    # Bax4 + Mito_I <-> Bax4:Mito_I --> Mito_A
    mito_binding = ReversibleSynthesis(
        A=Bax_M4,
        B=Mito_I,
        AB=Bax_M4_Mito_I,
        forward_rate=mito_binding_kf,
        reverse_rate=mito_binding_kr,
    )
    mito_activation = Conversion(A=Bax_M4_Mito_I, B=Mito_A, rate=mito_activation_kc)
    transport = PoreTransport(
        pore=Mito_A,
        Smac_M=Smac_M,
        Smac_C=Smac_C,
        CytoC_M=CytoC_M,
        CytoC_C=CytoC_C,
        smac_kf=2 * KF * MITO_RATE_SCALING,
        cytoc_kf=2 * KF * MITO_RATE_SCALING,
    )


MOMP11F_PARAMETERS = {
    "bax_translocation_kf": 1e-4,
    "bax_translocation_kr": 1e-4,
    "bax_dimerization_kf": KF * MITO_RATE_SCALING / 100,
    "bax_tetramerization_kf": KF * MITO_RATE_SCALING / 10,
}
"""Rate constants of :class:`MOMP11e` for cooperative pore assembly (Figure 11f)."""


# Pore to PARP
# ============


class CaspaseCascade(System):
    """Effector caspase activation and PARP cleavage.

    Released Smac and CytoC are activated in the cytosol,
    and Smac sequesters XIAP:

    - Smac_C <-> Smac_A
    - CytoC_C <-> CytoC_A
    - Smac_A + XIAP <-> Smac_A:XIAP

    Caspase reactions, including effectors, inhibitors and feedback initiators:

    - C8_A + C3_pro <-> C8_A:C3_pro --> C8_A + C3_A
    - XIAP + C3_A <-> XIAP:C3_A --> XIAP + C3_ub
    - C3_A + PARP_U <-> C3_A:PARP_U --> C3_A + PARP_C
    - C3_A + C6_pro <-> C3_A:C6_pro --> C3_A + C6_A
    - C6_A + C8_pro <-> C6_A:C8_pro --> C6_A + C8_A
    """

    Smac_C: Variable = initial(default=0)
    Smac_A: Variable = initial(default=0)
    CytoC_C: Variable = initial(default=0)
    CytoC_A: Variable = initial(default=0)
    C8_pro: Variable = initial(default=0)
    C8_A: Variable = initial(default=0)
    C3_pro: Variable = initial(default=1e4)
    C3_A: Variable = initial(default=0)
    C3_ub: Variable = initial(default=0)
    C6_pro: Variable = initial(default=1e4)
    C6_A: Variable = initial(default=0)
    PARP_U: Variable = initial(default=1e6)
    PARP_C: Variable = initial(default=0)
    XIAP: Variable = initial(default=1e5)

    Smac_A_XIAP = initial(default=0)
    C8_A_C3_pro = initial(default=0)
    XIAP_C3_A = initial(default=0)
    C3_A_PARP_U = initial(default=0)
    C3_A_C6_pro = initial(default=0)
    C6_A_C8_pro = initial(default=0)

    translocation_kf: Parameter = assign(default=1e-2)
    translocation_kr: Parameter = assign(default=1e-2)
    smac_xiap_kf: Parameter = assign(default=7e-6)
    smac_xiap_kr: Parameter = assign(default=KR)
    c3_activation_kf: Parameter = assign(default=1e-7)
    c3_activation_kr: Parameter = assign(default=KR)
    c3_activation_kc: Parameter = assign(default=KC)
    c3_ubiquitination_kf: Parameter = assign(default=2e-6)
    c3_ubiquitination_kr: Parameter = assign(default=KR)
    c3_ubiquitination_kc: Parameter = assign(default=1e-1)
    parp_cleavage_kf: Parameter = assign(default=KF)
    parp_cleavage_kr: Parameter = assign(default=1e-2)
    parp_cleavage_kc: Parameter = assign(default=KC)
    c6_activation_kf: Parameter = assign(default=KF)
    c6_activation_kr: Parameter = assign(default=KR)
    c6_activation_kc: Parameter = assign(default=KC)
    c8_activation_kf: Parameter = assign(default=3e-8)
    c8_activation_kr: Parameter = assign(default=KR)
    c8_activation_kc: Parameter = assign(default=KC)

    smac_activation = Equilibration(
        A=Smac_C,
        B=Smac_A,
        forward_rate=translocation_kf,
        reverse_rate=translocation_kr,
    )
    cytoc_activation = Equilibration(
        A=CytoC_C,
        B=CytoC_A,
        forward_rate=translocation_kf,
        reverse_rate=translocation_kr,
    )
    smac_xiap = ReversibleSynthesis(
        A=Smac_A,
        B=XIAP,
        AB=Smac_A_XIAP,
        forward_rate=smac_xiap_kf,
        reverse_rate=smac_xiap_kr,
    )
    c3_activation = MichaelisMenten(
        E=C8_A,
        S=C3_pro,
        ES=C8_A_C3_pro,
        P=C3_A,
        forward_rate=c3_activation_kf,
        reverse_rate=c3_activation_kr,
        catalytic_rate=c3_activation_kc,
    )
    c3_ubiquitination = MichaelisMenten(
        E=XIAP,
        S=C3_A,
        ES=XIAP_C3_A,
        P=C3_ub,
        forward_rate=c3_ubiquitination_kf,
        reverse_rate=c3_ubiquitination_kr,
        catalytic_rate=c3_ubiquitination_kc,
    )
    parp_cleavage = MichaelisMenten(
        E=C3_A,
        S=PARP_U,
        ES=C3_A_PARP_U,
        P=PARP_C,
        forward_rate=parp_cleavage_kf,
        reverse_rate=parp_cleavage_kr,
        catalytic_rate=parp_cleavage_kc,
    )
    c6_activation = MichaelisMenten(
        E=C3_A,
        S=C6_pro,
        ES=C3_A_C6_pro,
        P=C6_A,
        forward_rate=c6_activation_kf,
        reverse_rate=c6_activation_kr,
        catalytic_rate=c6_activation_kc,
    )
    c8_activation = MichaelisMenten(
        E=C6_A,
        S=C8_pro,
        ES=C6_A_C8_pro,
        P=C8_A,
        forward_rate=c8_activation_kf,
        reverse_rate=c8_activation_kr,
        catalytic_rate=c8_activation_kc,
    )


class Apoptosome(System):
    """Apoptosome formation and inhibition (Albeck et al. 2008).

    - CytoC_A + Apaf_I <-> CytoC_A:Apaf_I --> CytoC_A + Apaf_A
    - Apaf_A + C9 <-> Apop
    - Apop + C3_pro <-> Apop:C3_pro --> Apop + C3_A
    - Apop + XIAP <-> Apop:XIAP
    """

    CytoC_A: Variable = initial(default=0)
    Apaf_I: Variable = initial(default=1e5)
    Apaf_A: Variable = initial(default=0)
    C9: Variable = initial(default=1e5)
    Apop: Variable = initial(default=0)
    C3_pro: Variable = initial(default=0)
    C3_A: Variable = initial(default=0)
    XIAP: Variable = initial(default=0)

    CytoC_A_Apaf_I = initial(default=0)
    Apop_C3_pro = initial(default=0)
    Apop_XIAP = initial(default=0)

    apaf_activation_kf: Parameter = assign(default=5e-7)
    apaf_activation_kr: Parameter = assign(default=KR)
    apaf_activation_kc: Parameter = assign(default=KC)
    apop_formation_kf: Parameter = assign(default=5e-8)
    apop_formation_kr: Parameter = assign(default=KR)
    c3_activation_kf: Parameter = assign(default=5e-9)
    c3_activation_kr: Parameter = assign(default=KR)
    c3_activation_kc: Parameter = assign(default=KC)
    xiap_binding_kf: Parameter = assign(default=2e-6)
    xiap_binding_kr: Parameter = assign(default=KR)

    apaf_activation = MichaelisMenten(
        E=CytoC_A,
        S=Apaf_I,
        ES=CytoC_A_Apaf_I,
        P=Apaf_A,
        forward_rate=apaf_activation_kf,
        reverse_rate=apaf_activation_kr,
        catalytic_rate=apaf_activation_kc,
    )
    apop_formation = ReversibleSynthesis(
        A=Apaf_A,
        B=C9,
        AB=Apop,
        forward_rate=apop_formation_kf,
        reverse_rate=apop_formation_kr,
    )
    c3_activation = MichaelisMenten(
        E=Apop,
        S=C3_pro,
        ES=Apop_C3_pro,
        P=C3_A,
        forward_rate=c3_activation_kf,
        reverse_rate=c3_activation_kr,
        catalytic_rate=c3_activation_kc,
    )
    xiap_binding = ReversibleSynthesis(
        A=Apop,
        B=XIAP,
        AB=Apop_XIAP,
        forward_rate=xiap_binding_kf,
        reverse_rate=xiap_binding_kr,
    )


class DoubleApoptosome(System):
    """Apoptosome formation with two-step activation (ARM).

    Active Apaf already activates caspase-3,
    and active caspase-3 further activates Apaf into the apoptosome.

    - CytoC_A + Apaf_I <-> CytoC_A:Apaf_I --> CytoC_A + Apaf_A
    - Apaf_A + C3_pro <-> Apaf_A:C3_pro --> Apaf_A + C3_A
    - C3_A + Apaf_A <-> C3_A:Apaf_A --> C3_A + Apop
    - Apop + C3_pro <-> Apop:C3_pro --> Apop + C3_A
    - Apaf_A + XIAP <-> Apaf_A:XIAP
    """

    CytoC_A: Variable = initial(default=0)
    Apaf_I: Variable = initial(default=1e3)
    Apaf_A: Variable = initial(default=0)
    Apop: Variable = initial(default=0)
    C3_pro: Variable = initial(default=0)
    C3_A: Variable = initial(default=0)
    XIAP: Variable = initial(default=0)

    CytoC_A_Apaf_I = initial(default=0)
    Apaf_A_C3_pro = initial(default=0)
    C3_A_Apaf_A = initial(default=0)
    Apop_C3_pro = initial(default=0)
    Apaf_A_XIAP = initial(default=0)

    apaf_activation_kf: Parameter = assign(default=5e-7)
    apaf_activation_kr: Parameter = assign(default=KR)
    apaf_activation_kc: Parameter = assign(default=KC)
    apaf_c3_activation_kf: Parameter = assign(default=5e-9)
    apaf_c3_activation_kr: Parameter = assign(default=KR)
    apaf_c3_activation_kc: Parameter = assign(default=KC)
    apop_formation_kf: Parameter = assign(default=1.3e-6)
    apop_formation_kr: Parameter = assign(default=KR)
    apop_formation_kc: Parameter = assign(default=KC)
    c3_activation_kf: Parameter = assign(default=5e-9)
    c3_activation_kr: Parameter = assign(default=KR)
    c3_activation_kc: Parameter = assign(default=KC)
    xiap_binding_kf: Parameter = assign(default=2e-6)
    xiap_binding_kr: Parameter = assign(default=KR)

    apaf_activation = MichaelisMenten(
        E=CytoC_A,
        S=Apaf_I,
        ES=CytoC_A_Apaf_I,
        P=Apaf_A,
        forward_rate=apaf_activation_kf,
        reverse_rate=apaf_activation_kr,
        catalytic_rate=apaf_activation_kc,
    )
    apaf_c3_activation = MichaelisMenten(
        E=Apaf_A,
        S=C3_pro,
        ES=Apaf_A_C3_pro,
        P=C3_A,
        forward_rate=apaf_c3_activation_kf,
        reverse_rate=apaf_c3_activation_kr,
        catalytic_rate=apaf_c3_activation_kc,
    )
    apop_formation = MichaelisMenten(
        E=C3_A,
        S=Apaf_A,
        ES=C3_A_Apaf_A,
        P=Apop,
        forward_rate=apop_formation_kf,
        reverse_rate=apop_formation_kr,
        catalytic_rate=apop_formation_kc,
    )
    c3_activation = MichaelisMenten(
        E=Apop,
        S=C3_pro,
        ES=Apop_C3_pro,
        P=C3_A,
        forward_rate=c3_activation_kf,
        reverse_rate=c3_activation_kr,
        catalytic_rate=c3_activation_kc,
    )
    xiap_binding = ReversibleSynthesis(
        A=Apaf_A,
        B=XIAP,
        AB=Apaf_A_XIAP,
        forward_rate=xiap_binding_kf,
        reverse_rate=xiap_binding_kr,
    )


# Additions by Corbat et al. (2018)
# =================================


class CytosolicBcl2(System):
    """Cytosolic Bcl2 sequesters tBid, as in the original Albeck et al. (2008) model.

    - Bid_T + Bcl2c <-> Bid_T:Bcl2c
    """

    Bid_T: Variable = initial(default=0)
    Bcl2c: Variable = initial(default=2e4)

    Bid_T_Bcl2c = initial(default=0)

    kf: Parameter = assign(default=KF)
    kr: Parameter = assign(default=KR)

    binding = ReversibleSynthesis(
        A=Bid_T, B=Bcl2c, AB=Bid_T_Bcl2c, forward_rate=kf, reverse_rate=kr
    )


class IntrinsicStimulus(System):
    """Intrinsic stimulus that truncates Bid.

    - IntrinsicStimuli + Bid_U <-> IntrinsicStimuli:Bid_U --> IntrinsicStimuli + Bid_T
    """

    IntrinsicStimuli: Variable = initial(default=1e2)
    Bid_U: Variable = initial(default=0)
    Bid_T: Variable = initial(default=0)

    IntrinsicStimuli_Bid_U = initial(default=0)

    kf: Parameter = assign(default=KF)
    kr: Parameter = assign(default=KR)
    kc: Parameter = assign(default=KC)

    bid_cleavage = MichaelisMenten(
        E=IntrinsicStimuli,
        S=Bid_U,
        ES=IntrinsicStimuli_Bid_U,
        P=Bid_T,
        forward_rate=kf,
        reverse_rate=kr,
        catalytic_rate=kc,
    )


class DimerCleavage(System):
    """An enzyme cleaves a homodimeric biosensor into its monomers.

    - enzyme + dimer <-> enzyme:dimer --> enzyme + 2 monomer

    The dimer has two cleavage sites,
    so the effective association rate is twice the per-site rate ``kf``.
    """

    enzyme: Variable = initial(default=0)
    dimer: Variable = initial(default=0)
    monomer: Variable = initial(default=0)

    enzyme_dimer = initial(default=0)

    kf: Parameter = assign()
    kr: Parameter = assign()
    kc: Parameter = assign()

    cleavage = MichaelisMenten(
        E=enzyme,
        S=dimer,
        ES=enzyme_dimer,
        P=2 * monomer,
        forward_rate=2 * kf,
        reverse_rate=kr,
        catalytic_rate=kc,
    )


class Biosensor(System):
    """A homo-FRET biosensor (Corbat et al. 2018).

    A homodimer cleaved into monomers by its target enzyme,
    see :class:`DimerCleavage`.
    """

    enzyme: Variable = initial(default=0)
    dimer: Variable = initial(default=1e7)
    monomer: Variable = initial(default=0)

    kf: Parameter = assign()
    kr: Parameter = assign()
    kc: Parameter = assign(default=KC)

    cleavage = DimerCleavage(
        enzyme=enzyme, dimer=dimer, monomer=monomer, kf=kf, kr=kr, kc=kc
    )
