"""Apoptotic network models, composed from :mod:`.modules`.

SimBio counterparts of the PySB models:

- :class:`Albeck11b` to :class:`Albeck11f`: EARM 1.0 with each MOMP variant
  (:mod:`earm.albeck_modules`, with pore transport).
- :class:`AlbeckAsMatlab`: :func:`caspase_model.models.albeck_as_matlab`.
- :class:`Corbat2018`: :func:`caspase_model.models.corbat_2018`.
- :class:`ARM`: :func:`caspase_model.models.arm`.

The ligand ``L`` and the ``IntrinsicStimuli`` are top-level species,
so that the stimulus can be chosen on instantiation.
The models are extrinsic by default,
and :func:`intrinsic` returns the intrinsically stimulated model:

>>> from simbio import Simulator
>>> sim = Simulator(intrinsic(ARM))
"""

from simbio import System, Variable, initial

from .modules import (
    MOMP11F_PARAMETERS,
    Apoptosome,
    Biosensor,
    CaspaseCascade,
    CytosolicBcl2,
    DimerCleavage,
    DoubleApoptosome,
    IntrinsicStimulus,
    MOMP11b,
    MOMP11c,
    MOMP11d,
    MOMP11e,
    ReceptorToBid,
)


def intrinsic(model: type[System], stimulus: float = 1e2) -> System:
    """Replace the extrinsic stimulus (ligand) by an intrinsic one."""
    return model(L=0, IntrinsicStimuli=stimulus)


# EARM 1.0 (Albeck et al. 2008)
# =============================
#
#   ReceptorToBid -> MOMP -> CaspaseCascade <-> Apoptosome


class Albeck11b(System):
    """EARM 1.0 with the minimal MOMP model (Figure 11b)."""

    L: Variable = initial(default=3000)

    momp = MOMP11b()
    rec_to_bid = ReceptorToBid(L=L, Bid_U=momp.Bid_U, Bid_T=momp.Bid_T)
    caspases = CaspaseCascade(
        Smac_C=momp.Smac_C,
        CytoC_C=momp.CytoC_C,
        C8_pro=rec_to_bid.C8_pro,
        C8_A=rec_to_bid.C8_A,
    )
    apoptosome = Apoptosome(
        CytoC_A=caspases.CytoC_A,
        C3_pro=caspases.C3_pro,
        C3_A=caspases.C3_A,
        XIAP=caspases.XIAP,
    )


class Albeck11c(System):
    """EARM 1.0 with Bax oligomerization (Figure 11c)."""

    L: Variable = initial(default=3000)

    momp = MOMP11c()
    rec_to_bid = ReceptorToBid(L=L, Bid_U=momp.Bid_U, Bid_T=momp.Bid_T)
    caspases = CaspaseCascade(
        Smac_C=momp.Smac_C,
        CytoC_C=momp.CytoC_C,
        C8_pro=rec_to_bid.C8_pro,
        C8_A=rec_to_bid.C8_A,
    )
    apoptosome = Apoptosome(
        CytoC_A=caspases.CytoC_A,
        C3_pro=caspases.C3_pro,
        C3_A=caspases.C3_A,
        XIAP=caspases.XIAP,
    )


class Albeck11d(System):
    """EARM 1.0 with mitochondrial transport (Figure 11d)."""

    L: Variable = initial(default=3000)

    momp = MOMP11d()
    rec_to_bid = ReceptorToBid(L=L, Bid_U=momp.Bid_U, Bid_T=momp.Bid_T)
    caspases = CaspaseCascade(
        Smac_C=momp.Smac_C,
        CytoC_C=momp.CytoC_C,
        C8_pro=rec_to_bid.C8_pro,
        C8_A=rec_to_bid.C8_A,
    )
    apoptosome = Apoptosome(
        CytoC_A=caspases.CytoC_A,
        C3_pro=caspases.C3_pro,
        C3_A=caspases.C3_A,
        XIAP=caspases.XIAP,
    )


class Albeck11e(System):
    """EARM 1.0 with mitochondrial pore insertion (Figure 11e)."""

    L: Variable = initial(default=3000)

    momp = MOMP11e()
    rec_to_bid = ReceptorToBid(L=L, Bid_U=momp.Bid_U, Bid_T=momp.Bid_T)
    caspases = CaspaseCascade(
        Smac_C=momp.Smac_C,
        CytoC_C=momp.CytoC_C,
        C8_pro=rec_to_bid.C8_pro,
        C8_A=rec_to_bid.C8_A,
    )
    apoptosome = Apoptosome(
        CytoC_A=caspases.CytoC_A,
        C3_pro=caspases.C3_pro,
        C3_A=caspases.C3_A,
        XIAP=caspases.XIAP,
    )


class Albeck11f(System):
    """EARM 1.0 with cooperative pore assembly (Figure 11f)."""

    L: Variable = initial(default=3000)

    momp = MOMP11e(**MOMP11F_PARAMETERS)
    rec_to_bid = ReceptorToBid(L=L, Bid_U=momp.Bid_U, Bid_T=momp.Bid_T)
    caspases = CaspaseCascade(
        Smac_C=momp.Smac_C,
        CytoC_C=momp.CytoC_C,
        C8_pro=rec_to_bid.C8_pro,
        C8_A=rec_to_bid.C8_A,
    )
    apoptosome = Apoptosome(
        CytoC_A=caspases.CytoC_A,
        C3_pro=caspases.C3_pro,
        C3_A=caspases.C3_A,
        XIAP=caspases.XIAP,
    )


# Corbat et al. (2018)
# ====================


class AlbeckAsMatlab(System):
    """EARM 1.0 (Figure 11e) as in Albeck's original MATLAB implementation.

    Adds cytosolic Bcl2, which is absent in the EARM implementation,
    and the original C8-Bid association rate.
    """

    L: Variable = initial(default=3000)

    momp = MOMP11e()
    rec_to_bid = ReceptorToBid(
        L=L, Bid_U=momp.Bid_U, Bid_T=momp.Bid_T, bid_cleavage_kf=1e-7
    )
    caspases = CaspaseCascade(
        Smac_C=momp.Smac_C,
        CytoC_C=momp.CytoC_C,
        C8_pro=rec_to_bid.C8_pro,
        C8_A=rec_to_bid.C8_A,
    )
    apoptosome = Apoptosome(
        CytoC_A=caspases.CytoC_A,
        C3_pro=caspases.C3_pro,
        C3_A=caspases.C3_A,
        XIAP=caspases.XIAP,
    )
    bcl2c = CytosolicBcl2(Bid_T=momp.Bid_T)


class Corbat2018(System):
    """:class:`AlbeckAsMatlab` with caspase biosensors (Corbat et al. 2018)."""

    L: Variable = initial(default=1e3)
    IntrinsicStimuli: Variable = initial(default=0)

    momp = MOMP11e()
    rec_to_bid = ReceptorToBid(
        L=L, R=1e3, Bid_U=momp.Bid_U, Bid_T=momp.Bid_T, bid_cleavage_kf=1e-7
    )
    caspases = CaspaseCascade(
        Smac_C=momp.Smac_C,
        CytoC_C=momp.CytoC_C,
        C8_pro=rec_to_bid.C8_pro,
        C8_A=rec_to_bid.C8_A,
        XIAP=1e2,
    )
    apoptosome = Apoptosome(
        CytoC_A=caspases.CytoC_A,
        C3_pro=caspases.C3_pro,
        C3_A=caspases.C3_A,
        XIAP=caspases.XIAP,
    )
    bcl2c = CytosolicBcl2(Bid_T=momp.Bid_T)
    stimulus = IntrinsicStimulus(
        IntrinsicStimuli=IntrinsicStimuli, Bid_U=momp.Bid_U, Bid_T=momp.Bid_T
    )

    sCas3 = Biosensor(enzyme=caspases.C3_A, kf=1e-6 / 2, kr=1e-2)
    sCas8 = Biosensor(enzyme=rec_to_bid.C8_A, kf=1e-7 / 2, kr=1e-3)
    sCas9 = Biosensor(enzyme=apoptosome.Apop, kf=5e-9 / 2, kr=1e-3)


class ARM(System):
    """Apoptotic Reaction Model.

    :class:`Corbat2018` with a two-step apoptosome activation,
    where active Apaf already activates caspase-3 and cleaves the caspase-9 biosensor.
    """

    L: Variable = initial(default=1e3)
    IntrinsicStimuli: Variable = initial(default=0)

    momp = MOMP11e(Smac_M=1e5, CytoC_M=1e5)
    rec_to_bid = ReceptorToBid(
        L=L, Bid_U=momp.Bid_U, Bid_T=momp.Bid_T, bid_cleavage_kf=1e-7
    )
    caspases = CaspaseCascade(
        Smac_C=momp.Smac_C,
        CytoC_C=momp.CytoC_C,
        C8_pro=rec_to_bid.C8_pro,
        C8_A=rec_to_bid.C8_A,
        XIAP=1e4,
    )
    apoptosome = DoubleApoptosome(
        CytoC_A=caspases.CytoC_A,
        C3_pro=caspases.C3_pro,
        C3_A=caspases.C3_A,
        XIAP=caspases.XIAP,
    )
    bcl2c = CytosolicBcl2(Bid_T=momp.Bid_T)
    stimulus = IntrinsicStimulus(
        IntrinsicStimuli=IntrinsicStimuli, Bid_U=momp.Bid_U, Bid_T=momp.Bid_T
    )

    sCas3 = Biosensor(enzyme=caspases.C3_A, dimer=7.5e5, kf=2.8e-7, kr=1e-2)
    sCas8 = Biosensor(enzyme=rec_to_bid.C8_A, dimer=7.5e5, kf=5.4e-8, kr=1e-3)
    sCas9 = Biosensor(enzyme=apoptosome.Apop, dimer=7.5e5, kf=2.8e-7, kr=1e-3)
    # Active Apaf cleaves the caspase-9 biosensor at a slower rate
    sCas9_apaf = DimerCleavage(
        enzyme=apoptosome.Apaf_A,
        dimer=sCas9.dimer,
        monomer=sCas9.monomer,
        kf=2e-10,
        kr=1e-3,
        kc=1,
    )
