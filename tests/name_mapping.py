"""Name mapping between PySB species and SimBio variables.

A PySB species can map to different SimBio variables depending on the model,
so each one lists its candidates.
"""

BAX_A4 = (
    "Bax(bf={}, s1=1, s2=2, state='A') % Bax(bf=None, s1=3, s2=1, state='A') % "
    "Bax(bf=None, s1=4, s2=3, state='A') % Bax(bf=None, s1=2, s2=4, state='A')"
)
BAX_M4 = BAX_A4.replace("'A'", "'M'")
BAX_A4_BOUND = (
    "Bax(bf=1, s1=2, s2=3, state='A') % Bax(bf=None, s1=4, s2=2, state='A') % "
    "Bax(bf=None, s1=5, s2=4, state='A') % Bax(bf=None, s1=3, s2=5, state='A') % {}"
)
BAX_M4_BOUND = BAX_A4_BOUND.replace("'A'", "'M'")

name_mapping: dict[str, tuple[str, ...]] = {
    # Stimuli
    "L(bf=None)": ("L",),
    "IntrinsicStimuli(bf=None)": ("IntrinsicStimuli",),
    "Bid(bf=1, state='U') % IntrinsicStimuli(bf=1)": (
        "stimulus.IntrinsicStimuli_Bid_U",
    ),
    # ReceptorToBid
    "R(bf=None)": ("rec_to_bid.R",),
    "L(bf=1) % R(bf=1)": ("rec_to_bid.L_R",),
    "DISC(bf=None)": ("rec_to_bid.DISC",),
    "flip(bf=None)": ("rec_to_bid.flip",),
    "DISC(bf=1) % flip(bf=1)": ("rec_to_bid.DISC_flip",),
    "C8(bf=None, state='pro')": ("rec_to_bid.C8_pro",),
    "C8(bf=None, state='A')": ("rec_to_bid.C8_A",),
    "C8(bf=1, state='pro') % DISC(bf=1)": ("rec_to_bid.DISC_C8_pro",),
    "BAR(bf=None)": ("rec_to_bid.BAR",),
    "BAR(bf=1) % C8(bf=1, state='A')": ("rec_to_bid.BAR_C8_A",),
    "Bid(bf=1, state='U') % C8(bf=1, state='A')": ("rec_to_bid.C8_A_Bid_U",),
    # MOMP
    "Bid(bf=None, state='U')": ("momp.Bid_U",),
    "Bid(bf=None, state='T')": ("momp.Bid_T",),
    "Bax(bf=None, s1=None, s2=None, state='C')": ("momp.Bax_C",),
    "Bax(bf=1, s1=None, s2=None, state='C') % Bid(bf=1, state='T')": (
        "momp.Bid_T_Bax_C",
    ),
    "Bax(bf=None, s1=None, s2=None, state='A')": ("momp.Bax_A",),
    "Bax(bf=None, s1=None, s2=None, state='M')": ("momp.Bax_M",),
    "Bax(bf=None, s1=1, s2=None, state='A') % Bax(bf=None, s1=None, s2=1, state='A')": (
        "momp.Bax_A2",
    ),
    "Bax(bf=None, s1=1, s2=None, state='M') % Bax(bf=None, s1=None, s2=1, state='M')": (
        "momp.Bax_M2",
    ),
    BAX_A4.format("None"): ("momp.Bax_A4",),
    BAX_M4.format("None"): ("momp.Bax_M4",),
    "Bcl2(bf=None)": ("momp.Bcl2",),
    "Bax(bf=1, s1=None, s2=None, state='A') % Bcl2(bf=1)": (
        "momp.Bax_A_Bcl2",
        "momp.bcl2.monomer_Bcl2",
    ),
    "Bax(bf=1, s1=None, s2=None, state='M') % Bcl2(bf=1)": ("momp.bcl2.monomer_Bcl2",),
    "Bax(bf=1, s1=None, s2=2, state='A') % Bax(bf=None, s1=2, s2=None, state='A') % Bcl2(bf=1)": (
        "momp.bcl2.dimer_Bcl2",
    ),
    "Bax(bf=1, s1=None, s2=2, state='M') % Bax(bf=None, s1=2, s2=None, state='M') % Bcl2(bf=1)": (
        "momp.bcl2.dimer_Bcl2",
    ),
    BAX_A4_BOUND.format("Bcl2(bf=1)"): ("momp.bcl2.tetramer_Bcl2",),
    BAX_M4_BOUND.format("Bcl2(bf=1)"): ("momp.bcl2.tetramer_Bcl2",),
    "Mito(bf=None, state='I')": ("momp.Mito_I",),
    "Mito(bf=None, state='A')": ("momp.Mito_A",),
    BAX_M4_BOUND.format("Mito(bf=1, state='I')"): ("momp.Bax_M4_Mito_I",),
    # MOMP pore transport
    "Smac(bf=None, state='M')": ("momp.Smac_M",),
    "Smac(bf=None, state='C')": ("momp.Smac_C",),
    "CytoC(bf=None, state='M')": ("momp.CytoC_M",),
    "CytoC(bf=None, state='C')": ("momp.CytoC_C",),
    "Bax(bf=1, s1=None, s2=None, state='A') % Smac(bf=1, state='M')": (
        "momp.transport.pore_Smac_M",
    ),
    "Bax(bf=1, s1=None, s2=None, state='A') % CytoC(bf=1, state='M')": (
        "momp.transport.pore_CytoC_M",
    ),
    BAX_A4_BOUND.format("Smac(bf=1, state='M')"): ("momp.transport.pore_Smac_M",),
    BAX_A4_BOUND.format("CytoC(bf=1, state='M')"): ("momp.transport.pore_CytoC_M",),
    BAX_M4_BOUND.format("Smac(bf=1, state='M')"): ("momp.transport.pore_Smac_M",),
    BAX_M4_BOUND.format("CytoC(bf=1, state='M')"): ("momp.transport.pore_CytoC_M",),
    "Mito(bf=1, state='A') % Smac(bf=1, state='M')": ("momp.transport.pore_Smac_M",),
    "CytoC(bf=1, state='M') % Mito(bf=1, state='A')": ("momp.transport.pore_CytoC_M",),
    # CaspaseCascade
    "Smac(bf=None, state='A')": ("caspases.Smac_A",),
    "CytoC(bf=None, state='A')": ("caspases.CytoC_A",),
    "XIAP(bf=None)": ("caspases.XIAP",),
    "Smac(bf=1, state='A') % XIAP(bf=1)": ("caspases.Smac_A_XIAP",),
    "C3(bf=None, state='pro')": ("caspases.C3_pro",),
    "C3(bf=None, state='A')": ("caspases.C3_A",),
    "C3(bf=None, state='ub')": ("caspases.C3_ub",),
    "C3(bf=1, state='pro') % C8(bf=1, state='A')": ("caspases.C8_A_C3_pro",),
    "C3(bf=1, state='A') % XIAP(bf=1)": ("caspases.XIAP_C3_A",),
    "PARP(bf=None, state='U')": ("caspases.PARP_U",),
    "PARP(bf=None, state='C')": ("caspases.PARP_C",),
    "C3(bf=1, state='A') % PARP(bf=1, state='U')": ("caspases.C3_A_PARP_U",),
    "C6(bf=None, state='pro')": ("caspases.C6_pro",),
    "C6(bf=None, state='A')": ("caspases.C6_A",),
    "C3(bf=1, state='A') % C6(bf=1, state='pro')": ("caspases.C3_A_C6_pro",),
    "C6(bf=1, state='A') % C8(bf=1, state='pro')": ("caspases.C6_A_C8_pro",),
    # Apoptosome and DoubleApoptosome
    "Apaf(bf=None, state='I')": ("apoptosome.Apaf_I",),
    "Apaf(bf=None, state='A')": ("apoptosome.Apaf_A",),
    "Apaf(bf=1, state='I') % CytoC(bf=1, state='A')": ("apoptosome.CytoC_A_Apaf_I",),
    "C9(bf=None)": ("apoptosome.C9",),
    "Apop(bf=None)": ("apoptosome.Apop",),
    "Apop(bf=1) % C3(bf=1, state='pro')": ("apoptosome.Apop_C3_pro",),
    "Apop(bf=1) % XIAP(bf=1)": ("apoptosome.Apop_XIAP",),
    "Apaf(bf=1, state='A') % C3(bf=1, state='pro')": ("apoptosome.Apaf_A_C3_pro",),
    "Apaf(bf=1, state='A') % C3(bf=1, state='A')": ("apoptosome.C3_A_Apaf_A",),
    "Apaf(bf=1, state='A') % XIAP(bf=1)": ("apoptosome.Apaf_A_XIAP",),
    # CytosolicBcl2
    "Bcl2c(b=None)": ("bcl2c.Bcl2c",),
    "Bcl2c(b=1) % Bid(bf=1, state='T')": ("bcl2c.Bid_T_Bcl2c",),
    # Biosensors
    "sCas3(sl=1, bf=None) % sCas3(sl=1, bf=None)": ("sCas3.dimer",),
    "sCas3(sl=None, bf=None)": ("sCas3.monomer",),
    "C3(bf=1, state='A') % sCas3(sl=2, bf=1) % sCas3(sl=2, bf=None)": (
        "sCas3.cleavage.enzyme_dimer",
    ),
    "sCas8(sl=1, bf=None) % sCas8(sl=1, bf=None)": ("sCas8.dimer",),
    "sCas8(sl=None, bf=None)": ("sCas8.monomer",),
    "C8(bf=1, state='A') % sCas8(sl=2, bf=1) % sCas8(sl=2, bf=None)": (
        "sCas8.cleavage.enzyme_dimer",
    ),
    "sCas9(sl=1, bf=None) % sCas9(sl=1, bf=None)": ("sCas9.dimer",),
    "sCas9(sl=None, bf=None)": ("sCas9.monomer",),
    "Apop(bf=1) % sCas9(sl=2, bf=1) % sCas9(sl=2, bf=None)": (
        "sCas9.cleavage.enzyme_dimer",
    ),
    "Apaf(bf=1, state='A') % sCas9(sl=2, bf=1) % sCas9(sl=2, bf=None)": (
        "sCas9_apaf.enzyme_dimer",
    ),
}
