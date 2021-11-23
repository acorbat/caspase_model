from simbio import Compartment, Parameter, Species
from simbio.reactions import MichaelisMenten, ReversibleSynthesis

from .albeck import Albeck11ePoreTransport


class AlbeckAsMatlab(Albeck11ePoreTransport):
    Bcl2c = Species(2e4)
    k_Bid_U_C8_A = Parameter(1e-7)

    def add_reactions(self):
        yield ReversibleSynthesis(
            self.Bid.T, self.Bcl2c, self.Bid.T & self.Bcl2c, self.KF, self.KR
        )

    def override_reactions(self):
        yield MichaelisMenten(
            self.C8.A,
            self.Bid.U,
            self.C8.A & self.Bid.U,
            self.Bid.T,
            self.k_Bid_U_C8_A,
            self.KR,
            self.KC,
        )


class Corbat2018_extrinsic(AlbeckAsMatlab):
    XIAP = Species(1e2, override=True)
    L = Species(1e3, override=True)
    R = Species(1e3, override=True)

    class Sensor(Compartment):
        sCas3 = Species(0)
        sCas9 = Species(0)
        sCas8 = Species(0)

        sCas3_dimer = Species(1e7)
        sCas9_dimer = Species(1e7)
        sCas8_dimer = Species(1e7)

        kf_sCas3 = Parameter(1e-6)
        kr_sCas3 = Parameter(1e-2)
        kc_sCas3 = Parameter(1)

        kf_sCas9 = Parameter(5e-9)
        kr_sCas9 = Parameter(1e-3)
        kc_sCas9 = Parameter(1)

        kf_sCas8 = Parameter(1e-7)
        kr_sCas8 = Parameter(1e-3)
        kc_sCas8 = Parameter(1)

    def add_reactions(self):
        yield MichaelisMenten(
            self.C3.A,
            self.Sensor.sCas3_dimer,
            self.C3.A & self.Sensor.sCas3_dimer,
            2 * self.Sensor.sCas3,
            self.Sensor.kf_sCas3,
            self.Sensor.kr_sCas3,
            self.Sensor.kc_sCas3,
        )
        yield MichaelisMenten(
            self.Apop,
            self.Sensor.sCas9_dimer,
            self.Apop & self.Sensor.sCas9_dimer,
            2 * self.Sensor.sCas9,
            self.Sensor.kf_sCas9,
            self.Sensor.kr_sCas9,
            self.Sensor.kc_sCas9,
        )
        yield MichaelisMenten(
            self.C8.A,
            self.Sensor.sCas8_dimer,
            self.C8.A & self.Sensor.sCas8_dimer,
            2 * self.Sensor.sCas8,
            self.Sensor.kf_sCas8,
            self.Sensor.kr_sCas8,
            self.Sensor.kc_sCas8,
        )


class Corbat2018_intrinsic(Corbat2018_extrinsic):
    L = Species(0, override=True)  # Remove extrinsic stimuli
    IntrinsicStimuli = Species(1e2)

    def add_reactions(self):
        yield MichaelisMenten(
            self.IntrinsicStimuli,
            self.Bid.U,
            self.IntrinsicStimuli & self.Bid.U,
            self.Bid.T,
            self.KF,
            self.KR,
            self.KC,
        )


class ARM_extrinsic(AlbeckAsMatlab):
    L = Species(1e3, override=True)

    class Smac(AlbeckAsMatlab.Smac):
        M = Species(1e5, override=True)

    class CytoC(AlbeckAsMatlab.CytoC):
        M = Species(1e5, override=True)

    class Apaf(AlbeckAsMatlab.Apaf):
        I = Species(1e3, override=True)  # noqa: E741

    C9 = Species(0, override=True)
    XIAP = Species(1e4, override=True)
    k_Apaf_A_C3_pro = Parameter(5e-9)  # _Pore_to_PARP.k_Apop_C3_pro
    k_Apaf_A_C3_A = Parameter(1.3e-6)
    k_Apaf_A_XIAP = Parameter(2e-6)  # _Pore_to_PARP.k_C3_A_XIAP

    class Sensor(Compartment):
        sCas3 = Species(0)
        sCas9 = Species(0)
        sCas8 = Species(0)

        sCas3_dimer = Species(7.5e5)
        sCas9_dimer = Species(7.5e5)
        sCas8_dimer = Species(7.5e5)

        kf_sCas3 = Parameter(2 * 2.8e-7)
        kr_sCas3 = Parameter(1e-2)
        kc_sCas3 = Parameter(1)

        kf_sCas9 = Parameter(2 * 2.8e-7)
        kf_sCas92 = Parameter(2 * 2e-10)
        kr_sCas9 = Parameter(1e-3)
        kc_sCas9 = Parameter(1)

        kf_sCas8 = Parameter(2 * 5.4e-8)
        kr_sCas8 = Parameter(1e-3)
        kc_sCas8 = Parameter(1)

    def add_reactions(self):
        # Apoptosome formation
        #   aApaf + pC3 <-->  aApaf:pC3 --> aApaf + C3
        #   C3 + aApaf <-->  C3:aApaf --> C3 + Apop
        yield MichaelisMenten(
            self.Apaf.A,
            self.C3.pro,
            self.Apaf.A & self.C3.pro,
            self.C3.A,
            self.k_Apaf_A_C3_pro,
            self.KR,
            self.KC,
        )
        yield MichaelisMenten(
            self.C3.A,
            self.Apaf.A,
            self.Apaf.A & self.C3.A,
            self.Apop,
            self.k_Apaf_A_C3_A,
            self.KR,
            self.KC,
        )
        # Apoptosome-related inhibitors
        #   Apaf + XIAP <-->  Apaf:XIAP
        yield ReversibleSynthesis(
            self.Apaf.A, self.XIAP, self.Apaf.A & self.XIAP, self.k_Apaf_A_XIAP, self.KR
        )

        # Sensors
        yield MichaelisMenten(
            self.C3.A,
            self.Sensor.sCas3_dimer,
            self.C3.A & self.Sensor.sCas3_dimer,
            2 * self.Sensor.sCas3,
            self.Sensor.kf_sCas3,
            self.Sensor.kr_sCas3,
            self.Sensor.kc_sCas3,
        )
        yield MichaelisMenten(
            self.Apop,
            self.Sensor.sCas9_dimer,
            self.Apop & self.Sensor.sCas9_dimer,
            2 * self.Sensor.sCas9,
            self.Sensor.kf_sCas9,
            self.Sensor.kr_sCas9,
            self.Sensor.kc_sCas9,
        )
        yield MichaelisMenten(
            self.C8.A,
            self.Sensor.sCas8_dimer,
            self.C8.A & self.Sensor.sCas8_dimer,
            2 * self.Sensor.sCas8,
            self.Sensor.kf_sCas8,
            self.Sensor.kr_sCas8,
            self.Sensor.kc_sCas8,
        )

        # Add interaction between single activation
        yield MichaelisMenten(
            self.Apaf.A,
            self.Sensor.sCas9_dimer,
            self.Apaf.A & self.Sensor.sCas9_dimer,
            2 * self.Sensor.sCas9,
            self.Sensor.kf_sCas92,
            self.Sensor.kr_sCas9,
            self.Sensor.kc_sCas9,
        )

    @classmethod
    def _reactions_to_remove(self):
        #   aApaf + pC9 <-->  Apop
        yield ReversibleSynthesis(self.Apaf.A, self.C9, self.Apop, self.k_Apop, self.KR)
        #   Apop + XIAP <-->  Apop:XIAP
        yield ReversibleSynthesis(
            self.Apop, self.XIAP, self.Apop & self.XIAP, self.k_Apop_XIAP, self.KR
        )


# Remove some ARM reactions
for reaction in ARM_extrinsic._reactions_to_remove():
    for single_reaction in reaction.single_reactions():
        del ARM_extrinsic._reactions[single_reaction.reaction_balance()]


class ARM_intrinsic(ARM_extrinsic):
    L = Species(0, override=True)  # Remove extrinsic stimuli
    IntrinsicStimuli = Species(1e2)

    def add_reactions(self):
        yield MichaelisMenten(
            self.IntrinsicStimuli,
            self.Bid.U,
            self.IntrinsicStimuli & self.Bid.U,
            self.Bid.T,
            self.KF,
            self.KR,
            self.KC,
        )
