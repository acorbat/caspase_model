from simbio import Compartment, Parameter, Species
from simbio.reactions import CatalyzeConvert, MichaelisMenten, ReversibleSynthesis


class Base(Compartment):
    """Implements a base apoptosis model with 4 modules and their connections.

    Modules:
        - Extrinsic
            - Input: L
            - Output: C8.A
        - Mitochondria
            - Input: Bid.U
            - Output: CytoC.A and Smac.A
        - Intrinsic
            - Input: Apaf.I and XIAP
            - Output: Apop
        - Effector
            - Input: C3.pro

    Reactions:
        - Extrinsic to Mitochondria
            - C8.A + Bid.U <> C8.A : Bid.U > C8.A + Bid.T
        - Mitocondria to Intrinsic
            - CytoC.A + Apaf.I <> CytoC.A : Apaf.I > CytoC.A + Apaf.A
            - Smac.A + XIAP <> Smac.A : XIAP
        - Intrinsic to Effector
            - Apop + C3.pro <> Apop : C3.pro > Apop + C3.A
        - Extrinsic to Effector
            - C8.A + C3.pro <> C8.A : C3.pro > C8.A + C3.A

    """

    KF = Parameter(1e-6)
    KR = Parameter(1e-3)
    KC = Parameter(1)

    k_Apaf_I_Cyto_C = Parameter(5e-7)
    k_Apop_C3_pro = Parameter(5e-9)
    k_Smac_A_XIAP = Parameter(7e-6)
    k_C3_pro_C8_A = Parameter(1e-7)

    class Extrinsic(Compartment):
        KF = Parameter(1e-6)
        KR = Parameter(1e-3)
        KC = Parameter(1)

        L = Species(0)

        class C8(Compartment):
            pro = Species(2e4)
            A = Species(0)

    class Intrinsic(Compartment):
        KF = Parameter(1e-6)
        KR = Parameter(1e-3)
        KC = Parameter(1)

        class Apaf(Compartment):
            I = Species(0)  # noqa: #741
            A = Species(0)

        XIAP = Species(0)
        Apop = Species(0)

    class Effector(Compartment):
        KF = Parameter(1e-6)
        KR = Parameter(1e-3)
        KC = Parameter(1)

        class C3(Compartment):
            pro = Species(0)
            A = Species(0)

    class Mitochondria(Compartment):
        KF = Parameter(1e-6)
        KR = Parameter(1e-3)
        KC = Parameter(1)

        class Bid(Compartment):
            U = Species(0)
            T = Species(0)

        class CytoC(Compartment):
            A = Species(0)

        class Smac(Compartment):
            A = Species(0)

    def add_reactions(self):
        # Extrinsic to Mitocondria
        yield MichaelisMenten(
            self.Extrinsic.C8.A,
            self.Mitochondria.Bid.U,
            self.Extrinsic.C8.A & self.Mitochondria.Bid.U,
            self.Mitochondria.Bid.T,
            self.KF,
            self.KR,
            self.KC,
        )

        # Mitocondria to Intrinsic
        yield MichaelisMenten(
            self.Mitochondria.CytoC.A,
            self.Intrinsic.Apaf.I,
            self.Mitochondria.CytoC.A & self.Intrinsic.Apaf.I,
            self.Intrinsic.Apaf.A,
            self.k_Apaf_I_Cyto_C,
            self.KR,
            self.KC,
        )
        yield ReversibleSynthesis(
            self.Mitochondria.Smac.A,
            self.Intrinsic.XIAP,
            self.Mitochondria.Smac.A & self.Intrinsic.XIAP,
            self.k_Smac_A_XIAP,
            self.KR,
        )

        # Intrinsic to Effector
        yield MichaelisMenten(
            self.Intrinsic.Apop,
            self.Effector.C3.pro,
            self.Intrinsic.Apop & self.Effector.C3.pro,
            self.Effector.C3.A,
            self.k_Apop_C3_pro,
            self.KR,
            self.KC,
        )

        # Extrinsic to Effector
        yield MichaelisMenten(
            self.Extrinsic.C8.A,
            self.Effector.C3.pro,
            self.Extrinsic.C8.A & self.Effector.C3.pro,
            self.Effector.C3.A,
            self.k_C3_pro_C8_A,
            self.KR,
            self.KC,
        )


class Extrinsic(Base.Extrinsic):
    """Extrinsic module.

    Reactions:
        - L + R <--> L:R --> DISC
        - pC8 + DISC <--> DISC:pC8 --> C8 + DISC
        - flip + DISC <-->  flip:DISC
        - C8 + BAR <--> BAR:C8
    """

    k_L_R = Parameter(4e-7)
    k_DISC = Parameter(1e-5)

    R = Species(200)
    DISC = Species(0)
    flip = Species(100)
    BAR = Species(1e3)

    def add_reactions(self):
        yield CatalyzeConvert(
            self.L,
            self.R,
            self.L & self.R,
            self.DISC,
            self.k_L_R,
            self.KR,
            self.k_DISC,
        )
        yield MichaelisMenten(
            self.DISC,
            self.C8.pro,
            self.DISC & self.C8.pro,
            self.C8.A,
            self.KF,
            self.KR,
            self.KC,
        )
        yield ReversibleSynthesis(
            self.DISC, self.flip, self.DISC & self.flip, self.KF, self.KR
        )
        yield ReversibleSynthesis(
            self.BAR, self.C8.A, self.BAR & self.C8.A, self.KF, self.KR
        )
