from simbio import Compartment, Parameter, Species
from simbio.reactions import (
    CatalyzeConvert,
    Conversion,
    Equilibration,
    MichaelisMenten,
    ReversibleSynthesis,
)


class Base(Compartment):
    KF = Parameter(1e-6)
    KR = Parameter(1e-3)
    KC = Parameter(1)

    L = Species(3000)
    R = Species(200)
    DISC = Species(0)
    flip = Species(100)

    class C8(Compartment):
        pro = Species(2e4)
        A = Species(0)

    BAR = Species(1e3)

    class Bid(Compartment):
        U = Species(0)
        T = Species(0)
        M = Species(0)

    class Bax(Compartment):
        C = Species(0)
        M = Species(0)
        A = Species(0)

    Bcl2 = Species(0)

    class CytoC(Compartment):
        C = Species(0)
        M = Species(0)
        A = Species(0)

        transloc_rates = Parameter(1e-2)

        def add_reactions(self):
            yield Equilibration(
                self.C, self.A, self.transloc_rates, self.transloc_rates
            )

    class Smac(Compartment):
        C = Species(0)
        M = Species(0)
        A = Species(0)

        transloc_rates = Parameter(1e-2)

        def add_reactions(self):
            yield Equilibration(
                self.C, self.A, self.transloc_rates, self.transloc_rates
            )

    class Apaf(Compartment):
        I = Species(1e5)  # noqa: #741
        A = Species(0)

    Apop = Species(0)

    class C3(Compartment):
        pro = Species(1e4)
        A = Species(0)
        ub = Species(0)

    class C6(Compartment):
        pro = Species(1e4)
        A = Species(0)

    C9 = Species(1e5)

    class PARP(Compartment):
        U = Species(1e6)
        C = Species(0)

    XIAP = Species(1e5)

    k_L_R = Parameter(4e-7)
    k_DISC = Parameter(1e-5)
    k_Apaf_I_Cyto_C = Parameter(5e-7)
    k_Apop = Parameter(5e-8)
    k_Apop_C3_pro = Parameter(5e-9)
    k_Apop_XIAP = Parameter(2e-6)
    k_Smac_A_XIAP = Parameter(7e-6)
    k_C3_pro_C8_A = Parameter(1e-7)
    k_C3_A_XIAP = Parameter(2e-6)
    k_C3_ub = Parameter(1e-1)
    k_PARP_C = Parameter(1e-2)
    k_C6_A_C8_pro = Parameter(3e-8)

    # rec_to_bid

    # =====================
    # tBID Activation Rules
    # ---------------------
    #        L + R <--> L:R --> DISC
    #        pC8 + DISC <--> DISC:pC8 --> C8 + DISC
    #        Bid + C8 <--> Bid:C8 --> tBid + C8
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

        yield MichaelisMenten(
            self.C8.A,
            self.Bid.U,
            self.C8.A & self.Bid.U,
            self.Bid.T,
            self.KF,
            self.KR,
            self.KC,
        )

        # ---------------------
        # Inhibition Rules
        # ---------------------
        #        flip + DISC <-->  flip:DISC
        #        C8 + BAR <--> BAR:C8
        # ---------------------
        yield ReversibleSynthesis(
            self.DISC, self.flip, self.DISC & self.flip, self.KF, self.KR
        )
        yield ReversibleSynthesis(
            self.BAR, self.C8.A, self.BAR & self.C8.A, self.KF, self.KR
        )

        # Apoptosome formation
        # --------------------
        #   Apaf + cCytoC <-->  Apaf:cCytoC --> aApaf + cCytoC
        #   aApaf + pC9 <-->  Apop
        #   Apop + pC3 <-->  Apop:pC3 --> Apop + C3

        yield MichaelisMenten(
            self.CytoC.A,
            self.Apaf.I,
            self.CytoC.A & self.Apaf.I,
            self.Apaf.A,
            self.k_Apaf_I_Cyto_C,
            self.KR,
            self.KC,
        )
        yield ReversibleSynthesis(self.Apaf.A, self.C9, self.Apop, self.k_Apop, self.KR)
        yield MichaelisMenten(
            self.Apop,
            self.C3.pro,
            self.Apop & self.C3.pro,
            self.C3.A,
            self.k_Apop_C3_pro,
            self.KR,
            self.KC,
        )
        # Apoptosome-related inhibitors
        # -----------------------------
        #   Apop + XIAP <-->  Apop:XIAP
        #   cSmac + XIAP <-->  cSmac:XIAP

        yield ReversibleSynthesis(
            self.Apop, self.XIAP, self.Apop & self.XIAP, self.k_Apop_XIAP, self.KR
        )
        yield ReversibleSynthesis(
            self.Smac.A, self.XIAP, self.Smac.A & self.XIAP, self.k_Smac_A_XIAP, self.KR
        )

        # Caspase reactions
        # -----------------
        # Includes effectors, inhibitors, and feedback initiators:
        #
        #   pC3 + C8 <--> pC3:C8 --> C3 + C8 CSPS
        #   pC6 + C3 <--> pC6:C3 --> C6 + C3 CSPS
        #   XIAP + C3 <--> XIAP:C3 --> XIAP + C3_U CSPS
        #   PARP + C3 <--> PARP:C3 --> CPARP + C3 CSPS
        #   pC8 + C6 <--> pC8:C6 --> C8 + C6 CSPS

        yield MichaelisMenten(
            self.C8.A,
            self.C3.pro,
            self.C8.A & self.C3.pro,
            self.C3.A,
            self.k_C3_pro_C8_A,
            self.KR,
            self.KC,
        )
        yield MichaelisMenten(
            self.XIAP,
            self.C3.A,
            self.XIAP & self.C3.A,
            self.C3.ub,
            self.k_C3_A_XIAP,
            self.KR,
            self.k_C3_ub,
        )
        yield MichaelisMenten(
            self.C3.A,
            self.PARP.U,
            self.C3.A & self.PARP.U,
            self.PARP.C,
            self.KF,
            self.k_PARP_C,
            self.KC,
        )
        yield MichaelisMenten(
            self.C3.A,
            self.C6.pro,
            self.C3.A & self.C6.pro,
            self.C6.A,
            self.KF,
            self.KR,
            self.KC,
        )
        yield MichaelisMenten(
            self.C6.A,
            self.C8.pro,
            self.C6.A & self.C8.pro,
            self.C8.A,
            self.k_C6_A_C8_pro,
            self.KR,
            self.KC,
        )


class Albeck11b(Base):
    k_Bax_C_Bid_T = Parameter(1e-7)

    class Bid(Base.Bid):
        U = Species(1e5, override=True)

    class Bax(Base.Bax):
        C = Species(1e5, override=True)

    Bcl2 = Species(2e4, override=True)

    def add_reactions(self):
        yield MichaelisMenten(
            self.Bid.T,
            self.Bax.C,
            self.Bid.T & self.Bax.C,
            self.Bax.A,
            self.k_Bax_C_Bid_T,
            self.KR,
            self.KC,
        )
        yield ReversibleSynthesis(
            self.Bax.A,
            self.Bcl2,
            self.Bax.A & self.Bcl2,
            self.KF,
            self.KR,
        )


class Albeck11bPoreTransport(Albeck11b):
    k_pore = Parameter(10)

    class Smac(Albeck11b.Smac):
        M = Species(1e6, override=True)

    class CytoC(Albeck11b.CytoC):
        M = Species(1e6, override=True)

    def add_reactions(self):
        yield MichaelisMenten(
            self.Bax.A,
            self.Smac.M,
            self.Bax.A & self.Smac.M,
            self.Smac.C,
            self.KF,
            self.KR,
            self.k_pore,
        )
        yield MichaelisMenten(
            self.Bax.A,
            self.CytoC.M,
            self.Bax.A & self.CytoC.M,
            self.CytoC.C,
            self.KF,
            self.KR,
            self.k_pore,
        )


class Albeck11c(Base):
    k_Bax_C_Bid_T = Parameter(1e-7)

    class Bid(Base.Bid):
        U = Species(4e4, override=True)

    class Bax(Base.Bax):
        C = Species(1e5, override=True)

        A2 = Species(0)
        A4 = Species(0)

        KF = Parameter(1e-6)
        KR = Parameter(1e-3)

        def add_reactions(self):
            # Bax dimerizes/tetramerizes
            yield Equilibration(2 * self.A, self.A2, self.KF, self.KR)
            yield Equilibration(2 * self.A2, self.A4, self.KF, self.KR)

    Bcl2 = Species(2e4, override=True)

    def add_reactions(self):
        yield MichaelisMenten(
            self.Bid.T,
            self.Bax.C,
            self.Bid.T & self.Bax.C,
            self.Bax.A,
            self.k_Bax_C_Bid_T,
            self.KR,
            self.KC,
        )
        # Bcl2 inhibits Bax, Bax2, and Bax4
        for x in (self.Bax.A, self.Bax.A2, self.Bax.A4):
            yield ReversibleSynthesis(
                x,
                self.Bcl2,
                x & self.Bcl2,
                self.KF,
                self.KR,
            )


class Albeck11cPoreTransport(Albeck11c):
    k_pore = Parameter(10)

    class Smac(Albeck11c.Smac):
        M = Species(1e6, override=True)
        KF = Parameter(2 * 1e-6)

    class CytoC(Albeck11c.CytoC):
        M = Species(1e6, override=True)

    def add_reactions(self):
        yield MichaelisMenten(
            self.Bax.A4,
            self.Smac.M,
            self.Bax.A4 & self.Smac.M,
            self.Smac.C,
            self.Smac.KF,
            self.KR,
            self.k_pore,
        )
        yield MichaelisMenten(
            self.Bax.A4,
            self.CytoC.M,
            self.Bax.A4 & self.CytoC.M,
            self.CytoC.C,
            self.KF,
            self.KR,
            self.k_pore,
        )


class Albeck11d(Base):
    k_Bax_C_Bid_T = Parameter(1e-7)

    class Bid(Base.Bid):
        U = Species(4e4, override=True)

    class Bax(Base.Bax):
        C = Species(1e5, override=True)

        M2 = Species(0)
        M4 = Species(0)

        # Normalized by fractional volume of the mitochondrial membrane compartment
        KF = Parameter(1e-6 / 0.07)
        KR = Parameter(1e-3)
        transloc_rates = Parameter(1e-2)

        def add_reactions(self):
            # Active Bax translocates to the mitochondria
            yield Equilibration(
                self.A,
                self.M,
                self.transloc_rates,
                self.transloc_rates,
            )
            # Bax dimerizes/tetramerizes
            yield Equilibration(2 * self.M, self.M2, self.KF, self.KR)
            yield Equilibration(2 * self.M2, self.M4, self.KF, self.KR)

    Bcl2 = Species(2e4, override=True)

    def add_reactions(self):
        yield MichaelisMenten(
            self.Bid.T,
            self.Bax.C,
            self.Bid.T & self.Bax.C,
            self.Bax.A,
            self.k_Bax_C_Bid_T,
            self.KR,
            self.KC,
        )
        # Bcl2 inhibits Bax, Bax2, and Bax4
        for x in (self.Bax.M, self.Bax.M2, self.Bax.M4):
            yield ReversibleSynthesis(
                x,
                self.Bcl2,
                x & self.Bcl2,
                self.Bax.KF,
                self.KR,
            )


class Albeck11dPoreTransport(Albeck11d):
    k_pore = Parameter(10)

    class Smac(Albeck11d.Smac):
        M = Species(1e6, override=True)
        KF = Parameter(2 * 1e-6 / 0.07)

    class CytoC(Albeck11d.CytoC):
        M = Species(1e6, override=True)

    def add_reactions(self):
        yield MichaelisMenten(
            self.Bax.M4,
            self.Smac.M,
            self.Bax.M4 & self.Smac.M,
            self.Smac.C,
            self.Smac.KF,
            self.KR,
            self.k_pore,
        )
        yield MichaelisMenten(
            self.Bax.M4,
            self.CytoC.M,
            self.Bax.M4 & self.CytoC.M,
            self.CytoC.C,
            self.KF,
            self.KR,
            self.k_pore,
        )


class Albeck11e(Albeck11d):
    class Mito(Compartment):
        I = Species(5e5)  # noqa: #741
        A = Species(0)

    def add_reactions(self):
        yield ReversibleSynthesis(
            self.Bax.M4,
            self.Mito.I,
            self.Bax.M4 & self.Mito.I,
            self.Bax.KF,
            self.KR,
        )
        yield Conversion(self.Bax.M4 & self.Mito.I, self.Mito.A, self.KC)


class Albeck11ePoreTransport(Albeck11e):
    k_pore = Parameter(10)

    class Smac(Albeck11e.Smac):
        M = Species(1e6, override=True)
        KF = Parameter(2 * 1e-6 / 0.07)

    class CytoC(Albeck11e.CytoC):
        M = Species(1e6, override=True)
        KF = Parameter(2 * 1e-6 / 0.07)

    def add_reactions(self):
        yield MichaelisMenten(
            self.Mito.A,
            self.Smac.M,
            self.Mito.A & self.Smac.M,
            self.Smac.C,
            self.Smac.KF,
            self.KR,
            self.k_pore,
        )
        yield MichaelisMenten(
            self.Mito.A,
            self.CytoC.M,
            self.Mito.A & self.CytoC.M,
            self.CytoC.C,
            self.CytoC.KF,
            self.KR,
            self.k_pore,
        )


class Albeck11f(Albeck11e):
    class Bax(Albeck11e.Bax):
        transloc_rates = Parameter(1e-4, override=True)

        KF2 = Parameter(Albeck11e.Bax.KF.value / 100)
        KF4 = Parameter(Albeck11e.Bax.KF.value / 10)

        def override_reactions(self):
            yield Equilibration(
                self.A,
                self.M,
                self.transloc_rates,
                self.transloc_rates,
            )
            yield Equilibration(
                2 * self.M,
                self.M2,
                self.KF2,
                self.KR,
            )
            yield Equilibration(
                2 * self.M2,
                self.M4,
                self.KF4,
                self.KR,
            )


class Albeck11fPoreTransport(Albeck11f):
    k_pore = Parameter(10)

    class Smac(Albeck11f.Smac):
        M = Species(1e6, override=True)
        KF = Parameter(2 * 1e-6 / 0.07)

    class CytoC(Albeck11f.CytoC):
        M = Species(1e6, override=True)
        KF = Parameter(2 * 1e-6 / 0.07)

    def add_reactions(self):
        yield MichaelisMenten(
            self.Mito.A,
            self.Smac.M,
            self.Mito.A & self.Smac.M,
            self.Smac.C,
            self.Smac.KF,
            self.KR,
            self.k_pore,
        )
        yield MichaelisMenten(
            self.Mito.A,
            self.CytoC.M,
            self.Mito.A & self.CytoC.M,
            self.CytoC.C,
            self.CytoC.KF,
            self.KR,
            self.k_pore,
        )
