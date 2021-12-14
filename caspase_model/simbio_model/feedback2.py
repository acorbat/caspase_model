from simbio.reactions import MichaelisMenten

from .corbat import ARM


class NoFeedback(ARM):
    @classmethod
    def _reactions_to_remove(self):
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
            self.C3.A,
            self.Apaf.A,
            self.Apaf.A & self.C3.A,
            self.Apop,
            self.k_Apaf_A_C3_A,
            self.KR,
            self.KC,
        )


# Delete reactions: not yet implementend in SimBio
for reaction in NoFeedback._reactions_to_remove():
    for single_reaction in reaction.single_reactions():
        del NoFeedback._reactions[single_reaction.reaction_balance()]


class C3_to_C6(NoFeedback):
    def add_reactions(self):
        yield MichaelisMenten(
            self.C3.A,
            self.C6.pro,
            self.C3.A & self.C6.pro,
            self.C6.A,
            self.KF,
            self.KR,
            self.KC,
        )


class C3_to_C9(NoFeedback):
    def add_reactions(self):
        yield MichaelisMenten(
            self.C3.A,
            self.Apaf.A,
            self.Apaf.A & self.C3.A,
            self.Apop,
            self.k_Apaf_A_C3_A,
            self.KR,
            self.KC,
        )


class C3_to_Bid(NoFeedback):
    def add_reactions(self):
        yield MichaelisMenten(
            self.C3.A,
            self.Bid.U,
            self.C3.A & self.Bid.U,
            self.Bid.T,
            self.KF,
            self.KR,
            self.KC,
        )
