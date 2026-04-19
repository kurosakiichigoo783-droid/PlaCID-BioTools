from dataclasses import dataclass

@dataclass
class ScoringWeights:
    tm_deviation: float = 1.5
    gc_deviation: float = 0.5
    tm_difference: float = 2.0
    hairpin_dg: float = 0.5 # Penalty per kcal/mol for hairpins
    self_dimer_dg: float = 0.8 # Penalty per kcal/mol for self-dimers
    cross_dimer_dg: float = 1.0 # Penalty per kcal/mol for cross-dimers

class PrimerScorer:
    def __init__(self, params):
        self.params = params
        self.weights = ScoringWeights()

    def calculate_score(self, pair) -> float:
        penalty = 0.0
        # Forward primer
        penalty += abs(pair['PRIMER_LEFT_TM'] - self.params['PRIMER_OPT_TM']) * self.weights.tm_deviation
        penalty += abs(pair['PRIMER_LEFT_GC_PERCENT'] - self.params['PRIMER_OPT_GC_PERCENT']) * self.weights.gc_deviation
        # Reverse primer
        penalty += abs(pair['PRIMER_RIGHT_TM'] - self.params['PRIMER_OPT_TM']) * self.weights.tm_deviation
        penalty += abs(pair['PRIMER_RIGHT_GC_PERCENT'] - self.params['PRIMER_OPT_GC_PERCENT']) * self.weights.gc_deviation
        # Pair
        penalty += abs(pair['PRIMER_LEFT_TM'] - pair['PRIMER_RIGHT_TM']) * self.weights.tm_difference

        # Dimer and hairpin penalties (more negative dG is worse)
        if pair['PRIMER_LEFT_HAIRPIN_DG'] < -2:
            penalty += abs(pair['PRIMER_LEFT_HAIRPIN_DG']) * self.weights.hairpin_dg
        if pair['PRIMER_RIGHT_HAIRPIN_DG'] < -2:
            penalty += abs(pair['PRIMER_RIGHT_HAIRPIN_DG']) * self.weights.hairpin_dg
        if pair['PRIMER_LEFT_SELF_DIMER_DG'] < -6:
            penalty += abs(pair['PRIMER_LEFT_SELF_DIMER_DG']) * self.weights.self_dimer_dg
        if pair['PRIMER_RIGHT_SELF_DIMER_DG'] < -6:
            penalty += abs(pair['PRIMER_RIGHT_SELF_DIMER_DG']) * self.weights.self_dimer_dg
        if pair['PRIMER_PAIR_CROSS_DIMER_DG'] < -6:
            penalty += abs(pair['PRIMER_PAIR_CROSS_DIMER_DG']) * self.weights.cross_dimer_dg

        return max(0, 100 - penalty)
