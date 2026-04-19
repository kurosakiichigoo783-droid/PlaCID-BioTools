class DimerChecker:
    def __init__(self, dimer_threshold, hairpin_threshold):
        pass

    def check_self_dimer(self, seq):
        return self

    def check_hairpin(self, seq):
        return self
    
    def check_cross_dimer(self, seq1, seq2):
        return self

    @property
    def dG(self):
        return 0.0
