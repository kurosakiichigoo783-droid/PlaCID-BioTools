
import math

class ThermodynamicsCalculator:
    # Complement mapping
    COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

    # Nearest-Neighbor Parameters (SantaLucia 1998)
    NN_PARAMS = {
        'AA/TT': {'dH': -7.9, 'dS': -22.2},
        'AT/TA': {'dH': -7.2, 'dS': -20.4},
        'TA/AT': {'dH': -7.2, 'dS': -21.3},
        'CA/GT': {'dH': -8.5, 'dS': -22.7},
        'GT/CA': {'dH': -8.4, 'dS': -22.4},
        'CT/GA': {'dH': -7.8, 'dS': -21.0},
        'GA/CT': {'dH': -8.2, 'dS': -22.2},
        'CG/GC': {'dH': -10.6, 'dS': -27.2},
        'GC/CG': {'dH': -9.8, 'dS': -24.4},
        'GG/CC': {'dH': -8.0, 'dS': -19.9},
    }

    def __init__(self, primer_conc=250e-9, salt_conc=50e-3, mg_conc=1.5e-3, dntp_conc=0.2e-3):
        self.primer_conc = primer_conc
        self.salt_conc = salt_conc
        self.mg_conc = mg_conc
        self.dntp_conc = dntp_conc # This parameter is not directly used in the current Tm calculation but is passed from main.py

    def _get_complement(self, sequence):
        """Return complement of a DNA sequence"""
        return ''.join(self.COMPLEMENT.get(base, 'N') for base in sequence.upper())

    def calculate_tm(self, sequence):
        """
        Accurate Tm using nearest-neighbor thermodynamics
        """
        dH_total = 0  # Enthalpy
        dS_total = 0  # Entropy

        # Sum all nearest-neighbor pairs
        for i in range(len(sequence) - 1):
            pair = sequence[i:i+2]
            complement_pair = self._get_complement(pair) # Use _get_complement here
            key = f"{pair}/{complement_pair}"

            # Check both orientations for NN_PARAMS
            if key in self.NN_PARAMS:
                dH_total += self.NN_PARAMS[key]['dH']
                dS_total += self.NN_PARAMS[key]['dS']
            elif f"{complement_pair}/{pair}" in self.NN_PARAMS: # Check reverse complement pair
                 dH_total += self.NN_PARAMS[f"{complement_pair}/{pair}"]['dH']
                 dS_total += self.NN_PARAMS[f"{complement_pair}/{pair}"]['dS']
            else:
                # Handle cases where a pair might not be in NN_PARAMS (e.g., Ns or unusual pairs)
                # For now, we'll just skip them, but a more robust solution might add penalties
                pass

        # Salt correction (Owczarzy et al. 2008)
        # Assuming monovalent salt correction for now
        # More complex models can incorporate Mg2+ if needed
        dS_corrected = dS_total + 0.368 * (len(sequence) - 1) * math.log(self.salt_conc)

        # Final Tm calculation
        R = 1.987  # Gas constant (cal/(mol*K))
        # Convert dH to cal/mol and dS to cal/(mol*K)
        Tm_kelvin = (dH_total * 1000) / (dS_corrected + R * math.log(self.primer_conc / 4))
        Tm_celsius = Tm_kelvin - 273.15

        return round(Tm_celsius, 1)

    def calculate_gc_percent(self, sequence):
        """Calculate GC percentage of a DNA sequence"""
        gc_count = sequence.upper().count('G') + sequence.upper().count('C')
        return (gc_count / len(sequence)) * 100 if sequence else 0

    def calculate_end_stability(self, sequence):
        """Calculate ΔG of a short sequence (for 3' end stability)"""
        dG = 0
        for i in range(len(sequence) - 1):
            pair = sequence[i:i+2]
            comp = self._get_complement(pair)
            key = f"{pair}/{comp}"
            if key in self.NN_PARAMS:
                dH = self.NN_PARAMS[key]['dH']
                dS = self.NN_PARAMS[key]['dS'] / 1000 # Convert to kcal/mol/K
                dG += dH - 310.15 * dS # Assuming 37 C or 310.15 K
        return dG

    @staticmethod
    def estimate_product_tm(gc_percent: float, length: int) -> float:
        """Estimate amplicon melting temperature"""
        # Simple estimation for qPCR melting curve
        return 81.5 + 16.6 * (gc_percent / 100) + 41 * (gc_percent / 100) - 500 / length


# Removed standalone functions as they are now methods of ThermodynamicsCalculator
# COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
# NN_PARAMS = { ... }
