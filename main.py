"""
Main Primer Design Engine
Integrates all components for comprehensive primer design
"""

import subprocess
import tempfile
import os
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass, field
from concurrent.futures import ThreadPoolExecutor

from thermodynamics import ThermodynamicsCalculator
from dimer_checker import DimerChecker
from specificity import SpecificityChecker
from scorer import PrimerScorer


@dataclass
class PrimerCandidate:
    """Single primer candidate"""
    sequence: str
    start: int
    end: int
    length: int
    tm: float
    gc_percent: float
    self_dimer_dG: float
    hairpin_dG: float
    end_stability: Dict
    specificity_score: float = 1.0
    penalty: float = 0.0


@dataclass
class PrimerPair:
    """A pair of primers (forward + reverse)"""
    forward: PrimerCandidate
    reverse: PrimerCandidate
    product_size: int
    cross_dimer_dG: float
    tm_difference: float
    product_gc: float
    product_tm: float
    overall_score: float = 0.0
    rank: int = 0
    quality_badge: str = ""  # GOLD, SILVER, BRONZE


@dataclass
class DesignParameters:
    """Parameters for primer design"""
    # Tm settings
    tm_optimal: float = 60.0
    tm_min: float = 58.0
    tm_max: float = 62.0
    tm_max_diff: float = 2.0
    
    # Size settings
    primer_size_optimal: int = 20
    primer_size_min: int = 18
    primer_size_max: int = 25
    
    # Product settings
    product_size_min: int = 100
    product_size_max: int = 300
    product_size_optimal: int = 150
    
    # GC settings
    gc_optimal: float = 50.0
    gc_min: float = 40.0
    gc_max: float = 60.0
    gc_clamp: bool = True
    
    # Chemistry
    na_conc: float = 50e-3
    mg_conc: float = 1.5e-3
    dntp_conc: float = 0.2e-3
    primer_conc: float = 250e-9
    
    # Thresholds
    max_self_dimer_dG: float = -6.0
    max_hairpin_dG: float = -2.0
    max_cross_dimer_dG: float = -6.0
    
    # Options
    check_specificity: bool = True
    mask_snps: bool = True
    mask_repeats: bool = True
    num_primers_return: int = 5


class PrimerDesigner:
    """Main primer design engine"""
    
    def __init__(self, params: Optional[DesignParameters] = None):
        self.params = params or DesignParameters()
        
        self.thermo = ThermodynamicsCalculator(
            salt_conc=self.params.na_conc,
            mg_conc=self.params.mg_conc,
            dntp_conc=self.params.dntp_conc,
            primer_conc=self.params.primer_conc
        )
        
        self.dimer_checker = DimerChecker(
            dimer_threshold=self.params.max_self_dimer_dG,
            hairpin_threshold=self.params.max_hairpin_dG
        )
        
        self.scorer = PrimerScorer(self.params)
        
        if self.params.check_specificity:
            self.specificity_checker = SpecificityChecker()
    
    def design_primers(self, 
                       sequence: str, 
                       target_start: Optional[int] = None,
                       target_end: Optional[int] = None,
                       excluded_regions: Optional[List[Tuple[int, int]]] = None
                       ) -> List[PrimerPair]:
        """
        Main entry point for primer design
        
        Args:
            sequence: Template sequence
            target_start: Start of target region (optional)
            target_end: End of target region (optional)
            excluded_regions: List of (start, end) regions to avoid
            
        Returns:
            List of PrimerPair objects, ranked by score
        """
        # Step 1: Validate and preprocess sequence
        sequence = self._validate_sequence(sequence)
        
        # Step 2: Define target region
        if target_start is None:
            target_start = 0
        if target_end is None:
            target_end = len(sequence)
        
        # Step 3: Get candidate primers from Primer3
        primer3_results = self._run_primer3(
            sequence, target_start, target_end, excluded_regions
        )
        
        # Step 4: Evaluate each candidate pair
        primer_pairs = []
        
        for result in primer3_results:
            try:
                pair = self._evaluate_primer_pair(sequence, result)
                if pair:
                    primer_pairs.append(pair)
            except Exception as e:
                continue  # Skip problematic pairs
        
        # Step 5: Check specificity (if enabled)
        if self.params.check_specificity and primer_pairs:
            primer_pairs = self._check_specificity_batch(primer_pairs)
        
        # Step 6: Score and rank
        primer_pairs = self._score_and_rank(primer_pairs)
        
        # Step 7: Return top N
        return primer_pairs[:self.params.num_primers_return]
    
    def _validate_sequence(self, sequence: str) -> str:
        """Validate and clean input sequence"""
        sequence = sequence.upper().strip()
        sequence = ''.join(c for c in sequence if c in 'ATGCN')
        
        if len(sequence) < 50:
            raise ValueError("Sequence too short (minimum 50 bp)")
        
        # Check for too many Ns
        n_count = sequence.count('N')
        if n_count / len(sequence) > 0.1:
            raise ValueError("Too many ambiguous bases (>10%)")
        
        return sequence
    
    def _run_primer3(self, 
                     sequence: str, 
                     target_start: int, 
                     target_end: int,
                     excluded_regions: Optional[List[Tuple[int, int]]] = None
                     ) -> List[Dict]:
        """
        Run Primer3 to generate candidate primers
        """
        # Build Primer3 input
        primer3_input_lines = [
            'SEQUENCE_ID=target',
            f'SEQUENCE_TEMPLATE={sequence}',
            'PRIMER_TASK=generic',
            'PRIMER_PICK_LEFT_PRIMER=1',
            'PRIMER_PICK_RIGHT_PRIMER=1',
            f'PRIMER_PRODUCT_SIZE_RANGE={self.params.product_size_min}-{self.params.product_size_max}',
            '='
        ]
        primer3_input = '\n'.join(primer3_input_lines) + '\n'

        try:
            process = subprocess.Popen(
                ['primer3_core'],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            process.stdin.write(primer3_input)
            process.stdin.close()

            output = ""
            while True:
                line = process.stdout.readline()
                if not line:
                    break
                output += line
            
            process.wait(timeout=30)

            error = process.stderr.read()
            if error:
                print("\n--- Primer3 Error ---")
                print(error)
                print("---------------------\n")

        except Exception as e:
            print(f"An exception occurred while running primer3_core: {e}")
            output = ""

        # Parse Primer3 output
        return self._parse_primer3_output(output)
    
    def _parse_primer3_output(self, output: str) -> List[Dict]:
        """Parse Primer3 output into list of primer candidates"""
        results = []
        lines = output.strip().split('\n')
        data = {}
        
        for line in lines:
            if '=' in line and line != '=':
                key, value = line.split('=', 1)
                data[key] = value
        
        # Extract primer pairs
        i = 0
        while True:
            left_key = f'PRIMER_LEFT_{i}_SEQUENCE'
            right_key = f'PRIMER_RIGHT_{i}_SEQUENCE'
            
            if left_key not in data or right_key not in data:
                break
            
            left_pos = data.get(f'PRIMER_LEFT_{i}', '0,20').split(',')
            right_pos = data.get(f'PRIMER_RIGHT_{i}', '0,20').split(',')
            
            results.append({
                'left_sequence': data[left_key],
                'right_sequence': data[right_key],
                'left_start': int(left_pos[0]),
                'left_length': int(left_pos[1]),
                'right_start': int(right_pos[0]),
                'right_length': int(right_pos[1]),
                'left_tm': float(data.get(f'PRIMER_LEFT_{i}_TM', 60)),
                'right_tm': float(data.get(f'PRIMER_RIGHT_{i}_TM', 60)),
                'left_gc': float(data.get(f'PRIMER_LEFT_{i}_GC_PERCENT', 50)),
                'right_gc': float(data.get(f'PRIMER_RIGHT_{i}_GC_PERCENT', 50)),
                'product_size': int(data.get(f'PRIMER_PAIR_{i}_PRODUCT_SIZE', 150)),
            })
            
            i += 1
        
        return results
    
    def _evaluate_primer_pair(self, sequence: str, primer3_result: Dict) -> Optional[PrimerPair]:
        """Evaluate a single primer pair with all checks"""
        
        # Forward primer
        fwd_seq = primer3_result['left_sequence']
        fwd_self_dimer = self.dimer_checker.check_self_dimer(fwd_seq)
        fwd_hairpin = self.dimer_checker.check_hairpin(fwd_seq)
        fwd_end_stability = self.thermo.calculate_end_stability(fwd_seq)
        
        forward = PrimerCandidate(
            sequence=fwd_seq,
            start=primer3_result['left_start'],
            end=primer3_result['left_start'] + primer3_result['left_length'],
            length=primer3_result['left_length'],
            tm=self.thermo.calculate_tm(fwd_seq),  # Recalculate with our method
            gc_percent=self.thermo.calculate_gc_percent(fwd_seq),
            self_dimer_dG=fwd_self_dimer.dG,
            hairpin_dG=fwd_hairpin.dG,
            end_stability=fwd_end_stability
        )
        
        # Reverse primer
        rev_seq = primer3_result['right_sequence']
        rev_self_dimer = self.dimer_checker.check_self_dimer(rev_seq)
        rev_hairpin = self.dimer_checker.check_hairpin(rev_seq)
        rev_end_stability = self.thermo.calculate_end_stability(rev_seq)
        
        reverse = PrimerCandidate(
            sequence=rev_seq,
            start=primer3_result['right_start'] - primer3_result['right_length'] + 1,
            end=primer3_result['right_start'] + 1,
            length=primer3_result['right_length'],
            tm=self.thermo.calculate_tm(rev_seq),
            gc_percent=self.thermo.calculate_gc_percent(rev_seq),
            self_dimer_dG=rev_self_dimer.dG,
            hairpin_dG=rev_hairpin.dG,
            end_stability=rev_end_stability
        )
        
        # Cross-dimer check
        cross_dimer = self.dimer_checker.check_cross_dimer(fwd_seq, rev_seq)
        
        # Product analysis
        product_start = forward.end
        product_end = reverse.start
        product_seq = sequence[product_start:product_end]
        product_gc = self.thermo.calculate_gc_percent(product_seq) if product_seq else 50.0
        
        # Create pair
        pair = PrimerPair(
            forward=forward,
            reverse=reverse,
            product_size=primer3_result['product_size'],
            cross_dimer_dG=cross_dimer.dG,
            tm_difference=abs(forward.tm - reverse.tm),
            product_gc=product_gc,
            product_tm=self._estimate_product_tm(product_gc, len(product_seq))
        )
        
        return pair
    

    
    def _check_specificity_batch(self, pairs: List[PrimerPair]) -> List[PrimerPair]:
        """Check specificity for all primer pairs"""
        # Would integrate with BLAST here
        # For now, return all pairs
        return pairs
    
    def _score_and_rank(self, pairs: List[PrimerPair]) -> List[PrimerPair]:
        """Score and rank primer pairs"""
        for pair in pairs:
            pair.overall_score = self.scorer.calculate_score(pair)
            pair.quality_badge = self._assign_badge(pair.overall_score)
        
        # Sort by score (descending)
        pairs.sort(key=lambda x: x.overall_score, reverse=True)
        
        # Assign ranks
        for i, pair in enumerate(pairs):
            pair.rank = i + 1
        
        return pairs
    
    def _assign_badge(self, score: float) -> str:
        """Assign quality badge based on score"""
        if score >= 90:
            return "🥇 GOLD"
        elif score >= 75:
            return "🥈 SILVER"
        elif score >= 60:
            return "🥉 BRONZE"
        else:
            return "⚠️ LOW"

    def _estimate_product_tm(self, gc_percent: float, length: int) -> float:
        """Estimate product Tm based on GC content and length"""
        if length == 0:
            return 0.0
        return 64.9 + 41 * (gc_percent * length / 100 - 16.4) / length
