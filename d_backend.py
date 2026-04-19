
import primer3
import math
import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Tuple
from scorer import PrimerScorer
from thermodynamics import ThermodynamicsCalculator
from main import PrimerCandidate, PrimerPair, DesignParameters
from primer3_adapter import Primer3Adapter


# --- Main Backend Class ---

class PrimerDesignerBackend:
    def __init__(self, params: Dict):
        self.params = params
        self.scorer = PrimerScorer(params)
        self.thermo_calc = ThermodynamicsCalculator(
            salt_conc=self.params['NA_CONC'] / 1000.0,  # mM to M
            mg_conc=self.params['MG_CONC'] / 1000.0,  # mM to M
            dntp_conc=self.params['DNTP_CONC'] / 1000.0, # mM to M
            primer_conc=self.params['PRIMER_CONC'] / 1e9 # nM to M
        )
        self.primer3_adapter = Primer3Adapter(params)
        # Use primer3's own thermo calculations for consistency in dimer checks
        # These will be initialized with the parameters provided by the GUI
        self.thermo_params = primer3.thermoanalysis.ThermoAnalysis(
            mv_conc=self.params['NA_CONC'],
            dv_conc=self.params['MG_CONC'],
            dntp_conc=self.params['DNTP_CONC'],
            dna_conc=self.params['PRIMER_CONC'],
        )

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

    def _calculate_summary_statistics(self, primer_pairs: List[Dict], metrics: List[str]) -> Dict:
        """Calculate summary statistics for a list of primer pairs and specified metrics."""
        stats = {}
        for metric in metrics:
            values = [pair.get(metric) for pair in primer_pairs if pair.get(metric) is not None]
            if values:
                stats[metric] = {
                    "mean": np.mean(values),
                    "median": np.median(values),
                    "std_dev": np.std(values),
                    "min": np.min(values),
                    "max": np.max(values),
                }
        return stats

    def design_primers(self, sequence: str) -> Dict:
        if not sequence:
            raise ValueError("Sequence is empty.")

        seq_args = {
            'SEQUENCE_ID': 'target',
            'SEQUENCE_TEMPLATE': sequence,
        }
        
        target_start_param = None
        target_end_param = None

        target_start = self.params['TARGET_START']
        target_end = self.params['TARGET_END']

        if target_start is not None and target_end is not None:
            try:
                start = int(target_start)
                end = int(target_end)
                if not (0 <= start < end <= len(sequence)):
                    raise ValueError("Invalid target region coordinates.")
                target_start_param = start
                target_end_param = end
            except ValueError:
                raise ValueError("Target Start and End must be valid integers.")
        
        results = self.primer3_adapter.design_primers_with_primer3(
            sequence=sequence,
            target_start=target_start_param,
            target_end=target_end_param
        )
        
        if 'PRIMER_LEFT_EXPLAIN' in results and 'err' in results['PRIMER_LEFT_EXPLAIN']:
            raise ValueError(f"Primer3 error: {results['PRIMER_LEFT_EXPLAIN']}")
        
        num_returned = results.get('PRIMER_PAIR_NUM_RETURNED', 0)
        if num_returned == 0:
            return {'PRIMER_PAIR_NUM_RETURNED': 0, 'all_candidate_metrics': [], 'summary_statistics': {}}

        all_candidate_pairs = []
        for i in range(num_returned):
            fwd_seq = results[f'PRIMER_LEFT_{i}_SEQUENCE']
            rev_seq = results[f'PRIMER_RIGHT_{i}_SEQUENCE']

            # Create PrimerCandidate objects for forward and reverse primers
            forward_primer = PrimerCandidate(
                sequence=fwd_seq,
                start=results[f'PRIMER_LEFT_{i}'][0],
                end=results[f'PRIMER_LEFT_{i}'][0] + results[f'PRIMER_LEFT_{i}'][1],
                length=results[f'PRIMER_LEFT_{i}'][1],
                tm=self.thermo_params.calc_tm(fwd_seq),
                gc_percent=self.thermo_calc.calculate_gc_percent(fwd_seq),
                self_dimer_dG=self.thermo_params.calc_homodimer(fwd_seq).dg,
                hairpin_dG=self.thermo_params.calc_hairpin(fwd_seq).dg,
                end_stability={} # This needs to be calculated by our thermo_calc if needed
            )

            reverse_primer = PrimerCandidate(
                sequence=rev_seq,
                start=results[f'PRIMER_RIGHT_{i}'][0],
                end=results[f'PRIMER_RIGHT_{i}'][0] + results[f'PRIMER_RIGHT_{i}'][1],
                length=results[f'PRIMER_RIGHT_{i}'][1],
                tm=self.thermo_params.calc_tm(rev_seq),
                gc_percent=self.thermo_calc.calculate_gc_percent(rev_seq),
                self_dimer_dG=self.thermo_params.calc_homodimer(rev_seq).dg,
                hairpin_dG=self.thermo_params.calc_hairpin(rev_seq).dg,
                end_stability={} # This needs to be calculated by our thermo_calc if needed
            )
            
            cross_dimer_dG = self.thermo_params.calc_heterodimer(fwd_seq, rev_seq).dg

            # Create PrimerPair object
            primer_pair = PrimerPair(
                forward=forward_primer,
                reverse=reverse_primer,
                product_size=results[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'],
                cross_dimer_dG=cross_dimer_dG,
                tm_difference=abs(forward_primer.tm - reverse_primer.tm),
                product_gc=self.thermo_calc.calculate_gc_percent(sequence[forward_primer.end:reverse_primer.start]),
                product_tm=self.thermo_calc.estimate_product_tm(
                    gc_percent=self.thermo_calc.calculate_gc_percent(sequence[forward_primer.end:reverse_primer.start]),
                    length=len(sequence[forward_primer.end:reverse_primer.start])
                )
            )
            
            # Calculate score and assign badge
            pair_dict = { # Temporarily create a dict for scoring, could refactor scorer to use PrimerPair
                'PRIMER_LEFT_TM': forward_primer.tm,
                'PRIMER_RIGHT_TM': reverse_primer.tm,
                'PRIMER_LEFT_GC_PERCENT': forward_primer.gc_percent,
                'PRIMER_RIGHT_GC_PERCENT': reverse_primer.gc_percent,
                'PRIMER_LEFT_HAIRPIN_DG': forward_primer.hairpin_dG,
                'PRIMER_RIGHT_HAIRPIN_DG': reverse_primer.hairpin_dG,
                'PRIMER_LEFT_SELF_DIMER_DG': forward_primer.self_dimer_dG,
                'PRIMER_RIGHT_SELF_DIMER_DG': reverse_primer.self_dimer_dG,
                'PRIMER_PAIR_CROSS_DIMER_DG': primer_pair.cross_dimer_dG,
                'PRIMER_PAIR_PRODUCT_SIZE': primer_pair.product_size,
            }
            primer_pair.overall_score = self.scorer.calculate_score(pair_dict)
            primer_pair.quality_badge = self._assign_badge(primer_pair.overall_score)

            all_candidate_pairs.append(primer_pair)
            
        all_candidate_pairs.sort(key=lambda p: p.overall_score, reverse=True)

        # Select top N primers for the table display
        # The GUI expects a list of dictionaries, so we need to convert back
        top_primers_dicts = []
        for i, pair in enumerate(all_candidate_pairs[:self.params['PRIMER_NUM_RETURN']]):
            # Reconstruct primer3-like keys for GUI display
            pair_dict = {
                f'PRIMER_LEFT_{i}_SEQUENCE': pair.forward.sequence,
                f'PRIMER_RIGHT_{i}_SEQUENCE': pair.reverse.sequence,
                f'PRIMER_LEFT_{i}': [pair.forward.start, pair.forward.length],
                f'PRIMER_RIGHT_{i}': [pair.reverse.start, pair.reverse.length],
                f'PRIMER_LEFT_{i}_TM': pair.forward.tm,
                f'PRIMER_RIGHT_{i}_TM': pair.reverse.tm,
                f'PRIMER_LEFT_{i}_GC_PERCENT': pair.forward.gc_percent,
                f'PRIMER_RIGHT_{i}_GC_PERCENT': pair.reverse.gc_percent,
                f'PRIMER_LEFT_{i}_HAIRPIN_DG': pair.forward.hairpin_dG,
                f'PRIMER_RIGHT_{i}_HAIRPIN_DG': pair.reverse.hairpin_dG,
                f'PRIMER_LEFT_{i}_SELF_DIMER_DG': pair.forward.self_dimer_dG,
                f'PRIMER_RIGHT_{i}_SELF_DIMER_DG': pair.reverse.self_dimer_dG,
                f'PRIMER_PAIR_{i}_CROSS_DIMER_DG': pair.cross_dimer_dG,
                f'PRIMER_PAIR_{i}_PRODUCT_SIZE': pair.product_size,
                f'CUSTOM_SCORE_{i}': pair.overall_score,
                f'QUALITY_BADGE_{i}': pair.quality_badge
            }
            top_primers_dicts.append(pair_dict)

        final_results = {'PRIMER_PAIR_NUM_RETURNED': len(top_primers_dicts)}
        # Merge all dictionaries into final_results
        for d in top_primers_dicts:
            final_results.update(d)

        # Calculate summary statistics for all candidate primers (still using dicts)
        # Convert all_candidate_pairs (list of PrimerPair objects) to list of dicts for stats calculation
        all_candidate_dicts = []
        for pair_obj in all_candidate_pairs:
            all_candidate_dicts.append({
                'CUSTOM_SCORE': pair_obj.overall_score,
                'PRIMER_LEFT_TM': pair_obj.forward.tm,
                'PRIMER_RIGHT_TM': pair_obj.reverse.tm,
                'PRIMER_LEFT_GC_PERCENT': pair_obj.forward.gc_percent,
                'PRIMER_RIGHT_GC_PERCENT': pair_obj.reverse.gc_percent,
                'PRIMER_LEFT_HAIRPIN_DG': pair_obj.forward.hairpin_dG,
                'PRIMER_RIGHT_HAIRPIN_DG': pair_obj.reverse.hairpin_dG,
                'PRIMER_LEFT_SELF_DIMER_DG': pair_obj.forward.self_dimer_dG,
                'PRIMER_RIGHT_SELF_DIMER_DG': pair_obj.reverse.self_dimer_dG,
                'PRIMER_PAIR_CROSS_DIMER_DG': pair_obj.cross_dimer_dG,
                'PRIMER_PAIR_PRODUCT_SIZE': pair_obj.product_size,
            })
            
        metrics_to_analyze = [
            'CUSTOM_SCORE',
            'PRIMER_LEFT_TM', 'PRIMER_RIGHT_TM',
            'PRIMER_LEFT_GC_PERCENT', 'PRIMER_RIGHT_GC_PERCENT',
            'PRIMER_LEFT_HAIRPIN_DG', 'PRIMER_RIGHT_HAIRPIN_DG',
            'PRIMER_LEFT_SELF_DIMER_DG', 'PRIMER_RIGHT_SELF_DIMER_DG',
            'PRIMER_PAIR_CROSS_DIMER_DG',
            'PRIMER_PAIR_PRODUCT_SIZE',
        ]
        summary_statistics = self._calculate_summary_statistics(all_candidate_dicts, metrics_to_analyze)
        
        final_results['all_candidate_metrics'] = all_candidate_dicts
        final_results['summary_statistics'] = summary_statistics

        return final_results
