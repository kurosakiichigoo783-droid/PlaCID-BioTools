import primer3
from typing import List, Dict, Optional, Tuple
from main import PrimerCandidate, PrimerPair, DesignParameters # Assuming these are needed for type hinting

class Primer3Adapter:
    def __init__(self, params: Dict):
        self.params = params

    def design_primers_with_primer3(self, sequence: str, target_start: Optional[int], target_end: Optional[int]) -> Dict:
        seq_args = {
            'SEQUENCE_ID': 'target',
            'SEQUENCE_TEMPLATE': sequence,
        }
        
        if target_start is not None and target_end is not None:
            if not (0 <= target_start < target_end <= len(sequence)):
                raise ValueError("Invalid target region coordinates.")
            seq_args['SEQUENCE_TARGET'] = [target_start, target_end - target_start]
        
        p3_params = {
            'PRIMER_OPT_SIZE': self.params['PRIMER_OPT_SIZE'],
            'PRIMER_MIN_SIZE': self.params['PRIMER_MIN_SIZE'],
            'PRIMER_MAX_SIZE': self.params['PRIMER_MAX_SIZE'],
            'PRIMER_OPT_TM': self.params['PRIMER_OPT_TM'],
            'PRIMER_MIN_TM': self.params['PRIMER_MIN_TM'],
            'PRIMER_MAX_TM': self.params['PRIMER_MAX_TM'],
            'PRIMER_MAX_DIFF_TM': self.params['PRIMER_MAX_DIFF_TM'],
            'PRIMER_OPT_GC_PERCENT': self.params['PRIMER_OPT_GC_PERCENT'],
            'PRIMER_MIN_GC': self.params['PRIMER_MIN_GC'],
            'PRIMER_MAX_GC': self.params['PRIMER_MAX_GC'],
            'PRIMER_PRODUCT_SIZE_RANGE': [int(x) for x in self.params['PRIMER_PRODUCT_SIZE_RANGE'].split('-')],
            'PRIMER_NUM_RETURN': self.params['PRIMER_NUM_RETURN'] * 5, # Request more for internal filtering
            
            'PRIMER_SALT_MONOVALENT': self.params['NA_CONC'],
            'PRIMER_SALT_DIVALENT': self.params['MG_CONC'],
            'PRIMER_DNTP_CONC': self.params['DNTP_CONC'],
            'PRIMER_DNA_CONC': self.params['PRIMER_CONC'],
        }

        # Clean up None values or empty strings from p3_params
        p3_params = {k: v for k, v in p3_params.items() if v is not None and v != ''}

        results = primer3.bindings.designPrimers(seq_args, p3_params)
        
        if 'PRIMER_LEFT_EXPLAIN' in results and 'err' in results['PRIMER_LEFT_EXPLAIN']:
            raise ValueError(f"Primer3 error: {results['PRIMER_LEFT_EXPLAIN']}")
            
        return results
