# PlaCID-BioTools
Professional Primer Design Tool - Complete Blueprint 🧬🎯
Part 1: What Should It Have? (Complete Feature List)
Core Features
Feature	Description	Why Critical
Tm Calculation	Nearest-neighbor thermodynamic method (not basic formula)	Basic formulas give ±5°C error
GC% Analysis	Overall + sliding window + 3' end GC clamp	Ensures binding stability
Dimer Check	Self-dimer + cross-dimer + heterodimer	Prevents primer-primer binding
Hairpin Detection	Secondary structure prediction	Prevents self-folding
Specificity Check	BLAST against genome/transcriptome	No off-target amplification
SNP Avoidance	Mask SNPs in primer binding sites	Prevents allele dropout
Repeat Masking	Avoid low-complexity regions	Higher specificity
3' End Stability	ΔG calculation of last 5 bases	Critical for extension
Product Analysis	Secondary structure of amplicon	Ensures amplification
Input Options
text

✅ Single sequence (paste/type)
✅ FASTA file (single/multiple)
✅ GenBank file
✅ Gene ID/Symbol (auto-fetch from NCBI)
✅ Genomic coordinates (chr:start-end)
✅ BED file (batch targets)
✅ Excel/CSV with sequences
✅ SNP list (design flanking primers)
Output Options
text

✅ Interactive results table
✅ Primer visualization on sequence
✅ Detailed HTML report
✅ Excel export (IDT/Sigma format ready)
✅ CSV export
✅ FASTA of primers
✅ Plate map for ordering
✅ PDF report for publication
Part 2: Processing Pipeline for 100% Accuracy
text

┌─────────────────────────────────────────────────────────────────────────────┐
│                    PRIMER DESIGN PIPELINE (100% Accuracy)                   │
└─────────────────────────────────────────────────────────────────────────────┘

                              ┌──────────────┐
                              │   INPUT      │
                              │  Sequence    │
                              └──────┬───────┘
                                     │
                                     ▼
                    ┌────────────────────────────────┐
                    │  STEP 1: SEQUENCE VALIDATION   │
                    │  • Check valid nucleotides     │
                    │  • Check minimum length        │
                    │  • Identify ambiguous bases    │
                    └────────────────┬───────────────┘
                                     │
                                     ▼
                    ┌────────────────────────────────┐
                    │  STEP 2: SEQUENCE PREPROCESSING│
                    │  • Mask repeats (RepeatMasker) │
                    │  • Mask known SNPs (dbSNP)     │
                    │  • Mark low-complexity regions │
                    │  • Identify exon junctions     │
                    └────────────────┬───────────────┘
                                     │
                                     ▼
                    ┌────────────────────────────────┐
                    │  STEP 3: TARGET REGION DEFINE  │
                    │  • User-defined OR             │
                    │  • Auto-detect best region     │
                    │  • Exclude masked regions      │
                    └────────────────┬───────────────┘
                                     │
                                     ▼
                    ┌────────────────────────────────┐
                    │  STEP 4: PRIMER3 CORE ENGINE   │
                    │  • Generate candidate primers  │
                    │  • Apply Tm/GC/Size filters    │
                    │  • Nearest-neighbor Tm calc    │
                    │  • Generate 50-100 candidates  │
                    └────────────────┬───────────────┘
                                     │
                                     ▼
                    ┌────────────────────────────────┐
                    │  STEP 5: THERMODYNAMIC FILTER  │
                    │  • Self-dimer ΔG check         │
                    │  • Hairpin ΔG check            │
                    │  • 3' end stability (ΔG)       │
                    │  • Cross-dimer between F & R   │
                    └────────────────┬───────────────┘
                                     │
                                     ▼
                    ┌────────────────────────────────┐
                    │  STEP 6: SPECIFICITY CHECK     │
                    │  • Local BLAST vs genome       │
                    │  • In-silico PCR simulation    │
                    │  • Check for pseudogenes       │
                    │  • Off-target binding sites    │
                    └────────────────┬───────────────┘
                                     │
                                     ▼
                    ┌────────────────────────────────┐
                    │  STEP 7: AMPLICON ANALYSIS     │
                    │  • Product secondary structure │
                    │  • GC% distribution            │
                    │  • Melting curve prediction    │
                    └────────────────┬───────────────┘
                                     │
                                     ▼
                    ┌────────────────────────────────┐
                    │  STEP 8: SCORING & RANKING     │
                    │  • Weighted penalty score      │
                    │  • ML-based success prediction │
                    │  • Rank all primer pairs       │
                    └────────────────┬───────────────┘
                                     │
                                     ▼
                    ┌────────────────────────────────┐
                    │  STEP 9: FINAL VALIDATION      │
                    │  • Human review flags          │
                    │  • Quality badges (Gold/Silver)│
                    │  • Confidence score            │
                    └────────────────┬───────────────┘
                                     │
                                     ▼
                              ┌──────────────┐
                              │   OUTPUT     │
                              │  Best Primers│
                              └──────────────┘
