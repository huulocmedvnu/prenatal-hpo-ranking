# Prenatal HPO ranking

Code and data to reproduce the study *Fetal Human Phenotype Ontology terms did not improve gene ranking in prenatal
exome sequencing* (T.H.L. Tran and N.H.T. Tran).

The study takes the prenatal ultrasound findings and causal genes of 335 fetuses with a monogenic diagnosis from four
published prenatal exome sequencing cohorts, codes the findings as Human Phenotype Ontology (HPO) terms, and ranks
all annotated genes for each fetus with the 2020 and 2026 HPO releases, using identical queries. It then measures
which prenatal terms lack gene annotation and what supplementing gene annotations with fetal findings would change.
All data are public and de-identified.

## Quick start

```bash
git clone https://github.com/huulocmedvnu/prenatal-hpo-ranking.git
cd prenatal-hpo-ranking
pip install -r requirements.txt
jupyter lab reproduce.ipynb      # then Kernel > Restart Kernel and Run All Cells
```

The notebook downloads the two HPO releases and the supplementary table of the primary cohort, runs every analysis
step, redraws the figures in `figures/`, and compares every regenerated number with the committed `outputs/`.
It prints `REPRODUCED` when all agree. A full run takes about 20 minutes on a laptop (the permutation control is the
slowest step) and needs about 2 GB of disk.

## Pipeline

| Script | What it does |
|---|---|
| `download.sh` | HPO release v2026-09-01, the 2020-10-12 release (PyHPO 2.5.0 package data), Fu et al. 2022 Table S3 |
| `fu_cases.py` | Findings and genes of the primary cohort and the further cohorts, mapped to HPO terms |
| `rank.py` | Resnik best-match average ranking with both releases, primary outcome |
| `gaps.py` | Annotation of each finding to the causal gene, leakage check |
| `explore.py` | Post hoc analyses X1 to X4, including supplementation of gene annotations (X3) |
| `agreement.py` | Agreement between the two independent extractions |
| `sensitivity.py` | Further cohorts without fetuses with an uncertain diagnosis |
| `paper_analyses.py` | Rank distributions, the prenatal branch of the 2026 release, candidate annotations |
| `review_analyses.py`, `review_r6.py` | Post hoc analyses R1 to R7 (re-coding, alternative scores, mechanism, controls) |
| `paper_figures.py`, `flowchart.py` | Figures |

The analysis plan, fixed before any ranking was computed, and its dated amendments are in `PLAN.md`. Analyses added
after the primary result are marked post hoc there.

## Data

| Path | Content |
|---|---|
| `inputs_manual/term_map.csv` | Hand mapping of the 23 finding labels of the primary cohort without an exact HPO match |
| `inputs_manual/external_cases.csv` | Findings and genes of the three further cohorts, final version after adjudication |
| `inputs_manual/second_extraction/` | Both independent extractions before adjudication, and the adjudication record |
| `inputs_manual/prenatal_recoding.csv` | Fetal counterparts of the same meaning used in the re-coding analysis |
| `outputs/` | Every reported number |

Source studies: Fu et al. 2022 (Genome Med 14:123), Graham et al. 2026 (J Med Genet 63:328-336), Smogavec et al.
2022 (Eur J Hum Genet, PMC8991249) and Yao et al. 2026 (Acta Obstet Gynecol Scand 105:1950-1964). Only per-fetus
findings and gene symbols from their published tables are included here.

## Contact

Tran Thai Huu Loc, MD, Department of Obstetrics and Gynaecology, Hong Hung Hospital, Tay Ninh City, Viet Nam.
ORCID 0009-0009-1810-010X. tranthaihuuloc2@gmail.com
