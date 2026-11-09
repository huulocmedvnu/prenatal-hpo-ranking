# Comparing Human Phenotype Ontology releases highlights persistent gaps in prenatal gene-phenotype annotation

Code and data that reproduce every number and figure of the study by T.H.L. Tran and N.H.T. Tran (Department of
Obstetrics and Gynaecology, Hong Hung Hospital, Tay Ninh, Viet Nam).

The study takes the prenatal imaging findings and the source-reported causal genes of 335 fetuses with a monogenic
diagnosis from four published prenatal exome sequencing cohorts, codes the findings as Human Phenotype Ontology (HPO)
terms, and ranks all annotated genes for each fetus with the HPO releases of 12 October 2020 and 1 September 2026,
using identical queries. It then describes which prenatal terms lack gene annotation or a link to organ-system terms,
and what supplementing gene annotations with fetal findings would change. All data are public and de-identified.

Archived version: Zenodo, DOI [10.5281/zenodo.23385908](https://doi.org/10.5281/zenodo.23385908) (all versions; each release also has its
own DOI). Version 1.1.

## Contents

```
reproduce.ipynb          one-click reproduction (Run All)
scripts/                 analysis pipeline, run in the order listed below
inputs_manual/           data extracted by hand from the source articles, with both independent extractions
reference/               committed results: every reported number (JSON and tab-separated tables)
data/                    archived copies of the public inputs and their SHA-256 checksums
PLAN.md                  analysis plan, fixed before the first ranking, with dated amendments
requirements.txt         Python packages with the versions used
requirements-lock.txt    the complete environment used for the committed results (pip freeze)
```

`inputs/`, `outputs/` and `figures/` are created by a run and are not part of the repository.

## Requirements

- Python 3.12. The pinned numpy and scipy versions have no wheels for Python 3.13 or later.
- macOS or Linux with `bash`, `curl`, `unzip`, `tar` with xz support, and `sha256sum` or `shasum`. On Windows, use the
  Windows Subsystem for Linux (WSL).
- About 1 GB of free disk space. Internet access is used to download the inputs from their original sources. Without
  it, the archived copies in `data/` are used.

## Installation

```bash
git clone https://github.com/huulocmedvnu/prenatal-hpo-ranking.git
cd prenatal-hpo-ranking
python3.12 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt            # or requirements-lock.txt for the exact environment
python -c "import numpy, pandas, scipy, statsmodels, plotly, docx; print('environment ready')"
```

## Reproducing the study

**Option A, one click.** Start `jupyter lab`, open `reproduce.ipynb` and choose Kernel > Restart Kernel and Run All
Cells.

**Option B, command line.**

```bash
python -m nbconvert --to notebook --execute reproduce.ipynb --output executed.ipynb
```

`executed.ipynb` then holds all results and figures.

**Option C, step by step.** Run the scripts from the `scripts/` folder in this order. Each script reads the files
written by the scripts before it.

| Step | Script | Reads | Writes (in `outputs/`) |
|---|---|---|---|
| 1 | `bash download.sh` | internet or `data/` | `inputs/` (checked against `data/SHA256SUMS`) |
| 2 | `python fu_cases.py` | `inputs/hpo2026/hp.obo`, `inputs/fu2022/`, `inputs/hgnc/`, `inputs_manual/term_map.csv`, `inputs_manual/external_cases.csv` | `cases.tsv`, `unmatched_terms.tsv` |
| 3 | `python gene_audit.py` | `inputs/`, `cases.tsv` | `gene_audit.json`, `gene_audit.tsv` |
| 4 | `python rank.py` | `inputs/`, `cases.tsv` | `ranks.tsv`, `rank_summary.json` |
| 5 | `python gaps.py` | `inputs/`, `cases.tsv`, `ranks.tsv` | `gaps.json`, `gaps.tsv` |
| 6 | `python agreement.py` | `inputs/hpo2026/hp.obo`, `inputs_manual/` | `agreement.json`, `agreement_discrepancies.tsv` |
| 7 | `python sensitivity.py` | `ranks.tsv` | `sensitivity.json` |
| 8 | `python explore.py` | `inputs/`, `cases.tsv` | `explore.json`, `explore_x3.tsv` |
| 9 | `python paper_analyses.py` | `inputs/`, `cases.tsv`, `ranks.tsv`, `rank_summary.json`, `explore_x3.tsv` | `paper.json`, `prenatal_branch_terms.tsv`, `candidate_annotations.tsv` |
| 10 | `python robustness.py` | `inputs/`, `inputs_manual/prenatal_recoding.csv`, `cases.tsv`, `ranks.tsv` | `robustness.json`, `robustness_r4.tsv`, `fetus_ranks.tsv` |
| 11 | `python robustness_r6.py` | `inputs/`, `inputs_manual/prenatal_recoding.csv`, `cases.tsv`, `robustness_r4.tsv` | `robustness_r6.json`, `robustness_r6.tsv` |
| 12 | `python paired_ci.py` | | nothing (checks the interval function against R PropCIs values) |
| 13 | `python paired_differences.py` | `ranks.tsv`, `fetus_ranks.tsv` | `paired_ci.json`, `paired_ci.tsv` |
| 14 | `python annotation_groups.py` | `inputs/`, `prenatal_branch_terms.tsv`, `cases.tsv` | `annotation_groups.json`, `fetal_term_linkage.tsv` |
| 15 | `python flowchart.py`, `python paper_figures.py` | outputs above | `figures/` |

## Expected result

The last cell of the notebook compares every file in `outputs/` with the file of the same name in `reference/`
(numbers to a relative tolerance of 1e-9, text exactly) and prints

```
25 committed result files, 25 regenerated
missing: [] extra: [] different: []
REPRODUCED
```

A full run, including the downloads, takes about 3 minutes on a laptop. The permutation control in step 10 is the
slowest step. The pipeline is deterministic: the permutation control uses a fixed random generator seed, and repeated
runs give identical files.

If the notebook prints `NOT REPRODUCED`, the lists name the files that differ. Please open an issue with the output of
`pip freeze` and the operating system.

## Where each result is

| Result | Script | File in `reference/` | Key or column |
|---|---|---|---|
| Cohorts, findings and genes | `fu_cases.py` | `cases.tsv` | one row per fetus and finding |
| Gene symbols in sources and releases | `gene_audit.py` | `gene_audit.json`, `gene_audit.tsv` | `renamed_in_sources`, `harmonised` |
| Primary outcome, top-10 with 2020 and 2026 | `rank.py` | `rank_summary.json` | `fu2022.top10_2020`, `fu2022.top10_2026`, `fu2022.mcnemar_top10` |
| Rank of every fetus | `rank.py` | `ranks.tsv` | `rank_2020`, `rank_2026` |
| Paired differences with 95% CI | `paired_differences.py` | `paired_ci.tsv` | `diff_pp`, `ci` |
| Further cohorts, with and without uncertain diagnoses | `rank.py`, `sensitivity.py` | `rank_summary.json`, `sensitivity.json` | `further` |
| Further cohorts without the Yao cohort (fetal MRI) | `paired_differences.py` | `paired_ci.tsv` | rows `further_without_yao` |
| Rank distribution (top 1, 5, 10, 20, 50, median) | `paper_analyses.py` | `paper.json` | `p1` |
| Prenatal branch of the 2026 release | `paper_analyses.py` | `paper.json`, `prenatal_branch_terms.tsv` | `p2` |
| Annotation by clinical group, fetal term linkage | `annotation_groups.py` | `annotation_groups.json`, `fetal_term_linkage.tsv` | `groups`, `fetal_anomaly_linkage` |
| Findings annotated to the causal gene, by organ system, by inheritance | `gaps.py` | `gaps.json` | `finding_exact_2026`, `by_system_fetus_none_2026`, `by_inheritance` |
| Terms created after 2020 used by the cohorts | `explore.py`, `paper_analyses.py` | `explore.json`, `paper.json` | `x2`, `p5` |
| Supplementation of the causal gene (noncompetitive upper bound) | `explore.py` | `explore.json`, `explore_x3.tsv` | `x3` |
| Supplementation of all genes of the primary cohort (competitive), permutation control | `robustness.py` | `robustness.json`, `fetus_ranks.tsv` | `r5c`, `r5b` |
| Candidate gene and fetal phenotype annotations | `paper_analyses.py` | `candidate_annotations.tsv` | one row per gene and term |
| Alternative scores, own gene sets, mechanism, re-coding | `robustness.py`, `robustness_r6.py` | `robustness.json`, `robustness_r6.json` | `r1` to `r4`, `r6`, `r7` |
| Agreement of the two extractions | `agreement.py` | `agreement.json`, `agreement_discrepancies.tsv` | per cohort |

Labels: X1 to X4 are the exploratory analyses of amendment 2026-10-14, P1 to P5 the descriptive analyses of amendment
2026-10-15c, and R1 to R7 the robustness analyses of amendments 2026-10-15d and 2026-10-15e in `PLAN.md`. The
"external cohorts" of the earlier plan entries are the further cohorts. The cohort key `aogs2026` is the cohort of
Yao et al. 2026.

## Methods in brief

- **Gene annotation set.** For each release, all HPO terms linked to a gene in the gene annotation file, closed under
  their ancestors and restricted to Phenotypic abnormality (HP:0000118). Genes are named by their current HGNC
  symbol through their NCBI Gene ID (`genes.py`), so that a gene renamed between releases or written with a previous
  symbol in a source is the same gene. Loci without an approved HGNC record (regulatory regions such as H19-ICR) are
  left out and listed in `gene_audit.json`.
- **Information content** of a term $t$: $IC(t) = -\ln(n_t / N)$, with $n_t$ the number of annotated genes that carry
  $t$ and $N$ the number of annotated genes in the release.
- **Score** of gene $g$ for a fetus with query terms $Q$ (Resnik best-match average, query to gene):
  $S(g) = \frac{1}{|Q|} \sum_{q \in Q} \max_{t \in A(q) \cap G(g)} IC(t)$, with $A(q)$ the ancestors of $q$
  (including $q$) and $G(g)$ the annotation set of $g$.
- **Rank.** Genes annotated in both releases are ranked by score, ties take their mean rank, and a causal gene without
  annotation takes the last rank. With two causal genes, the better rank is used.
- **Identical queries.** A term absent from the 2020 release is replaced by its nearest ancestor present in that
  release (shortest is_a path, ties broken by the lowest term ID).
- **Paired differences** (second minus first, percentage points) use the Tango score 95% confidence interval
  (`paired_ci.py`). Paired proportions use the exact McNemar test.
- **Supplementation** adds the findings of the other fetuses with the same causal gene, and their ancestors, to the
  annotation of that gene. The fetus under evaluation never contributes to its own gene. X3 supplements only the
  causal gene (noncompetitive upper bound), R5c supplements all genes of the primary cohort at once (competitive).

The full plan, with the date of every change and the analyses marked post hoc, is in `PLAN.md`.

## Data

| Path | Content |
|---|---|
| `inputs_manual/term_map.csv` | Mapping of the 23 finding labels of the primary cohort without an exact HPO match |
| `inputs_manual/external_cases.csv` | Findings and genes of the three further cohorts after adjudication, with source table and row |
| `inputs_manual/second_extraction/` | Both independent extractions before adjudication and the adjudication record |
| `inputs_manual/prenatal_recoding.csv` | Fetal terms with the same meaning as postnatal terms, used in the re-coding analysis (R1) |
| `data/inputs.tar.xz` | Unmodified copies of the downloaded public inputs (see `data/README.md`) |

The findings and genes of the primary cohort are read directly from its supplementary table (Fu et al. 2022, Genome
Med 14:123, Table S3). Further cohorts: Graham et al. 2026 (J Med Genet 63:328-336), Smogavec et al. 2022 (Eur J Hum
Genet, PMC8991249) and Yao et al. 2026 (Acta Obstet Gynecol Scand 105:1950-1964). Only per-fetus findings and gene
symbols from their published tables are included. Findings of the Yao cohort combine neurosonography and fetal MRI,
as reported in the source table.

## Use of an AI assistant

Claude Opus 5.5 (Anthropic), a generative AI assistant, was used under the authors' direction for three tasks.

1. **Analysis code.** The assistant wrote the scripts and the notebook in this repository. T.H.L.T. reviewed the code
   and its outputs. The scoring was implemented twice, as a loop over genes (`rank.py`) and as a similarity matrix
   (`robustness.py`), and the two give identical ranks for every fetus with both releases (checked at every run). The
   confidence interval function is checked against values from the R package PropCIs at every run.
2. **Data extraction.** The findings and genes of the three further cohorts, and the HPO mapping of the hand-mapped
   labels of the primary cohort, were each extracted twice, in two separate assistant sessions without access to each
   other (amendment 2026-10-15 of `PLAN.md`, where the sessions are called agents, and the note of 2026-11-09). Both versions are in `inputs_manual/second_extraction/`. Every discrepancy was resolved against the source
   tables and recorded (`adjudication.csv`). T.H.L.T. reviewed the extracted data and all decisions, and N.H.T.T.
   checked the mapping of the imaging findings to HPO terms.
3. **Text.** The assistant drafted and edited text, including this file. The authors take full responsibility for
   the content.

Every step from the public inputs to the reported numbers is in this repository and can be rerun without the
assistant.

## Troubleshooting

- `pip install` fails on numpy or scipy: the Python version is not 3.12. Create the environment with `python3.12`.
- A figure step fails with a kaleido error: reinstall the pinned version (`pip install kaleido==0.2.1`). The numbers
  do not depend on the figures.
- `download failed, using the archived copy`: a source link is unavailable. The archived copy is identical and the
  run continues.
- `sha256sum: WARNING` or a checksum error at the end of step 1: delete `inputs/` and run again.

## License and citation

Code (`scripts/`, `reproduce.ipynb`): MIT License (`LICENSE`). Extracted data and results (`inputs_manual/`,
`reference/`): Creative Commons Attribution 4.0 (`LICENSE-data`). The archived public inputs in `data/` keep their own
licenses (`data/README.md`). Please cite the Zenodo record (`CITATION.cff`).

## Contact

Thai Huu Loc Tran, MD, Department of Obstetrics and Gynaecology, Hong Hung Hospital, Tay Ninh City, Viet Nam.
ORCID 0009-0009-1810-010X. tranthaihuuloc2@gmail.com
