# PLAN: prenatal ultrasound phenotypes and HPO-based gene prioritisation

Fixed 2026-10-14, before any ranking was computed. Changes go below as dated amendments.

## Question

Prenatal exome analysis ranks candidate genes by the similarity between the fetal ultrasound findings, coded as
Human Phenotype Ontology (HPO) terms, and the phenotypes annotated to each gene. HPO annotations derive mostly from
postnatal descriptions, and prenatal terms were expanded in 2022 (Dhombres 2022, PMID 35872606). If fetal findings
are poorly represented, the causal gene should rank low when only ultrasound findings are used, and the expansion
should have improved its rank. This study measures (1) how highly the causal gene ranks from ultrasound findings
alone, (2) whether the 2026 HPO release ranks it higher than a 2020 release on the same fetuses, and (3) which
diagnoses still lack annotation of the observed fetal findings.

## Novelty check (2026-10-14)

Closest work: Berrylyzer preprint 2026 (Exomiser and others with sequence data, one HPO version, closed data),
Ardiles-Ruesjas 2026 Prenat Diagn (PMID 41486082, 79 diagnoses, gene linked to any finding, no ranking), commercial
tools on in-house trios. No study found that ranks all genes from published per-case ultrasound findings and
compares HPO releases before and after the prenatal expansion.

## Data

- Primary cohort: Fu et al. 2022 Genome Med (PMID 36307859), Additional file 2, Table S3: 229 diagnosed fetuses,
  ultrasound findings written as HPO term labels, gene, inheritance. Labels are matched to HPO terms by exact
  label or synonym (2026 release), unmatched labels mapped by hand in inputs_manual/ with the rule "nearest HPO term
  that is equal or more general".
- External cohorts (secondary, pooled): Graham 2026 J Med Genet (PMC13151435, Table S1), Smogavec 2022 Eur J Hum
  Genet (PMC8991249, Table 1), prenatal CNS cohort AOGS 2026 (PMC13394886, Table 1). Free text mapped by hand with
  the same rule. Only prenatal ultrasound findings are used (no postnatal or parental features).
- HPO releases: 2020 (hp.obo 2020-10-12, phenotype_to_genes.txt 2020-10, from the pyhpo 2.5.0 package) and
  2026 (v2026-09-01, hp.obo and genes_to_phenotype.txt).

## Ranking (identical code for both releases)

- Gene annotation set: all HPO terms linked to the gene in that release, closed under ancestors.
- Information content (IC) of a term = -log(share of annotated genes carrying the term) in that release.
- Score of gene g for a fetus with terms Q = mean over q in Q of the maximum IC of a common ancestor of q and any
  term of g (query-to-gene best-match average, Resnik).
- Rank of the causal gene among all annotated genes of the release (ties take the mean rank). A causal gene without
  annotation takes the last rank. Main rank scale: rank among genes annotated in both releases, so the gene
  universe is the same. Fetuses with two diagnostic genes use the better rank.
- Case terms present in 2020 are used in both releases (terms created after 2020 are mapped to their nearest 2020
  ancestor for the 2020 run), so the query is identical and only the knowledge base changes.

## Outcomes

- Primary: top-10 rate of the causal gene with the 2026 release, Fu cohort, with 95% Wilson interval. Paired
  comparison of 2020 and 2026: Wilcoxon signed-rank on log10 rank and McNemar test on top-10.
- Secondary: top-1 and top-50 rates, the same in the external cohorts. Decomposition: 2026 ontology with 2020
  annotations against 2020 ontology with 2020 annotations. Share of fetuses whose causal gene has no annotated term
  in the organ system of any finding (top-level HPO branch) in 2026, by organ system. Per inheritance.
- Leakage check: genes whose 2026 annotations cite PMID 36307859 in phenotype.hpoa are listed and the primary
  analysis is repeated without them.

## Limits stated in advance

Findings are those reported by the authors, not standardised. Ranking uses phenotype only, with no variant
evidence. Fu findings reflect a median of 25 weeks, not first-trimester findings.

## Amendments

## Result note 2026-10-14: Fu cohort

All 513 findings of 229 fetuses were mapped (490 exact, 23 by hand in inputs_manual/term_map.csv). Common gene
universe 4,354 genes. Top-10 rate 13.1% (2020) and 14.4% (2026, 95% CI 10.4 to 19.5), McNemar p 0.51 (6 gained, 3
lost). Median rank 90 and 91.5. Log rank changed significantly (Wilcoxon p 0.0003) with more fetuses ranked lower in
2026 (165) than higher (56). Exact annotation of the finding to the causal gene: 53.2% of findings (2020) and 56.1%
(2026). In 2026, 41 fetuses (17.9%) had none of their findings annotated to the causal gene. Findings in the
prenatal development branch were annotated in 34.4%. No 2026 annotation cites the Fu paper (no leakage).

## Amendment 2026-10-14 (exploratory, after the Fu result, before these analyses were run)

Labelled post hoc in any report.
- X1 Dilution: fetuses with both a finding in the prenatal development branch (HP:0001197, for example hydrops,
  increased nuchal translucency, oligohydramnios, polyhydramnios, growth restriction) and at least one structural
  finding outside it. Rank with structural findings only against all findings (2026), paired Wilcoxon on log rank
  and top-10 McNemar.
- X2 Vocabulary against annotation: number of case findings whose term did not exist in 2020. For the 2026 release,
  mean number of annotated genes per term (direct and closed under descendants) in the prenatal development branch
  against other branches.
- X3 Annotation rescue (leave one fetus out): for genes diagnosed in at least 2 Fu fetuses, the gene's 2026
  annotation is extended with the findings of the other fetuses with the same gene (closed under ancestors, IC
  unchanged), and the held-out fetus is ranked again. When the external cohorts are available, Fu findings are added
  for all Fu genes and external fetuses with a Fu gene are ranked.
- X4 Metric robustness: primary comparison repeated with the maximum IC match (single best finding) and with the
  symmetric best-match average.
- 2026-10-14: the symmetric best-match average in X4 was dropped before any of its results were seen, because the
  gene-to-query direction over all genes did not finish in reasonable time. X4 uses the maximum IC match only.

## Result note 2026-10-14: external cohorts and exploratory analyses

External cohorts (106 fetuses: Graham 39, Smogavec 31, AOGS 36): top-10 5.7% (2020) and 7.5% (2026), McNemar p 0.50,
median rank 170 and 150.
X1: in 47 fetuses with both prenatal-branch and structural findings, adding the prenatal findings did not lower the
causal gene's rank (top-10 3 with structural findings only, 8 with all findings, McNemar p 0.06, median rank 98.5
and 81.5, Wilcoxon p 0.07).
X2: only 4 of 513 Fu findings use terms created after 2020. In 2026 the prenatal branch has 267 terms, 142 (53%) with
no annotated gene, and a mean of 7.2 directly annotated genes per term against 13.6 elsewhere.
X3: for 131 fetuses with a recurrent gene (37 genes), adding the findings of the other fetuses with the same gene
raised the top-10 rate from 26 (19.8%) to 35 (26.7%), McNemar p 0.004 (9 gained, 0 lost), median rank 59 to 46.5.
X4: maximum IC match gave top-10 6.6% (2020) and 6.1% (2026), McNemar p 1.0.

## Amendment 2026-10-15: second independent extraction (before any comparison)

The external cohorts and the 23 hand-mapped Fu labels were extracted and mapped once. A second,
independent extraction repeats both from the sources without access to the first files, writing inputs_manual/second_extraction/. A script
(scripts/agreement.py) reports agreement at three levels: fetus and gene, finding (present in both), and HPO term
(identical, or one an ancestor of the other). Discrepancies are listed for the clinician (T.H.L.T.) to resolve
against the source table. The adjudicated file replaces the first extraction and all analyses are rerun. Reported
in Methods as double extraction with clinician adjudication only after the adjudication is done.

## Amendment 2026-10-15b: adjudication and sensitivity analysis (before rerunning)

Agreement (outputs/agreement.json): all 106 external fetuses found by both extractions, causal genes identical in
105, HPO terms identical in 94 to 100% per cohort, identical or related in 97 to 100%. Fu hand-mapped labels 19 of 23
identical. The 17 discrepancies were resolved against the source tables (second_extraction/adjudication.csv, later
reviewed by the clinician). Rules fixed in the process: findings reported as suspected are kept
as reported, a finding naming two entities gets both terms, terms follow the "equal or more general" rule (Left
ventricular pseudoaneurysm becomes Abnormal left ventricle morphology, so one of the post-2020 Fu terms goes).
Sensitivity analysis: the external comparison is repeated without 12 fetuses whose diagnostic status is uncertain in
the source (a single heterozygous variant in a recessive gene, a susceptibility or candidate gene, or a variant
discordant with the stated disorder): AOGS cases 7, 13, 15, 17, 19, 20, 21, 29, 32 and Smogavec cases 17, 18, 22.

## Result note 2026-10-15: after adjudication

Fu cohort unchanged (top-10 13.1% and 14.4%, McNemar p 0.51, 17.9% with no exact annotation, X3 19.8% to 26.7%).
Post-2020 terms among the 513 Fu findings: 3 (Left ventricular pseudoaneurysm now maps to an older term).
External cohorts: top-10 7.5% (2020) and 9.4% (2026), McNemar p 0.50 (before adjudication 5.7% and 7.5%).
Sensitivity without the 12 uncertain diagnoses: 94 fetuses, 7.4% and 9.6%, McNemar p 0.50.

## Amendment 2026-10-15c: descriptive analyses for the report (before running)

Post hoc and descriptive, reported as such. P1: top-1, 5, 10, 20 and 50 shares and median rank (IQR) per release and
cohort, cumulative rank curves. P2: the 2026 terms under Abnormality of prenatal development or birth (HP:0001197)
by sub-branch (and, under Fetal anomaly, by organ sub-branch), whether each term is new since the 2020 release, and
whether any gene is annotated to it or a descendant. P3: the X3 fetuses with the largest rank change, as worked
examples. P4: a resource table of candidate gene to fetal phenotype annotations from all four cohorts (finding terms
whose causal gene carries neither the term nor a descendant in the 2026 release), with the number of supporting
fetuses. P5: findings coded with terms created after 2020, per cohort, and the further cohorts ranked with the
unmapped 2026 terms. No new hypothesis tests.

## Note 2026-10-15: clinician review of the adjudication

The clinician (T.H.L.T.) reviewed the extracted data and the 17 adjudication decisions
(second_extraction/adjudication.csv) and made no changes.


## Amendment 2026-10-15d: robustness analyses (before running, all post hoc)

R1 prenatal re-coding: findings whose term has a fetal counterpart of the same meaning in the 2026 prenatal branch
(inputs_manual/prenatal_recoding.csv, 11 pairs, e.g. Microcephaly to Decreased fetal head circumference) are
re-coded, and the 2026 ranking is repeated in all cohorts. R2 alternative scores: symmetric Resnik best-match
average (mean of query-to-gene and gene-to-query, over the direct annotations of the gene) and Lin best-match
average, 2020 against 2026, primary cohort (terms without annotated genes take the information content of a term
annotated to one gene). R3 unrestricted gene space: each release ranks all of its own annotated genes, reported
as top-10 share and percentile rank, with the fetuses whose causal gene is annotated in only one release. R4
mechanism of the lower 2026 ranks: change in the score of the causal gene, in the information content of its best
matches and in the number of genes scoring above it. R5 cross-cohort supplementation: genes of the further cohorts
that are also diagnosed in the primary cohort receive the findings of the primary-cohort fetuses, and the further
fetuses are ranked again; and a permutation control for X3, in which the added findings are replaced by the same
number of findings drawn at random from primary-cohort fetuses with a different causal gene (1,000 permutations).
Added before running: R5c, a realistic variant of X3 in which the annotations of all genes of the primary cohort
(not only the causal gene of the fetus under evaluation) are supplemented with the findings of the other fetuses,
leaving the fetus under evaluation out, so that competing genes also gain annotations.

## Result note 2026-10-15: R1 to R5

Matrix scores reproduced rank.py exactly (max difference 0). R1: 55 findings re-coded in 49 fetuses, top-10 8 to 4
(McNemar p 0.12), median rank 64 to 87. R2: symmetric Resnik 37 and 37, Lin 29 and 34 (p 0.18). R3: own gene space
31 and 28 (p 0.55), median percentile 2.01% and 2.02%. R4: of 165 fetuses ranked lower in 2026, the causal-gene score
fell in 152 and more genes scored above it in 89. Median query-term information content 3.69 to 3.47. R5a: 51 further
fetuses with a gene shared with the primary cohort, top-10 7 to 12 (p 0.06), 14 better, 0 worse. R5b: random findings
gave a mean of 28.6 in the top 10 (95th percentile 31) against 35 observed, permutation p 0.001. R5c: all genes
supplemented, 26 to 31 (p 0.23). The supplementation gain is specific but smaller when competing genes are also
annotated, and the paper reports both.

## Amendment 2026-10-15e: dual coding and a worked example (before running, post hoc)

R6: for the fetuses of R1, the re-coded findings are queried with both the original and the fetal term (dual
coding), and the 2026 ranking is compared with the original and the fetal-only coding. R7: for the primary-cohort
fetus whose causal gene fell most between 2020 and 2026 (largest increase in log rank), the best-matching term of
each query term for the causal gene and its information content are reported for both releases, as a worked example.
Result 2026-10-15: R6 dual coding, top-10 7 against 8 with the original terms and 4 with fetal terms only (49
fetuses, McNemar p 1.0). R7: fetus 167 (SETD5, increased nuchal translucency only) fell from rank 16 to 714.5 because
the 2026 release no longer annotates SETD5 with Increased nuchal translucency (2020: via ORPHA:199), and 66 genes
carried that term in 2026 against 31 in 2020. Verified in the raw annotation files.

## Note 2026-10-15: co-author role

N.H.T.T. checked the mapping of the ultrasound findings to HPO terms (stated by T.H.L.T.).
