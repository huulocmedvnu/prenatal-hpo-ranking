# Archived inputs

`inputs.tar.xz` holds unmodified copies of the public files that `scripts/download.sh` downloads. The script uses a
copy only when the download from the original source fails or the downloaded file differs from `SHA256SUMS`, and it
checks every file against `SHA256SUMS` at the end.

| File | Source | Terms |
|---|---|---|
| `hpo2026/hp.obo`, `hpo2026/genes_to_phenotype.txt`, `hpo2026/phenotype.hpoa` | Human Phenotype Ontology release v2026-09-01, github.com/obophenotype/human-phenotype-ontology | HPO license (hpo.jax.org/license), unmodified |
| `hpo2020/hp.obo`, `hpo2020/phenotype_to_genes.txt` | HPO release of 12 October 2020 as distributed in PyHPO 2.5.0 (pypi.org/project/pyhpo) | HPO license, unmodified |
| `hgnc/hgnc_complete_set.txt` | HGNC complete set of 2026-09-04, genenames.org | CC0 |
| `fu2022/additional_file_2.docx` | Fu et al. 2022, Genome Med 14:123, Additional file 2 | CC BY 4.0 |

Please cite the original sources when using these files.
