#!/usr/bin/env bash
# Public inputs into inputs/ (gitignored). Records each source in inputs/RELEASE.txt.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p inputs/hpo2026 inputs/hpo2020 inputs/fu2022
get() { [ -s "$2" ] || { curl -fsSL --retry 5 -o "$2.part" "$1" && mv "$2.part" "$2" && echo "$2 <- $1 $(date +%F)" >> inputs/RELEASE.txt; }; }
H=https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2026-09-01
for f in hp.obo genes_to_phenotype.txt phenotype.hpoa; do get $H/$f inputs/hpo2026/$f; done
get https://files.pythonhosted.org/packages/bf/1e/12da46f40260a8802f840ab5984072f7e4db54d7df8abcf68f82c27e8e6e/pyhpo-2.5.0-py3-none-any.whl inputs/hpo2020/pyhpo-2.5.0.whl
[ -s inputs/hpo2020/hp.obo ] || (cd inputs/hpo2020 && unzip -o -j -q pyhpo-2.5.0.whl 'pyhpo/data/*' -d .)
get "https://static-content.springer.com/esm/art%3A10.1186%2Fs13073-022-01130-x/MediaObjects/13073_2022_1130_MOESM2_ESM.docx" inputs/fu2022/additional_file_2.docx
