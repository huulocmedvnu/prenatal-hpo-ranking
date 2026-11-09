#!/usr/bin/env bash
# Public inputs into inputs/ (gitignored). Each file is downloaded from its original source. When a download fails or
# the downloaded file differs from data/SHA256SUMS, the file is taken from the archived copy in data/inputs.tar.xz
# (unmodified copies, archived with the code on Zenodo). A file already present and verified is not downloaded again.
# Every file is checked against data/SHA256SUMS at the end. Sources are recorded in inputs/RELEASE.txt.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p inputs/hpo2026 inputs/hpo2020 inputs/fu2022 inputs/hgnc
ARCHIVE=data/inputs.tar.xz
if command -v sha256sum >/dev/null; then SHA="sha256sum"; else SHA="shasum -a 256"; fi
# true when the file matches its line in data/SHA256SUMS
verified() { grep " $1\$" data/SHA256SUMS | $SHA -c --status - 2>/dev/null; }
archived() {
    [ -f "$ARCHIVE" ] || { echo "no archived copy available for $1" >&2; exit 1; }
    tar -xJf "$ARCHIVE" -C inputs "${1#inputs/}"
    echo "$1 <- $ARCHIVE (archived copy) $(date +%F)" >> inputs/RELEASE.txt
}
get() {
    [ -s "$2" ] && verified "$2" && return 0
    if curl -fsSL --retry 5 -o "$2.part" "$1"; then
        mv "$2.part" "$2"
        if verified "$2"; then
            echo "$2 <- $1 $(date +%F)" >> inputs/RELEASE.txt
        else
            echo "checksum differs from data/SHA256SUMS, using the archived copy: $1" >&2; archived "$2"
        fi
    else
        rm -f "$2.part"; echo "download failed, using the archived copy: $1" >&2; archived "$2"
    fi
}
H=https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2026-09-01
for f in hp.obo genes_to_phenotype.txt phenotype.hpoa; do get $H/$f inputs/hpo2026/$f; done
# HPO release of 12 October 2020, as distributed in the PyHPO 2.5.0 package
if ! verified inputs/hpo2020/hp.obo || ! verified inputs/hpo2020/phenotype_to_genes.txt; then
    W=https://files.pythonhosted.org/packages/bf/1e/12da46f40260a8802f840ab5984072f7e4db54d7df8abcf68f82c27e8e6e/pyhpo-2.5.0-py3-none-any.whl
    if curl -fsSL --retry 5 -o inputs/hpo2020/pyhpo-2.5.0.whl "$W" \
        && (cd inputs/hpo2020 && unzip -o -j -q pyhpo-2.5.0.whl 'pyhpo/data/*' -d .); then
        echo "inputs/hpo2020 <- $W $(date +%F)" >> inputs/RELEASE.txt
        for f in hp.obo phenotype_to_genes.txt; do
            verified inputs/hpo2020/$f || { echo "checksum differs, using the archived copy: inputs/hpo2020/$f" >&2; archived inputs/hpo2020/$f; }
        done
    else
        echo "download or unpacking failed, using the archived copy: $W" >&2
        archived inputs/hpo2020/hp.obo; archived inputs/hpo2020/phenotype_to_genes.txt
    fi
fi
get "https://static-content.springer.com/esm/art%3A10.1186%2Fs13073-022-01130-x/MediaObjects/13073_2022_1130_MOESM2_ESM.docx" inputs/fu2022/additional_file_2.docx
# HGNC complete set closest to HPO v2026-09-01, for gene symbol harmonisation
get https://storage.googleapis.com/public-download-files/hgnc/archive/archive/monthly/tsv/hgnc_complete_set_2026-09-04.txt inputs/hgnc/hgnc_complete_set.txt
$SHA -c --quiet data/SHA256SUMS
echo "inputs verified against data/SHA256SUMS"
