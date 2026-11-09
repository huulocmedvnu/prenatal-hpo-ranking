"""Audit of gene identifiers. For each source-reported causal gene: its current HGNC symbol and NCBI Gene ID, and
whether the gene is annotated (with a term under Phenotypic abnormality) in the 2020 and 2026 releases, first when
genes are matched by the symbol written in each file and then after harmonisation through NCBI Gene IDs
(scripts/genes.py). Also lists the annotated loci without an approved HGNC record, which are left out of the ranking. Writes outputs/gene_audit.tsv
(one row per causal gene) and outputs/gene_audit.json."""
import os, json
import pandas as pd
from rank import KB, IN, OUT
from genes import hgnc, symbol_from_source


def release_ids(release):
    """NCBI Gene ID and symbol pairs as written in the gene annotation file of a release."""
    if release == "2020":
        A = pd.read_csv(os.path.join(IN, "hpo2020", "phenotype_to_genes.txt"), sep="\t", comment="#", header=None,
                        usecols=[2, 3], names=["ncbi", "gene"], dtype=str)
    else:
        A = pd.read_csv(os.path.join(IN, "hpo2026", "genes_to_phenotype.txt"), sep="\t",
                        usecols=["ncbi_gene_id", "gene_symbol"], dtype=str).rename(columns={"ncbi_gene_id": "ncbi", "gene_symbol": "gene"})
    return A.drop_duplicates()


if __name__ == "__main__":
    by_id, approved, prev, alias = hgnc()
    id_of = {s: i for i, s in by_id.items()}
    raw = {r: KB(r, r, harmonise=False) for r in ("2020", "2026")}
    har = {r: KB(r, r) for r in ("2020", "2026")}
    C = pd.read_csv(os.path.join(OUT, "cases.tsv"), sep="\t", dtype={"case": str})
    F = C.drop_duplicates(["cohort", "case"])
    fetuses = {}
    for _, f in F.iterrows():
        for g in f.source_genes.split(";"):
            fetuses.setdefault(g, []).append(f"{f.cohort}|{f.case}")
    rows = []
    for g, fs in sorted(fetuses.items()):
        h = symbol_from_source(g)
        if g in approved:
            how = "approved"
        elif g in prev:
            how = "previous symbol" if len(prev[g]) == 1 else "ambiguous previous symbol (kept)"
        elif g in alias:
            how = "alias" if len(alias[g]) == 1 else "ambiguous alias (kept)"
        else:
            how = "unresolved (kept)"
        row = dict(source_symbol=g, hgnc_symbol=h, ncbi_gene_id=id_of.get(h, ""), resolved_by=how, fetuses=len(fs),
                   cases=";".join(fs))
        for r in ("2020", "2026"):
            row[f"by_symbol_{r}"] = g in raw[r].gene_terms
            row[f"harmonised_{r}"] = h in har[r].gene_terms
        rows.append(row)
    T = pd.DataFrame(rows)
    T.to_csv(os.path.join(OUT, "gene_audit.tsv"), sep="\t", index=False)

    def absent(prefix):
        a20, a26 = ~T[f"{prefix}_2020"], ~T[f"{prefix}_2026"]
        return dict(absent_2020=[int(a20.sum()), int(T[a20].fetuses.sum())],
                    absent_2026=[int(a26.sum()), int(T[a26].fetuses.sum())],
                    absent_either=[int((a20 | a26).sum()), int(T[a20 | a26].fetuses.sum())],
                    genes=T[a20 | a26].source_symbol.tolist())
    names = {r: {g: sorted(set(d.gene)) for g, d in release_ids(r).assign(h=lambda d: d.ncbi.map(lambda i: by_id.get(i))).dropna(subset=["h"]).groupby("h")}
             for r in ("2020", "2026")}
    # annotated loci without an approved HGNC record, left out of the ranking (symbol as written in the release)
    no_hgnc = {}
    for r in ("2020", "2026"):
        ids = release_ids(r)
        ids = ids[~ids.ncbi.isin(set(by_id)) & ids.gene.isin(set(raw[r].gene_terms))]
        no_hgnc[r] = sorted(f"NCBIGene:{i} ({g})" for i, g in zip(ids.ncbi, ids.gene))
    S = dict(
        causal_genes=len(T), fetuses=len(F), resolved_by=T.resolved_by.value_counts().to_dict(),
        renamed_in_sources=T[T.source_symbol != T.hgnc_symbol][["source_symbol", "hgnc_symbol", "resolved_by"]].to_dict("records"),
        causal_genes_named_differently_in_releases=[
            dict(hgnc_symbol=h, symbol_2020=names["2020"].get(h, []), symbol_2026=names["2026"].get(h, []))
            for h in T.hgnc_symbol if names["2020"].get(h) and names["2026"].get(h) and names["2020"][h] != names["2026"][h]],
        universe_by_symbol=dict(genes_2020=len(raw["2020"].gene_terms), genes_2026=len(raw["2026"].gene_terms),
                                common=len(set(raw["2020"].gene_terms) & set(raw["2026"].gene_terms))),
        universe_harmonised=dict(genes_2020=len(har["2020"].gene_terms), genes_2026=len(har["2026"].gene_terms),
                                 common=len(set(har["2020"].gene_terms) & set(har["2026"].gene_terms))),
        loci_without_hgnc_record=no_hgnc,
        by_symbol=absent("by_symbol"), harmonised=absent("harmonised"))
    with open(os.path.join(OUT, "gene_audit.json"), "w", encoding="utf-8") as f:
        json.dump(S, f, indent=1)
    print(json.dumps(S, indent=1))
