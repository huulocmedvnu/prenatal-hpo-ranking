"""R6 dual coding and R7 worked example (PLAN.md, amendment 2026-10-15e, post hoc). Writes outputs/review_r6.json."""
import os, json
import numpy as np, pandas as pd
from scipy.stats import binomtest, wilcoxon
from rank import KB, nearest_in
from review_analyses import Mat, mcnemar, rank_in

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "outputs")

if __name__ == "__main__":
    C = pd.read_csv(os.path.join(OUT, "cases.tsv"), sep="\t", dtype={"case": str})
    k20, k26 = KB("2020", "2020"), KB("2026", "2026")
    common = set(k20.gene_terms) & set(k26.gene_terms)
    m20, m26 = Mat(k20, "2020", common), Mat(k26, "2026", common)
    cases = {(c, i): (d.genes.iloc[0].split(";"), list(dict.fromkeys(d.hpo_id)), d.finding.tolist())
             for (c, i), d in C.groupby(["cohort", "case"], sort=False)}
    rec = pd.read_csv(os.path.join(ROOT, "inputs_manual", "prenatal_recoding.csv")).set_index("term").fetal_term.to_dict()
    rows = []
    for key, (genes, q26, _) in cases.items():
        if not any(q in rec for q in q26):
            continue
        fetal = list(dict.fromkeys(rec.get(q, q) for q in q26))
        dual = list(dict.fromkeys(q26 + [rec[q] for q in q26 if q in rec]))
        rows.append(dict(cohort=key[0], case=key[1], original=m26.rank(q26, genes), fetal=m26.rank(fetal, genes), dual=m26.rank(dual, genes)))
    D = pd.DataFrame(rows)
    D.to_csv(os.path.join(OUT, "review_r6.tsv"), sep="\t", index=False)
    out = {"r6": dict(fetuses=len(D), **{f"top10_{c}": int((D[c] <= 10).sum()) for c in ("original", "fetal", "dual")},
                      **{f"median_{c}": float(D[c].median()) for c in ("original", "fetal", "dual")},
                      dual_vs_original=mcnemar(D.original <= 10, D.dual <= 10),
                      dual_worse=int((D.dual > D.original).sum()), dual_better=int((D.dual < D.original).sum()),
                      dual_wilcoxon_p=float(wilcoxon(np.log10(D.original), np.log10(D.dual)).pvalue))}
    print(out, flush=True)
    # R7 worked example
    R4 = pd.read_csv(os.path.join(OUT, "review_r4.tsv"), sep="\t", dtype={"case": str})
    R4["fall"] = np.log10(R4.rank_2026) - np.log10(R4.rank_2020)
    r = R4.sort_values("fall", ascending=False).iloc[0]
    key = ("fu2022", r.case); genes, q26, findings = cases[key]
    Q = list(dict.fromkeys(filter(None, (nearest_in(k20.terms, k26.anc, k26.terms, q) for q in q26))))
    g = next(x for x in genes if x in m26.gi)
    ex = dict(case=r.case, gene=g, findings=findings, rank_2020=float(r.rank_2020), rank_2026=float(r.rank_2026), terms=[])
    for q in Q:
        item = dict(query=q, query_name=k26.terms[q]["name"])
        for lab, k, m in (("2020", k20, m20), ("2026", k26, m26)):
            gt = k.gene_terms[g]
            best = max(((k.ic.get(a, 0.0), a) for a in (k.anc[q] & k.anc_desc_pa() & gt)), default=(0.0, None))
            item[f"best_{lab}"] = best[1]; item[f"best_name_{lab}"] = k.terms[best[1]]["name"] if best[1] else None
            item[f"ic_{lab}"] = round(best[0], 2)
            item[f"genes_with_query_{lab}"] = int(sum(q in s for s in k.gene_terms.values()))
        ex["terms"].append(item)
    out["r7"] = ex
    print(json.dumps(ex, indent=1))
    json.dump(out, open(os.path.join(OUT, "review_r6.json"), "w"), indent=1)
