"""Descriptive analyses (PLAN.md, amendment 2026-10-15c).
P1 rank distribution per release and cohort, P2 the prenatal branch of HPO 2026 by sub-branch, novelty and gene
annotation, P3 worked examples from X3, P4 candidate gene to fetal phenotype annotations from all cohorts, P5 use of post-2020 terms per cohort.
Writes outputs/paper.json, outputs/prenatal_branch_terms.tsv, outputs/candidate_annotations.tsv."""
import os, json
import numpy as np, pandas as pd
from rank import KB

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "outputs")
PRENATAL, FETAL = "HP:0001197", "HP:0034057"

if __name__ == "__main__":
    R = pd.read_csv(os.path.join(OUT, "ranks.tsv"), sep="\t", dtype={"case": str})
    C = pd.read_csv(os.path.join(OUT, "cases.tsv"), sep="\t", dtype={"case": str})
    k20, k26 = KB("2020", "2020"), KB("2026", "2026")
    t26 = k26.terms
    out = {}
    # P1
    p1 = {}
    for name, d in [("fu2022", R[R.cohort == "fu2022"]), ("further", R[R.cohort != "fu2022"])]:
        s = {"fetuses": len(d)}
        for lab in ("2020", "2026"):
            r = d[f"rank_{lab}"]
            s[lab] = {f"top{k}": [int((r <= k).sum()), round(100 * (r <= k).mean(), 1)] for k in (1, 5, 10, 20, 50)}
            s[lab]["median"], s[lab]["q1"], s[lab]["q3"] = (float(np.percentile(r, q)) for q in (50, 25, 75))
        p1[name] = s
    out["p1"] = p1
    # P2
    nd = lambda t: not t26[t].get("obsolete")
    branch = [t for t, a in k26.anc.items() if PRENATAL in a and t != PRENATAL and nd(t)]
    closed = pd.Series([t for s in k26.gene_terms.values() for t in s]).value_counts()
    top = {c for c in branch if PRENATAL in t26[c]["parents"]}
    fetal_kids = {c for c in branch if any(p in t26[c]["parents"] for p in (FETAL,))}
    fetal_sub = {c for c in branch if set(t26[c]["parents"]) & fetal_kids}
    rows = []
    for t in branch:
        a = k26.anc[t]
        sub = next(iter(sorted(a & top)), "")
        org = next(iter(sorted(a & fetal_sub)), "") if FETAL in a else ""
        rows.append(dict(term=t, name=t26[t]["name"], sub_branch=t26[sub]["name"] if sub else "",
                         fetal_sub_branch=t26[org]["name"] if org else "",
                         new_since_2020=t not in k20.terms or bool(k20.terms[t].get("obsolete")),
                         genes=int(closed.get(t, 0))))
    B = pd.DataFrame(rows).sort_values(["sub_branch", "fetal_sub_branch", "name"])
    B.to_csv(os.path.join(OUT, "prenatal_branch_terms.tsv"), sep="\t", index=False)
    B["no_gene"] = B.genes == 0
    out["p2"] = dict(terms=len(B), no_gene=int(B.no_gene.sum()), new=int(B.new_since_2020.sum()),
                     new_no_gene=int((B.new_since_2020 & B.no_gene).sum()),
                     fetal_anomaly_terms=int((B.sub_branch == "Fetal anomaly").sum()),
                     fetal_anomaly_no_gene=int(((B.sub_branch == "Fetal anomaly") & B.no_gene).sum()),
                     by_sub_branch=B.groupby("sub_branch").agg(terms=("term", "size"), no_gene=("no_gene", "sum"),
                                                               new=("new_since_2020", "sum")).astype(int).to_dict("index"),
                     by_fetal_sub_branch=B[B.fetal_sub_branch != ""].groupby("fetal_sub_branch").agg(
                         terms=("term", "size"), no_gene=("no_gene", "sum")).astype(int).to_dict("index"))
    out["p2"]["examples_genes"] = {n: int(closed.get(t, 0)) for n, t in (("Microcephaly", "HP:0000252"),
                                   ("Decreased fetal head circumference", "HP:0020306"), ("Decreased biparietal diameter", "HP:0020259"),
                                   ("Abnormality of ductus venosus blood flow", "HP:0010947"), ("Fixed fetal head flexion", "HP:6001473"))}
    # P3
    X3 = pd.read_csv(os.path.join(OUT, "explore_x3.tsv"), sep="\t", dtype={"case": str})
    X3["gain"] = np.log10(X3.before) - np.log10(X3.after)
    ex = X3.sort_values("gain", ascending=False).head(5)
    fu = C[C.cohort == "fu2022"]
    out["p3"] = [dict(gene=r.gene, case=r.case, before=float(r.before), after=float(r.after),
                      findings=fu[fu.case == r.case].finding.tolist(),
                      other_fetuses=int((X3.gene == r.gene).sum() - 1)) for r in ex.itertuples()]
    # P4
    rows = []
    for (coh, case), d in C.groupby(["cohort", "case"]):
        for g in d.genes.iloc[0].split(";"):
            gt = k26.gene_terms.get(g, set())
            for t in dict.fromkeys(d.hpo_id):
                if t not in t26 or t26[t].get("obsolete"):
                    continue
                desc = any(t in k26.anc[x] for x in gt)
                if not desc:
                    rows.append(dict(gene=g, term=t, name=t26[t]["name"], cohort=coh, case=case,
                                     gene_in_hpo_2026=g in k26.gene_terms))
    P = pd.DataFrame(rows)
    T = (P.groupby(["gene", "term", "name", "gene_in_hpo_2026"])
         .agg(fetuses=("case", "size"), cohorts=("cohort", lambda s: ";".join(sorted(set(s)))),
              cases=("case", lambda s: ";".join(map(str, s))))
         .reset_index().sort_values(["fetuses", "gene"], ascending=[False, True]))
    T.to_csv(os.path.join(OUT, "candidate_annotations.tsv"), sep="\t", index=False)
    out["p4"] = dict(pairs=len(T), genes=int(T.gene.nunique()), pairs_two_or_more=int((T.fetuses >= 2).sum()),
                     genes_absent_from_hpo=int(T[~T.gene_in_hpo_2026].gene.nunique()))
    # P5 use of post-2020 terms per cohort (amendment 2026-10-15c, descriptive)
    C["new"] = C.hpo_id.map(lambda q: q not in k20.terms or bool(k20.terms[q].get("obsolete")))
    S = json.load(open(os.path.join(OUT, "rank_summary.json")))
    out["p5"] = {coh: dict(findings=int(len(d)), findings_new=int(d.new.sum()), fetuses_new=int(d[d.new].case.nunique()))
                 for coh, d in C.groupby("cohort")}
    out["p5"]["further_fullterms_top10"] = S["external"]["top10_2026_fullterms"]
    out["p5"]["further_fullterms_median"] = S["external"]["median_rank_2026_fullterms"]
    json.dump(out, open(os.path.join(OUT, "paper.json"), "w"), indent=1)
    print(json.dumps(out, indent=1)[:5000])
