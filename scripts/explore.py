"""Exploratory analyses X1 to X4 (PLAN amendment). Writes outputs/explore.json and outputs/explore_x3.tsv."""
import os, json
import numpy as np, pandas as pd
from scipy.stats import wilcoxon, binomtest
from rank import KB, nearest_in, rank_of

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "outputs")
PRENATAL = "HP:0001197"


def mcnemar(a, b):
    g, l = int((~a & b).sum()), int((a & ~b).sum())
    return dict(gained=g, lost=l, p=float(binomtest(g, g + l).pvalue) if g + l else 1.0)


def scores_metric(k, query, metric):
    genes = list(k.gene_terms); pa = k.anc_desc_pa()
    qa = [k.anc[q] & pa for q in query if q in k.anc]
    M = np.zeros((len(qa), len(genes)))
    for i, a in enumerate(qa):
        ics = sorted(((k.ic.get(t, 0.0), t) for t in a), reverse=True)
        for j, g in enumerate(genes):
            M[i, j] = next((v for v, t in ics if t in k.gene_terms[g]), 0.0)
    if metric == "max":
        s = M.max(axis=0)
    elif metric == "symmetric":
        rev = np.zeros(len(genes))
        for j, g in enumerate(genes):
            gt = k.gene_terms[g]
            leaves = [t for t in gt if not any((t in k.anc[u]) and u != t for u in gt)][:200]
            vals = []
            for t in leaves:
                best = max((max((k.ic.get(x, 0.0) for x in (k.anc[t] & k.anc[q] & pa)), default=0.0) for q in query if q in k.anc), default=0.0)
                vals.append(best)
            rev[j] = np.mean(vals) if vals else 0.0
        s = (M.mean(axis=0) + rev) / 2
    else:
        s = M.mean(axis=0)
    return pd.Series(s, index=genes)


if __name__ == "__main__":
    C = pd.read_csv(os.path.join(OUT, "cases.tsv"), sep="\t")
    C = C[C.cohort == "fu2022"]
    kb = {"2020": KB("2020", "2020"), "2026": KB("2026", "2026")}
    k26 = kb["2026"]; t20, t26 = kb["2020"].terms, k26.terms
    common = set(kb["2020"].gene_terms) & set(k26.gene_terms)
    cases = {case: (d.genes.iloc[0].split(";"), list(dict.fromkeys(d.hpo_id))) for case, d in C.groupby("case", sort=False)}
    q20 = {c: list(dict.fromkeys(filter(None, (nearest_in(t20, k26.anc, t26, q) for q in qs)))) for c, (g, qs) in cases.items()}
    out = {}
    # X1 dilution (2026 terms as recorded)
    is_pre = lambda q: PRENATAL in k26.anc.get(q, set())
    rows = []
    for c, (genes, qs) in cases.items():
        pre, struct = [q for q in qs if is_pre(q)], [q for q in qs if not is_pre(q)]
        if pre and struct:
            rows.append(dict(case=c, all=rank_of(k26.scores(qs), genes, common), structural=rank_of(k26.scores(struct), genes, common)))
    X1 = pd.DataFrame(rows)
    out["x1"] = dict(fetuses=len(X1), median_all=float(X1["all"].median()), median_structural=float(X1.structural.median()),
                     top10_all=int((X1["all"] <= 10).sum()), top10_structural=int((X1.structural <= 10).sum()),
                     worse_with_prenatal=int((X1["all"] > X1.structural).sum()), better_with_prenatal=int((X1["all"] < X1.structural).sum()),
                     wilcoxon_p=float(wilcoxon(np.log10(X1["all"]), np.log10(X1.structural)).pvalue),
                     mcnemar=mcnemar(X1.structural <= 10, X1["all"] <= 10))
    print(json.dumps(out, default=str)[-600:], flush=True)
    # X2 vocabulary against annotation
    allq = [q for _, qs in cases.values() for q in qs]
    new_terms = [q for q in allq if q not in t20 or t20[q].get("obsolete")]
    direct = pd.read_csv(os.path.join(ROOT, "inputs", "hpo2026", "genes_to_phenotype.txt"), sep="\t", usecols=["gene_symbol", "hpo_id"])
    dcount = direct.groupby("hpo_id").gene_symbol.nunique()
    closed = pd.Series([t for s in k26.gene_terms.values() for t in s]).value_counts()
    pa = k26.anc_desc_pa()
    br = pd.DataFrame({"term": sorted(pa)}); br["prenatal"] = br.term.map(is_pre) & (br.term != PRENATAL)
    br["direct"] = br.term.map(dcount).fillna(0); br["closed"] = br.term.map(closed).fillna(0)
    out["x2"] = dict(case_findings=len(allq), findings_new_after_2020=len(new_terms), new_terms=sorted(set(new_terms)),
                     prenatal_terms=int(br.prenatal.sum()), other_terms=int((~br.prenatal).sum()),
                     prenatal_terms_with_no_gene=int(((br.closed == 0) & br.prenatal).sum()),
                     other_terms_with_no_gene=int(((br.closed == 0) & ~br.prenatal).sum()),
                     median_genes_closed_prenatal=float(br.closed[br.prenatal].median()), median_genes_closed_other=float(br.closed[~br.prenatal].median()),
                     mean_genes_direct_prenatal=round(float(br.direct[br.prenatal].mean()), 2), mean_genes_direct_other=round(float(br.direct[~br.prenatal].mean()), 2))
    print(json.dumps(out, default=str)[-600:], flush=True)
    # X3 rescue, leave one fetus out
    by_gene = {}
    for c, (genes, qs) in cases.items():
        for g in genes:
            by_gene.setdefault(g, []).append(c)
    rows = []
    for g, cs in by_gene.items():
        if len(cs) < 2 or g not in k26.gene_terms:
            continue
        for c in cs:
            extra = set().union(*[k26.anc[q] for o in cs if o != c for q in cases[o][1] if q in k26.anc]) & pa
            before = rank_of(k26.scores(cases[c][1]), [g], common)
            saved = k26.gene_terms[g]; k26.gene_terms[g] = saved | extra
            after = rank_of(k26.scores(cases[c][1]), [g], common)
            k26.gene_terms[g] = saved
            rows.append(dict(gene=g, case=c, before=before, after=after))
    X3 = pd.DataFrame(rows); X3.to_csv(os.path.join(OUT, "explore_x3.tsv"), sep="\t", index=False)
    out["x3"] = dict(fetuses=len(X3), genes=int(X3.gene.nunique()), median_before=float(X3.before.median()), median_after=float(X3.after.median()),
                     top10_before=int((X3.before <= 10).sum()), top10_after=int((X3.after <= 10).sum()),
                     wilcoxon_p=float(wilcoxon(np.log10(X3.before), np.log10(X3.after)).pvalue), mcnemar=mcnemar(X3.before <= 10, X3.after <= 10))
    print(json.dumps(out, default=str)[-600:], flush=True)
    # X4 metric robustness
    rr = {}
    for m in ["max"]:
        r20 = np.array([rank_of(scores_metric(kb["2020"], q20[c], m), cases[c][0], common) for c in cases])
        r26 = np.array([rank_of(scores_metric(k26, q20[c], m), cases[c][0], common) for c in cases])
        rr[m] = dict(top10_2020=round(100 * float((r20 <= 10).mean()), 1), top10_2026=round(100 * float((r26 <= 10).mean()), 1),
                     median_2020=float(np.median(r20)), median_2026=float(np.median(r26)), mcnemar=mcnemar(r20 <= 10, r26 <= 10))
        print(m, rr[m], flush=True)
    out["x4"] = rr
    json.dump(out, open(os.path.join(OUT, "explore.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "x2"}, indent=1)); print({k: v for k, v in out["x2"].items() if k != "new_terms"})
