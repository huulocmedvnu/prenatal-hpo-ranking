"""Gene ranking from ultrasound findings (PLAN, ranking). Resnik best-match average query-to-gene. Releases: 2020
(pyhpo 2.5.0 data) and 2026 (v2026-09-01), plus the decomposition 2026 ontology with 2020 annotations. Writes
outputs/ranks.tsv and outputs/rank_summary.json."""
import os, json
import numpy as np, pandas as pd
from scipy.stats import wilcoxon, binomtest
from statsmodels.stats.proportion import proportion_confint
from hpo import load_obo, ancestors_map

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IN, OUT = os.path.join(ROOT, "inputs"), os.path.join(ROOT, "outputs")
PA = "HP:0000118"


def resolve(terms):
    alt = {a: t for t, d in terms.items() for a in d["alt"]}
    def f(t):
        t = alt.get(t, t)
        d = terms.get(t)
        if d and d.get("obsolete") and d.get("replaced_by"):
            t = d["replaced_by"]
        return t if t in terms and not terms[t].get("obsolete") else None
    return f


def annotations(release):
    if release == "2020":
        A = pd.read_csv(os.path.join(IN, "hpo2020", "phenotype_to_genes.txt"), sep="\t", comment="#", header=None,
                        usecols=[0, 3], names=["hpo_id", "gene"])
    else:
        A = pd.read_csv(os.path.join(IN, "hpo2026", "genes_to_phenotype.txt"), sep="\t", usecols=["hpo_id", "gene_symbol"])[["hpo_id", "gene_symbol"]]
        A.columns = ["hpo_id", "gene"]
    return A.drop_duplicates()


class KB:
    def __init__(self, obo_release, ann_release):
        self.terms = load_obo(os.path.join(IN, f"hpo{obo_release}", "hp.obo"))
        self.anc = ancestors_map(self.terms)
        res = resolve(self.terms)
        A = annotations(ann_release)
        A["hpo_id"] = A.hpo_id.map(res)
        A = A.dropna()
        self.gene_terms = {}
        for g, ts in A.groupby("gene").hpo_id:
            closed = set().union(*[self.anc[t] for t in ts])
            if PA in closed:
                self.gene_terms[g] = closed & self.anc_desc_pa()
        n = len(self.gene_terms)
        cnt = pd.Series([t for s in self.gene_terms.values() for t in s]).value_counts()
        self.ic = (-np.log(cnt / n)).to_dict()
        self.res = res

    def anc_desc_pa(self):
        if not hasattr(self, "_pa"):
            self._pa = {t for t, a in self.anc.items() if PA in a}
        return self._pa

    def scores(self, query):
        q_anc = [self.anc[q] & self.anc_desc_pa() for q in query if q in self.anc]
        genes = list(self.gene_terms)
        out = np.zeros(len(genes))
        for qa in q_anc:
            ics = sorted(((self.ic.get(t, 0.0), t) for t in qa), reverse=True)
            for i, g in enumerate(genes):
                s = self.gene_terms[g]
                out[i] += next((v for v, t in ics if t in s), 0.0)
        return pd.Series(out / max(len(q_anc), 1), index=genes)


def nearest_in(terms_old, anc_new, terms_new, t):
    """Nearest ancestor of t (new ontology) present and active in the old ontology, by shortest is_a path."""
    frontier, seen = [t], {t}
    while frontier:
        hit = [x for x in frontier if x in terms_old and not terms_old[x].get("obsolete")]
        if hit:
            return sorted(hit)[0]
        nxt = []
        for x in frontier:
            for p in terms_new.get(x, {}).get("parents", []):
                if p not in seen:
                    seen.add(p); nxt.append(p)
        frontier = nxt
    return None


def rank_of(sc, genes, universe):
    s = sc.reindex(sorted(universe)).fillna(-1.0)
    r = s.rank(ascending=False, method="average")
    best = [r[g] for g in genes if g in r.index and sc.get(g) is not None and g in sc.index]
    return min(best) if best else float(len(s))


if __name__ == "__main__":
    C = pd.read_csv(os.path.join(OUT, "cases.tsv"), sep="\t")
    kb = {"2020": KB("2020", "2020"), "2026": KB("2026", "2026"), "2026onto_2020ann": KB("2026", "2020")}
    common = set(kb["2020"].gene_terms) & set(kb["2026"].gene_terms)
    t20, t26 = kb["2020"].terms, kb["2026"].terms
    rows = []
    for (coh, case), d in C.groupby(["cohort", "case"], sort=False):
        genes = d.genes.iloc[0].split(";")
        q26 = list(dict.fromkeys(d.hpo_id))
        q20 = list(dict.fromkeys(filter(None, (nearest_in(t20, kb["2026"].anc, t26, q) for q in q26))))
        r = dict(cohort=coh, case=case, section=d.section.iloc[0], genes=";".join(genes), inheritance=d.inheritance.iloc[0],
                 n_terms=len(q26), n_terms_changed=int(len(set(q26) - set(q20))),
                 annotated_2020=any(g in kb["2020"].gene_terms for g in genes),
                 annotated_2026=any(g in kb["2026"].gene_terms for g in genes))
        for lab, k in kb.items():
            sc = k.scores(q20)
            r[f"rank_{lab}"] = rank_of(sc, genes, common)
            r[f"rank_all_{lab}"] = rank_of(sc, genes, set(k.gene_terms))
        r["rank_2026_fullterms"] = rank_of(kb["2026"].scores(q26), genes, common)
        rows.append(r)
    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(OUT, "ranks.tsv"), sep="\t", index=False)
    S = {"common_genes": len(common), "genes_2020": len(kb["2020"].gene_terms), "genes_2026": len(kb["2026"].gene_terms)}
    for coh, d in list(R.groupby("cohort")) + [("external", R[R.cohort != "fu2022"])]:
        n = len(d); s = {"fetuses": n}
        for lab in ["2020", "2026", "2026onto_2020ann", "2026_fullterms"]:
            col = f"rank_{lab}"
            for k in [1, 10, 50]:
                x = int((d[col] <= k).sum()); lo, hi = proportion_confint(x, n, method="wilson")
                s[f"top{k}_{lab}"] = [x, round(100 * x / n, 1), round(100 * lo, 1), round(100 * hi, 1)]
            s[f"median_rank_{lab}"] = float(d[col].median())
        a, b = np.log10(d.rank_2020), np.log10(d.rank_2026)
        s["wilcoxon_p"] = float(wilcoxon(a, b).pvalue)
        s["improved"], s["worsened"] = int((b < a).sum()), int((b > a).sum())
        t20, t26 = d.rank_2020 <= 10, d.rank_2026 <= 10
        n01, n10 = int((~t20 & t26).sum()), int((t20 & ~t26).sum())
        s["mcnemar_top10"] = dict(gained=n01, lost=n10, p=float(binomtest(n01, n01 + n10).pvalue) if n01 + n10 else 1.0)
        s["unannotated_2020"], s["unannotated_2026"] = int((~d.annotated_2020).sum()), int((~d.annotated_2026).sum())
        S[coh] = s
    json.dump(S, open(os.path.join(OUT, "rank_summary.json"), "w"), indent=1)
    print(json.dumps(S, indent=1))
