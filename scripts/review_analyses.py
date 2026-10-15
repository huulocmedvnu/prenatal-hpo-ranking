"""Robustness analyses R1 to R5 (PLAN.md, amendment 2026-10-15d, all post hoc). Scores are computed from a
matrix of Resnik similarities between query terms and the direct annotation terms of the genes, which reproduces the
query-to-gene score of rank.py (checked below).
R1 prenatal re-coding, R2 symmetric Resnik and Lin best-match average, R3 unrestricted gene space, R4 mechanism of the
lower 2026 ranks, R5 cross-cohort supplementation, permutation control for X3 and all-gene supplementation (R5c).
Writes outputs/review.json."""
import os, json
import numpy as np, pandas as pd
from scipy.stats import binomtest, wilcoxon
from rank import KB, nearest_in, annotations

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "outputs")
rng = np.random.default_rng(2026)


def mcnemar(before, after):
    g, l = int((~before & after).sum()), int((before & ~after).sum())
    return dict(gained=g, lost=l, p=float(binomtest(g, g + l).pvalue) if g + l else 1.0)


def rank_in(scores, i):
    s = scores[i]
    return float((scores > s).sum() + ((scores == s).sum() + 1) / 2)


class Mat:
    """Resnik similarity between any query term and the direct annotation terms of the genes of one release."""
    def __init__(self, k, ann_release, universe):
        self.k, pa = k, k.anc_desc_pa()
        A = annotations(ann_release); A["hpo_id"] = A.hpo_id.map(k.res); A = A.dropna()
        A = A[A.hpo_id.isin(pa)]
        direct = {g: sorted(set(ts)) for g, ts in A.groupby("gene").hpo_id}
        self.genes = sorted(universe)
        self.terms = sorted({t for g in self.genes for t in direct[g]})
        ti = {t: i for i, t in enumerate(self.terms)}
        self.cols = [np.array([ti[t] for t in direct[g]]) for g in self.genes]
        self.gi = {g: i for i, g in enumerate(self.genes)}
        self.icmax = float(-np.log(1 / len(k.gene_terms)))
        self.ict = np.array([self.ic(t) for t in self.terms])
        self.cache = {}

    def ic(self, t):
        return self.k.ic.get(t, self.icmax)

    def resnik(self, q, terms):
        aq = self.k.anc.get(q, set()) & self.k.anc_desc_pa()
        icq = sorted(((self.k.ic.get(a, 0.0), a) for a in aq), reverse=True)
        return np.array([next((v for v, a in icq if a in self.k.anc[t]), 0.0) for t in terms])

    def pair(self, q, e):
        key = (q, e)
        if key not in self.cache:
            self.cache[key] = float(self.resnik(q, [e])[0])
        return self.cache[key]

    def extra_best(self, Q, E):
        return np.array([max(self.pair(q, e) for e in E) for q in Q])

    def score_one(self, Q, g, E=None):
        Q = [q for q in Q if q in self.k.anc]
        qg = np.vstack([self.row(q) for q in Q])[:, self.cols[self.gi[g]]].max(axis=1)
        if E:
            qg = np.maximum(qg, self.extra_best(Q, E))
        return qg.mean()

    def row(self, q):
        if q not in self.cache:
            self.cache[q] = self.resnik(q, self.terms)
        return self.cache[q]

    def scores(self, Q, metric="resnik", extra=None):
        """Score of every gene for the query terms Q. extra: {gene: [terms]} added to the annotation of those genes."""
        Q = [q for q in Q if q in self.k.anc]
        R = np.vstack([self.row(q) for q in Q])
        if metric == "lin":
            R = 2 * R / (np.array([self.ic(q) for q in Q])[:, None] + self.ict[None, :])
        out = np.empty(len(self.genes))
        for j, c in enumerate(self.cols):
            M = R[:, c]
            qg = M.max(axis=1)
            if extra and self.genes[j] in extra:
                E = extra[self.genes[j]]
                if E:
                    qg = np.maximum(qg, self.extra_best(Q, E))
            if metric == "symmetric":
                out[j] = (qg.mean() + M.max(axis=0).mean()) / 2
            else:
                out[j] = qg.mean()
        return out

    def rank(self, Q, genes, **kw):
        s = self.scores(Q, **kw)
        return min((rank_in(s, self.gi[g]) for g in genes if g in self.gi), default=float(len(s)))


if __name__ == "__main__":
    C = pd.read_csv(os.path.join(OUT, "cases.tsv"), sep="\t", dtype={"case": str})
    RK = pd.read_csv(os.path.join(OUT, "ranks.tsv"), sep="\t", dtype={"case": str})
    k20, k26 = KB("2020", "2020"), KB("2026", "2026")
    common = set(k20.gene_terms) & set(k26.gene_terms)
    m20, m26 = Mat(k20, "2020", common), Mat(k26, "2026", common)
    t20, t26 = k20.terms, k26.terms
    cases = {(c, i): (d.genes.iloc[0].split(";"), list(dict.fromkeys(d.hpo_id)), d.finding.tolist())
             for (c, i), d in C.groupby(["cohort", "case"], sort=False)}
    q20 = {key: list(dict.fromkeys(filter(None, (nearest_in(t20, k26.anc, t26, q) for q in v[1])))) for key, v in cases.items()}
    fu = [key for key in cases if key[0] == "fu2022"]
    out = {}
    # check: matrix score reproduces rank.py
    chk = [(m26.rank(q20[key], cases[key][0]), RK[(RK.cohort == key[0]) & (RK.case == key[1])].rank_2026.iloc[0]) for key in fu[:25]]
    out["check_max_abs_diff"] = float(max(abs(a - b) for a, b in chk))
    print("check", out["check_max_abs_diff"], flush=True)

    # R1 prenatal re-coding (2026, all cohorts, unmapped 2026 terms as the baseline)
    rec = pd.read_csv(os.path.join(ROOT, "inputs_manual", "prenatal_recoding.csv")).set_index("term").fetal_term.to_dict()
    rows = []
    for key, (genes, q26, _) in cases.items():
        if not any(q in rec for q in q26):
            continue
        qr = list(dict.fromkeys(rec.get(q, q) for q in q26))
        rows.append(dict(cohort=key[0], case=key[1], n_recoded=sum(q in rec for q in q26),
                         before=m26.rank(q26, genes), after=m26.rank(qr, genes)))
    R1 = pd.DataFrame(rows)
    b, a = R1.before <= 10, R1.after <= 10
    out["r1"] = dict(fetuses=len(R1), findings_recoded=int(R1.n_recoded.sum()), top10_before=int(b.sum()), top10_after=int(a.sum()),
                     mcnemar=mcnemar(b, a), median_before=float(R1.before.median()), median_after=float(R1.after.median()),
                     worse=int((R1.after > R1.before).sum()), better=int((R1.after < R1.before).sum()),
                     wilcoxon_p=float(wilcoxon(np.log10(R1.before), np.log10(R1.after)).pvalue))
    print("r1", out["r1"], flush=True)

    # R2 alternative scores, primary cohort
    r2 = {}
    for metric in ("symmetric", "lin"):
        r_20 = np.array([m20.rank(q20[key], cases[key][0], metric=metric) for key in fu])
        r_26 = np.array([m26.rank(q20[key], cases[key][0], metric=metric) for key in fu])
        r2[metric] = dict(top10_2020=int((r_20 <= 10).sum()), top10_2026=int((r_26 <= 10).sum()),
                          median_2020=float(np.median(r_20)), median_2026=float(np.median(r_26)),
                          mcnemar=mcnemar(r_20 <= 10, r_26 <= 10))
        print("r2", metric, r2[metric], flush=True)
    out["r2"] = r2

    # R3 unrestricted gene space
    n20, n26 = len(k20.gene_terms), len(k26.gene_terms)
    r3 = {}
    for name, d in [("fu2022", RK[RK.cohort == "fu2022"]), ("further", RK[RK.cohort != "fu2022"])]:
        p20, p26 = 100 * d.rank_all_2020 / n20, 100 * d.rank_all_2026 / n26
        r3[name] = dict(fetuses=len(d), genes_2020=n20, genes_2026=n26,
                        top10_2020=int((d.rank_all_2020 <= 10).sum()), top10_2026=int((d.rank_all_2026 <= 10).sum()),
                        mcnemar=mcnemar(d.rank_all_2020 <= 10, d.rank_all_2026 <= 10),
                        median_pct_2020=round(float(p20.median()), 2), median_pct_2026=round(float(p26.median()), 2),
                        only_2026=int((~d.annotated_2020 & d.annotated_2026).sum()), only_2020=int((d.annotated_2020 & ~d.annotated_2026).sum()),
                        neither=int((~d.annotated_2020 & ~d.annotated_2026).sum()))
    out["r3"] = r3
    print("r3", r3, flush=True)

    # R4 mechanism, primary cohort
    rows = []
    for key in fu:
        genes, Q = cases[key][0], q20[key]
        s20, s26 = m20.scores(Q), m26.scores(Q)
        g = min((x for x in genes if x in m26.gi), key=lambda x: rank_in(s26, m26.gi[x]), default=None)
        if g is None:
            continue
        i = m26.gi[g]
        rows.append(dict(case=key[1], score_2020=s20[i], score_2026=s26[i], above_2020=int((s20 > s20[i]).sum()),
                         above_2026=int((s26 > s26[i]).sum()), tied_2020=int((s20 == s20[i]).sum() - 1), tied_2026=int((s26 == s26[i]).sum() - 1),
                         qic_2020=float(np.mean([k20.ic.get(q, m20.icmax) for q in Q])), qic_2026=float(np.mean([k26.ic.get(q, m26.icmax) for q in Q])),
                         rank_2020=rank_in(s20, i), rank_2026=rank_in(s26, i)))
    R4 = pd.DataFrame(rows)
    R4.to_csv(os.path.join(OUT, "review_r4.tsv"), sep="\t", index=False)
    worse = R4.rank_2026 > R4.rank_2020
    out["r4"] = dict(fetuses=len(R4), worse=int(worse.sum()),
                     score_lower_2026=int((R4.score_2026 < R4.score_2020).sum()), score_higher_2026=int((R4.score_2026 > R4.score_2020).sum()),
                     median_score_2020=round(float(R4.score_2020.median()), 3), median_score_2026=round(float(R4.score_2026.median()), 3),
                     median_above_2020=float(R4.above_2020.median()), median_above_2026=float(R4.above_2026.median()),
                     median_query_ic_2020=round(float(R4.qic_2020.median()), 3), median_query_ic_2026=round(float(R4.qic_2026.median()), 3),
                     worse_score_lower=int((worse & (R4.score_2026 < R4.score_2020)).sum()),
                     worse_more_above=int((worse & (R4.above_2026 > R4.above_2020)).sum()),
                     wilcoxon_score_p=float(wilcoxon(R4.score_2020, R4.score_2026).pvalue))
    print("r4", out["r4"], flush=True)

    # R5 supplementation
    by_gene = {}
    for key in fu:
        for g in cases[key][0]:
            by_gene.setdefault(g, []).append(key)
    terms_of = lambda keys: sorted({q for k_ in keys for q in cases[k_][1] if q in k26.anc})
    # R5a cross-cohort: further fetuses whose gene is also diagnosed in the primary cohort
    rows = []
    for key in [k_ for k_ in cases if k_[0] != "fu2022"]:
        genes, q26 = cases[key][0], cases[key][1]
        shared = [g for g in genes if g in by_gene and g in m26.gi]
        if not shared:
            continue
        extra = {g: terms_of(by_gene[g]) for g in shared}
        rows.append(dict(cohort=key[0], case=key[1], before=m26.rank(q26, genes), after=m26.rank(q26, genes, extra=extra)))
    R5 = pd.DataFrame(rows)
    out["r5a"] = dict(fetuses=len(R5), top10_before=int((R5.before <= 10).sum()), top10_after=int((R5.after <= 10).sum()),
                      mcnemar=mcnemar(R5.before <= 10, R5.after <= 10), median_before=float(R5.before.median()), median_after=float(R5.after.median()),
                      better=int((R5.after < R5.before).sum()), worse=int((R5.after > R5.before).sum()))
    print("r5a", out["r5a"], flush=True)
    # X3 reproduction, permutation control and all-gene supplementation (primary cohort, 2026 unmapped terms as in X3)
    x3 = [(key, g) for g, ks in by_gene.items() if len(ks) >= 2 and g in m26.gi for key in ks]
    base = {key: m26.scores(cases[key][1]) for key, _ in x3}
    def rank_with(key, g, E):
        s = base[key].copy(); i = m26.gi[g]
        s[i] = m26.score_one(cases[key][1], g, E) if E else s[i]
        return rank_in(s, i)
    obs = np.array([rank_with(key, g, terms_of([o for o in by_gene[g] if o != key])) for key, g in x3])
    before = np.array([rank_in(base[key], m26.gi[g]) for key, g in x3])
    pool = {key: [q for q in cases[key][1] if q in k26.anc] for key in fu}
    sizes = [len(terms_of([o for o in by_gene[g] if o != key])) for key, g in x3]
    others = {(key, g): [q for o in fu if g not in cases[o][0] for q in pool[o]] for key, g in x3}
    n_perm, perm = 1000, []
    # the rank of a gene with extra terms depends only on the best Resnik match of each query term, so cache rows
    for _ in range(n_perm):
        top = 0
        for (key, g), n in zip(x3, sizes):
            E = sorted(set(rng.choice(others[(key, g)], size=n, replace=False))) if n else []
            top += rank_with(key, g, E) <= 10
        perm.append(top)
    perm = np.array(perm)
    out["r5b"] = dict(fetuses=len(x3), top10_before=int((before <= 10).sum()), top10_observed=int((obs <= 10).sum()),
                      perm_mean=round(float(perm.mean()), 1), perm_q95=float(np.percentile(perm, 95)),
                      perm_p=float((1 + (perm >= (obs <= 10).sum()).sum()) / (n_perm + 1)), n_perm=n_perm)
    print("r5b", out["r5b"], flush=True)
    # R5c all genes supplemented, leave the fetus under evaluation out
    rows = []
    for key, g in x3:
        extra = {h: terms_of([o for o in ks if o != key]) for h, ks in by_gene.items() if h in m26.gi}
        extra = {h: E for h, E in extra.items() if E}
        s = m26.scores(cases[key][1], extra=extra)
        rows.append(rank_in(s, m26.gi[g]))
    r5c = np.array(rows)
    out["r5c"] = dict(fetuses=len(x3), top10_before=int((before <= 10).sum()), top10_after=int((r5c <= 10).sum()),
                      mcnemar=mcnemar(before <= 10, r5c <= 10), median_before=float(np.median(before)), median_after=float(np.median(r5c)))
    print("r5c", out["r5c"], flush=True)
    json.dump(out, open(os.path.join(OUT, "review.json"), "w"), indent=1)
