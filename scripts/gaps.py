"""Secondary outcomes (PLAN), primary (Fu) cohort: annotation gaps, ranking by mode of inheritance and leakage check.
For each finding (query term as used in the 2020 and 2026 runs): does the causal gene carry the term or a descendant
(exact match), and does it carry any term in the same organ system (top-level branch under Phenotypic abnormality)?
By organ system: among fetuses with a finding in the system, the share whose causal gene carries no term of that
system in 2026. By inheritance: top-10 counts with both releases for autosomal dominant, autosomal recessive and
X-linked disorders. Leakage: genes whose 2026 phenotype.hpoa annotations cite PMID:36307859, and the top-10 shares of
the primary cohort without the fetuses diagnosed with them. Writes outputs/gaps.tsv, outputs/gaps.json."""
import os, json
import pandas as pd
from rank import KB, nearest_in, PA
from genes import symbol_from_id
from statsmodels.stats.proportion import proportion_confint
from scipy.stats import binomtest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IN, OUT = os.path.join(ROOT, "inputs"), os.path.join(ROOT, "outputs")


def read_hpoa(path):
    """phenotype.hpoa with its leading comment lines skipped ('#' also occurs inside some disease names)."""
    with open(path, encoding="utf-8") as f:
        skip = 0
        while f.readline().startswith("#"):
            skip += 1
    return pd.read_csv(path, sep="\t", skiprows=skip, dtype=str)


if __name__ == "__main__":
    C = pd.read_csv(os.path.join(OUT, "cases.tsv"), sep="\t", dtype={"case": str})
    C = C[C.cohort == "fu2022"]
    R = pd.read_csv(os.path.join(OUT, "ranks.tsv"), sep="\t", dtype={"case": str})
    R = R[R.cohort == "fu2022"]
    kb = {"2020": KB("2020", "2020"), "2026": KB("2026", "2026")}
    top = {c for c, d in kb["2026"].terms.items() if PA in d["parents"]}
    t20, t26 = kb["2020"].terms, kb["2026"].terms
    rows, per_system = [], []
    for (coh, case), d in C.groupby(["cohort", "case"], sort=False):
        genes = d.genes.iloc[0].split(";")
        for q26 in dict.fromkeys(d.hpo_id):
            q = nearest_in(t20, kb["2026"].anc, t26, q26)
            branch = kb["2026"].anc[q26] & top
            r = dict(cohort=coh, case=case, term=q26, branch=";".join(sorted(t26[b]["name"] for b in branch)))
            for lab, k in kb.items():
                gt = set().union(*[k.gene_terms.get(g, set()) for g in genes])
                desc = {t for t, a in k.anc.items() if q in a}
                r[f"exact_{lab}"] = bool(gt & desc)
                r[f"system_{lab}"] = bool(gt & {t for t in gt if k.anc[t] & branch})
            # each organ system of the finding separately: does the causal gene carry any term of that system?
            gt26 = set().union(*[kb["2026"].gene_terms.get(g, set()) for g in genes])
            per_system += [dict(case=case, system=t26[b]["name"], covered=any(b in kb["2026"].anc[t] for t in gt26)) for b in branch]
            rows.append(r)
    G = pd.DataFrame(rows); G.to_csv(os.path.join(OUT, "gaps.tsv"), sep="\t", index=False)
    F = G.groupby(["cohort", "case"]).agg(**{f"{m}_{l}": (f"{m}_{l}", "any") for m in ["exact", "system"] for l in ["2020", "2026"]})
    out = {"findings": len(G), "fetuses": len(F)}
    for col in ["exact_2020", "exact_2026", "system_2020", "system_2026"]:
        x = int(G[col].sum()); out[f"finding_{col}"] = [x, round(100 * x / len(G), 1)]
        y = int((~F[col]).sum()); lo, hi = proportion_confint(y, len(F), method="wilson")
        out[f"fetus_none_{col}"] = [y, round(100 * y / len(F), 1), round(100 * lo, 1), round(100 * hi, 1)]
    out["by_branch_exact_2026"] = G.assign(b=G.branch.str.split(";").str[0]).groupby("b").exact_2026.agg(["size", "mean"]).round(3).sort_values("size", ascending=False).head(10).to_dict("index")
    # fetuses with at least one finding in the system, and among them those whose causal gene carries no term of it
    sysf = pd.DataFrame(per_system).groupby(["system", "case"]).covered.any().reset_index()
    out["by_system_fetus_none_2026"] = {b: dict(fetuses=len(d), none=int((~d.covered).sum()),
                                                pct=round(100 * float((~d.covered).mean()), 1))
                                        for b, d in sysf.groupby("system") if len(d) >= 5}
    inh = C.drop_duplicates(["cohort", "case"]).set_index("case").inheritance.astype(str)
    Ri = R.assign(group=R.case.astype(str).map(lambda c: "X-linked" if inh[c].startswith("X") else inh[c]))
    out["by_inheritance"] = {}
    for gname, d in Ri.groupby("group"):
        a, b = d.rank_2020 <= 10, d.rank_2026 <= 10
        g_, l_ = int((~a & b).sum()), int((a & ~b).sum())
        out["by_inheritance"][gname] = dict(fetuses=len(d), top10_2020=int(a.sum()), top10_2026=int(b.sum()),
                                            gained=g_, lost=l_, p=float(binomtest(g_, g_ + l_).pvalue) if g_ + l_ else 1.0)
    H = read_hpoa(os.path.join(IN, "hpo2026", "phenotype.hpoa"))
    dis = set(H.database_id[H.reference.fillna("").str.contains("36307859")])
    g2p = pd.read_csv(os.path.join(IN, "hpo2026", "genes_to_phenotype.txt"), sep="\t", usecols=["ncbi_gene_id", "disease_id"], dtype=str)
    leak = sorted(set(g2p.ncbi_gene_id[g2p.disease_id.isin(dis)].map(symbol_from_id)))
    keep = ~R.genes.str.split(";").apply(lambda gs: any(g in leak for g in gs))
    Rk = R[keep]
    out["leakage"] = dict(diseases=len(dis), genes=leak, fetuses_removed=int((~keep).sum()),
                          top10_2020=round(100 * float((Rk.rank_2020 <= 10).mean()), 1), top10_2026=round(100 * float((Rk.rank_2026 <= 10).mean()), 1))
    with open(os.path.join(OUT, "gaps.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))
