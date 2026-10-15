"""Secondary outcomes (PLAN), primary (Fu) cohort: annotation gaps and leakage check. For each finding (query term as used in the 2020
and 2026 runs): does the causal gene carry the term or a descendant (exact match), and does it carry any term in the
same organ system (top-level branch under Phenotypic abnormality)? Leakage: genes whose 2026 phenotype.hpoa
annotations cite PMID:36307859, primary ranks repeated without them. Writes outputs/gaps.tsv, outputs/gaps.json."""
import os, json
import numpy as np, pandas as pd
from rank import KB, nearest_in, PA
from statsmodels.stats.proportion import proportion_confint

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IN, OUT = os.path.join(ROOT, "inputs"), os.path.join(ROOT, "outputs")

if __name__ == "__main__":
    C = pd.read_csv(os.path.join(OUT, "cases.tsv"), sep="\t")
    C = C[C.cohort == "fu2022"]
    R = pd.read_csv(os.path.join(OUT, "ranks.tsv"), sep="\t")
    R = R[R.cohort == "fu2022"]
    kb = {"2020": KB("2020", "2020"), "2026": KB("2026", "2026")}
    top = {c for c, d in kb["2026"].terms.items() if PA in d["parents"]}
    t20, t26 = kb["2020"].terms, kb["2026"].terms
    rows = []
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
            rows.append(r)
    G = pd.DataFrame(rows); G.to_csv(os.path.join(OUT, "gaps.tsv"), sep="\t", index=False)
    F = G.groupby(["cohort", "case"]).agg(**{f"{m}_{l}": (f"{m}_{l}", "any") for m in ["exact", "system"] for l in ["2020", "2026"]})
    out = {"findings": len(G), "fetuses": len(F)}
    for col in ["exact_2020", "exact_2026", "system_2020", "system_2026"]:
        x = int(G[col].sum()); out[f"finding_{col}"] = [x, round(100 * x / len(G), 1)]
        y = int((~F[col]).sum()); lo, hi = proportion_confint(y, len(F), method="wilson")
        out[f"fetus_none_{col}"] = [y, round(100 * y / len(F), 1), round(100 * lo, 1), round(100 * hi, 1)]
    out["by_branch_exact_2026"] = G.assign(b=G.branch.str.split(";").str[0]).groupby("b").exact_2026.agg(["size", "mean"]).round(3).sort_values("size", ascending=False).head(10).to_dict("index")
    H = pd.read_csv(os.path.join(IN, "hpo2026", "phenotype.hpoa"), sep="\t", comment="#", dtype=str)
    dis = set(H.database_id[H.reference.fillna("").str.contains("36307859")])
    g2p = pd.read_csv(os.path.join(IN, "hpo2026", "genes_to_phenotype.txt"), sep="\t", usecols=["gene_symbol", "disease_id"])
    leak = sorted(set(g2p.gene_symbol[g2p.disease_id.isin(dis)]))
    keep = ~R.genes.str.split(";").apply(lambda gs: any(g in leak for g in gs))
    Rk = R[keep]
    out["leakage"] = dict(diseases=len(dis), genes=leak, fetuses_removed=int((~keep).sum()),
                          top10_2020=round(100 * float((Rk.rank_2020 <= 10).mean()), 1), top10_2026=round(100 * float((Rk.rank_2026 <= 10).mean()), 1))
    json.dump(out, open(os.path.join(OUT, "gaps.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))
