"""Sensitivity analysis (PLAN.md, amendment 2026-10-15b): external cohorts without the 12 fetuses whose
diagnostic status is uncertain in the source. Reads outputs/ranks.tsv, writes outputs/sensitivity.json."""
import os, json
import pandas as pd
from scipy.stats import binomtest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
UNCERTAIN = {"aogs2026": {"7", "13", "15", "17", "19", "20", "21", "29", "32"}, "smogavec2022": {"17", "18", "22"}}

if __name__ == "__main__":
    R = pd.read_csv(os.path.join(ROOT, "outputs", "ranks.tsv"), sep="\t", dtype={"case": str})
    E = R[R.cohort != "fu2022"]
    keep = E[~E.apply(lambda r: r.case in UNCERTAIN.get(r.cohort, set()), axis=1)]
    t20, t26 = keep.rank_2020 <= 10, keep.rank_2026 <= 10
    g, l = int((~t20 & t26).sum()), int((t20 & ~t26).sum())
    out = dict(excluded=len(E) - len(keep), fetuses=len(keep), top10_2020=int(t20.sum()), top10_2026=int(t26.sum()),
               pct_2020=round(100 * t20.mean(), 1), pct_2026=round(100 * t26.mean(), 1),
               mcnemar=dict(gained=g, lost=l, p=float(binomtest(g, g + l).pvalue) if g + l else 1.0))
    json.dump(out, open(os.path.join(ROOT, "outputs", "sensitivity.json"), "w"), indent=1)
    print(out)
