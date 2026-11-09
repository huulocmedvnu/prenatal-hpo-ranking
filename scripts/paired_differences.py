"""Paired absolute differences in top-k proportions with Tango score 95% confidence intervals.
Second minus first, in percentage points. Reads outputs/ranks.tsv and outputs/fetus_ranks.tsv (written by robustness.py), writes
outputs/paired_ci.json and outputs/paired_ci.tsv."""
import os, json
import pandas as pd
from scipy.stats import binomtest
from paired_ci import paired
from sensitivity import UNCERTAIN

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "outputs")
KS = (1, 5, 10, 20, 50)


def row(name, k, a, b):
    r = paired(a <= k, b <= k)
    g, l = r["gained"], r["lost"]
    return dict(comparison=name, top=k, **r, mcnemar_p=float(binomtest(g, g + l).pvalue) if g + l else 1.0)


if __name__ == "__main__":
    R = pd.read_csv(os.path.join(OUT, "ranks.tsv"), sep="\t", dtype={"case": str})
    V = pd.read_csv(os.path.join(OUT, "fetus_ranks.tsv"), sep="\t", dtype={"case": str})
    fu, ext = R[R.cohort == "fu2022"], R[R.cohort != "fu2022"]
    ext94 = ext[~ext.apply(lambda r: r.case in UNCERTAIN.get(r.cohort, set()), axis=1)]
    rows = []
    # the Yao cohort used fetal MRI as well as ultrasound, and the source does not separate the two
    ext_us = ext[ext.cohort != "aogs2026"]
    for name, d in (("primary", fu), ("further", ext), ("further_certain", ext94), ("further_without_yao", ext_us)):
        rows += [row(f"{name}: 2026 vs 2020", k, d.rank_2020, d.rank_2026) for k in KS]
    rows.append(row("primary: 2026 original terms vs 2020", 10, fu.rank_2020, fu.rank_2026_fullterms))
    rows.append(row("primary: 2026 ontology with 2020 annotations vs 2020", 10, fu.rank_2020, fu.rank_2026onto_2020ann))
    rows.append(row("primary: all genes of each release, 2026 vs 2020", 10, fu.rank_all_2020, fu.rank_all_2026))
    rows.append(row("further: all genes of each release, 2026 vs 2020", 10, ext.rank_all_2020, ext.rank_all_2026))
    for a, name in (("score_symmetric", "primary: symmetric Resnik, 2026 vs 2020"), ("score_lin", "primary: Lin, 2026 vs 2020"),
                    ("supp_causal_only", "supplementation, causal gene only (noncompetitive)"),
                    ("supp_all_genes", "supplementation, all genes of the primary cohort (competitive)"),
                    ("supp_further", "supplementation from primary cohort, further cohorts")):
        d = V[V.analysis == a]
        rows.append(row(name, 10, d["first"], d["second"]))
    T = pd.DataFrame(rows)
    T.to_csv(os.path.join(OUT, "paired_ci.tsv"), sep="\t", index=False)
    with open(os.path.join(OUT, "paired_ci.json"), "w", encoding="utf-8") as f:
        json.dump(rows, f, indent=1)
    pd.set_option("display.width", 200)
    print(T.to_string(index=False))
