"""Agreement between the first and the second independent extraction (PLAN.md, amendment 2026-10-15).
External cohorts (first extraction as saved before adjudication, second_extraction/external_cases_1_before_adjudication.csv, against second_extraction/external_cases_2.csv):
1. fetus: case present in both, and causal gene set identical;
2. HPO term, within fetuses present in both: each term of one extraction is identical to a term of the other,
   related (one an ancestor of the other), or unmatched.
Fu hand-mapped labels (term_map.csv against second_extraction/term_map_2.csv): identical, related or different term.
Writes outputs/agreement.json and outputs/agreement_discrepancies.tsv (the list for clinician adjudication)."""
import os, re, json
import pandas as pd
from hpo import load_obo, ancestors_map

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAN, OUT = os.path.join(ROOT, "inputs_manual"), os.path.join(ROOT, "outputs")


def norm_case(x):
    """Case label without words or footnote marks (Smogavec prints 12*, 14a, 19#)."""
    m = re.search(r"\d+", str(x))
    return m.group(0) if m else str(x).strip().lower()


def genes(s):
    return frozenset(g.strip().upper() for g in str(s).replace(",", ";").split(";") if g.strip())


def relation(t, others, anc):
    if t in others:
        return "identical"
    if any(o in anc.get(t, ()) or t in anc.get(o, ()) for o in others):
        return "related"
    return "unmatched"


if __name__ == "__main__":
    terms = load_obo(os.path.join(ROOT, "inputs", "hpo2026", "hp.obo"))
    anc = ancestors_map(terms)
    name = lambda t: terms.get(t, {}).get("name", "")
    A = pd.read_csv(os.path.join(MAN, "second_extraction", "external_cases_1_before_adjudication.csv"))
    B = pd.read_csv(os.path.join(MAN, "second_extraction", "external_cases_2.csv"))
    for d in (A, B):
        d["key"] = d.cohort + "|" + d.case.map(norm_case)
    res, disc = {}, []
    for coh in sorted(set(A.cohort) | set(B.cohort)):
        a, b = A[A.cohort == coh], B[B.cohort == coh]
        ka, kb = set(a.key), set(b.key)
        both = sorted(ka & kb)
        for k in sorted(ka - kb):
            disc.append(dict(level="fetus", cohort=coh, case=k.split("|")[1], first="present", second="absent"))
        for k in sorted(kb - ka):
            disc.append(dict(level="fetus", cohort=coh, case=k.split("|")[1], first="absent", second="present"))
        gene_same, cnt = 0, {"identical": 0, "related": 0, "unmatched": 0}
        for k in both:
            ga, gb = genes(a[a.key == k].genes.iloc[0]), genes(b[b.key == k].genes.iloc[0])
            gene_same += ga == gb
            if ga != gb:
                disc.append(dict(level="gene", cohort=coh, case=k.split("|")[1], first=";".join(sorted(ga)), second=";".join(sorted(gb))))
            ta, tb = set(a[a.key == k].hpo_id), set(b[b.key == k].hpo_id)
            for x, y, side in [(ta, tb, "first"), (tb, ta, "second")]:
                for t in sorted(x):
                    r = relation(t, y, anc)
                    cnt[r] += 1
                    if r != "identical":
                        src = (a if side == "first" else b)
                        f = src[(src.key == k) & (src.hpo_id == t)].finding.iloc[0]
                        disc.append(dict(level=f"term ({r})", cohort=coh, case=k.split("|")[1], only_in=side, finding=f,
                                         hpo_id=t, hpo_name=name(t),
                                         other_terms="; ".join(f"{o} {name(o)}" for o in sorted(y))))
        n = sum(cnt.values())
        res[coh] = dict(fetuses_first=len(ka), fetuses_second=len(kb), fetuses_both=len(both),
                        gene_identical=int(gene_same), term_mentions=n,
                        **{f"terms_{k}": v for k, v in cnt.items()},
                        pct_terms_identical=round(100 * cnt["identical"] / n, 1) if n else None,
                        pct_terms_identical_or_related=round(100 * (cnt["identical"] + cnt["related"]) / n, 1) if n else None)
    F, G = pd.read_csv(os.path.join(MAN, "second_extraction", "term_map_1_before_adjudication.csv")), pd.read_csv(os.path.join(MAN, "second_extraction", "term_map_2.csv"))
    m = F.merge(G, on="finding", how="outer", suffixes=("_first", "_second"))
    fu = {"identical": 0, "related": 0, "unmatched": 0}
    for _, r in m.iterrows():
        rel = relation(r.hpo_id_first, {r.hpo_id_second}, anc) if pd.notna(r.hpo_id_first) and pd.notna(r.hpo_id_second) else "unmatched"
        fu[rel] += 1
        if rel != "identical":
            disc.append(dict(level=f"fu label ({rel})", cohort="fu2022", finding=r.finding,
                             first=f"{r.hpo_id_first} {name(r.hpo_id_first)}", second=f"{r.hpo_id_second} {name(r.hpo_id_second)}"))
    res["fu2022_labels"] = dict(labels=len(m), **fu)
    json.dump(res, open(os.path.join(OUT, "agreement.json"), "w"), indent=1)
    D = pd.DataFrame(disc)
    D.to_csv(os.path.join(OUT, "agreement_discrepancies.tsv"), sep="\t", index=False)
    print(json.dumps(res, indent=1)); print(D.level.value_counts().to_string() if len(D) else "no discrepancies")
