"""Annotation of the prenatal branch by clinical group and linkage of fetal terms to organ-system terms.
Group: terms under Fetal anomaly (structural and physiological abnormalities of the fetus) against the other
sub-branches of Abnormality of prenatal development or birth (placenta and cord, amniotic fluid, delivery, APGAR score,
gestational age at birth, embryonic development).
Linkage: a term is linked when at least one of its ancestors lies under Phenotypic abnormality outside the prenatal
branch, so that a query with the term also matches genes annotated with the corresponding organ-system term (for
example Fetal pleural effusion and Pleural effusion). Unlinked terms reach organ-system annotations only through the
prenatal branch. Reads outputs/prenatal_branch_terms.tsv and outputs/cases.tsv, writes outputs/annotation_groups.json
and outputs/fetal_term_linkage.tsv."""
import os, json
import pandas as pd
from rank import KB, PA, OUT

PRENATAL = "HP:0001197"


if __name__ == "__main__":
    k = KB("2026", "2026")
    B = pd.read_csv(os.path.join(OUT, "prenatal_branch_terms.tsv"), sep="\t")
    branch = {t for t, a in k.anc.items() if PRENATAL in a and not k.terms[t].get("obsolete")}
    assert set(B.term) <= branch
    organ = k.anc_desc_pa() - branch - {PA}
    B["no_gene"] = B.genes == 0
    B["group"] = (B.sub_branch == "Fetal anomaly").map({True: "Fetal anomaly", False: "Other prenatal sub-branches"})
    B["organ_ancestors"] = B.term.map(lambda t: sorted(k.anc[t] & organ))
    B["linked"] = B.organ_ancestors.map(bool)
    # nearest organ-system ancestor by name, for the table (several are possible, the most specific are listed)
    def most_specific(ts):
        ts = set(ts)
        return "; ".join(sorted(k.terms[t]["name"] for t in ts if not any(t in k.anc[o] and o != t for o in ts)))
    B["organ_terms"] = B.organ_ancestors.map(most_specific)

    out = {"groups": {}}
    for g, d in B.groupby("group"):
        out["groups"][g] = dict(terms=len(d), no_gene=int(d.no_gene.sum()), new_since_2020=int(d.new_since_2020.sum()),
                                new_no_gene=int((d.new_since_2020 & d.no_gene).sum()))
    F = B[B.group == "Fetal anomaly"].copy()
    out["fetal_anomaly_linkage"] = dict(
        terms=len(F), linked=int(F.linked.sum()), unlinked=int((~F.linked).sum()),
        unlinked_no_gene=int((~F.linked & F.no_gene).sum()), linked_no_gene=int((F.linked & F.no_gene).sum()),
        unlinked_by_organ=F[~F.linked].groupby(F.fetal_sub_branch.fillna("(none)")).size().sort_values(ascending=False).to_dict())
    # use of unlinked fetal terms in the coded findings of the four cohorts
    C = pd.read_csv(os.path.join(OUT, "cases.tsv"), sep="\t", dtype={"case": str})
    unl = set(F[~F.linked].term)
    used = C[C.hpo_id.isin(unl)]
    out["findings_with_unlinked_fetal_terms"] = dict(findings=len(used), fetuses=int(used[["cohort", "case"]].drop_duplicates().shape[0]),
                                                     terms=sorted(used.hpo_id.unique()))
    F[["term", "name", "fetal_sub_branch", "new_since_2020", "genes", "linked", "organ_terms"]] \
        .sort_values(["linked", "fetal_sub_branch", "name"]).to_csv(os.path.join(OUT, "fetal_term_linkage.tsv"), sep="\t", index=False)
    with open(os.path.join(OUT, "annotation_groups.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))
