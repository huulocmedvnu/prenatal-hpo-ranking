"""Fu 2022 Table S3 (Additional file 2): one row per diagnosed fetus with ultrasound findings and genes. Findings
are split on commas and matched to HPO 2026 by exact label or synonym (case-insensitive). Hand mappings come from
inputs_manual/term_map.csv. Writes outputs/cases.tsv (one row per fetus and finding) and outputs/unmatched_terms.tsv."""
import os, re
import pandas as pd, docx
from hpo import load_obo, label_index

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IN, OUT = os.path.join(ROOT, "inputs"), os.path.join(ROOT, "outputs")


def table_s3():
    t = docx.Document(os.path.join(IN, "fu2022", "additional_file_2.docx")).tables[2]
    hdr = [c.text.strip() for c in t.rows[0].cells]
    rows, sect, last = [], "", None
    for r in t.rows[1:]:
        cells = [c.text.strip() for c in r.cells]
        if len(set(cells)) == 1 or (cells[0] and all(c in ("", cells[0]) for c in cells)):
            sect = cells[0]; continue
        d = dict(zip(hdr, cells)); d["section"] = sect
        if not d["Case ID"]:
            d["Case ID"] = last["Case ID"]; d["Ultrasound findings"] = d["Ultrasound findings"] or last["Ultrasound findings"]
        last = d; rows.append(d)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    T = table_s3()
    T["gene"] = T["Gene(OMIM ID)"].str.replace(r"\(.*", "", regex=True).str.replace("*", "", regex=False).str.strip()
    T = T[T.gene != ""]
    terms = load_obo(os.path.join(IN, "hpo2026", "hp.obo")); idx = label_index(terms)
    hand = {}
    p = os.path.join(ROOT, "inputs_manual", "term_map.csv")
    if os.path.exists(p):
        H = pd.read_csv(p); hand = dict(zip(H.finding.str.lower(), H.hpo_id))
    rows, un = [], []
    for case, d in T.groupby("Case ID", sort=False):
        genes = ";".join(sorted(set(d.gene)))
        finds = [f.strip(" .") for f in re.split(r"[,;，]", d["Ultrasound findings"].iloc[0]) if f.strip(" .")]
        for f in finds:
            hid = idx.get(f.lower()) or hand.get(f.lower())
            if hid is None:
                un.append(dict(case=case, finding=f))
            rows.append(dict(cohort="fu2022", case=case, section=d.section.iloc[0], genes=genes,
                             inheritance=d.Inheritance.iloc[0], finding=f, hpo_id=hid))
    C = pd.DataFrame(rows)
    E = pd.read_csv(os.path.join(ROOT, "inputs_manual", "external_cases.csv"))
    E = E.rename(columns={"genes": "genes"}).assign(section="", case=E.case.astype(str))
    C = pd.concat([C, E[["cohort", "case", "section", "genes", "inheritance", "finding", "hpo_id"]]], ignore_index=True)
    C.to_csv(os.path.join(OUT, "cases.tsv"), sep="\t", index=False)
    U = pd.DataFrame(un)
    (U.groupby("finding").size().rename("n").reset_index() if len(U) else pd.DataFrame(columns=["finding", "n"])).to_csv(
        os.path.join(OUT, "unmatched_terms.tsv"), sep="\t", index=False)
    print(C.groupby("cohort").case.nunique().to_dict(), "fetuses", C.case.nunique(), "findings", len(C), "matched", int(C.hpo_id.notna().sum()), "genes", T.gene.nunique())
    print(C.section.value_counts().to_string()); print(U.finding.value_counts().to_string() if len(U) else "all matched")
