"""Figure S1 of Supplementary Methods S1: one row per method step (annotation set and information content, score of a
gene, identical queries, gene symbol harmonisation, supplementation, paired comparison). The worked example is
primary-cohort fetus 96 (PKHD1, enlarged and hyperechogenic kidneys). All numbers are computed here or read from
outputs/. Drawn with plotly shapes in the style of the study flow (flowchart.py): one layout pixel prints at 0.5 pt.
Writes figures/s1/figure_s1.{html,pdf,png}."""
import os, ast, json, math
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from rank import KB, nearest_in, rank_of, PA

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT, FIG = os.path.join(ROOT, "outputs"), os.path.join(ROOT, "figures", "s1")
C = pio.templates["plotly_white"].layout.colorway
FS, SMALL = 16, 14
INK, GREY, DIV = "#111111", "#8C96A8", "#DCE1EA"
R20, R26, SUPP, PRE, COMP, GENE = C[2], C[0], C[1], C[3], C[4], C[5]
rgba = lambda h, a: f"rgba({int(h[1:3], 16)},{int(h[3:5], 16)},{int(h[5:7], 16)},{a})"

# ---------- numbers ----------
kb = {y: KB(y, y) for y in ("2020", "2026")}
common = set(kb["2020"].gene_terms) & set(kb["2026"].gene_terms)
CASE, GENE_EX = "96", "PKHD1"
cases = pd.read_csv(os.path.join(OUT, "cases.tsv"), sep="\t")
d = cases[(cases.cohort == "fu2022") & (cases.case == CASE)]
assert d.genes.iloc[0] == GENE_EX
q26 = list(dict.fromkeys(d.hpo_id))
q20 = [nearest_in(kb["2020"].terms, kb["2026"].anc, kb["2026"].terms, q) for q in q26]
ex = {}
for y, q in (("2020", q20), ("2026", q26)):
    k = kb[y]
    best = []
    for t in q:
        qa = k.anc[t] & k.anc_desc_pa()
        best.append(max(((k.ic.get(a, 0.0), a) for a in qa if a in k.gene_terms[GENE_EX]), default=(0.0, None)))
    sc = k.scores(q)
    s = sc.reindex(sorted(common))
    ex[y] = dict(best=best, score=sc[GENE_EX], rank=rank_of(sc, [GENE_EX], common),
                 tied=int((s.round(9) == round(sc[GENE_EX], 9)).sum()))
R = json.load(open(os.path.join(OUT, "rank_summary.json")))
assert (ex["2020"]["rank"], ex["2026"]["rank"]) == (28.0, 4.0)
k26 = kb["2026"]
# chain from the first query term to the root, through the parent with most ancestors (the deepest path)
chain = [q26[0]]
while chain[-1] != PA:
    chain.append(sorted((p for p in k26.terms[chain[-1]]["parents"] if p in k26.anc_desc_pa()), key=lambda p: -len(k26.anc[p]))[0])
chain = chain[::-1]
ngenes = lambda t: sum(t in s for s in k26.gene_terms.values())
want = ["Phenotypic abnormality", "Abnormality of the kidney", "Abnormal renal morphology", "Enlarged kidney"]
show = [t for t in chain if k26.terms[t]["name"] in want]
assert k26.terms[show[0]]["name"] == "Phenotypic abnormality"

P = {(p["comparison"], p["top"]): p for p in json.load(open(os.path.join(OUT, "paired_ci.json")))}
prim = P[("primary: 2026 vs 2020", 10)]
nc = P[("supplementation, causal gene only (noncompetitive)", 10)]
co = P[("supplementation, all genes of the primary cohort (competitive)", 10)]
GA = json.load(open(os.path.join(OUT, "gene_audit.json")))
audit = pd.read_csv(os.path.join(OUT, "gene_audit.tsv"), sep="\t")
renamed = audit[audit.source_symbol != audit.hgnc_symbol]
gba = audit[audit.source_symbol == "GBA"].iloc[0]
n_find = len(cases)
t20 = kb["2020"].terms
new_terms = cases[~cases.hpo_id.map(lambda h: h in t20 and not t20[h].get("obsolete"))]
EDEMA_NEW = "HP:0025672"   # Fetal skin edema, created after 2020
EDEMA_OLD = nearest_in(t20, k26.anc, k26.terms, EDEMA_NEW)
assert k26.terms[EDEMA_OLD]["name"] == "Edema"

# ---------- drawing helpers ----------
W_, TOP = 960, 64
HS = [190, 250, 182, 182, 206, 190]
ROWS = len(HS)
H_ = TOP + sum(HS) + 8
fig = go.Figure()
fig.update_xaxes(range=[0, W_], visible=False); fig.update_yaxes(range=[H_, 0], visible=False)


def text(x, y, t, size=FS, bold=False, col=INK, anchor="center"):
    fig.add_annotation(x=x, y=y, text=f"<b>{t}</b>" if bold else t, showarrow=False, xanchor=anchor, yanchor="middle",
                       font=dict(size=size, color=col))


def line(x0, y0, x1, y1, col, width=2, dash=None):
    fig.add_shape(type="line", x0=x0, y0=y0, x1=x1, y1=y1, line=dict(color=col, width=width, dash=dash))


def circle(x, y, r, col, fill=None, width=2, dash=None):
    fig.add_shape(type="circle", x0=x - r, x1=x + r, y0=y - r, y1=y + r, line=dict(color=col, width=width, dash=dash),
                  fillcolor=fill or "rgba(0,0,0,0)")


def rect(x0, y0, x1, y1, col, fill=None, width=1.5, dash=None):
    fig.add_shape(type="rect", x0=x0, x1=x1, y0=y0, y1=y1, line=dict(color=col, width=width, dash=dash),
                  fillcolor=fill or "rgba(0,0,0,0)")


def chip(x, y, t, col, w=None, size=SMALL, bold=False):
    w = w or 9 + 7.6 * len(t)
    fig.add_shape(type="rect", x0=x - w / 2, x1=x + w / 2, y0=y - 13, y1=y + 13, line=dict(color=col, width=1.5),
                  fillcolor=rgba(col, 0.12))
    text(x, y, t, size, bold)
    return w


def arrow(x0, y0, x1, y1, col, width=2, dash=None):
    fig.add_annotation(x=x1, y=y1, ax=x0, ay=y0, xref="x", yref="y", axref="x", ayref="y", text="", showarrow=True,
                       arrowhead=2, arrowsize=1, arrowwidth=width, arrowcolor=col)


def hbar(x0, y, val, vmax, length, col, label, h=13):
    w = length * val / vmax
    fig.add_shape(type="rect", x0=x0, x1=x0 + max(w, 1.5), y0=y - h / 2, y1=y + h / 2, line=dict(width=0), fillcolor=col)
    text(x0 + max(w, 1.5) + 6, y, label, SMALL, anchor="left")


# fetus icon shared with the study flow (flowchart.py), taken from its source so the drawing is identical
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "flowchart.py"), encoding="utf-8").read()
ns = {"math": math}
for node in ast.parse(src).body:
    tgt = node.targets[0] if isinstance(node, ast.Assign) else None
    if isinstance(tgt, ast.Tuple) and any(getattr(t, "id", "") in ("WALL", "SKIN") for t in tgt.elts):
        exec(compile(ast.Module([node], []), "flowchart", "exec"), ns)
    if isinstance(node, ast.FunctionDef) and node.name in ("ell", "fetus_icon"):
        exec(compile(ast.Module([node], []), "flowchart", "exec"), ns)
fetus = lambda x, y, s: ns["fetus_icon"](fig, x, y, s)

# ---------- layout: header and rows ----------
X1, X2, X3 = 14, 250, 642
text(X1, 26, "Step", 18, True, anchor="left")
text(X2, 26, "Illustration", 18, True, anchor="left")
text(X3, 26, "Example from the data", 18, True, anchor="left")
line(8, TOP - 12, W_ - 8, TOP - 12, DIV, 1.5)
rows = [("a", "Weight of a term", "information content", "Section 1", R26),
        ("b", "Match of a fetus", "with a gene", "Section 1", R26),
        ("c", "Identical queries", "", "Section 1", R20),
        ("d", "Gene symbol", "harmonisation", "Section 2", GENE),
        ("e", "Supplementation", "", "Section 3", SUPP),
        ("f", "Paired comparison", "", "Section 5", INK)]
Y = [TOP + sum(HS[:i]) for i in range(ROWS)]
for i, (k, t1, t2, sec, col) in enumerate(rows):
    y0, RH = Y[i], HS[i]
    fig.add_shape(type="rect", x0=X1, x1=X1 + 4, y0=y0 + 14, y1=y0 + RH - 18, line=dict(width=0), fillcolor=col)
    text(X1 + 16, y0 + 34, f"{k}  {t1}", FS, True, col, "left")
    if t2:
        text(X1 + 36, y0 + 56, t2, FS, True, col, "left")
    text(X1 + 36, y0 + (80 if t2 else 58), sec, SMALL, False, INK, "left")
    if i:
        line(8, y0, W_ - 8, y0, DIV, 1.5)
    line(X2 - 14, y0 + 12, X2 - 14, y0 + RH - 12, DIV, 1)
    line(X3 - 14, y0 + 12, X3 - 14, y0 + RH - 12, DIV, 1)

# a: weight of a term
y0 = Y[0]
labels = ["Phenotypic<br>abnormality", "Abnormality<br>of the<br>kidney", "Abnormal<br>renal<br>morphology", "Enlarged<br>kidney"]
assert [k26.terms[t]["name"] for t in show] == ["Phenotypic abnormality", "Abnormality of the kidney", "Abnormal renal morphology", "Enlarged kidney"]
xs = [290, 392, 494, 590]
yn = y0 + 100
for j in range(3):
    line(xs[j], yn, xs[j + 1], yn, R26, 2.5)
for j, (x, t, lab) in enumerate(zip(xs, show, labels)):
    circle(x, yn, 9, R26, R26, 0)
    text(x, yn - 36 - (8 if lab.count("<br>") == 2 else 0), lab, SMALL)
    text(x, yn + 30, f"{ngenes(t):,}", SMALL, True)
text(440, yn + 50, "genes carrying the term", SMALL)
arrow(300, y0 + 172, 580, y0 + 172, GREY, 1.5)
text(440, y0 + 184, "general to specific", 13)
text(X3, y0 + 30, "Weight, HPO 2026", SMALL, True, anchor="left")
for j, t in enumerate(show):
    yy = y0 + 58 + 28 * j
    text(X3, yy, labels[j].replace("<br>", " "), SMALL, anchor="left")
    v = abs(k26.ic.get(t, 0.0))
    hbar(X3 + 206, yy, v, 5, 64, R26, f"{v:.2f}")
text(X3, y0 + 172, f"Enlarged kidney: −ln({ngenes(show[-1])} / {len(k26.gene_terms):,}) = {abs(k26.ic[show[-1]]):.2f}", SMALL, anchor="left")

# b: match of fetus 96 with PKHD1, 2020 and 2026
y0 = Y[1]
fetus(282, y0 + 68, 0.56)
names = [k26.terms[t]["name"].replace("kidneys", "kidney") for t in q26]
text(282, y0 + 114, f"Fetus {CASE}", SMALL, True)
text(282, y0 + 146, "enlarged and<br>hyperechogenic<br>kidneys", 12)
parent = ex["2020"]["best"][1][1]
assert k26.terms[parent]["name"] == "Abnormal renal morphology"
text(484, y0 + 50, "Abnormal renal morphology", 12)
for cx, y, col in ((414, "2020", R20), (554, "2026", R26)):
    k = kb[y]
    py_, cy_ = y0 + 74, y0 + 140
    kids = [(cx - 40, q26[0]), (cx + 40, q26[1])]
    text(cx, y0 + 24, f"HPO {y}", SMALL, True, col)
    for kx, t in kids:
        line(cx, py_ + 9, kx, cy_ - 9, col, 2)
    circle(cx, py_, 9, col, col, 0)
    for jj, (kx, t) in enumerate(kids):
        t_y = q20[jj] if y == "2020" else t
        on = t_y in k.gene_terms[GENE_EX]
        circle(kx, cy_, 15, PRE, None, 2.5)
        circle(kx, cy_, 9, col, col if on else "white", 2.5)
        best = ex[y]["best"][jj]
        if best[1] == parent:
            fig.add_shape(type="path", path=f"M{kx + 13},{cy_ - 12} Q{kx + 22},{py_ + 4} {cx + 11},{py_ + 2}",
                          line=dict(color=PRE, width=2, dash="dot"), fillcolor="rgba(0,0,0,0)")
            text(kx + 16, (py_ + cy_) / 2 - 4, f"{best[0]:.2f}", 13, True, INK, "left")
        else:
            text(kx, cy_ + 26, f"{best[0]:.2f}", 13, True)
    text(cx - 40, cy_ + 50, "Enlarged", 12); text(cx + 40, cy_ + 50, "Hyper-<br>echogenic", 12)
ly = y0 + 214
circle(262, ly, 11, PRE, None, 2.5); text(280, ly, "fetal finding", 13, anchor="left")
circle(452, ly, 7, GREY, GREY, 0); text(466, ly, f"annotated to <i>{GENE_EX}</i>", 13, anchor="left")
text(250, ly + 24, "0.00: weight of the match", 13, anchor="left")
circle(452, ly + 24, 7, GREY, "white", 2); text(466, ly + 24, "not annotated", 13, anchor="left")
text(X3, y0 + 26, f"Match with <i>{GENE_EX}</i>", SMALL, True, anchor="left")
text(X3 + 228, y0 + 26, "2020", SMALL, True, R20); text(X3 + 286, y0 + 26, "2026", SMALL, True, R26)
for j, nm in enumerate(names + ["Score (mean)"]):
    yy = y0 + 56 + 30 * j
    text(X3, yy, nm, SMALL, j == 2, anchor="left")
    v20 = ex["2020"]["best"][j][0] if j < 2 else ex["2020"]["score"]
    v26 = ex["2026"]["best"][j][0] if j < 2 else ex["2026"]["score"]
    text(X3 + 228, yy, f"{v20:.2f}", SMALL, j == 2); text(X3 + 286, yy, f"{v26:.2f}", SMALL, j == 2)
yy = y0 + 56 + 96
text(X3, yy, f"Rank among {len(common):,} genes", SMALL, True, anchor="left")
text(X3 + 228, yy, f"{ex['2020']['rank']:.0f}", SMALL, True, R20); text(X3 + 286, yy, f"{ex['2026']['rank']:.0f}", SMALL, True, R26)
text(X3, y0 + 186, f"2026: {ex['2026']['tied']} genes share the top score,", 13, anchor="left")
text(X3, y0 + 206, f"mean rank {ex['2026']['rank']:.0f}", 13, anchor="left")

# c: identical queries
y0 = Y[2]
tx, ty = 290, y0 + 52
text(tx - 30, y0 + 20, "Finding coded in 2026", SMALL, True, PRE, "left")
circle(tx, ty, 10, R26, R26, 0); text(tx + 18, ty, "Edema", SMALL, anchor="left")
line(tx, ty + 10, tx, ty + 70, PRE, 2.5)
circle(tx, ty + 80, 10, PRE, "white", 2.5); text(tx + 18, ty + 80, "Fetal skin edema<br>created after 2020", 13, anchor="left")
fig.add_shape(type="path", path=f"M{tx - 14},{ty + 74} Q{tx - 44},{ty + 40} {tx - 14},{ty + 6}",
              line=dict(color=R20, width=2, dash="dot"), fillcolor="rgba(0,0,0,0)")
text(tx - 30, ty + 118, "nearest term present in 2020", 13, False, INK, "left")
chip(540, y0 + 52, "Query: Edema", R20, w=118, bold=True)
arrow(tx + 70, ty, 480, y0 + 52, GREY, 1.5)
arrow(540, y0 + 66, 506, y0 + 104, R20, 1.8); arrow(540, y0 + 66, 576, y0 + 104, R26, 1.8)
text(500, y0 + 118, "HPO 2020", SMALL, True, R20); text(584, y0 + 118, "HPO 2026", SMALL, True, R26)
text(X3, y0 + 26, f"Findings of the four cohorts ({n_find})", SMALL, True, anchor="left")
n_new = len(new_terms)
fig.add_shape(type="rect", x0=X3, x1=X3 + 280 * (n_find - n_new) / n_find, y0=y0 + 50, y1=y0 + 68, line=dict(width=0), fillcolor=R20)
fig.add_shape(type="rect", x0=X3 + 280 * (n_find - n_new) / n_find, x1=X3 + 280, y0=y0 + 50, y1=y0 + 68, line=dict(width=0), fillcolor=PRE)
text(X3, y0 + 92, f"{n_find - n_new} with a term present in 2020", SMALL, anchor="left")
circle(X3 + 6, y0 + 92, 0.1, R20)
text(X3, y0 + 118, f"{n_new} with a term created after 2020,", SMALL, anchor="left")
text(X3, y0 + 140, "queried with the nearest 2020 ancestor", SMALL, anchor="left")
text(X3, y0 + 166, "Same query for both releases", SMALL, True, anchor="left")

# d: gene symbol harmonisation
y0 = Y[3]
srcs = [("Source table", "GBA", y0 + 46), ("HPO 2020", "GBA", y0 + 104), ("HPO 2026", "GBA1", y0 + 150)]
for where, sym, yy in srcs:
    text(250, yy, where, SMALL, anchor="left")
    chip(378, yy, sym, GREY, w=58, bold=True)
chip(572, y0 + 60, f"HGNC: {gba.hgnc_symbol}", GENE, w=104, bold=True)
chip(490, y0 + 150, f"NCBI Gene {gba.ncbi_gene_id}", GENE, w=116)
arrow(408, y0 + 46, 518, y0 + 58, GREY, 1.5)
text(462, y0 + 36, "previous symbol", 12)
arrow(408, y0 + 104, 446, y0 + 137, GREY, 1.5)
arrow(408, y0 + 150, 430, y0 + 150, GREY, 1.5)
arrow(516, y0 + 136, 566, y0 + 75, GENE, 1.8)
text(X3, y0 + 26, "Causal genes with a previous symbol", SMALL, True, anchor="left")
for j, (_, r) in enumerate(renamed.sort_values("source_symbol").iterrows()):
    text(X3, y0 + 54 + 26 * j, f"<i>{r.source_symbol}</i> to <i>{r.hgnc_symbol}</i>", SMALL, anchor="left")
loci = GA["loci_without_hgnc_record"]
text(X3, y0 + 142, "Loci without an HGNC record removed", SMALL, True, anchor="left")
text(X3, y0 + 166, f"{len(loci['2020'])} in 2020, {len(loci['2026'])} in 2026", SMALL, anchor="left")


# e: supplementation
y0 = Y[4]
def supp_panel(cx, title, compete):
    text(cx, y0 + 20, title, SMALL, True, COMP if compete else SUPP)
    gy = y0 + 138
    gx = [cx - 58, cx, cx + 58]
    # fetus under evaluation, contributes to no gene
    fetus(cx - 58, y0 + 72, 0.36)
    circle(cx - 58, y0 + 73, 22, GREY, None, 1.5, "dot")
    for dx in (0, 58):
        fetus(cx + dx, y0 + 72, 0.36)
    arrow(cx - 6, y0 + 94, cx - 50, gy - 18, SUPP, 2)
    arrow(cx + 52, y0 + 94, cx - 44, gy - 14, SUPP, 2)
    for j, x in enumerate(gx):
        causal = j == 0
        col = SUPP if causal else GREY
        circle(x, gy, 15, col, rgba(col, 0.2), 2)
        if compete and not causal:
            fetus(x + 22, gy - 28, 0.22)
            arrow(x + 16, gy - 20, x + 8, gy - 13, COMP, 1.8)
    text(gx[0], gy + 30, "causal gene", 13)
    text((gx[1] + gx[2]) / 2, gy + 30, "other genes", 13)
    text(cx - 58, y0 + 41, "fetus ranked", 12)
    return gx

supp_panel(340, "Noncompetitive", False)
supp_panel(520, "Competitive", True)
text(X3, y0 + 26, f"Causal gene in the top 10, {nc['n']} fetuses", SMALL, True, anchor="left")
for j, (lab, v, col) in enumerate((("HPO 2026 annotations", nc["first"], R26), ("Competitive", co["second"], COMP),
                                   ("Noncompetitive", nc["second"], SUPP))):
    yy = y0 + 62 + 34 * j
    text(X3, yy, lab, SMALL, anchor="left")
    hbar(X3 + 186, yy, v, 40, 100, col, str(v), 15)
assert nc["first"] == co["first"]

# f: paired comparison
y0 = Y[5]
both, lost, gained = prim["first"] - prim["lost"], prim["lost"], prim["gained"]
neither = prim["n"] - both - lost - gained
gx0, gy0, cw, ch = 418, y0 + 62, 88, 44
text(gx0 + cw, y0 + 24, "HPO 2026 top 10", SMALL, True, R26)
text(gx0 + cw / 2, y0 + 46, "yes", SMALL); text(gx0 + 1.5 * cw, y0 + 46, "no", SMALL)
text(gx0 - 112, gy0 + ch, "HPO 2020<br>top 10", SMALL, True, R20)
text(gx0 - 14, gy0 + ch / 2, "yes", SMALL, anchor="right"); text(gx0 - 14, gy0 + 1.5 * ch, "no", SMALL, anchor="right")
cells = [(0, 0, both, None), (1, 0, lost, R20), (0, 1, gained, R26), (1, 1, neither, None)]
for cx_, cy_, v, col in cells:
    rect(gx0 + cx_ * cw, gy0 + cy_ * ch, gx0 + (cx_ + 1) * cw, gy0 + (cy_ + 1) * ch, col or GREY, rgba(col, 0.18) if col else None, 2 if col else 1)
    lab = {R20: "only 2020", R26: "only 2026"}.get(col)
    text(gx0 + (cx_ + 0.5) * cw, gy0 + (cy_ + 0.5) * ch - (8 if lab else 0), str(v), FS, bool(col))
    if lab:
        text(gx0 + (cx_ + 0.5) * cw, gy0 + (cy_ + 0.5) * ch + 12, lab, 12)
text(gx0 + cw, gy0 + 2 * ch + 22, f"{prim['n']} fetuses", SMALL)
text(X3, y0 + 26, "Difference in top-10 proportion", SMALL, True, anchor="left")
text(X3, y0 + 50, f"({gained} − {lost}) / {prim['n']} = {prim['diff_pp']:.1f} percentage points", SMALL, anchor="left")
ax0, ax1, ay = X3 + 10, X3 + 290, y0 + 116
v2x = lambda v: ax0 + (v + 5) / 10 * (ax1 - ax0)
line(ax0, ay + 22, ax1, ay + 22, INK, 1.2)
for v in range(-5, 6, 5):
    line(v2x(v), ay + 22, v2x(v), ay + 28, INK, 1.2); text(v2x(v), ay + 40, (f"{v:+d}" if v else "0").replace("-", "−"), SMALL)
line(v2x(0), ay - 18, v2x(0), ay + 22, GREY, 1, "dot")
line(v2x(prim["ci"][0]), ay, v2x(prim["ci"][1]), ay, INK, 2)
for v in prim["ci"]:
    line(v2x(v), ay - 6, v2x(v), ay + 6, INK, 2)
fig.add_shape(type="path", path=f"M{v2x(prim['diff_pp'])},{ay - 9} L{v2x(prim['diff_pp']) + 9},{ay} L{v2x(prim['diff_pp'])},{ay + 9} "
              f"L{v2x(prim['diff_pp']) - 9},{ay} Z", line=dict(width=0), fillcolor=INK)
text(X3, y0 + 80, f"Tango 95% CI {prim['ci'][0]:.1f} to {prim['ci'][1]:.1f}".replace("-", "−"), SMALL, anchor="left")

fig.update_layout(showlegend=False, plot_bgcolor="white", paper_bgcolor="white", width=W_ + 4, height=H_ + 4,
                  margin=dict(l=2, r=2, t=2, b=2), font=dict(size=FS, color=INK))
os.makedirs(FIG, exist_ok=True)
fig.write_html(os.path.join(FIG, "figure_s1.html"), include_plotlyjs="cdn")
fig.write_image(os.path.join(FIG, "figure_s1.pdf"))
fig.write_image(os.path.join(FIG, "figure_s1.png"), scale=600 / 144)
json.dump({"case": CASE, "gene": GENE_EX, "rank_2020": ex["2020"]["rank"], "rank_2026": ex["2026"]["rank"],
           "score_2020": round(ex["2020"]["score"], 2), "score_2026": round(ex["2026"]["score"], 2),
           "tied_2026": ex["2026"]["tied"], "findings": n_find, "findings_new_terms": n_new},
          open(os.path.join(OUT, "figure_s1.json"), "w"), indent=1)
print("figure S1 written", ex["2020"]["rank"], ex["2026"]["rank"], ex["2026"]["tied"], n_new)
