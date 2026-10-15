"""Study flow schematic, drawn with plotly shapes: published cohorts, the two HPO releases with an identical query,
and the analyses. All counts from outputs/. One layout pixel prints at 0.5 pt, so 16 px fonts print at 8 pt.
Writes figures/study_flow/study_flow.{html,pdf,png}."""
import os, json, math
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT, FIG = os.path.join(ROOT, "outputs"), os.path.join(ROOT, "figures", "study_flow")
C = pio.templates["plotly_white"].layout.colorway
FS = 16
INK, TGREY, GREY, DIV = "#2A3F5F", "#5A6478", "#8C96A8", "#DCE1EA"
# colour = meaning, shared with figures/fig_main: HPO 2026 blue, supplemented annotations red, prenatal branch purple
R20, R26, SUPP, PRE, COH, RANK = C[2], C[0], C[1], C[3], C[4], C[5]

R = json.load(open(os.path.join(OUT, "rank_summary.json")))
X = json.load(open(os.path.join(OUT, "explore.json")))
cases = pd.read_csv(os.path.join(OUT, "cases.tsv"), sep="\t")
per = cases.groupby("cohort").agg(findings=("hpo_id", "size"), fetuses=("case", "nunique"))
fu_genes = cases[cases.cohort == "fu2022"].genes.str.split(";").explode().nunique()
n_all = int(per.fetuses.sum())

W_, H_ = 960, 552
fig = go.Figure()
fig.update_xaxes(range=[0, W_], visible=False); fig.update_yaxes(range=[H_, 0], visible=False)
rgba = lambda h, a: f"rgba({int(h[1:3], 16)},{int(h[3:5], 16)},{int(h[5:7], 16)},{a})"


def text(x, y, t, size=FS, bold=False, col=INK, anchor="center"):
    fig.add_annotation(x=x, y=y, text=f"<b>{t}</b>" if bold else t, showarrow=False, xanchor=anchor, yanchor="middle",
                       font=dict(size=size, color=col))


def path(d, line, fill=None, width=2):
    fig.add_shape(type="path", path=d, line=dict(color=line, width=width), fillcolor=fill or "rgba(0,0,0,0)")


def circle(x, y, r, line, fill=None, width=2):
    fig.add_shape(type="circle", x0=x - r, x1=x + r, y0=y - r, y1=y + r, line=dict(color=line, width=width),
                  fillcolor=fill or "rgba(0,0,0,0)")


def line(x0, y0, x1, y1, col, width=2, dash=None):
    fig.add_shape(type="line", x0=x0, y0=y0, x1=x1, y1=y1, line=dict(color=col, width=width, dash=dash))


def chevron(x, y):
    path(f"M{x},{y - 15} L{x + 18},{y} L{x},{y + 15}", GREY, None, 3)


WALL, WALL_D, FLUID, FLUID_SPOT = "#FCD9C0", "#F08A4B", "#8C3B1E", "#A8502C"
SKIN, SKIN_D, SKIN_L, CORD, LINE = "#FBD3B8", "#F3AE86", "#FFF1E6", "#E0603A", "#D9763F"


def ell(cx, cy, rx, ry, ang=0, n=40):
    a = math.radians(ang)
    pts = [(cx + rx * math.cos(t) * math.cos(a) - ry * math.sin(t) * math.sin(a),
            cy + rx * math.cos(t) * math.sin(a) + ry * math.sin(t) * math.cos(a)) for t in [2 * math.pi * k / n for k in range(n)]]
    return "M" + " L".join(f"{x:.2f},{y:.2f}" for x, y in pts) + " Z"


def fetus_icon(fig, cx, cy, s):
    """Uterus with amniotic fluid and a curled fetus facing right, hand at the mouth, legs folded, umbilical cord.
    Local units: about 100 wide; s scales to layout pixels."""
    q = lambda x, y: (cx + x * s, cy + y * s)
    T = lambda *pts: " ".join(f"{q(*p)[0]:.1f},{q(*p)[1]:.1f}" for p in pts)
    def shp(d, fill, line=None, w=0):
        fig.add_shape(type="path", path=d, fillcolor=fill, line=dict(color=line or fill, width=w))
    def E(x, y, rx, ry, fill, ang=0, line=None, w=0):
        X, Y = q(x, y); shp(ell(X, Y, rx * s, ry * s, ang), fill, line, w)
    def capsule(x0, y0, x1, y1, r, fill):
        X0, Y0 = q(x0, y0); X1, Y1 = q(x1, y1)
        fig.add_shape(type="line", x0=X0, y0=Y0, x1=X1, y1=Y1, line=dict(color=fill, width=2 * r * s))
        for X, Y in ((X0, Y0), (X1, Y1)):
            fig.add_shape(type="circle", x0=X - r * s, x1=X + r * s, y0=Y - r * s, y1=Y + r * s, fillcolor=fill, line=dict(width=0))
    # uterus wall, cervix, fluid
    E(0, 0, 50, 44, WALL)
    shp(f"M{T((24, 34))} C{T((34, 38), (40, 46), (35, 51))} C{T((30, 55), (22, 50), (16, 40))} Z", WALL)
    E(0, -1, 44, 38, WALL_D)
    E(0, -1, 41, 35, FLUID)
    for x, y, rx, ry, a in ((-14, -26, 9, 3.5, 0), (26, -14, 3, 6, 20), (-30, 10, 3, 6, -10), (14, 26, 7, 2.5, 0)):
        E(x, y, rx, ry, FLUID_SPOT, a)
    # cord from the belly to the wall
    fig.add_shape(type="path", path=f"M{T((6, 14))} C{T((14, 6), (6, -2), (16, -8))} C{T((24, -14), (18, -22), (28, -28))}",
                  line=dict(color=CORD, width=2.6 * s), fillcolor="rgba(0,0,0,0)")
    # body: back on the left, belly facing the limbs
    shp(f"M{T((-18, -2))} C{T((-28, 12), (-14, 30), (6, 27))} C{T((20, 24), (24, 12), (14, 4))} C{T((6, -2), (-8, -6), (-18, -2))} Z", SKIN, LINE, 1)
    E(7, 17, 12, 7, SKIN_D, -28)                   # thigh
    capsule(16, 10, 21, 21, 3.2, SKIN_D)           # lower leg
    E(24, 22, 4, 2.4, SKIN_D, 10)                  # foot
    # head
    E(-11, -11, 16, 15, SKIN, 0, LINE, 1)
    E(-15, -18, 6, 3, SKIN_L, -20)                 # highlight
    E(-12, -7, 3, 4, SKIN_D)                       # ear
    fig.add_shape(type="circle", x0=q(-1.2, -12)[0], x1=q(1.2, -12)[0], y0=q(0, -14.2)[1], y1=q(0, -9.8)[1], fillcolor=LINE, line=dict(width=0))
    # arm, hand at the mouth
    capsule(-4, 6, 2, 10, 3, SKIN_D)
    capsule(2, 10, 6, -1, 2.8, SKIN_D)
    E(6.5, -2, 3.2, 3, SKIN_D)


def scan(cx, cy, s, col):
    """Cohort icon: uterus with a fetus, (col kept for the call signature)."""
    fetus_icon(fig, cx, cy + 4 * s, 0.85 * s)


def tree(cx, cy, col, prenatal):
    """Ontology: root, three levels, filled node = term with genes. With prenatal=True an extra purple branch whose
    hollow nodes carry no gene annotation."""
    nodes = {"r": (cx, cy), "a": (cx - 45, cy + 40), "b": (cx, cy + 40), "c": (cx + 45, cy + 40),
             "a1": (cx - 60, cy + 80), "a2": (cx - 30, cy + 80), "b1": (cx, cy + 80), "c1": (cx + 45, cy + 80)}
    edges = [("r", "a"), ("r", "b"), ("r", "c"), ("a", "a1"), ("a", "a2"), ("b", "b1"), ("c", "c1")]
    for u, v in edges:
        line(*nodes[u], *nodes[v], col, 2)
    for k, (x, y) in nodes.items():
        circle(x, y, 7, col, col, 0)
    if prenatal:
        px, py = cx + 85, cy + 40
        line(cx, cy, px, py, PRE, 2.5)
        kids = [(px - 15, cy + 80), (px + 15, cy + 80), (px + 15, cy + 115)]
        line(px, py, *kids[0], PRE, 2); line(px, py, *kids[1], PRE, 2); line(*kids[1], *kids[2], PRE, 2)
        circle(px, py, 7, PRE, PRE, 0); circle(*kids[0], 7, PRE, PRE, 0)
        for k in kids[1:]:
            circle(*k, 7, PRE, "white", 2.5)


def ranklist(x0, y0, highlight):
    for i, w in enumerate((70, 60, 52, 44, 36)):
        col = RANK if i == highlight else DIV
        fig.add_shape(type="rect", x0=x0, x1=x0 + w, y0=y0 + 13 * i, y1=y0 + 13 * i + 9, line=dict(width=0), fillcolor=col)


# column 1: cohorts
text(12, 22, "Published cohorts", 18, True, anchor="left")
scan(56, 90, 0.95, COH)
line(100, 52, 100, 140, COH, 3)
text(112, 66, "Fu et al. 2022, primary", FS, True, anchor="left")
text(112, 90, f"{per.fetuses['fu2022']} fetuses, {fu_genes} genes", anchor="left")
text(112, 114, f"{per.findings['fu2022']} ultrasound findings", col=TGREY, anchor="left")
line(20, 160, 300, 160, DIV, 1.5)
further = [("graham2026", "Graham et al. 2026"), ("smogavec2022", "Smogavec et al. 2022"), ("aogs2026", "Yao et al. 2026")]
for i, (k, name) in enumerate(further):
    y = 200 + 62 * i
    scan(52, y - 4, 0.62, COH)
    line(100, y - 24, 100, y + 24, COH, 3)
    text(112, y - 11, name, FS, True, anchor="left")
    text(112, y + 12, f"{per.fetuses[k]} fetuses, {per.findings[k]} findings", col=TGREY, anchor="left")
text(160, 400, "Monogenic diagnosis,", FS, True); text(160, 422, "prenatal findings only", FS, True)
text(160, 448, "Findings coded as HPO terms", col=TGREY)

# column 2: two releases, one query
chevron(312, 440)
text(338, 22, "Two HPO releases", 18, True, anchor="left")
tree(405, 70, R20, False); tree(535, 70, R26, True)
text(405, 222, "HPO 2020", FS, True, R20); text(405, 245, f"{R['genes_2020']:,} genes", col=TGREY)
text(550, 222, "HPO 2026", FS, True, R26); text(550, 245, f"{R['genes_2026']:,} genes", col=TGREY)
line(335, 275, 645, 275, DIV, 1.5)
circle(352, 305, 7, PRE, PRE, 0)
text(368, 305, f"Prenatal branch, {X['x2']['prenatal_terms']} terms", anchor="left")
circle(352, 335, 7, PRE, "white", 2.5)
text(368, 335, "term without gene annotation", col=TGREY, anchor="left")
line(335, 365, 645, 365, DIV, 1.5)
text(490, 395, "Identical query terms", FS, True)
text(490, 419, "post-2020 terms replaced by", col=TGREY); text(490, 441, "their nearest 2020 ancestor", col=TGREY)
text(490, 470, f"{R['common_genes']:,} genes common to both", FS, True)

# column 3: analyses
chevron(652, 440)
text(685, 22, "Analyses", 18, True, anchor="left")
ranklist(685, 62, 2)
text(772, 66, "Rank of the", FS, True, anchor="left"); text(772, 88, "causal gene", FS, True, anchor="left")
text(772, 110, "Resnik best-match", col=TGREY, anchor="left"); text(772, 132, "top 10, 2020 vs 2026", col=TGREY, anchor="left")
line(685, 160, 955, 160, DIV, 1.5)
circle(720, 205, 9, PRE, PRE, 0)
for dx in (-22, 22):
    line(720, 214, 720 + dx, 246, PRE, 2); circle(720 + dx, 252, 8, PRE, "white", 2.5)
text(772, 205, "Annotation gap", FS, True, anchor="left")
text(772, 229, "prenatal terms", col=TGREY, anchor="left"); text(772, 251, "without genes", col=TGREY, anchor="left")
line(685, 290, 955, 290, DIV, 1.5)
gx, gy = 720, 360
for dx, dy in ((-36, -38), (32, -40), (0, 48)):
    scan(gx + dx, gy + dy - 4, 0.34, COH)
    fig.add_annotation(x=gx + 0.35 * dx, y=gy + 0.35 * dy, ax=gx + 0.75 * dx, ay=gy + 0.75 * dy, xref="x", yref="y", axref="x",
                       ayref="y", text="", showarrow=True, arrowhead=2, arrowsize=1, arrowwidth=2, arrowcolor=SUPP)
circle(gx, gy, 13, SUPP, rgba(SUPP, 0.25), 2.5)
text(gx, gy, "<i>G</i>", 14, True, SUPP)
text(772, 330, "Annotation", FS, True, anchor="left"); text(772, 352, "supplementation", FS, True, anchor="left")
text(772, 376, "post hoc, leave one", col=TGREY, anchor="left"); text(772, 398, "fetus out", col=TGREY, anchor="left")
text(772, 424, f"{X['x3']['fetuses']} fetuses, {X['x3']['genes']} genes", anchor="left")

line(10, 500, 950, 500, DIV, 1.5)
text(480, 528, f"{n_all} fetuses from four published cohorts, {R['common_genes']:,} genes ranked for each fetus in two HPO releases",
     FS, True)
fig.update_layout(showlegend=False, plot_bgcolor="white", template="plotly_white", width=W_ + 4, height=H_ + 4,
                  margin=dict(l=2, r=2, t=2, b=2), font=dict(size=FS))
os.makedirs(FIG, exist_ok=True)
fig.write_html(os.path.join(FIG, "study_flow.html"), include_plotlyjs="cdn")
fig.write_image(os.path.join(FIG, "study_flow.pdf"))
fig.write_image(os.path.join(FIG, "study_flow.png"), scale=600 / 144)
print("study flow written")
