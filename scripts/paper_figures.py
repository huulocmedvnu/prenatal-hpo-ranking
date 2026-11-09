"""Figures 2 to 4 (Figure 1 is scripts/flowchart.py). One layout pixel prints at 0.5 pt, so 16 px
fonts print at 8 pt and 960 px at 169 mm. Colour = meaning, as in the study flow: HPO 2020 green, HPO 2026 blue,
supplemented annotations red (causal gene only) and orange (all genes of the primary cohort), prenatal branch
purple.
Figure 2: cumulative share of fetuses by rank of the causal gene, 2020 against 2026, primary and further cohorts.
Figure 3: terms of the prenatal branch of HPO 2026 with and without gene annotation, by sub-branch (sub-branches
with at least five terms) and, under Fetal anomaly, by organ sub-branch (at least four terms). All terms are listed
in outputs/prenatal_branch_terms.tsv.
Figure 4: rank with the 2026 annotations, with the annotations of all genes of the primary cohort supplemented and
with only the causal gene supplemented by the findings of other fetuses (post hoc, leave one fetus out).
Writes figures/paper/figure<k>.{html,pdf,png}."""
import os
import numpy as np, pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT, FIG = os.path.join(ROOT, "outputs"), os.path.join(ROOT, "figures", "paper")
C = pio.templates["plotly_white"].layout.colorway
FS = 16
R20, R26, SUPP, PRE, COMP, GREY = C[2], C[0], C[1], C[3], C[4], "#B0B7C3"


def save(fig, k, w, h, margin=None):
    os.makedirs(FIG, exist_ok=True)
    fig.update_layout(template="plotly_white", width=w, height=h, font=dict(size=FS), margin=margin or dict(l=80, r=20, t=50, b=60))
    fig.update_annotations(font_size=FS)
    fig.write_html(os.path.join(FIG, f"figure{k}.html"), include_plotlyjs="cdn")
    fig.write_image(os.path.join(FIG, f"figure{k}.pdf")); fig.write_image(os.path.join(FIG, f"figure{k}.png"), scale=600 / 144)


def ecdf(fig, r, col, ranks, name, color, show, dash=None):
    grid = np.unique(np.concatenate([[1], np.logspace(0, np.log10(max(ranks.max(), 2)), 300)]))
    fig.add_trace(go.Scatter(x=grid, y=[100 * (ranks <= g).mean() for g in grid], mode="lines", name=name, legendgroup=name,
                             showlegend=show, line=dict(color=color, width=2.5, shape="hv", dash=dash)), r, col)


def figure2():
    R = pd.read_csv(os.path.join(OUT, "ranks.tsv"), sep="\t")
    sets = [("Primary cohort", R[R.cohort == "fu2022"]), ("Further cohorts", R[R.cohort != "fu2022"])]
    fig = make_subplots(1, 2, shared_yaxes=True, horizontal_spacing=0.06,
                        subplot_titles=[f"{n} (n = {len(d)} fetuses)" for n, d in sets])
    for j, (n, d) in enumerate(sets, 1):
        ecdf(fig, 1, j, d.rank_2020, "HPO 2020", R20, j == 1)
        ecdf(fig, 1, j, d.rank_2026, "HPO 2026", R26, j == 1)
        fig.add_vline(x=10, line=dict(color="grey", dash="dot", width=1), row=1, col=j)
        fig.update_xaxes(type="log", title_text="Rank of the causal gene", row=1, col=j)
    fig.update_yaxes(title_text="Fetuses with the causal gene<br>at or above the rank (%)", range=[0, 100], row=1, col=1)
    fig.update_layout(legend=dict(orientation="h", x=0.5, xanchor="center", y=-0.25))
    save(fig, 2, 960, 470, dict(l=90, r=20, t=50, b=110))


def figure3():
    B = pd.read_csv(os.path.join(OUT, "prenatal_branch_terms.tsv"), sep="\t")
    B["no_gene"] = B.genes == 0
    top = B[B.sub_branch != "Fetal anomaly"].groupby("sub_branch").agg(n=("term", "size"), none=("no_gene", "sum"))
    top = top[top.n >= 5].sort_values("n")
    fet = B[B.fetal_sub_branch.notna() & (B.fetal_sub_branch != "")].groupby("fetal_sub_branch").agg(n=("term", "size"), none=("no_gene", "sum"))
    fet = fet[fet.n >= 4].sort_values("n")
    fig = make_subplots(1, 2, horizontal_spacing=0.42, column_widths=[0.5, 0.5],
                        subplot_titles=["Outside Fetal anomaly", "Fetal anomaly"])
    for j, T in ((1, top), (2, fet)):
        lab = [s.replace("Abnormal fetal ", "").replace("Abnormalities of ", "").replace("Abnormality of the ", "")
               .replace(" morphology", "") for s in T.index]
        lab = [s[0].upper() + s[1:] for s in lab]
        fig.add_trace(go.Bar(y=lab, x=T.n - T.none, orientation="h", name="With gene annotation", marker_color=PRE,
                             showlegend=j == 1, legendgroup="a"), 1, j)
        fig.add_trace(go.Bar(y=lab, x=T.none, orientation="h", name="Without gene annotation", marker_color="white",
                             marker_line=dict(color=PRE, width=2), showlegend=j == 1, legendgroup="b"), 1, j)
        fig.update_xaxes(title_text="HPO terms", row=1, col=j)
    fig.update_layout(barmode="stack", legend=dict(orientation="h", x=0.5, xanchor="center", y=-0.18))
    save(fig, 3, 960, 520, dict(l=230, r=20, t=50, b=110))


def figure4():
    """Fetuses whose causal gene is diagnosed in at least one other primary-cohort fetus, queried with their original
    2026 terms: 2026 annotations, annotations of all genes of the primary cohort supplemented (competitive), and only
    the causal gene supplemented (noncompetitive upper bound)."""
    V = pd.read_csv(os.path.join(OUT, "fetus_ranks.tsv"), sep="\t", dtype={"case": str})
    comp, upper = V[V.analysis == "supp_all_genes"], V[V.analysis == "supp_causal_only"]
    assert (comp["first"].values == upper["first"].values).all()
    fig = make_subplots(1, 1)
    ecdf(fig, 1, 1, upper["first"], "HPO 2026 annotations", R26, True)
    ecdf(fig, 1, 1, comp["second"], "All genes of the primary cohort supplemented (competitive)", COMP, True)
    ecdf(fig, 1, 1, upper["second"], "Causal gene only supplemented (noncompetitive upper bound)", SUPP, True, dash="dash")
    fig.add_vline(x=10, line=dict(color="grey", dash="dot", width=1))
    fig.update_xaxes(type="log", title_text="Rank of the causal gene")
    fig.update_yaxes(title_text="Fetuses with the causal gene<br>at or above the rank (%)", range=[0, 100])
    fig.update_layout(legend=dict(orientation="v", x=0, xanchor="left", y=-0.2, yanchor="top"))
    save(fig, 4, 760, 560, dict(l=90, r=20, t=30, b=190))


if __name__ == "__main__":
    figure2(); figure3(); figure4()
    print("paper figures written")
