# -*- coding: utf-8 -*-
"""
22_ext_figure_panel.py
----------------------
Compose the four external-validation panels of **GSE65858** into a single
main-text figure (Figure 10 of the manuscript), plus a standalone summary
table (Table 6).

Panels
~~~~~~
  A  gene set score convergence (Spearman rho matrix)
  B  correlation of the consensus TLS score with microenvironmental programmes
  C  consensus TLS score by HPV16 RNA status
  D  Kaplan-Meier: consensus score (null) vs the pre-specified positive control

Everything is read back from `results/21_ext_*`, so this script does not
re-run the statistics — it only draws them.
"""
from __future__ import annotations

import io
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RES = os.path.join(ROOT, "results")
DISCOVERY_CUT = 0.2489
CTRL = "G7_Tcell_inflamed_cytolytic"


def km(t, e):
    ts = np.unique(t[e == 1])
    surv, pts = 1.0, []
    for tt in ts:
        n = (t >= tt).sum()
        d = ((t == tt) & (e == 1)).sum()
        if n > 0:
            surv *= (1 - d / n)
        pts.append((tt, surv))
    return (np.array([0.0] + [a for a, _ in pts]),
            np.array([1.0] + [b for _, b in pts]))


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams["font.sans-serif"] = ["DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False

    d = pd.read_csv(os.path.join(RES, "21_ext_cts_scores.csv"))
    summ = json.load(io.open(os.path.join(RES, "21_ext_summary.json"),
                            encoding="utf-8"))

    sets_all = [c for c in d.columns
                if c.startswith(("G1_", "G2_", "G3_", "G4_", "G5_", "G6_", "G7_"))]
    nice = {s: s.split("_", 1)[1] for s in sets_all}

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 9.2), dpi=300)

    # ---- A: convergence heatmap ----------------------------------------
    ax = axes[0][0]
    M = np.ones((len(sets_all), len(sets_all)))
    for i, a in enumerate(sets_all):
        for j, b in enumerate(sets_all):
            M[i, j] = stats.spearmanr(d[a], d[b]).statistic
    im = ax.imshow(M, cmap="RdYlBu_r", vmin=0.4, vmax=1.0)
    lab = [nice[s] for s in sets_all]
    ax.set_xticks(range(len(lab)))
    ax.set_yticks(range(len(lab)))
    ax.set_xticklabels(lab, rotation=40, ha="right", fontsize=6.4)
    ax.set_yticklabels(lab, fontsize=6.4)
    for i in range(len(lab)):
        for j in range(len(lab)):
            ax.text(j, i, "%.2f" % M[i, j], ha="center", va="center",
                    fontsize=5.4, color="white" if M[i, j] > 0.88 else "black")
    fig.colorbar(im, ax=ax, shrink=0.82, label="Spearman rho")
    ax.set_title("A  Gene set score convergence", fontsize=9.5, loc="left")

    # ---- B: programme correlation with CTS -----------------------------
    ax = axes[0][1]
    pc = pd.DataFrame(summ["programme_correlation"])
    pc = pc.sort_values("spearman_rho_with_CTS")
    cols = ["#A32D2D" if v > 0 else "#0F6E56"
            for v in pc["spearman_rho_with_CTS"]]
    ax.barh(pc["programme"], pc["spearman_rho_with_CTS"], color=cols, alpha=0.85)
    for y, v in enumerate(pc["spearman_rho_with_CTS"]):
        ax.text(v + (0.02 if v >= 0 else -0.02), y, "%.2f" % v,
                va="center", ha="left" if v >= 0 else "right", fontsize=6)
    ax.axvline(0, color="black", lw=0.7)
    ax.set_xlim(-0.30, 1.0)
    ax.set_xlabel("Spearman rho with consensus TLS score", fontsize=8)
    ax.tick_params(labelsize=6.8)
    ax.set_title("B  What the consensus score tracks", fontsize=9.5, loc="left")

    # ---- C: CTS by HPV16 RNA status ------------------------------------
    ax = axes[1][0]
    order = [k for k in ["DNA-", "DNA+RNA-", "DNA+RNA+"] if k in set(d["hpv_rna"])]
    data = [d.loc[d["hpv_rna"] == k, "CTS"].dropna().values for k in order]
    bp = ax.boxplot(data, patch_artist=True, showfliers=False, widths=0.55)
    ax.set_xticks(range(1, len(order) + 1))
    ax.set_xticklabels(order, fontsize=8)
    for patch, col in zip(bp["boxes"], ["#4A6FA5", "#C58C3A", "#A32D2D"]):
        patch.set_facecolor(col)
        patch.set_alpha(0.75)
    rng = np.random.default_rng(0)
    for i, y in enumerate(data):
        ax.scatter(rng.normal(i + 1, 0.06, len(y)), y, s=5,
                   color="#333333", alpha=0.35, zorder=3)
    ax.axhline(DISCOVERY_CUT, color="grey", ls="--", lw=1)
    ax.text(0.015, DISCOVERY_CUT, " discovery cut-off",
            transform=ax.get_yaxis_transform(), fontsize=6.5, va="bottom",
            color="grey")
    ax.set_ylabel("Consensus TLS score", fontsize=8)
    ax.set_xlabel("HPV16 status", fontsize=8)
    pv = summ["hpv"].get("mannwhitney_DNAposRNApos_vs_DNAneg", {}).get("p")
    ax.set_title("C  CTS by HPV16 status (DNA+RNA+ vs DNA-  p = %s)"
                 % ("%.4g" % pv if pv else "n.a."), fontsize=9.5, loc="left")

    # ---- D: KM, CTS vs positive control --------------------------------
    ax = axes[1][1]
    sub = d.dropna(subset=["os_days", "os_event", "CTS"])
    sub = sub[sub["os_days"] > 0]
    t = sub["os_days"].to_numpy(float)
    e = sub["os_event"].to_numpy(int)
    ctrl_row = next((r for r in summ["per_signature_survival"]
                     if r["score"] == CTRL), None)
    cut_c = ctrl_row["optimal_cut"] if ctrl_row else None

    def logrank_p(vec, cutv):
        m = vec > cutv
        t1, e1, t2, e2 = t[m], e[m], t[~m], e[~m]
        times = np.unique(np.concatenate([t1, t2]))
        O1 = E1 = V = 0.0
        for tt in times:
            n1, n2 = (t1 >= tt).sum(), (t2 >= tt).sum()
            n = n1 + n2
            d1 = ((t1 == tt) & (e1 == 1)).sum()
            d2 = ((t2 == tt) & (e2 == 1)).sum()
            dd = d1 + d2
            if n <= 1 or dd == 0:
                continue
            O1 += d1
            E1 += dd * n1 / n
            V += dd * (n1 / n) * (1 - n1 / n) * (n - dd) / (n - 1)
        if V <= 0:
            return 1.0, m.sum()
        return float(stats.chi2.sf((O1 - E1) ** 2 / V, 1)), int(m.sum())

    vecs = [("Consensus TLS score", sub["CTS"].to_numpy(float), DISCOVERY_CUT),
            ("T cell-inflamed (control)",
             sub[CTRL].to_numpy(float), cut_c)]
    styles = [("-", "#1f3b73"), ("--", "#1f3b73"),
              ("-", "#A32D2D"), ("--", "#A32D2D")]
    k = 0
    for title, vec, cutv in vecs:
        m = vec > cutv
        pv2, nh = logrank_p(vec, cutv)
        for mask, lab, ls, col in [(m, "high", styles[k][0], styles[k][1]),
                                   (~m, "low", styles[k + 1][0], styles[k + 1][1])]:
            xs, ys = km(t[mask], e[mask])
            ax.step(xs, ys, where="post", ls=ls, color=col, lw=1.6,
                    label="%s, %s (n=%d)" % (title, lab, mask.sum()))
        ax.text(0.98, 0.94 - 0.08 * (k // 2),
                "%s: p = %.3f" % (title.split(" (")[0], pv2),
                transform=ax.transAxes, ha="right", fontsize=7, color=styles[k][1])
        k += 2
    ax.set_xlabel("Overall survival (days)", fontsize=8)
    ax.set_ylabel("Survival probability", fontsize=8)
    ax.set_ylim(0, 1.03)
    ax.tick_params(labelsize=7)
    ax.legend(frameon=False, fontsize=5.8, loc="lower left")
    ax.set_title("D  Survival in the validation cohort", fontsize=9.5, loc="left")

    fig.suptitle("External validation in GSE65858 (n = %d, %d deaths; "
                 "EBI ArrayExpress mirror E-GEOD-65858)"
                 % (len(sub), int(e.sum())), fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    out = os.path.join(RES, "fig_ext_validation.png")
    fig.savefig(out)
    plt.close(fig)
    print("[fig] written", out)

    # ---- Table 6 --------------------------------------------------------
    rows = []
    def a(key, lab, fmt="%s"):
        rows.append([lab, fmt % summ[key]])
    rows = [
        ["Cohort", "GSE65858 (E-GEOD-65858)"],
        ["Platform", "Illumina HumanHT-12 v4"],
        ["Patients analysed", str(summ["n_samples"])],
        ["Deaths", str(summ["n_deaths"])],
        ["Genes scored", str(summ["n_genes_scored"])],
        ["HPV16 RNA annotation available", str(summ["hpv"]["n_annotated"])],
        ["Chemokine signature pairs, rho (min-max)",
         "%.3f-%.3f" % (summ["convergence"]["chemokine_pairs_min_rho"],
                        summ["convergence"]["chemokine_pairs_max_rho"])],
        ["CTS vs generic T cell-inflamed, rho",
         "%.3f-%.3f" % (
             min(x["rho"] for x in summ["convergence"]["chemokine_vs_Tcell_rho"]
                 if x["set"].startswith(("G1_", "G2_", "G3_"))),
             max(x["rho"] for x in summ["convergence"]["chemokine_vs_Tcell_rho"]
                 if x["set"].startswith(("G1_", "G2_", "G3_"))))],
        ["CTS vs checkpoint programme, rho",
         "%.3f" % next(p["spearman_rho_with_CTS"] for p in summ["programme_correlation"]
                       if p["programme"] == "Checkpoint")],
        ["CTS vs GC B cell programme, rho",
         "%.3f" % next(p["spearman_rho_with_CTS"] for p in summ["programme_correlation"]
                       if p["programme"] == "GC B cells")],
        ["CTS vs plasma cell programme, rho",
         "%.3f" % next(p["spearman_rho_with_CTS"] for p in summ["programme_correlation"]
                       if p["programme"] == "Plasma cells")],
        ["CTS, univariable HR per SD",
         "%.3f (p = %.3f)" % (next(r["hr"] for r in summ["survival"]
                                   if r["analysis"].startswith("continuous CTS (univariable")),
                              next(r["p"] for r in summ["survival"]
                                   if r["analysis"].startswith("continuous CTS (univariable")))],
        ["CTS, adjusted HR per SD (stage + age)",
         "%.3f (p = %.3f)" % (next(r["hr"] for r in summ["survival"]
                                   if r["analysis"].startswith("continuous CTS (adjusted")),
                              next(r["p"] for r in summ["survival"]
                                   if r["analysis"].startswith("continuous CTS (adjusted")))],
        ["Positive control: optimal cut-off p",
         "%.4f" % next(r["p_at_cut"] for r in summ["per_signature_survival"]
                       if r["score"] == CTRL and r.get("p_at_cut") is not None)],
        ["Median CTS, HPV16 DNA+RNA+ vs DNA-",
         "%.3f vs %.3f (p = %.4g)" % (
             summ["hpv"]["groups"]["DNA+RNA+"]["median_cts"],
             summ["hpv"]["groups"]["DNA-"]["median_cts"],
             summ["hpv"]["mannwhitney_DNAposRNApos_vs_DNAneg"]["p"])],
        ["Highest-CTS expression cluster",
         "%s (median %.3f; Kruskal-Wallis p = %.3g)" % (
             summ["cluster_kruskal"]["highest_median_cluster"],
             next(c["median"] for c in summ["cluster_cts"]
                  if c["cluster"] == summ["cluster_kruskal"]["highest_median_cluster"]),
             summ["cluster_kruskal"]["p"])],
    ]
    tbl = pd.DataFrame(rows, columns=["Item", "Value"])
    p6 = os.path.join(RES, "table6_external_validation.csv")
    tbl.to_csv(p6, index=False, encoding="utf-8-sig")
    print("[tbl] written", p6)
    for r in rows:
        print("   %-46s %s" % (r[0], r[1]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
