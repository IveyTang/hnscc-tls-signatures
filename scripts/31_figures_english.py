# -*- coding: utf-8 -*-
"""选题四 / Step 31 —— 图 1–8 英文重绘 + DCA 补 clinical-only 曲线

两个目的：
  1) 原 fig_km_cts / fig_forest / fig_roc / fig_convergence / fig_microenv /
     fig_nomogram / fig_calibration / fig_dca 的标题与坐标轴为中文
     （"总生存率""阈概率""净获益"等），不能直接投英文刊 → 全部重绘为英文。
  2) 终审 S-3：06 的 DCA 只有 CTS+临床 / CTS 单独 / treat-all，缺 clinical-only，
     正文却写成"加评分带来增量获益"。本脚本补第四条曲线并量化差值。

图 9–11 与 Figure S1（脚本 21/24/26 产出）本就是英文，不动。

输出（写入 results/fig_en/，不改动旧 PNG）：
  fig_km_cts.png, fig_forest_multivariate.png, fig_roc_timedependent.png,
  fig_signature_convergence.png, fig_microenv_correlation.png,
  fig_nomogram.png, fig_calibration.png, fig_dca.png
  results/31_dca_summary.json
"""
import json
import os
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter

warnings.filterwarnings("ignore")

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7,
    "axes.linewidth": 0.7,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
})
RED, TEAL, BLUE, GREY = "#A32D2D", "#0F6E56", "#2E6DA4", "#8A8985"

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROC = os.path.join(ROOT, "data", "processed")
RES = os.path.join(ROOT, "results")
OUT = os.path.join(RES, "fig_en")
os.makedirs(OUT, exist_ok=True)

SETS = {
    "G1": "G1  Coppola 12-chemokine",
    "G2": "G2  11-chemokine (HCC)",
    "G3": "G3  HNSCC 13-chemokine",
    "G4": "G4  Maturation / germinal centre",
    "G5": "G5  Th1 / IFN-γ",
    "G6": "G6  HEV / lymphoid stroma",
    "G7": "G7  T cell-inflamed (control)",
}
PROG = {
    "检查点": "Checkpoint",
    "Treg": "Regulatory T cells",
    "Tfh": "T follicular helper",
    "CD8 T": "CD8 T cells",
    "IFN-γ 程序": "IFN-γ programme",
    "B 细胞": "B cells",
    "巨噬细胞": "Macrophages",
    "内皮/HEV": "Endothelium / HEV",
    "浆细胞": "Plasma cells",
    "生发中心 B": "Germinal-centre B cells",
    "CAF": "CAFs",
    "上皮/肿瘤": "Epithelium (control)",
}

df = pd.read_csv(os.path.join(PROC, "hnsc_tls_scores.csv"), index_col=0)


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("[ok] %s" % name)


# ==================================================== Figure 1  KM
fig, ax = plt.subplots(figsize=(4.0, 3.4))
for grp, clr in [("High", RED), ("Low", TEAL)]:
    sub = df[df["CTS_group_opt"] == grp]
    kmf = KaplanMeierFitter().fit(sub["os_time_days"] / 30.44, sub["os_event"],
                                  label="CTS-high (n = %d)" % len(sub))
    kmf.plot_survival_function(ax=ax, color=clr, ci_show=True,
                               ci_alpha=0.16, linewidth=1.3)
ax.set_xlabel("Time since diagnosis (months)")
ax.set_ylabel("Overall survival")
ax.set_xlim(0, 180)
ax.set_ylim(0, 1.02)
ax.text(0.04, 0.06, "Log-rank $\\chi^2$ = 5.72, p = 0.0168",
        transform=ax.transAxes, fontsize=7)
ax.legend(loc="upper right", frameon=False)
ax.grid(alpha=0.2, linewidth=0.4)
save(fig, "fig_km_cts.png")

# ==================================================== Figure 2  forest
mv = pd.read_csv(os.path.join(RES, "04_cox_multivariate.csv"))
labels = {
    "CTS": "CTS (per unit)",
    "age": "Age (per year)",
    "male": "Male sex",
    "stage3_II": "Stage II vs I",
    "stage3_III": "Stage III vs I",
    "stage3_IVA": "Stage IVA vs I",
    "stage3_IVB": "Stage IVB vs I",
    "stage3_IVC": "Stage IVC vs I",
}
order = ["CTS", "age", "male", "stage3_II", "stage3_III", "stage3_IVA",
         "stage3_IVB", "stage3_IVC"]
keep = mv[mv["var"].isin(order)].copy()
keep["lab"] = keep["var"].map(labels)
keep = keep.set_index("var").loc[[v for v in order if v in set(keep["var"])]].reset_index()

fig, ax = plt.subplots(figsize=(4.4, 3.2))
y = np.arange(len(keep))[::-1]
ax.errorbar(keep["HR"], y,
            xerr=[keep["HR"] - keep["lower"], keep["upper"] - keep["HR"]],
            fmt="o", color=RED, ecolor=GREY, capsize=2.5, markersize=4,
            linewidth=0.9, elinewidth=0.9)
ax.axvline(1, color=GREY, linestyle="--", linewidth=0.7)
ax.set_yticks(y)
ax.set_yticklabels(keep["lab"])
ax.set_xscale("log")
ticks = [0.3, 0.5, 1.0, 2.0, 4.0, 8.0]
ax.set_xticks(ticks)
ax.set_xticklabels(["0.3", "0.5", "1.0", "2.0", "4.0", "8.0"])
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
ax.set_xlabel("Hazard ratio for overall survival (95% CI, log scale)")
ax.grid(axis="x", alpha=0.2, linewidth=0.4)
for yy, (_, r) in zip(y, keep.iterrows()):
    ax.text(9.5, yy, "%.3f" % r["p"], va="center", fontsize=6.5, color="#333333")
ax.text(9.5, len(keep) - 0.35, "p", va="center", fontsize=6.5, style="italic")
ax.set_xlim(0.25, 20)
ax.set_ylim(-0.7, len(keep) + 0.2)
save(fig, "fig_forest_multivariate.png")

# ==================================================== Figure 3  time-dependent AUC
auc = pd.read_csv(os.path.join(RES, "04_time_dependent_auc.csv"))
name = {"CTS": "Consensus TLS score", "G1_Coppola_12chemokine": "Coppola 12-chemokine",
        "G3_HNSCC_13chemokine": "HNSCC 13-chemokine",
        "G7_Tcell_inflamed_cytolytic": "T cell-inflamed (control)"}
style = {"CTS": (RED, "o", "-"), "G1_Coppola_12chemokine": (BLUE, "s", "--"),
         "G3_HNSCC_13chemokine": (TEAL, "^", "-."),
         "G7_Tcell_inflamed_cytolytic": (GREY, "D", "-")}
fig, ax = plt.subplots(figsize=(4.0, 3.2))
for s, (clr, mk, ls) in style.items():
    sub = auc[auc["score"] == s].sort_values("time_months")
    if sub.empty:
        continue
    ax.plot(sub["time_months"], sub["AUC"], marker=mk, markersize=4, linewidth=1.2,
            color=clr, linestyle=ls, label=name[s])
ax.axhline(0.5, color=GREY, linestyle=":", linewidth=0.7)
ax.set_xticks([12, 36, 60])
ax.set_xlabel("Time after diagnosis (months)")
ax.set_ylabel("Time-dependent AUC")
ax.set_ylim(0.45, 0.62)
ax.legend(loc="upper left", frameon=False, fontsize=6.5)
ax.grid(alpha=0.2, linewidth=0.4)
save(fig, "fig_roc_timedependent.png")

# ==================================================== Figure 4  heatmap
cols = [c for c in df.columns if c.startswith("G") and "__pca1" not in c]
cm = df[cols].corr(method="spearman")
fig, ax = plt.subplots(figsize=(4.6, 3.9))
im = ax.imshow(cm.values, cmap="RdYlBu_r", vmin=0, vmax=1)
ax.set_xticks(range(len(cm)))
ax.set_xticklabels([c.split("_")[0] for c in cm.columns])
ax.set_yticks(range(len(cm)))
ax.set_yticklabels([SETS[c.split("_")[0]] for c in cm.index])
for i in range(len(cm)):
    for j in range(len(cm)):
        v = cm.values[i, j]
        ax.text(j, i, "%.2f" % v, ha="center", va="center", fontsize=6.5,
                color="white" if v > 0.75 or v < 0.35 else "black")
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)
cb = fig.colorbar(im, ax=ax, shrink=0.72, pad=0.02)
cb.set_label("Spearman ρ", fontsize=7)
cb.ax.tick_params(labelsize=6.5)
save(fig, "fig_signature_convergence.png")

# ==================================================== Figure 5  marker correlations
corr = pd.read_csv(os.path.join(RES, "06_cts_marker_correlation.csv"))
corr["lab"] = corr["marker"].map(lambda x: PROG.get(x, x))
corr = corr.sort_values("spearman_rho")
fig, ax = plt.subplots(figsize=(4.3, 3.4))
ax.barh(corr["lab"], corr["spearman_rho"],
        color=[RED if v > 0 else TEAL for v in corr["spearman_rho"]], height=0.62)
ax.axvline(0, color="#4A4A48", linewidth=0.7)
ax.set_xlabel("Spearman ρ with the consensus TLS score")
ax.set_xlim(-0.25, 1.0)
ax.grid(axis="x", alpha=0.2, linewidth=0.4)
ax.tick_params(axis="y", length=0)
save(fig, "fig_microenv_correlation.png")

# ==================================================== Figure 6  nomogram (points)
pts = pd.read_csv(os.path.join(RES, "06_nomogram_points.csv"))
pretty = {"CTS": "Consensus TLS score", "age": "Age", "male": "Male sex",
          "stage_iv": "Stage IV"}
pts["lab"] = pts["var"].map(lambda x: pretty.get(x, x))
pts = pts.sort_values("points_per_unit")

fig, ax = plt.subplots(figsize=(5.0, 2.6))
ax.axis("off")
n = len(pts)
ax.set_xlim(0.0, 1.0)
ax.set_ylim(-0.55, n + 0.95)
# 顶端 0–100 分值轴
ax.plot([0.34, 0.96], [n + 0.35] * 2, color="#2C2C2A", linewidth=0.7)
for i in range(0, 101, 10):
    x = 0.34 + 0.62 * i / 100
    ax.plot([x, x], [n + 0.35, n + 0.22], color="#2C2C2A", linewidth=0.5)
    ax.text(x, n + 0.48, str(i), ha="center", fontsize=6.5)
ax.text(0.30, n + 0.35, "Points", ha="right", va="center", fontsize=7.5)
for yy, (_, r) in zip(np.arange(n)[::-1], pts.iterrows()):
    frac = r["points_per_unit"] / 100.0
    ax.plot([0.34, 0.34 + 0.62 * frac], [yy, yy], color=BLUE, linewidth=1.4)
    ax.text(0.325, yy, r["lab"], ha="right", va="center", fontsize=7.5)
    ax.text(0.34 + 0.62 * frac + 0.012, yy, "%.1f" % r["points_per_unit"],
            va="center", fontsize=6.5, color="#4A4A48")
save(fig, "fig_nomogram.png")

# ==================================================== Figure 7  calibration
cal = pd.read_csv(os.path.join(RES, "06_calibration.csv"))
fig, ax = plt.subplots(figsize=(3.4, 3.2))
for t, clr, mk in zip((12, 36, 60), (BLUE, RED, TEAL), ("o", "s", "^")):
    sub = cal[cal["time"] == t]
    if sub.empty:
        continue
    ax.plot(1 - sub["pred_risk_mean"], sub["observed_surv"], marker=mk, markersize=4,
            linewidth=1.1, color=clr, label="%d months" % t)
ax.plot([0.2, 1.0], [0.2, 1.0], linestyle="--", color=GREY, linewidth=0.7)
ax.set_xlabel("Predicted overall survival")
ax.set_ylabel("Observed overall survival (Kaplan–Meier)")
ax.set_xlim(0.15, 1.0)
ax.set_ylim(0.15, 1.05)
ax.legend(loc="upper left", frameon=False)
ax.grid(alpha=0.2, linewidth=0.4)
save(fig, "fig_calibration.png")

# ==================================================== Figure 8  DCA（含 clinical-only）
def prep(d0):
    d = d0.copy()
    d = d[d["os_time_days"].notna() & d["os_event"].notna()]
    d["os_event"] = d["os_event"].astype(int)
    d["age"] = pd.to_numeric(d["age_at_index"], errors="coerce")
    d["male"] = (d["gender"] == "male").astype(int)
    d["stage_iv"] = d["stage"].astype(str).str.startswith("IV").astype(int)
    return d


d = prep(pd.read_csv(os.path.join(PROC, "hnsc_tls_scores.csv"), index_col=0))
d = d[["CTS", "age", "male", "stage_iv", "os_time_days", "os_event"]].dropna().copy()


def fit(cols):
    c = CoxPHFitter(penalizer=0.05)
    c.fit(d[cols + ["os_time_days", "os_event"]], "os_time_days", "os_event")
    return c


def risk36(cph):
    bh = cph.baseline_cumulative_hazard_.copy()
    bh.columns = ["H0"]
    ev = np.sort(d.loc[d["os_event"] == 1, "os_time_days"].unique())
    h0 = np.interp(ev, bh.index.values, bh["H0"].values)
    lp = cph.predict_log_partial_hazard(d).values
    s0 = np.exp(-h0)
    idx = min(np.searchsorted(ev, 36 * 30.44), len(s0) - 1)
    return 1.0 - s0[idx] ** np.exp(lp)


m_comb = fit(["CTS", "age", "male", "stage_iv"])
m_clin = fit(["age", "male", "stage_iv"])
r_comb, r_clin = risk36(m_comb), risk36(m_clin)
print("[model] combined CTS p=%.3f  stageIV HR=%.3f"
      % (m_comb.summary.loc["CTS", "p"], np.exp(m_comb.summary.loc["stage_iv", "coef"])))


def dca(pred, thr):
    n = len(d)
    ev = d["os_event"].values
    nb = []
    for pt in thr:
        yh = pred >= pt
        tp = int(((yh == 1) & (ev == 1)).sum())
        fp = int(((yh == 1) & (ev == 0)).sum())
        nb.append(tp / n - fp / n * (pt / (1 - pt)))
    e = float(ev.mean())
    return np.array(nb), e - (1 - e) * (thr / (1 - thr))


thr = np.arange(0.05, 0.81, 0.01)
nb_comb, nb_all = dca(r_comb, thr)
nb_clin, _ = dca(r_clin, thr)
nb_cts, _ = dca(d["CTS"].rank(pct=True).values, thr)

pd.DataFrame({"threshold": thr, "nb_combined": nb_comb, "nb_clinical_only": nb_clin,
              "nb_cts_only": nb_cts, "nb_treat_all": nb_all}).to_csv(
    os.path.join(RES, "31_dca_curves.csv"), index=False)

fig, ax = plt.subplots(figsize=(4.4, 3.4))
ax.plot(thr, nb_comb, color=RED, linewidth=1.5, label="CTS + clinical")
ax.plot(thr, nb_clin, color=TEAL, linewidth=1.4, linestyle="--", label="Clinical only")
ax.plot(thr, nb_cts, color=BLUE, linewidth=1.2, linestyle="-.", label="CTS alone")
ax.plot(thr, nb_all, color=GREY, linewidth=1.0, linestyle=":", label="Treat all")
ax.axhline(0, color="#4A4A48", linewidth=0.7)
ax.set_xlim(0.05, 0.80)
ax.set_ylim(-0.05, max(0.35, float(np.nanmax(nb_comb)) * 1.12))
ax.set_xlabel("Threshold probability")
ax.set_ylabel("Net benefit")
ax.legend(loc="upper right", frameon=False)
ax.grid(alpha=0.2, linewidth=0.4)
save(fig, "fig_dca.png")

band = (thr >= 0.10) & (thr <= 0.65)
summ = {
    "n": int(len(d)),
    "band": [0.10, 0.65],
    "mean_nb_combined": round(float(nb_comb[band].mean()), 4),
    "mean_nb_clinical_only": round(float(nb_clin[band].mean()), 4),
    "mean_nb_cts_only": round(float(nb_cts[band].mean()), 4),
    "mean_delta_combined_minus_clinical": round(float((nb_comb - nb_clin)[band].mean()), 4),
    "delta_range": [round(float((nb_comb - nb_clin)[band].min()), 4),
                    round(float((nb_comb - nb_clin)[band].max()), 4)],
    "frac_thresholds_combined_above_clinical": round(
        float(((nb_comb - nb_clin)[band] > 0).mean()), 3),
    "cox_combined_CTS_p": round(float(m_comb.summary.loc["CTS", "p"]), 4),
    "cox_combined_stageIV_HR": round(float(np.exp(m_comb.summary.loc["stage_iv", "coef"])), 3),
}
json.dump(summ, open(os.path.join(RES, "31_dca_summary.json"), "w"), indent=1)
print(json.dumps(summ, indent=1))
print("\n[ok] 英文图件写入 %s" % OUT)
