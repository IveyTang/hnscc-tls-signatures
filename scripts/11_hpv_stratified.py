# -*- coding: utf-8 -*-
"""选题四 / Step 11 —— HPV 分层与混杂分析

目的：检验 CTS 的预后价值是否被 HPV 状态混杂。

分析内容：
  A. CTS 在不同 HPV 层内的分布（预测 HPV+ / HPV-）
  B. 分层 KM：HPV+ 与 HPV- 内部各自做 CTS 高低分组生存比较
  C. 分层 Cox：各层内 CTS（连续）的单因素 HR
  D. 全队列多因素 Cox：CTS + 预测 HPV + 年龄 + 性别 + 分期；比较纳入 HPV 前后 CTS 的 HR
  E. 交互检验：CTS × HPV（似然比检验）
  F. 敏感性分析：仅在直接标注子集（ISH/p16/报告，n=84）中做 E
  G. 亚部位分层：仅口咽（HPV 相关）与仅非口咽

输出：
  results/11_hpv_cts_summary.json
  results/11_hpv_stratified_cox.csv
  results/11_hpv_interaction.csv
  results/fig_km_cts_by_hpv.png
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
from lifelines.statistics import logrank_test

warnings.filterwarnings("ignore")
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROC = os.path.join(ROOT, "data", "processed")
RES = os.path.join(ROOT, "results")
CLR_HI, CLR_LO = "#A32D2D", "#0F6E56"

d = pd.read_csv(os.path.join(PROC, "hnsc_with_hpv.csv"))
d["time_m"] = d["os_time_days"] / 30.44
d["age"] = d["age_at_index"].astype(float)
d["male"] = (d["sex_at_birth"] == "male").astype(int)
d["stage4"] = d["stage"].astype(str).isin(["IVA", "IVB", "IVC"]).astype(int)
d["hpv_pred_bin"] = (d["hpv_pred"] == "Positive").astype(int)
d["hpv_final_bin"] = (d["hpv_final"] == "Positive").astype(int)
d["hpv_label_bin"] = np.where(d["hpv_label"] == "Positive", 1,
                              np.where(d["hpv_label"] == "Negative", 0, np.nan))
d["cts_hi"] = (d["CTS"] >= 0.249).astype(int)

out = {}

# ---------------------------------------------------------------- A. 分布
g_pos = d.loc[d["hpv_pred_bin"] == 1, "CTS"]
g_neg = d.loc[d["hpv_pred_bin"] == 0, "CTS"]
from scipy.stats import mannwhitneyu
u, p = mannwhitneyu(g_pos, g_neg, alternative="two-sided")
out["A_distribution"] = {
    "n_hpv_pred_pos": int(len(g_pos)), "n_hpv_pred_neg": int(len(g_neg)),
    "cts_median_pos": round(float(g_pos.median()), 4),
    "cts_median_neg": round(float(g_neg.median()), 4),
    "mwu_p": float(p),
    "cts_high_fraction_pos": round(float((g_pos >= 0.249).mean()), 4),
    "cts_high_fraction_neg": round(float((g_neg >= 0.249).mean()), 4),
}

# ---------------------------------------------------------------- B. 分层 KM
km_rows = []
strata_def = {
    "HPV-positive (predicted)": d["hpv_pred_bin"] == 1,
    "HPV-negative (predicted)": d["hpv_pred_bin"] == 0,
    "Oropharynx": d["subsite"].astype(str) == "口咽",
    "Non-oropharynx": d["subsite"].astype(str) != "口咽",
}
for name, mask in strata_def.items():
    sub = d[mask]
    if len(sub) < 20:
        continue
    hi = sub["cts_hi"] == 1
    lo = sub["cts_hi"] == 0
    if hi.sum() < 5 or lo.sum() < 5:
        km_rows.append({"stratum": name, "n": int(len(sub)),
                        "n_high": int(hi.sum()), "n_low": int(lo.sum()),
                        "median_os_high": None, "median_os_low": None,
                        "logrank_p": None, "note": "组内样本不足以分层"})
        continue
    lr = logrank_test(sub.loc[hi, "time_m"], sub.loc[lo, "time_m"],
                      sub.loc[hi, "os_event"], sub.loc[lo, "os_event"])
    kmf_h = KaplanMeierFitter().fit(sub.loc[hi, "time_m"], sub.loc[hi, "os_event"])
    kmf_l = KaplanMeierFitter().fit(sub.loc[lo, "time_m"], sub.loc[lo, "os_event"])
    km_rows.append({
        "stratum": name, "n": int(len(sub)),
        "n_high": int(hi.sum()), "n_low": int(lo.sum()),
        "n_events": int(sub["os_event"].sum()),
        "median_os_high": round(float(kmf_h.median_survival_time_), 1),
        "median_os_low": round(float(kmf_l.median_survival_time_), 1),
        "logrank_p": float(lr.p_value),
    })
km_df = pd.DataFrame(km_rows)
km_df.to_csv(os.path.join(RES, "11_hpv_stratified_km.csv"), index=False,
             encoding="utf-8-sig")
out["B_stratified_km"] = km_rows

# ---------------------------------------------------------------- C/D. Cox
def fit_cox(df, cols, penalizer=0.01):
    c = CoxPHFitter(penalizer=penalizer)
    c.fit(df[cols + ["time_m", "os_event"]], "time_m", "os_event")
    return c

cox_rows = []
# 各层内 CTS 单因素（连续）
for name, mask in strata_def.items():
    sub = d[mask].dropna(subset=["CTS", "time_m", "os_event"])
    if sub["os_event"].sum() < 10:
        continue
    c = fit_cox(sub, ["CTS"])
    s = c.summary.loc["CTS"]
    cox_rows.append({
        "model": f"CTS (continuous), {name}", "n": int(len(sub)),
        "events": int(sub["os_event"].sum()),
        "HR": round(float(s["exp(coef)"]), 3),
        "CI_low": round(float(s["exp(coef) lower 95%"]), 3),
        "CI_high": round(float(s["exp(coef) upper 95%"]), 3),
        "p": float(s["p"]),
    })

# 全队列：纳入 HPV 前后
base = d.dropna(subset=["CTS", "age", "time_m", "os_event", "stage4", "male"])
c1 = fit_cox(base, ["CTS", "age", "male", "stage4"])
c2 = fit_cox(base, ["CTS", "age", "male", "stage4", "hpv_pred_bin"])
for tag, c in [("without HPV", c1), ("with HPV", c2)]:
    s = c.summary.loc["CTS"]
    cox_rows.append({
        "model": f"Multivariable CTS ({tag})", "n": int(len(base)),
        "events": int(base["os_event"].sum()),
        "HR": round(float(s["exp(coef)"]), 3),
        "CI_low": round(float(s["exp(coef) lower 95%"]), 3),
        "CI_high": round(float(s["exp(coef) upper 95%"]), 3),
        "p": float(s["p"]),
    })
s = c2.summary.loc["hpv_pred_bin"]
cox_rows.append({
    "model": "Multivariable HPV-positive (predicted)", "n": int(len(base)),
    "events": int(base["os_event"].sum()),
    "HR": round(float(s["exp(coef)"]), 3),
    "CI_low": round(float(s["exp(coef) lower 95%"]), 3),
    "CI_high": round(float(s["exp(coef) upper 95%"]), 3),
    "p": float(s["p"]),
})
cox_df = pd.DataFrame(cox_rows)
cox_df.to_csv(os.path.join(RES, "11_hpv_stratified_cox.csv"), index=False,
              encoding="utf-8-sig")
out["CD_cox"] = cox_rows

# ---------------------------------------------------------------- E. 交互
def interaction(df, tag):
    df = df.dropna(subset=["CTS", "age", "male", "stage4", "time_m", "os_event",
                           "hpv_pred_bin"]).copy()
    df["CTS_x_HPV"] = df["CTS"] * df["hpv_pred_bin"]
    c_main = CoxPHFitter(penalizer=0.01).fit(
        df[["CTS", "hpv_pred_bin", "age", "male", "stage4", "time_m", "os_event"]],
        "time_m", "os_event")
    c_int = CoxPHFitter(penalizer=0.01).fit(
        df[["CTS", "hpv_pred_bin", "CTS_x_HPV", "age", "male", "stage4",
            "time_m", "os_event"]], "time_m", "os_event")
    from lifelines.statistics import proportional_hazard_test  # noqa
    lr = 2 * (c_int.log_likelihood_ - c_main.log_likelihood_)
    from scipy.stats import chi2
    p_int = float(chi2.sf(lr, 1))
    s = c_int.summary.loc["CTS_x_HPV"]
    return {
        "analysis": tag, "n": int(len(df)),
        "LR_chi2": round(float(lr), 3), "df": 1, "p_interaction": p_int,
        "HR_interaction": round(float(s["exp(coef)"]), 3),
        "CI_low": round(float(s["exp(coef) lower 95%"]), 3),
        "CI_high": round(float(s["exp(coef) upper 95%"]), 3),
    }

int_rows = [interaction(d, "Full cohort (predicted HPV)")]
lab = d.dropna(subset=["hpv_label_bin"]).copy()
if lab["os_event"].sum() >= 10:
    lab["hpv_pred_bin"] = lab["hpv_label_bin"].astype(int)
    int_rows.append(interaction(lab, "Directly labelled subset (n=84)"))
int_df = pd.DataFrame(int_rows)
int_df.to_csv(os.path.join(RES, "11_hpv_interaction.csv"), index=False,
              encoding="utf-8-sig")
out["E_interaction"] = int_rows

# ---------------------------------------------------------------- F. 标注子集分层
lab_int = d.dropna(subset=["hpv_label_bin"])
sub_rows = []
for lbl, nm in [(1, "HPV-positive (histology/ISH)"), (0, "HPV-negative (histology/ISH)")]:
    sub = lab_int[lab_int["hpv_label_bin"] == lbl]
    if len(sub) < 15 or sub["os_event"].sum() < 5:
        sub_rows.append({"stratum": nm, "n": int(len(sub)),
                         "events": int(sub["os_event"].sum()), "note": "事件不足"})
        continue
    c = fit_cox(sub.dropna(subset=["CTS", "time_m", "os_event"]), ["CTS"], penalizer=0.1)
    s = c.summary.loc["CTS"]
    sub_rows.append({
        "stratum": nm, "n": int(len(sub)), "events": int(sub["os_event"].sum()),
        "cts_median": round(float(sub["CTS"].median()), 3),
        "HR": round(float(s["exp(coef)"]), 3),
        "CI_low": round(float(s["exp(coef) lower 95%"]), 3),
        "CI_high": round(float(s["exp(coef) upper 95%"]), 3),
        "p": float(s["p"]),
    })
out["F_labelled_subset"] = sub_rows

# ---------------------------------------------------------------- G. 图
fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
for ax, (name, mask) in zip(axes, [("HPV-positive (predicted)", d["hpv_pred_bin"] == 1),
                                   ("HPV-negative (predicted)", d["hpv_pred_bin"] == 0)]):
    sub = d[mask]
    for flag, clr, lbl in [(1, CLR_HI, "CTS-high"), (0, CLR_LO, "CTS-low")]:
        ss = sub[sub["cts_hi"] == flag]
        if len(ss) < 3:
            continue
        KaplanMeierFitter().fit(ss["time_m"], ss["os_event"],
                                label=f"{lbl} (n={len(ss)})").plot_survival_function(
            ax=ax, ci_show=False, color=clr, linewidth=1.8)
    hi, lo = sub["cts_hi"] == 1, sub["cts_hi"] == 0
    if hi.sum() > 3 and lo.sum() > 3:
        lr = logrank_test(sub.loc[hi, "time_m"], sub.loc[lo, "time_m"],
                          sub.loc[hi, "os_event"], sub.loc[lo, "os_event"])
        ax.set_title(f"{name}\nlog-rank p = {lr.p_value:.3f}   (n = {len(sub)})", fontsize=10)
    else:
        ax.set_title(name, fontsize=10)
    ax.set_xlabel("Time (months)")
    ax.set_ylabel("Overall survival")
    ax.set_ylim(0, 1.02)
    ax.set_xlim(0, 120)
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(alpha=0.25)
plt.tight_layout()
plt.savefig(os.path.join(RES, "fig_km_cts_by_hpv.png"), dpi=300)
plt.close()

with open(os.path.join(RES, "11_hpv_cts_summary.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)

print(json.dumps(out, ensure_ascii=False, indent=1))
