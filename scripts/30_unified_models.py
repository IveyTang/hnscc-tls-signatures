# -*- coding: utf-8 -*-
"""选题四 / Step 30 —— 统一多变量口径 + DCA 补 clinical-only 曲线 + 判别力 bootstrap CI

背景（2026-10-21 终审 S-2/S-3/S-10）：
  S-2 稿件中出现两套"多变量模型"：04 的 ridge 0.1 全协变量集（HR 0.757），
      11 的 penalizer 0.01 精简集（HR 0.660），两处都称 the multivariable model。
      本脚本把 HPV 调整改到主模型协变量集上，全文只保留一套多变量口径。
  S-3 06 的 DCA 只画了 CTS+临床 / CTS 单独 / treat-all，没有 clinical-only，
      正文却写成"加评分带来增量获益"。本脚本补 clinical-only 曲线并量化差值。
  S-10 C-index 与 36 月 AUC 的组间差异无 CI，措辞只能用"无优势证据"。
      本脚本给 1,000 次 bootstrap CI。

输出：
  results/30_unified_models.json
  results/30_dca_curves.csv
  results/30_bootstrap_discrimination.csv
"""
import json
import os
import warnings

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index

warnings.filterwarnings("ignore")
np.random.seed(20261016)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROC = os.path.join(ROOT, "data", "processed")
RES = os.path.join(ROOT, "results")

out = {}

# ============================================================ 读数据
d = pd.read_csv(os.path.join(PROC, "hnsc_with_hpv.csv"))
d = d[d["os_time_days"].notna() & d["os_event"].notna()].copy()
d["os_event"] = d["os_event"].astype(int)
d["age"] = pd.to_numeric(d["age_at_index"], errors="coerce")
d["male"] = (d["sex_at_birth"].astype(str).str.lower() == "male").astype(int)
d["stage_iv"] = d["stage"].astype(str).str.startswith("IV").astype(int)
d["hpv"] = (d["hpv_pred"] == "Positive").astype(int)
d["stage_c"] = d["stage"].astype(str)
d["subsite_c"] = d["subsite"].astype(str)
d["smoke_c"] = d["smoking_status"].fillna("Unknown").astype(str)
print("[info] n=%d  events=%d" % (len(d), int(d["os_event"].sum())))


def fit(df, cols, pen):
    c = CoxPHFitter(penalizer=pen)
    c.fit(df[cols + ["os_time_days", "os_event"]], "os_time_days", "os_event")
    return c


def hr_of(c, term):
    s = c.summary.loc[term]
    return {"HR": round(float(np.exp(s["coef"])), 3),
            "CI_low": round(float(np.exp(s["coef lower 95%"])), 3),
            "CI_high": round(float(np.exp(s["coef upper 95%"])), 3),
            "p": float(s["p"])}


# ============================================================ 1. 主模型（04 口径）
main_base = d[["CTS", "age", "male", "os_time_days", "os_event",
               "stage_c", "subsite_c", "smoke_c", "hpv"]].dropna().copy()
dm = pd.get_dummies(main_base, columns=["stage_c", "subsite_c", "smoke_c"],
                    drop_first=True, dtype=float)
cov_main = [c for c in dm.columns
            if c not in ("os_time_days", "os_event", "hpv")]
# 删除零方差列，与 04 行为一致
zero = [c for c in cov_main if dm[c].nunique() <= 1]
if zero:
    dm = dm.drop(columns=zero)
    cov_main = [c for c in cov_main if c not in zero]

m_main = fit(dm, cov_main, 0.1)
out["main_model_04"] = {
    "n": int(len(dm)), "events": int(dm["os_event"].sum()),
    "covariates": cov_main, "penalizer": 0.1, "CTS": hr_of(m_main, "CTS"),
}
print("[main] CTS %s" % out["main_model_04"]["CTS"])

# ============================================================ 2. 主模型 + HPV（统一口径）
m_main_hpv = fit(dm, cov_main + ["hpv"], 0.1)
out["main_model_plus_hpv"] = {
    "n": int(len(dm)), "events": int(dm["os_event"].sum()),
    "penalizer": 0.1,
    "CTS": hr_of(m_main_hpv, "CTS"),
    "HPV_positive": hr_of(m_main_hpv, "hpv"),
}
print("[main+HPV] CTS %s ; HPV %s"
      % (out["main_model_plus_hpv"]["CTS"], out["main_model_plus_hpv"]["HPV_positive"]))

# ============================================================ 3. 精简模型（11 口径，留作对照）
p_base = d.dropna(subset=["CTS", "age", "male", "stage_iv", "os_time_days", "os_event"])
m_par = fit(p_base, ["CTS", "age", "male", "stage_iv"], 0.01)
m_par_hpv = fit(p_base, ["CTS", "age", "male", "stage_iv", "hpv"], 0.01)
out["parsimonious_11"] = {
    "n": int(len(p_base)), "penalizer": 0.01,
    "CTS_without_hpv": hr_of(m_par, "CTS"),
    "CTS_with_hpv": hr_of(m_par_hpv, "CTS"),
}
print("[parsimonious] %s -> %s" % (out["parsimonious_11"]["CTS_without_hpv"],
                                   out["parsimonious_11"]["CTS_with_hpv"]))

# ============================================================ 4. DCA：补 clinical-only
dd = d.dropna(subset=["CTS", "age", "male", "stage_iv", "os_time_days", "os_event"]).copy()
m_clin = fit(dd, ["age", "male", "stage_iv"], 0.05)
m_comb = fit(dd, ["CTS", "age", "male", "stage_iv"], 0.05)


def surv_at(cph, df, t_month):
    bh = cph.baseline_cumulative_hazard_.copy()
    bh.columns = ["H0"]
    ev = np.sort(df.loc[df["os_event"] == 1, "os_time_days"].unique())
    h0 = np.interp(ev, bh.index.values, bh["H0"].values)
    lp = cph.predict_log_partial_hazard(df).values
    s0 = np.exp(-h0)
    idx = min(np.searchsorted(ev, t_month * 30.44), len(s0) - 1)
    return 1.0 - s0[idx] ** np.exp(lp)


def dca(time, event, pred, thr):
    n = len(time)
    nb = []
    for pt in thr:
        yh = pred >= pt
        tp = int(((yh == 1) & (event == 1)).sum())
        fp = int(((yh == 1) & (event == 0)).sum())
        nb.append(tp / n - fp / n * (pt / (1 - pt)))
    ev = float((event == 1).mean())
    nb_all = ev - (1 - ev) * (thr / (1 - thr))
    return np.array(nb), nb_all


dd["risk_clin"] = surv_at(m_clin, dd, 36)
dd["risk_comb"] = surv_at(m_comb, dd, 36)
thr = np.arange(0.05, 0.81, 0.01)
nb_comb, nb_all = dca(dd["os_time_days"].values, dd["os_event"].values, dd["risk_comb"].values, thr)
nb_clin, _ = dca(dd["os_time_days"].values, dd["os_event"].values, dd["risk_clin"].values, thr)
nb_cts, _ = dca(dd["os_time_days"].values, dd["os_event"].values,
                dd["CTS"].rank(pct=True).values, thr)

pd.DataFrame({"threshold": thr, "nb_combined": nb_comb, "nb_clinical_only": nb_clin,
              "nb_cts_only": nb_cts, "nb_treat_all": nb_all}).to_csv(
    os.path.join(RES, "30_dca_curves.csv"), index=False)

band = (thr >= 0.10) & (thr <= 0.65)
out["dca"] = {
    "n": int(len(dd)),
    "band": [0.10, 0.65],
    "mean_nb_combined": round(float(nb_comb[band].mean()), 4),
    "mean_nb_clinical_only": round(float(nb_clin[band].mean()), 4),
    "mean_nb_cts_only": round(float(nb_cts[band].mean()), 4),
    "mean_delta_combined_minus_clinical": round(float((nb_comb - nb_clin)[band].mean()), 4),
    "delta_min": round(float((nb_comb - nb_clin)[band].min()), 4),
    "delta_max": round(float((nb_comb - nb_clin)[band].max()), 4),
    "frac_thresholds_combined_above_clinical": round(
        float(((nb_comb - nb_clin)[band] > 0).mean()), 3),
    "cox_combined_CTS": hr_of(m_comb, "CTS"),
    "cox_clinical_stage_iv": hr_of(m_clin, "stage_iv"),
}
print("[dca] %s" % json.dumps(out["dca"], ensure_ascii=False))

# ============================================================ 5. 判别力 bootstrap CI
def td_auc(score, time, event, t):
    score = np.asarray(score, float)
    time = np.asarray(time, float)
    event = np.asarray(event)
    case = (event == 1) & (time <= t)
    ctrl = time > t
    y = np.concatenate([np.ones(int(case.sum())), np.zeros(int(ctrl.sum()))])
    s = np.concatenate([score[case], score[ctrl]])
    if y.sum() == 0 or (1 - y).sum() == 0:
        return np.nan
    o = np.argsort(s)
    r = np.empty(len(s), float)
    r[o] = np.arange(1, len(s) + 1)
    n1, n0 = y.sum(), (1 - y).sum()
    return float((r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


b = d.dropna(subset=["CTS", "G7_Tcell_inflamed_cytolytic", "os_time_days", "os_event"])
t36 = 36 * 30.44
rows = []
rng = np.random.default_rng(20261016)
for i in range(1000):
    idx = rng.integers(0, len(b), len(b))
    s = b.iloc[idx]
    c_cts = concordance_index(s["os_time_days"], s["CTS"], s["os_event"])
    c_g7 = concordance_index(s["os_time_days"], s["G7_Tcell_inflamed_cytolytic"], s["os_event"])
    a_cts = td_auc(-s["CTS"].values, s["os_time_days"].values, s["os_event"].values, t36)
    a_g7 = td_auc(-s["G7_Tcell_inflamed_cytolytic"].values, s["os_time_days"].values,
                  s["os_event"].values, t36)
    rows.append({"c_CTS": c_cts, "c_G7": c_g7, "c_diff": c_cts - c_g7,
                 "auc36_CTS": a_cts, "auc36_G7": a_g7, "auc36_diff": a_cts - a_g7})
bs = pd.DataFrame(rows).dropna()
bs.to_csv(os.path.join(RES, "30_bootstrap_discrimination.csv"), index=False)


def ci(x):
    return [round(float(np.percentile(x, 2.5)), 4), round(float(np.percentile(x, 97.5)), 4)]


out["discrimination"] = {
    "n": int(len(b)),
    "cindex_CTS": round(float(concordance_index(b["os_time_days"], b["CTS"], b["os_event"])), 4),
    "cindex_G7": round(float(concordance_index(b["os_time_days"],
                                               b["G7_Tcell_inflamed_cytolytic"],
                                               b["os_event"])), 4),
    "cindex_diff_CI": ci(bs["c_diff"]),
    "cindex_diff_p_two_sided": round(float(2 * min((bs["c_diff"] <= 0).mean(),
                                                   (bs["c_diff"] >= 0).mean())), 4),
    "auc36_CTS": round(float(td_auc(-b["CTS"].values, b["os_time_days"].values,
                                    b["os_event"].values, t36)), 4),
    "auc36_G7": round(float(td_auc(-b["G7_Tcell_inflamed_cytolytic"].values,
                                   b["os_time_days"].values, b["os_event"].values, t36)), 4),
    "auc36_diff_CI": ci(bs["auc36_diff"]),
    "auc36_diff_p_two_sided": round(float(2 * min((bs["auc36_diff"] <= 0).mean(),
                                                  (bs["auc36_diff"] >= 0).mean())), 4),
}
print("[disc] %s" % json.dumps(out["discrimination"], ensure_ascii=False))

json.dump(out, open(os.path.join(RES, "30_unified_models.json"), "w"),
          ensure_ascii=False, indent=1)
print("\n[ok] 写入 results/30_*")
