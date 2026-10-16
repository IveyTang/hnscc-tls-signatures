# -*- coding: utf-8 -*-
"""
选题四 / Step 5-6 —— 免疫微环境解析 + 临床转化模型（列线图 / 校准 / DCA）

内容：
  Step 5 微环境
    * 基于标记基因的免疫细胞/基质评分（B 细胞、浆细胞、Tfh、生发中心 B、
      CD8 T、Treg、CAF、内皮细胞、巨噬细胞、IFN-γ 程序）
    * CTS 与各评分的 Spearman 相关
    * CTS-high vs low 的差异表达基因数 + 富集方向概览
  Step 6 临床转化
    * 多因素模型（CTS + 临床变量）
    * 基于 Cox 系数的列线图（手动实现，不依赖 R 的 rms）
    * 12/36/60 月 校准曲线（KM 法）
    * 决策曲线分析（DCA）

输出：
  results/06_microenvironment.csv
  results/06_cts_marker_correlation.csv
  results/06_nomogram_points.csv
  results/fig_microenv_correlation.png
  results/fig_nomogram.png
  results/fig_calibration.png
  results/fig_dca.png
  results/06_model_summary.json
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
from scipy import stats

warnings.filterwarnings("ignore")
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROC = os.path.join(ROOT, "data", "processed")
RES = os.path.join(ROOT, "results")

MARKERS = {
    "B 细胞":        ["MS4A1", "CD19", "CD79A", "CD79B", "TCL1A", "BANK1"],
    "生发中心 B":     ["AICDA", "BCL6", "MEF2B", "RGS13", "MKI67"],
    "浆细胞":         ["MZB1", "JCHAIN", "IGHG1", "IGKC", "XBP1", "PRDM1", "SDC1"],
    "Tfh":           ["CXCL13", "IL21", "BCL6", "ICOS", "PDCD1", "SH2D1A"],
    "CD8 T":         ["CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "NKG7"],
    "Treg":          ["FOXP3", "IL2RA", "CTLA4", "IKZF2"],
    "CAF":           ["ACTA2", "FAP", "COL1A1", "COL1A2", "PDGFRB", "TAGLN", "THY1"],
    "内皮/HEV":       ["PECAM1", "VWF", "CLDN5", "PDPN", "CD34"],
    "巨噬细胞":       ["CD68", "CD163", "CSF1R", "MRC1"],
    "IFN-γ 程序":     ["IFNG", "STAT1", "IDO1", "CXCL9", "CXCL10", "GBP1", "IRF1"],
    "检查点":         ["CD274", "PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT"],
    "上皮/肿瘤":      ["EPCAM", "KRT5", "KRT14", "KRT13", "KRT4"],
}


def prep(df):
    d = df.copy()
    d = d[d["os_time_days"].notna() & d["os_event"].notna()]
    d["os_event"] = d["os_event"].astype(int)
    d["age"] = pd.to_numeric(d["age_at_index"], errors="coerce")
    d["male"] = (d["gender"] == "male").astype(int)
    d["stage_iv"] = d["stage"].astype(str).str.startswith("IV").astype(int)
    return d


def marker_scores(expr_log2, sample_index):
    out = {}
    for name, genes in MARKERS.items():
        present = [g for g in genes if g in expr_log2.index]
        if len(present) < 2:
            continue
        z = expr_log2.loc[present]
        z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1).replace(0, np.nan), axis=0)
        out[name] = z.mean(axis=0).reindex(sample_index)
    return pd.DataFrame(out)


def dca_curve(time, event, pred_prob, thresholds):
    """决策曲线分析：净获益 = TP/n - FP/n * (pt/(1-pt))"""
    n = len(time)
    nb_model, nb_all = [], []
    for pt in thresholds:
        yhat = pred_prob >= pt
        tp = int(((yhat == 1) & (event == 1)).sum())
        fp = int(((yhat == 1) & (event == 0)).sum())
        nb_model.append(tp / n - fp / n * (pt / (1 - pt)))
        nb_all.append(((event == 1).mean()) - (1 - (event == 1).mean()) * (pt / (1 - pt)))
    return np.array(nb_model), np.array(nb_all)


def calibration(df, col, times=(12, 36, 60), nbin=4):
    """简化校准曲线：按预测风险四分位分组，比较 KM 实测生存 vs 平均预测生存"""
    res = []
    pred = df[col]
    bins = pd.qcut(pred.rank(method="first"), nbin, labels=False)
    for t in times:
        tm = t * 30.44
        for b in range(nbin):
            sub = df[bins == b]
            if len(sub) < 5:
                continue
            kmf = KaplanMeierFitter().fit(sub["os_time_days"], sub["os_event"])
            obs = float(kmf.predict(tm))
            res.append({"time": t, "bin": b, "n": len(sub),
                        "pred_risk_mean": float(sub[col].mean()), "observed_surv": obs})
    return pd.DataFrame(res)


def main():
    df = pd.read_csv(os.path.join(PROC, "hnsc_tls_scores.csv"), index_col=0)
    expr = pd.read_csv(os.path.join(PROC, "hnsc_gene_matrix.csv"), index_col=0)
    df = prep(df)
    print("[info] 样本 %d，事件 %d" % (len(df), int(df["os_event"].sum())))

    # ---------- 表达矩阵按样本对齐 ----------
    prim_cols = [c for c in expr.columns if "|" in c and c.split("|")[1] == "Primary Tumor"]
    e = expr[prim_cols]
    e.columns = [c.split("|")[0][:12] for c in prim_cols]
    e = e.groupby(level=0, axis=1).mean()
    e = np.log2(e + 1.0)
    common = [b for b in df.index if b in e.columns]
    e = e[common]

    # ---------- Step 5 微环境 ----------
    ms = marker_scores(e, e.columns)
    ms = ms.reindex(df.index)
    ms["CTS"] = df["CTS"].values

    corr_rows = []
    for c in [x for x in ms.columns if x != "CTS"]:
        rho, p = stats.spearmanr(ms[c], ms["CTS"], nan_policy="omit")
        corr_rows.append({"marker": c, "spearman_rho": round(float(rho), 4),
                          "p": float(p)})
    corr = pd.DataFrame(corr_rows).sort_values("spearman_rho", ascending=False)
    corr.to_csv(os.path.join(RES, "06_cts_marker_correlation.csv"), index=False)
    ms.to_csv(os.path.join(RES, "06_microenvironment.csv"), encoding="utf-8-sig")
    print("\n--- CTS 与免疫/基质评分的相关性 ---")
    print(corr.to_string(index=False))

    fig, ax = plt.subplots(figsize=(6.0, 4.2), dpi=300)
    c2 = corr.sort_values("spearman_rho")
    colors = ["#A32D2D" if v > 0 else "#0F6E56" for v in c2["spearman_rho"]]
    ax.barh(c2["marker"], c2["spearman_rho"], color=colors, height=0.62)
    ax.axvline(0, color="#5F5E5A", linewidth=0.8)
    ax.set_xlabel("与 CTS 的 Spearman 相关系数")
    ax.set_title("CTS 的免疫微环境关联特征", fontsize=10)
    ax.grid(axis="x", alpha=0.25, linewidth=0.5)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "fig_microenv_correlation.png"))
    plt.close(fig)

    # ---------- Step 6 多因素模型 ----------
    d = df[["CTS", "age", "male", "stage_iv", "os_time_days", "os_event"]].dropna()
    cph = CoxPHFitter(penalizer=0.05)
    cph.fit(d, "os_time_days", "os_event")
    summ = cph.summary.reset_index().rename(columns={"covariate": "var"})
    summ["HR"] = np.exp(summ["coef"])

    base_haz = cph.baseline_cumulative_hazard_
    base_haz.columns = ["H0"]
    ev_times = np.sort(d.loc[d["os_event"] == 1, "os_time_days"].unique())
    h0 = np.interp(ev_times, base_haz.index.values, base_haz["H0"].values)
    lp = cph.predict_log_partial_hazard(d).values
    S0 = np.exp(-h0)

    def surv_at(t_m):
        idx = np.searchsorted(ev_times, t_m * 30.44)
        idx = min(idx, len(S0) - 1)
        return S0[idx] ** np.exp(lp)

    # 列线图点数表（以最大 |β| 变量为 100 分基准）
    pts = summ.copy()
    coefs = dict(zip(summ["var"], summ["coef"]))
    max_abs = max(abs(v) for v in coefs.values())
    pts["points_per_unit"] = [abs(v) / max_abs * 100 for v in coefs.values()]
    pts["beta"] = [coefs[v] for v in pts["var"]]
    pts_out = pts[["var", "beta", "HR", "points_per_unit", "p"]]
    pts_out.to_csv(os.path.join(RES, "06_nomogram_points.csv"), index=False)

    # 列线图渲染（示意刻度）
    fig, ax = plt.subplots(figsize=(6.2, 0.5 * len(pts_out) + 1.5), dpi=300)
    ax.axis("off")
    n_row = len(pts_out)
    ax.set_xlim(0, 1.0)
    ax.set_ylim(-0.8, n_row + 1.1)
    ypos = np.arange(n_row)[::-1]
    ax.plot([0.15, 0.95], [n_row + 0.4] * 2, color="#2C2C2A", linewidth=0.8)
    for i in range(0, 101, 10):
        x = 0.15 + 0.8 * i / 100
        ax.plot([x, x], [n_row + 0.4, n_row + 0.25], color="#2C2C2A", linewidth=0.6)
        ax.text(x, n_row + 0.58, str(i), ha="center", fontsize=6.5)
    ax.text(0.025, n_row + 0.4, "分值", fontsize=8, va="center")
    for y, (_, r) in zip(ypos, pts_out.iterrows()):
        frac = min(r["points_per_unit"], 100) / 100
        ax.plot([0.15, 0.15 + 0.8 * frac], [y, y], color="#378ADD", linewidth=1.6)
        ax.text(0.142, y, r["var"], ha="right", va="center", fontsize=8)
        ax.text(0.15 + 0.8 * frac + 0.012, y, "%.1f" % r["points_per_unit"],
                va="center", fontsize=7, color="#5F5E5A")
    ax.set_title("CTS + 临床变量 列线图（点数为相对权重）", fontsize=10, y=0.98)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "fig_nomogram.png"))
    plt.close(fig)

    # 12/36/60 月风险与校准
    for t in (12, 36, 60):
        d["risk_%dm" % t] = 1 - surv_at(t)
    cal = calibration(d, "risk_36m", times=(12, 36, 60))
    cal.to_csv(os.path.join(RES, "06_calibration.csv"), index=False)

    fig, ax = plt.subplots(figsize=(4.6, 4.2), dpi=300)
    for t, clr in zip((12, 36, 60), ["#378ADD", "#A32D2D", "#0F6E56"]):
        sub = cal[cal["time"] == t]
        if sub.empty:
            continue
        ax.plot(1 - sub["pred_risk_mean"], sub["observed_surv"], marker="o",
                markersize=4, linewidth=1.2, color=clr, label="%d 月" % t)
    ax.plot([0.2, 1.0], [0.2, 1.0], linestyle="--", color="#B4B2A9", linewidth=0.8)
    ax.set_xlabel("模型预测生存率")
    ax.set_ylabel("KM 实测生存率")
    ax.set_title("校准曲线", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    ax.grid(alpha=0.25, linewidth=0.5)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "fig_calibration.png"))
    plt.close(fig)

    # DCA
    thr = np.arange(0.05, 0.81, 0.01)
    nb_m, nb_all = dca_curve(d["os_time_days"].values, d["os_event"].values,
                             d["risk_36m"].values, thr)
    nb_cts, _ = dca_curve(d["os_time_days"].values, d["os_event"].values,
                          d["CTS"].rank(pct=True).values, thr)
    fig, ax = plt.subplots(figsize=(5.4, 4.0), dpi=300)
    ax.plot(thr, nb_m, color="#A32D2D", linewidth=1.4, label="CTS + 临床模型")
    ax.plot(thr, nb_cts, color="#378ADD", linewidth=1.2, label="CTS 单独")
    ax.plot(thr, nb_all, color="#888780", linestyle="--", linewidth=1.0, label="全部干预")
    ax.axhline(0, color="#B4B2A9", linewidth=0.8)
    ax.set_ylim(-0.05, max(0.35, float(np.nanmax(nb_m)) * 1.15))
    ax.set_xlabel("阈概率")
    ax.set_ylabel("净获益")
    ax.set_title("36 个月决策曲线分析（DCA）", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    ax.grid(alpha=0.25, linewidth=0.5)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "fig_dca.png"))
    plt.close(fig)

    summary = {
        "n_model": int(len(d)),
        "cox_multivariable": summ[["var", "HR", "p"]].round(4).to_dict("records"),
        "top_marker_correlations": corr.head(5).to_dict("records"),
        "top_negative_correlations": corr.tail(3).to_dict("records"),
        "nomogram_points": pts_out.round(3).to_dict("records"),
    }
    json.dump(summary, open(os.path.join(RES, "06_model_summary.json"), "w"),
              ensure_ascii=False, indent=1, default=str)
    print("\n[ok] Step5-6 输出完成")


if __name__ == "__main__":
    main()
