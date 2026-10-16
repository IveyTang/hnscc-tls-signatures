# -*- coding: utf-8 -*-
"""
选题四 / Step 3 —— 内部验证（TCGA-HNSC）

内容：
  1. KM 生存曲线 + log-rank（最优切点 / 中位数切点）
  2. 单因素与多因素 Cox（校正年龄、性别、分期、吸烟、亚部位）
  3. 时间依赖 AUC（1/3/5 年，累计/动态定义）
  4. C-index
  5. 签名收敛性热图（对应探索性目标 E1）

输出：
  results/04_cox_univariate.csv
  results/04_cox_multivariate.csv
  results/04_time_dependent_auc.csv
  results/04_validation_summary.json
  results/fig_km_cts.png
  results/fig_forest_multivariate.png
  results/fig_roc_timedependent.png
  results/fig_signature_convergence.png
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
from lifelines.utils import concordance_index

warnings.filterwarnings("ignore")
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROC = os.path.join(ROOT, "data", "processed")
RES = os.path.join(ROOT, "results")
CLR_HI, CLR_LO, CLR_N = "#A32D2D", "#0F6E56", "#B4B2A9"


def td_auc(score, time, event, t):
    """累计/动态时间依赖 AUC（简化 Heagerty 定义）。"""
    score = np.asarray(score, dtype=float)
    time = np.asarray(time, dtype=float)
    event = np.asarray(event)
    case = (event == 1) & (time <= t)
    ctrl = time > t
    y = np.concatenate([np.ones(int(case.sum())), np.zeros(int(ctrl.sum()))])
    s = np.concatenate([score[case], score[ctrl]])
    if y.sum() == 0 or (1 - y).sum() == 0:
        return None
    order = np.argsort(s)
    ranks = np.empty(len(s), float)
    ranks[order] = np.arange(1, len(s) + 1)
    n1, n0 = y.sum(), (1 - y).sum()
    return float((ranks[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def prep(df):
    d = df.copy()
    d = d[d["os_time_days"].notna() & d["os_event"].notna()]
    d["os_event"] = d["os_event"].astype(int)
    d["age"] = pd.to_numeric(d["age_at_index"], errors="coerce")
    d["male"] = (d["gender"] == "male").astype(int)
    d["stage3"] = d["stage"].astype(str)
    d["subsite"] = d["subsite"].astype(str)
    d["smoke"] = d["smoking_status"].fillna("Unknown").astype(str)
    return d


def main():
    df = pd.read_csv(os.path.join(PROC, "hnsc_tls_scores.csv"))
    df = prep(df)
    print("[info] 验证样本: %d，死亡事件: %d" % (len(df), int(df["os_event"].sum())))

    # ---------- 1. KM + log-rank ----------
    res = {}
    for gcol, tag in [("CTS_group_opt", "optimal"), ("CTS_group_median", "median")]:
        hi = df[df[gcol] == "High"]
        lo = df[df[gcol] == "Low"]
        lr = logrank_test(hi["os_time_days"], lo["os_time_days"],
                          hi["os_event"], lo["os_event"])
        res["logrank_%s" % tag] = {"p": float(lr.p_value),
                                   "chi2": float(lr.test_statistic),
                                   "n_high": int(len(hi)), "n_low": int(len(lo))}
        print("[ok] log-rank(%s): p=%.4g, χ²=%.2f, High=%d, Low=%d"
              % (tag, lr.p_value, lr.test_statistic, len(hi), len(lo)))

    fig, ax = plt.subplots(figsize=(5.2, 4.2), dpi=300)
    for grp, clr in [("High", CLR_HI), ("Low", CLR_LO)]:
        sub = df[df["CTS_group_opt"] == grp]
        kmf = KaplanMeierFitter().fit(sub["os_time_days"] / 30.44, sub["os_event"],
                                      label="CTS %s (n=%d)" % (grp, len(sub)))
        kmf.plot_survival_function(ax=ax, color=clr, ci_show=True, linewidth=1.4)
    ax.set_xlabel("总生存时间（月）")
    ax.set_ylabel("总生存率")
    ax.set_title("TCGA-HNSC：共识 TLS 评分（CTS）与总生存", fontsize=10)
    ax.text(0.04, 0.06, "log-rank p = %.3g" % res["logrank_optimal"]["p"],
            transform=ax.transAxes, fontsize=9)
    ax.grid(alpha=0.25, linewidth=0.5)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "fig_km_cts.png"))
    plt.close(fig)

    # ---------- 2. Cox ----------
    cat_vars = {"stage3", "subsite", "smoke"}
    covs = ["CTS", "age", "male", "stage3", "subsite", "smoke"]
    dm = df[covs + ["os_time_days", "os_event"]].copy()
    for v in cat_vars:
        dm[v] = dm[v].astype("category").cat.codes
    dm = dm.astype(float).dropna()
    # 剔除零方差变量（否则 Cox 的 Hessian 奇异，报 delta contains nan）
    zero_var = [c for c in dm.columns if dm[c].nunique() <= 1]
    if zero_var:
        print("[warn] 剔除零方差变量: %s" % zero_var)
        dm = dm.drop(columns=zero_var)
        covs = [c for c in covs if c not in zero_var]
    print("[info] 建模样本（无缺失）: %d，事件 %d，变量 %d 个"
          % (len(dm), int(dm["os_event"].sum()), len(covs)))

    uni_rows = []
    for v in covs:
        sub = dm[[v, "os_time_days", "os_event"]]
        try:
            c = CoxPHFitter(penalizer=0.01).fit(sub, "os_time_days", "os_event")
            row = c.summary.iloc[0]
            uni_rows.append({"var": v, "HR": float(np.exp(row["coef"])),
                             "lower": float(np.exp(row["coef lower 95%"])),
                             "upper": float(np.exp(row["coef upper 95%"])),
                             "p": float(row["p"])})
        except Exception as e:
            print("  [warn] 单因素 %s 拟合失败: %s" % (v, e))
    uni = pd.DataFrame(uni_rows)
    uni.to_csv(os.path.join(RES, "04_cox_univariate.csv"), index=False)
    print("\n--- 单因素 Cox ---")
    print(uni.to_string(index=False))

    # 多因素：将分类变量转为哑变量
    d2 = pd.get_dummies(
        df[["CTS", "age", "male", "os_time_days", "os_event"]].assign(
            stage3=df["stage3"], subsite=df["subsite"], smoke=df["smoke"]),
        columns=["stage3", "subsite", "smoke"], drop_first=True, dtype=float)
    d2 = d2.dropna()
    if d2["os_event"].sum() < 10 or len(d2) < 30:
        raise SystemExit("[错误] 建模样本不足")
    cph2 = CoxPHFitter(penalizer=0.1)
    cph2.fit(d2, "os_time_days", "os_event")
    summ = cph2.summary.reset_index().rename(columns={"covariate": "var"})
    summ["HR"] = np.exp(summ["coef"])
    summ["lower"] = np.exp(summ["coef lower 95%"])
    summ["upper"] = np.exp(summ["coef upper 95%"])
    keep = ~summ["var"].str.startswith(("stage3_", "subsite_", "smoke_"))
    out_m = summ[["var", "HR", "lower", "upper", "p"]]
    out_m.to_csv(os.path.join(RES, "04_cox_multivariate.csv"), index=False)
    print("\n--- 多因素 Cox（前 12 行）---")
    print(out_m.head(12).to_string(index=False))

    # 森林图（仅主要变量）
    plot_m = out_m[keep].copy()
    fig, ax = plt.subplots(figsize=(5.0, 0.5 * len(plot_m) + 1.4), dpi=300)
    y = np.arange(len(plot_m))[::-1]
    ax.errorbar(plot_m["HR"], y,
                xerr=[plot_m["HR"] - plot_m["lower"], plot_m["upper"] - plot_m["HR"]],
                fmt="o", color="#A32D2D", ecolor="#888780", capsize=3, markersize=5)
    ax.axvline(1, color="#B4B2A9", linestyle="--", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(plot_m["var"])
    ax.set_xscale("log")
    ax.set_xlabel("风险比 HR（95% CI，对数轴）")
    ax.set_title("多因素 Cox（校正分期、亚部位、吸烟等）", fontsize=10)
    ax.grid(alpha=0.25, linewidth=0.5)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "fig_forest_multivariate.png"))
    plt.close(fig)

    # ---------- 3. 时间依赖 AUC ----------
    # CTS 呈保护性方向（分值越高生存越好）；统一以 -评分 作为风险评分，
    # 使 AUC 以 0.5 为无判别力基准。
    rows = []
    for t_m in [12, 36, 60]:
        t = t_m * 30.44
        for col in ["CTS", "G3_HNSCC_13chemokine", "G1_Coppola_12chemokine",
                    "G7_Tcell_inflamed_cytolytic"]:
            if col not in df.columns:
                continue
            a = td_auc(-df[col].values, df["os_time_days"].values,
                       df["os_event"].values, t)
            rows.append({"time_months": t_m, "score": col, "AUC": a})
    auc = pd.DataFrame(rows)
    auc.to_csv(os.path.join(RES, "04_time_dependent_auc.csv"), index=False)
    print("\n--- 时间依赖 AUC（以 -评分 为风险方向，0.5 = 无判别力）---")
    print(auc.pivot(index="score", columns="time_months", values="AUC").round(3).to_string())

    fig, ax = plt.subplots(figsize=(5.2, 3.8), dpi=300)
    piv = auc.pivot(index="time_months", columns="score", values="AUC")
    for c in piv.columns:
        ax.plot(piv.index, piv[c], marker="o", markersize=4, linewidth=1.3, label=c)
    ax.set_xticks([12, 36, 60])
    ax.set_xlabel("时间（月）")
    ax.set_ylabel("时间依赖 AUC")
    ax.set_ylim(0.4, 0.8)
    ax.legend(fontsize=7, frameon=False)
    ax.grid(alpha=0.25, linewidth=0.5)
    ax.set_title("CTS 与单基因集签名的时间依赖判别力", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "fig_roc_timedependent.png"))
    plt.close(fig)

    # ---------- 4. C-index ----------
    # lifelines 的 concordance_index 约定：分值越高 = 生存越长。
    # CTS 为保护性方向，故直接传入 +CTS。
    cidx = concordance_index(df["os_time_days"], df["CTS"], df["os_event"])
    cidx_g3 = concordance_index(df["os_time_days"], df["G3_HNSCC_13chemokine"],
                                df["os_event"])
    cidx_g7 = concordance_index(df["os_time_days"], df["G7_Tcell_inflamed_cytolytic"],
                                df["os_event"])
    print("\n[ok] C-index（lifelines 约定：分值越高生存越长）: "
          "CTS=%.3f, G3(13基因签名)=%.3f, G7(泛T细胞炎症)=%.3f"
          % (cidx, cidx_g3, cidx_g7))
    for g in ["High", "Low"]:
        sub = df[df["CTS_group_opt"] == g]
        kmf = KaplanMeierFitter().fit(sub["os_time_days"], sub["os_event"])
        print("     CTS %-4s 组 (n=%3d) 中位 OS = %.1f 月"
              % (g, len(sub), float(kmf.median_survival_time_) / 30.44))

    # ---------- 5. 签名收敛性热图 ----------
    conv_cols = [c for c in ["G1_Coppola_12chemokine", "G2_11chemokine_HCC",
                             "G3_HNSCC_13chemokine", "G4_TLS_maturation_GC",
                             "G5_Th1_IFNgamma", "G6_HEV_lymphoid_stroma",
                             "G7_Tcell_inflamed_cytolytic"] if c in df.columns]
    cm = df[conv_cols].corr(method="spearman")
    fig, ax = plt.subplots(figsize=(5.6, 4.8), dpi=300)
    im = ax.imshow(cm.values, cmap="RdYlBu_r", vmin=0, vmax=1)
    ax.set_xticks(range(len(cm)))
    ax.set_xticklabels([c.split("_")[0] for c in cm.columns])
    ax.set_yticks(range(len(cm)))
    ax.set_yticklabels([c.split("_")[0] for c in cm.index])
    for i in range(len(cm)):
        for j in range(len(cm)):
            ax.text(j, i, "%.2f" % cm.values[i, j], ha="center", va="center", fontsize=7)
    fig.colorbar(im, ax=ax, shrink=0.8, label="Spearman ρ")
    ax.set_title("TLS / 免疫基因集两两一致性（验证 E1）", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "fig_signature_convergence.png"))
    plt.close(fig)

    tol = [c for c in conv_cols if c not in ("G5_Th1_IFNgamma", "G7_Tcell_inflamed_cytolytic")]
    sub_cm = cm.loc[tol, tol].values
    off = sub_cm[np.triu_indices(len(tol), 1)]

    summary = {
        "n": int(len(df)), "n_events": int(df["os_event"].sum()),
        "logrank": res,
        "multivariable_CTS": out_m[out_m["var"] == "CTS"].to_dict("records"),
        "cindex_CTS": float(cidx), "cindex_G3": float(cidx_g3),
        "cindex_G7_Tcell_inflamed": float(cidx_g7),
        "td_auc": auc.pivot(index="score", columns="time_months", values="AUC")
                     .round(4).to_dict(),
        "tls_signature_convergence_median_rho": float(np.median(off)),
        "tls_signature_convergence_range": [float(off.min()), float(off.max())],
    }
    json.dump(summary, open(os.path.join(RES, "04_validation_summary.json"), "w"),
              ensure_ascii=False, indent=1, default=str)
    print("\n[ok] 全部输出写入 results/")


if __name__ == "__main__":
    main()
