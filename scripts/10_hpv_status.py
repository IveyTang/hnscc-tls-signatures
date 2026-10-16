# -*- coding: utf-8 -*-
"""构建 TCGA-HNSC 的 HPV 状态：直接标注 + 转录组分类器。

标注来源（三源合并，优先级 ISH > p16 > 病理报告文本）：
  1. BCR XML  hnsc:hpv_status_by_ish_testing  —— 直接病毒核酸原位杂交（最可靠）
  2. BCR XML  hnsc:hpv_status_by_p16_testing  —— p16 IHC 替代指标（口咽部可靠）
  3. Pathology Report PDF 文本                —— 补充来源

分类器：107 基因宿主表达面板（含 CDKN2A/p16、RB 通路、鳞状分化轴、免疫轴）
  - 采用嵌套 LOO 交叉验证：每折内独立做特征选择，避免信息泄漏
  - 单基因 CDKN2A 作为参考基线
  - 输出 LOO-AUC / 灵敏度 / 特异度

生物学合理性检查（决定结果可否写入论文）：
  - 预测 HPV+ 应在口咽部富集
  - 预测 HPV+ 应有更好的 OS
  - 预测 HPV+ 应有更高的 CTS（TLS 评分）
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.preprocessing import StandardScaler

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw")
PROC = os.path.join(ROOT, "data", "processed")
RES = os.path.join(ROOT, "results")
os.makedirs(RES, exist_ok=True)

OROPHARYNX = {"口咽"}


# ---------------------------------------------------------------- 标签构建
def consensus_label(row):
    """返回 (hpv_call, source)。优先级 ISH > p16 > 报告文本。"""
    ish = row.get("hpv_ish")
    p16 = row.get("hpv_p16")
    rep = row.get("hpv_from_report")
    if isinstance(ish, str) and ish in ("Positive", "Negative"):
        return ish, "ISH"
    if isinstance(p16, str) and p16 in ("Positive", "Negative"):
        return p16, "p16"
    if isinstance(rep, str) and rep in ("Positive", "Negative"):
        return rep, "pathology report"
    return None, None


def build_labels(score_df):
    bcr = pd.read_csv(os.path.join(RAW, "tcga_hnsc_bcr_hpv.csv"))
    # CSV 同时含 XML 解析出的 barcode 与 submitter_id，二者等价，先删一列再统一列名
    if "barcode" in bcr.columns:
        bcr = bcr.drop(columns=["barcode"])
    bcr = bcr.rename(columns={"submitter_id": "barcode"})
    keep = ["barcode", "hpv_p16", "hpv_ish", "anatomic_subdivision",
            "pack_years", "tobacco_smoking_history", "margin_status", "lvi",
            "pni", "ln_examined", "ln_positive", "ene", "histologic_grade",
            "bcr_cancer_status"]
    bcr = bcr[[c for c in keep if c in bcr.columns]]

    rep_p = os.path.join(RAW, "tcga_hnsc_path_report_hpv.csv")
    if os.path.exists(rep_p):
        rep = pd.read_csv(rep_p).rename(columns={"submitter_id": "barcode"})
        rep = rep[["barcode", "hpv_from_report"]].dropna(subset=["hpv_from_report"])
        # 同一患者多份报告：任一 Positive 优先
        rep = rep.groupby("barcode")["hpv_from_report"].agg(
            lambda s: "Positive" if (s == "Positive").any() else "Negative").reset_index()
        bcr = bcr.merge(rep, on="barcode", how="left")
    else:
        bcr["hpv_from_report"] = None

    df = score_df.merge(bcr, on="barcode", how="left")
    lab = df.apply(consensus_label, axis=1, result_type="expand")
    df["hpv_label"], df["hpv_label_source"] = lab[0], lab[1]
    return df


# ---------------------------------------------------------------- 分类器
def nested_loo(X, y, k_list=(3, 5, 8, 12, 20), seed=0):
    """嵌套 LOO：每折内独立做单变量特征选择，再拟合 L2 logistic。

    返回 {k: {'auc':..., 'y_true':..., 'y_prob':...}}
    """
    n = len(y)
    out = {}
    for k in k_list:
        probs = np.zeros(n)
        for i in range(n):
            tr = np.ones(n, bool)
            tr[i] = False
            Xtr, ytr = X[tr], y[tr]
            # --- 折内特征选择：按 |标准化均值差| 排序 ---
            mu1, mu0 = Xtr[ytr == 1].mean(0), Xtr[ytr == 0].mean(0)
            sd = Xtr.std(0) + 1e-9
            score = np.abs(mu1 - mu0) / sd
            top = np.argsort(-score)[:k]
            sc = StandardScaler().fit(Xtr[:, top])
            clf = LogisticRegression(C=1.0, max_iter=5000, solver="lbfgs")
            clf.fit(sc.transform(Xtr[:, top]), ytr)
            probs[i] = clf.predict_proba(sc.transform(X[i:i + 1, top]))[0, 1]
        out[k] = {"auc": float(roc_auc_score(y, probs)),
                  "y_true": y.copy(), "y_prob": probs}
    return out


def youden_threshold(y, p):
    fpr, tpr, thr = roc_curve(y, p)
    j = tpr - fpr
    i = int(np.argmax(j))
    return float(thr[i]), float(tpr[i]), float(1 - fpr[i])


def main():
    score_df = pd.read_csv(os.path.join(PROC, "hnsc_tls_scores.csv"))
    expr = pd.read_csv(os.path.join(PROC, "hnsc_gene_matrix_hpv.csv"), index_col=0)

    # 样本列名 -> barcode（列形如 'TCGA-XX-XXXX-01A|Primary Tumor'）
    def to_barcode(c):
        return c.split("|")[0][:12]

    ecols = {}
    for c in expr.columns:
        bc = to_barcode(c)
        ecols.setdefault(bc, []).append(c)
    emat = pd.DataFrame(
        {bc: expr[cs].mean(axis=1) for bc, cs in ecols.items()})
    emat = np.log2(emat + 1.0)

    df = build_labels(score_df)
    df = df.set_index("barcode", drop=False)
    common = [b for b in df.index if b in emat.columns]
    E = emat[common].T                      # 样本 × 基因
    df = df.loc[common]

    print("[info] 分析集 %d 例，表达矩阵可用 %d 例，基因 %d 个"
          % (len(score_df), len(df), E.shape[1]))

    # ---------- 三源标注统计 ----------
    print("\n--- HPV 直接标注可用性 ---")
    for col in ["hpv_p16", "hpv_ish", "hpv_from_report"]:
        if col in df.columns:
            v = df[col].value_counts(dropna=True)
            print("  %-18s %s" % (col, dict(v)))

    print("\n--- 共识标签（优先 ISH > p16 > 报告） ---")
    print("  " + str(df["hpv_label_source"].value_counts(dropna=False).to_dict()))
    print("  标签分布:", df["hpv_label"].value_counts(dropna=False).to_dict())

    lab = df["hpv_label"].notna()
    print("\n  可用于训练/验证: %d 例（阳性 %d / 阴性 %d）"
          % (lab.sum(), int((df.loc[lab, "hpv_label"] == "Positive").sum()),
             int((df.loc[lab, "hpv_label"] == "Negative").sum())))

    # 标注 vs 部位的一致性（p16 已知口咽部可靠）
    print("\n--- 直接标注 × 解剖部位 ---")
    ct = pd.crosstab(df.loc[lab, "subsite"], df.loc[lab, "hpv_label"])
    print(ct.to_string())

    if lab.sum() < 30:
        print("\n[stop] 标注样本不足 30 例，无法建模。")
        return

    # ---------- 分类器 ----------
    Xlab = E.loc[df.index[lab]].values
    ylab = (df.loc[lab, "hpv_label"] == "Positive").astype(int).values

    # 单基因 CDKN2A 基线
    cidx = list(E.columns).index("CDKN2A") if "CDKN2A" in E.columns else None
    base_auc = None
    if cidx is not None:
        base_auc = float(roc_auc_score(ylab, Xlab[:, cidx]))
        print("\n[基线] 单基因 CDKN2A(p16) AUC = %.3f" % base_auc)

    print("\n[run] 嵌套 LOO 交叉验证（特征选择在折内进行）...")
    sys.stdout.flush()
    res = nested_loo(Xlab, ylab)
    print("\n--- 嵌套 LOO 性能 ---")
    print("  %-6s %-8s %-10s %-10s" % ("k", "AUC", "灵敏度", "特异度"))
    best_k, best_auc = None, -1
    for k, r in res.items():
        t, se, sp = youden_threshold(ylab, r["y_prob"])
        print("  %-6d %-8.3f %-10.3f %-10.3f" % (k, r["auc"], se, sp))
        if r["auc"] > best_auc:
            best_auc, best_k = r["auc"], k
    print("\n[ok] 最优 k=%d, LOO-AUC=%.3f" % (best_k, best_auc))

    # 最终模型：全标注集上重拟合
    mu1, mu0 = Xlab[ylab == 1].mean(0), Xlab[ylab == 0].mean(0)
    sd = Xlab.std(0) + 1e-9
    score = np.abs(mu1 - mu0) / sd
    top = np.argsort(-score)[:best_k]
    top_genes = list(np.array(E.columns)[top])
    sc = StandardScaler().fit(Xlab[:, top])
    clf = LogisticRegression(C=1.0, max_iter=5000)
    clf.fit(sc.transform(Xlab[:, top]), ylab)
    print("     入选基因: %s" % ", ".join(top_genes))
    print("     系数: %s" % ", ".join(
        "%s=%+.2f" % (g, c) for g, c in zip(top_genes, clf.coef_[0])))

    thr, se, sp = youden_threshold(ylab, res[best_k]["y_prob"])
    print("     LOO-Youden 阈值 = %.4f（灵敏度 %.3f / 特异度 %.3f）" % (thr, se, sp))

    # ---------- 应用于全队列 ----------
    Xall = E.values
    p_all = clf.predict_proba(sc.transform(Xall[:, top]))[:, 1]
    df["hpv_prob"] = p_all
    df["hpv_pred"] = np.where(p_all > thr, "Positive", "Negative")
    df["hpv_final"] = df["hpv_label"].where(df["hpv_label"].notna(), df["hpv_pred"])
    df["hpv_source"] = np.where(df["hpv_label"].notna(),
                                df["hpv_label_source"], "classifier")

    print("\n--- 全队列 HPV 状态 ---")
    print("  " + str(df["hpv_final"].value_counts().to_dict()))
    print("  " + str(df["hpv_source"].value_counts().to_dict()))

    # ---------- 生物学合理性检查 ----------
    print("\n=== 合理性检查 ===")
    print("1) 预测 HPV+ 是否在口咽富集")
    sub_pred = df[df["hpv_source"] == "classifier"]
    ct2 = pd.crosstab(sub_pred["subsite"], sub_pred["hpv_pred"])
    ct2["阳性率"] = (100 * ct2.get("Positive", 0) / ct2.sum(axis=1)).round(1)
    print(ct2.to_string())

    print("\n2) 预测 HPV+ 是否 OS 更优")
    from lifelines import CoxPHFitter
    from lifelines.statistics import logrank_test
    d = df.dropna(subset=["os_time_days", "os_event"]).copy()
    d["hpv_pos"] = (d["hpv_final"] == "Positive").astype(int)
    g = d["hpv_pos"] == 1
    lr = logrank_test(d.loc[g, "os_time_days"], d.loc[~g, "os_time_days"],
                      d.loc[g, "os_event"], d.loc[~g, "os_event"])
    cph = CoxPHFitter().fit(d[["hpv_pos", "os_time_days", "os_event"]],
                            "os_time_days", "os_event")
    hr = float(np.exp(cph.params_["hpv_pos"]))
    p = float(cph.summary.loc["hpv_pos", "p"])
    print("   log-rank p=%.4g | Cox HR=%.2f (p=%.4g)  [n=%d, 事件 %d]"
          % (lr.p_value, hr, p, len(d), int(d["os_event"].sum())))

    print("\n3) 预测 HPV+ 的 CTS 是否更高")
    from scipy import stats
    a = d.loc[d["hpv_final"] == "Positive", "CTS"]
    b = d.loc[d["hpv_final"] == "Negative", "CTS"]
    if len(a) > 2 and len(b) > 2:
        u = stats.mannwhitneyu(a, b)
        print("   CTS 中位数  HPV+ %.3f vs HPV- %.3f | MWU p=%.4g"
              % (a.median(), b.median(), u.pvalue))

    print("\n4) 直接标注病例上的分类器一致性（即 LOO 结果）")
    rt, rs, rp = youden_threshold(ylab, res[best_k]["y_prob"])
    print("   准确率 %.3f | AUC %.3f | 灵敏度 %.3f | 特异度 %.3f"
          % ((np.array(res[best_k]["y_prob"]) > rt).astype(int).__eq__(ylab).mean(),
             best_auc, rs, rp))

    # ---------- 输出 ----------
    outcols = ["barcode", "subsite", "stage", "t_stage", "n_stage", "age_at_index",
               "sex_at_birth", "pack_years", "os_time_days", "os_event", "CTS",
               "hpv_label", "hpv_label_source", "hpv_prob", "hpv_pred",
               "hpv_final", "hpv_source"]
    out = df[[c for c in outcols if c in df.columns]]
    out.to_csv(os.path.join(PROC, "hnsc_hpv_status.csv"),
               index=False, encoding="utf-8-sig")

    # 供后续合并用：展开在分析集上
    df.to_csv(os.path.join(PROC, "hnsc_with_hpv.csv"),
              index=False, encoding="utf-8-sig")

    summ = {
        "n_analysis": int(len(df)),
        "n_labeled": int(lab.sum()),
        "n_labeled_pos": int((df.loc[lab, "hpv_label"] == "Positive").sum()),
        "label_sources": {k: int(v) for k, v in
                          df["hpv_label_source"].value_counts().items()},
        "cdkn2a_auc": base_auc,
        "best_k": int(best_k),
        "loo_auc": float(best_auc),
        "loo_sens": float(se), "loo_spec": float(sp),
        "threshold": float(thr),
        "top_genes": top_genes,
        "coefficients": {g: float(c) for g, c in zip(top_genes, clf.coef_[0])},
        "n_pred_pos": int((df["hpv_final"] == "Positive").sum()),
        "n_pred_neg": int((df["hpv_final"] == "Negative").sum()),
        "os_logrank_p": float(lr.p_value), "os_hr": hr, "os_hr_p": p,
        "cts_pos_median": float(a.median()), "cts_neg_median": float(b.median()),
        "cts_mwu_p": float(u.pvalue),
        "per_k": {str(k): {"auc": v["auc"]} for k, v in res.items()},
    }
    with open(os.path.join(RES, "10_hpv_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summ, fh, ensure_ascii=False, indent=1)

    # LOO 概率存盘，供画 ROC 用
    np.save(os.path.join(RES, "10_loo_y.npy"), ylab)
    np.save(os.path.join(RES, "10_loo_p.npy"), res[best_k]["y_prob"])
    print("\n[ok] 输出 -> data/processed/hnsc_hpv_status.csv")


if __name__ == "__main__":
    main()
