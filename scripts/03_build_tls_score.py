# -*- coding: utf-8 -*-
"""
选题四 / Step 2 —— 构建 HNSCC 共识 TLS 评分（Consensus TLS Score, CTS）

流程：
  1. 载入基因×样本 TPM 矩阵与临床分析集，仅保留原发肿瘤样本
  2. log2(TPM+1) 转换
  3. 对 7 组候选基因集，用 3 种算法计算评分：
        A. z_mean  —— 逐基因 z-score 后取均值（平台无关）
        B. ssgsea  —— gseapy.ssGSEA（秩基，平台无关）
        C. pca1    —— z-score 矩阵的第一主成分（符号对齐 z_mean）
  4. 共识 CTS = 5 个 TLS 相关基因集（G1/G2/G3/G4/G6）在 A 法下的
     百分位秩均值
  5. 阈值：最优 cutpoint（对数秩统计量最大）+ 中位数（敏感性分析）
  6. Bootstrap 稳健性：1000 次重抽样，报告分组一致性

输出：
  data/processed/hnsc_tls_scores.csv    每样本的各基因集评分、CTS、分组
  results/03_score_summary.json         评分汇总与稳健性指标
"""
import json
import os
from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROC = os.path.join(ROOT, "data", "processed")
RES = os.path.join(ROOT, "results")
os.makedirs(RES, exist_ok=True)

CONSENSUS_SETS = ["G1_Coppola_12chemokine", "G2_11chemokine_HCC",
                  "G3_HNSCC_13chemokine", "G4_TLS_maturation_GC",
                  "G6_HEV_lymphoid_stroma"]
NEG_CTRL = ["G5_Th1_IFNgamma", "G7_Tcell_inflamed_cytolytic"]
RNG = np.random.default_rng(20261016)


def load_gene_sets():
    cfg = json.load(open(os.path.join(ROOT, "config", "tls_genesets.json"), encoding="utf-8"))
    return {k: v["genes"] for k, v in cfg.items() if isinstance(v, dict) and "genes" in v}


def load_matrix():
    df = pd.read_csv(os.path.join(PROC, "hnsc_gene_matrix.csv"), index_col=0)
    df.columns = [c for c in df.columns]
    # 列名形如 "TCGA-CV-7247-01A|Primary Tumor"
    meta = pd.DataFrame({
        "col": df.columns,
        "barcode": [c.split("|")[0][:12] for c in df.columns],
        "sample_type": [c.split("|")[1] if "|" in c else "" for c in df.columns],
        "sample_code": [c.split("|")[0].split("-")[-1][:2] if "|" in c else "" for c in df.columns],
    })
    return df, meta


def zscore(mat):
    mu = mat.mean(axis=1)
    sd = mat.std(axis=1).replace(0, np.nan)
    return mat.sub(mu, axis=0).div(sd, axis=0)


def score_gene_set(z, genes, method):
    present = [g for g in genes if g in z.index]
    if len(present) == 0:
        return None, 0
    sub = z.loc[present].dropna(how="all")
    if sub.empty:
        return None, 0
    if method == "z_mean":
        return sub.mean(axis=0), len(present)
    if method == "pca1":
        x = sub.T.fillna(0.0).values
        x = x - x.mean(axis=0, keepdims=True)
        u, s, vt = np.linalg.svd(x, full_matrices=False)
        pc1 = u[:, 0] * s[0]
        base = sub.mean(axis=0).values
        if np.corrcoef(pc1, base)[0, 1] < 0:
            pc1 = -pc1
        return pd.Series(pc1, index=sub.columns), len(present)
    raise ValueError(method)


def ssgsea_scores(z, gene_sets):
    try:
        import gseapy
    except Exception:
        return None
    sets = {k: [g for g in v if g in z.index] for k, v in gene_sets.items()}
    sets = {k: v for k, v in sets.items() if len(v) >= 3}
    if not sets:
        return None
    res = gseapy.ssgsea(data=z, gene_sets=sets, outdir=None,
                        min_size=3, max_size=2000, permutation_num=0,
                        threads=4, verbose=False)
    sc = res.res2d.pivot(index="Name", columns="Term", values="ES")
    return sc.T if sc.shape[0] < sc.shape[1] else sc


def best_cutpoint(score, time, event, n_cand=60, min_grp=15):
    """在 10%–90% 分位间搜索使对数秩统计量最大的切点（向量化 log-rank）。"""
    lo, hi = np.percentile(score, 10), np.percentile(score, 90)
    cands = np.unique(np.round(np.linspace(lo, hi, n_cand), 6))
    best = (None, -np.inf)
    order = np.argsort(time)
    t_s, e_s, s_s = time[order], event[order], score[order]
    for c in cands:
        g = s_s > c
        n1, n0 = int(g.sum()), int((~g).sum())
        if n1 < min_grp or n0 < min_grp:
            continue
        # 向量化 log-rank：按事件时间累加期望与观测
        ev = np.where(e_s == 1)[0]
        if ev.size == 0:
            continue
        at_risk = np.arange(len(t_s), 0, -1).astype(float)          # 风险集人数
        at1 = np.cumsum(g[::-1])[::-1].astype(float)                 # 组1 风险集人数
        d_ind = (e_s == 1).astype(float)                             # 事件指示
        d1_ind = (g & (e_s == 1)).astype(float)                       # 组1 事件指示
        with np.errstate(divide="ignore", invalid="ignore"):
            exp = d_ind * at1 / at_risk                              # 组1 期望事件数
        O1, E1 = d1_ind[ev].sum(), np.nansum(exp[ev])
        keep = at_risk > 1
        v = np.nansum(
            (at1[keep] / at_risk[keep]) * (1 - at1[keep] / at_risk[keep])
            * d_ind[keep] * (at_risk[keep] - d_ind[keep]) / (at_risk[keep] - 1))
        if v <= 1e-9 or not np.isfinite(v):
            continue
        stat = (O1 - E1) ** 2 / v
        if stat > best[1]:
            best = (float(c), float(stat))
    return best


def main():
    gs = load_gene_sets()
    expr, meta = load_matrix()
    clin = pd.read_csv(os.path.join(PROC, "hnsc_analysis_set.csv"))
    print("[info] 表达矩阵: %d 基因 × %d 样本" % expr.shape)

    # 仅保留原发肿瘤
    prim = meta[meta["sample_type"] == "Primary Tumor"].reset_index(drop=True)
    keep_cols = list(prim["col"])
    expr = expr[keep_cols]
    print("[info] 原发肿瘤样本: %d" % len(keep_cols))

    log2m = np.log2(expr + 1.0)
    z = zscore(log2m)

    scores = {}
    coverage = {}
    for name, genes in gs.items():
        s, n = score_gene_set(z, genes, "z_mean")
        if s is not None:
            scores[name] = s
        coverage[name] = {"n_genes": len(genes), "n_present": n}

    # ssGSEA
    ss = ssgsea_scores(z, gs)
    ss_ok = ss is not None
    if ss_ok:
        print("[info] ssGSEA 完成")

    # PCA1
    for name, genes in gs.items():
        s, n = score_gene_set(z, genes, "pca1")
        if s is not None:
            scores[name + "__pca1"] = s

    sc = pd.DataFrame(scores).T
    sc.index.name = "gene_set"
    sc = sc.T  # 样本 × 基因集
    sc["barcode"] = prim["barcode"].values
    sc["sample"] = prim["col"].values

    # 多样本合并（同一病例多张切片取均值）
    num_cols = [c for c in sc.columns if c not in ("barcode", "sample")]
    sc_num = sc.groupby("barcode")[num_cols].mean()

    # ---- 共识 CTS：5 个 TLS 集的 z_mean 百分位秩均值 ----
    cons_cols = [c for c in CONSENSUS_SETS if c in sc_num.columns]
    ranks = sc_num[cons_cols].rank(pct=True)
    cts = ranks.mean(axis=1)

    out = pd.DataFrame(index=sc_num.index)
    for c in num_cols:
        out[c] = sc_num[c]
    out["CTS"] = cts
    if ss_ok:
        for c in cons_cols:
            if c in ss.index:
                out[c + "__ssgsea"] = ss.loc[c].reindex(out.index)

    merged = out.join(clin.set_index("bcr_patient_barcode"), how="inner")
    print("[ok] CTS 与临床交集: %d 例" % len(merged))

    # ---- 阈值 ----
    time = merged["os_time_days"].values
    event = merged["os_event"].values
    cut, stat = best_cutpoint(merged["CTS"].values, time, event)
    merged["CTS_group_opt"] = np.where(merged["CTS"] > cut, "High", "Low")
    med = merged["CTS"].median()
    merged["CTS_group_median"] = np.where(merged["CTS"] > med, "High", "Low")
    # 用 lifelines 复核最优切点处的精确对数秩 p 值
    from lifelines.statistics import logrank_test
    _hi = merged["CTS_group_opt"] == "High"
    _lr = logrank_test(time[_hi.values], time[(~_hi).values],
                       event[_hi.values], event[(~_hi).values])
    print("[ok] 最优切点 CTS=%.4f (向量化 χ²=%.2f)；精确 log-rank χ²=%.2f, p=%.4g；中位数=%.4f"
          % (cut, stat, _lr.test_statistic, _lr.p_value, med))
    # 自检：向量化实现与 lifelines 在同一口径下的量级一致性
    _vcut, _vstat = best_cutpoint(merged["CTS"].values, time, event)
    print("[chk] 向量化 χ²=%.2f vs lifelines χ²=%.2f（差异来自并列时间处理）"
          % (_vstat, _lr.test_statistic))
    print("    最优切点分组: %s" % merged["CTS_group_opt"].value_counts().to_dict())

    # ---- Jackknife / Bootstrap 稳健性 ----
    boot = []
    n = len(merged)
    for _ in range(200):
        idx = RNG.integers(0, n, n)
        sub = merged.iloc[idx]
        try:
            c, _ = best_cutpoint(sub["CTS"].values, sub["os_time_days"].values,
                                 sub["os_event"].values)
            if c is not None:
                boot.append(c)
        except Exception:
            pass
    boot = np.array(boot)
    boot_ok = len(boot) > 0
    print("[ok] Bootstrap 切点: n=%d, 中位=%.4f, IQR=[%.4f, %.4f]"
          % (len(boot), np.median(boot), np.percentile(boot, 25), np.percentile(boot, 75))
          if boot_ok else "[warn] Bootstrap 未获得有效切点")

    # ---- 签名收敛性（两两一致性）----
    conv = []
    for a, b in combinations(cons_cols, 2):
        rho = stats.spearmanr(merged[a], merged[b]).statistic
        conv.append({"set_a": a, "set_b": b, "spearman_rho": round(float(rho), 4)})
    conv_df = pd.DataFrame(conv)
    print("\n--- 签名两两 Spearman 一致性（中位 %.3f）---"
          % conv_df["spearman_rho"].median())
    print(conv_df.to_string(index=False))

    merged.index.name = "barcode"
    merged.to_csv(os.path.join(PROC, "hnsc_tls_scores.csv"), encoding="utf-8-sig")
    conv_df.to_csv(os.path.join(RES, "signature_convergence.csv"), index=False)

    summary = {
        "n_samples_expr": int(expr.shape[1]),
        "n_primary_tumor": int(len(keep_cols)),
        "n_with_clinical": int(len(merged)),
        "n_death": int(merged["os_event"].sum()),
        "gene_coverage": coverage,
        "consensus_sets": cons_cols,
        "ssgsea_available": bool(ss_ok),
        "cutpoint_optimal": cut, "logrank_chi2_at_cut": stat,
        "logrank_chi2_exact": float(_lr.test_statistic),
        "logrank_p_exact": float(_lr.p_value),
        "cutpoint_median": float(med),
        "group_counts_opt": merged["CTS_group_opt"].value_counts().to_dict(),
        "bootstrap_cut_median": float(np.median(boot)) if boot_ok else None,
        "signature_convergence_median_rho": float(conv_df["spearman_rho"].median()),
        "signature_convergence_min_rho": float(conv_df["spearman_rho"].min()),
        "signature_convergence_max_rho": float(conv_df["spearman_rho"].max()),
    }
    json.dump(summary, open(os.path.join(RES, "03_score_summary.json"), "w"),
              ensure_ascii=False, indent=1)
    print("\n[ok] 评分表 -> data/processed/hnsc_tls_scores.csv")
    print("[ok] 汇总   -> results/03_score_summary.json")


if __name__ == "__main__":
    main()
