# -*- coding: utf-8 -*-
"""
选题四 / Step 4 —— 外部验证（GEO 队列）

⚠️ 前置条件（本机 GEO 不可达，需人工下载后放入 data/raw/geo/）：
    GSE65858/  GSE65858_series_matrix.txt.gz  +  GPL*（平台注释）
    GSE159067/ GSE159067_series_matrix.txt.gz +  GPL*

用法：
    python scripts/05_external_validation.py --dataset GSE65858
    python scripts/05_external_validation.py --dataset GSE159067 --ici

设计要点（对应方案约束 C4）：
    * 跨平台一致性依赖**秩/相对评分**，绝不直接比较原始表达值
    * 外部队列**不重新调参**，cutpoint 沿用 TCGA 训练集
    * 报告 HR / C-index / 校准表现

输出：
    results/05_external_<dataset>.json
    results/fig_external_<dataset>.png
"""
import argparse
import gzip
import json
import os
import re
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.utils import concordance_index

warnings.filterwarnings("ignore")
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw", "geo")
PROC = os.path.join(ROOT, "data", "processed")
RES = os.path.join(ROOT, "results")
CONSENSUS_SETS = ["G1_Coppola_12chemokine", "G2_11chemokine_HCC",
                  "G3_HNSCC_13chemokine", "G4_TLS_maturation_GC", "G6_HEV_lymphoid_stroma"]


def load_gene_sets():
    cfg = json.load(open(os.path.join(ROOT, "config", "tls_genesets.json"), encoding="utf-8"))
    return {k: v["genes"] for k, v in cfg.items() if isinstance(v, dict) and "genes" in v}


def read_series_matrix(path):
    """解析 GEO series matrix：返回 (表达矩阵 gene×sample, 样本注释, 特征行)"""
    opener = gzip.open if path.endswith(".gz") else open
    expr_lines, meta_rows, feat = [], {}, []
    with opener(path, "rt", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("!Sample_"):
                p = [x.strip('"') for x in line.split("\t")]
                meta_rows.setdefault(p[0], []).extend(p[1:])
            elif line.startswith("!series_matrix_table_begin"):
                for l2 in fh:
                    l2 = l2.rstrip("\n")
                    if l2.startswith("!series_matrix_table_end"):
                        break
                    expr_lines.append(l2)
    if not expr_lines:
        raise ValueError("未找到表达数据表")
    hdr = [x.strip('"') for x in expr_lines[0].split("\t")]
    data = [l.split("\t") for l in expr_lines[1:] if l.strip()]
    feats = [d[0].strip('"') for d in data]
    vals = np.array([[float(x) if x not in ("", "NA", "null") else np.nan for x in d[1:]]
                     for d in data])
    expr = pd.DataFrame(vals, index=feats, columns=hdr[1:])
    return expr, meta_rows, hdr[1:]


def map_probes(expr, annot_path=None):
    """探针 -> 基因符号。若提供 GPL 注释文件则用之，否则退回启发式（探针本身即基因符号）。"""
    if annot_path and os.path.exists(annot_path):
        mapping = {}
        opener = gzip.open if annot_path.endswith(".gz") else open
        with opener(annot_path, "rt", encoding="utf-8", errors="ignore") as fh:
            header, idx_id, idx_sym = None, None, None
            for line in fh:
                if line.startswith(("^", "!", "#")):
                    continue
                p = line.rstrip("\n").split("\t")
                if header is None:
                    header = p
                    for cand in ("ID", "ID_REF", "Probe Set ID", "SPOT_ID"):
                        if cand in header:
                            idx_id = header.index(cand)
                            break
                    for cand in ("Gene symbol", "GENE_SYMBOL", "Symbol",
                                 "gene_assignment", "Gene Symbol"):
                        if cand in header:
                            idx_sym = header.index(cand)
                            break
                    if idx_id is None or idx_sym is None:
                        return None
                    continue
                if len(p) <= max(idx_id, idx_sym):
                    continue
                sym = p[idx_sym]
                m = re.search(r"([A-Za-z0-9\-\.]+)\s*//", sym) or re.match(r"([A-Z0-9\-\.]+)", sym)
                if m:
                    mapping[p[idx_id]] = m.group(1).strip()
        if mapping:
            e = expr.copy()
            e["__sym"] = [mapping.get(i, None) for i in e.index]
            e = e.dropna(subset=["__sym"])
            e = e.groupby("__sym").mean()
            return e
    return None


def score_ranks(expr_log2, genes):
    """平台无关评分：基因表达百分位秩的均值。"""
    present = [g for g in genes if g in expr_log2.index]
    if len(present) < 3:
        return None, len(present)
    sub = expr_log2.loc[present]
    r = sub.rank(pct=True, axis=1)
    return r.mean(axis=0), len(present)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--ici", action="store_true", help="按 ICI 队列处理（疗效分组）")
    ap.add_argument("--time-col", default=None)
    ap.add_argument("--event-col", default=None)
    ap.add_argument("--cut-from", default=os.path.join(RES, "03_score_summary.json"))
    args = ap.parse_args()

    ds_dir = os.path.join(RAW, args.dataset)
    if not os.path.isdir(ds_dir):
        print("=" * 72)
        print("[阻塞] 未找到 %s" % ds_dir)
        print("本机 GEO（NCBI）不可达（实测 HTTP 000）。请在其他联网环境下载：")
        print("  1. %s_series_matrix.txt.gz" % args.dataset)
        print("  2. 平台注释文件（GPL*_family.soft.gz 或 GPL*.annot.gz）")
        print("放入：%s" % ds_dir)
        print("=" * 72)
        sys.exit(2)

    files = os.listdir(ds_dir)
    smat = [f for f in files if f.endswith("series_matrix.txt.gz") or f.endswith("series_matrix.txt")]
    if not smat:
        print("[错误] %s 目录下未找到 series matrix 文件" % ds_dir)
        sys.exit(2)
    annot = [f for f in files if f.upper().startswith("GPL")]
    annot_path = os.path.join(ds_dir, annot[0]) if annot else None

    print("[info] 读取 %s" % smat[0])
    expr, meta, samples = read_series_matrix(os.path.join(ds_dir, smat[0]))
    print("[info] 原始矩阵: %d 探针 × %d 样本" % expr.shape)

    mapped = map_probes(expr, annot_path)
    if mapped is None:
        print("[warn] 未找到平台注释，按启发式将行索引当作基因符号")
        mapped = expr
    print("[info] 映射后: %d 基因" % mapped.shape[0])

    log2m = np.log2(mapped.clip(lower=0) + 1.0)
    gs = load_gene_sets()
    scores, cov = {}, {}
    for name, genes in gs.items():
        s, n = score_ranks(log2m, genes)
        cov[name] = {"n_genes": len(genes), "n_present": n}
        if s is not None:
            scores[name] = s
    sc = pd.DataFrame(scores)
    cons = [c for c in CONSENSUS_SETS if c in sc.columns]
    sc["CTS"] = sc[cons].rank(pct=True).mean(axis=1)

    # 沿用训练集切点转换为百分位阈值
    cut = None
    if os.path.exists(args.cut_from):
        cut = json.load(open(args.cut_from)).get("cutpoint_optimal")
    if cut is not None:
        sc["CTS_group"] = np.where(sc["CTS"] > cut, "High", "Low")
    else:
        sc["CTS_group"] = np.where(sc["CTS"] > sc["CTS"].median(), "High", "Low")

    out = {"dataset": args.dataset, "n_samples": int(sc.shape[0]),
           "gene_coverage": cov, "cutpoint_used": cut}

    # 生存验证（若该 GEO 队列提供生存字段，需在 --time-col/--event-col 指定并另接临床表）
    clin_path = os.path.join(ds_dir, "clinical.csv")
    if os.path.exists(clin_path) and args.time_col and args.event_col:
        cl = pd.read_csv(clin_path)
        key = cl.columns[0]
        d = sc.join(cl.set_index(key), how="inner")
        d = d[d[args.time_col].notna() & d[args.event_col].notna()]
        hi, lo = d[d["CTS_group"] == "High"], d[d["CTS_group"] == "Low"]
        lr = logrank_test(hi[args.time_col], lo[args.time_col],
                          hi[args.event_col].astype(int), lo[args.event_col].astype(int))
        cidx = concordance_index(d[args.time_col], -d["CTS"], d[args.event_col].astype(int))
        out["logrank_p"] = float(lr.p_value)
        out["cindex"] = float(cidx)
        print("[ok] 外部验证 log-rank p=%.4g, C-index=%.3f" % (lr.p_value, cidx))

        fig, ax = plt.subplots(figsize=(5.2, 4.2), dpi=300)
        for grp, clr in [("High", "#A32D2D"), ("Low", "#0F6E56")]:
            sub = d[d["CTS_group"] == grp]
            KaplanMeierFitter().fit(sub[args.time_col], sub[args.event_col].astype(int),
                                    label="CTS %s (n=%d)" % (grp, len(sub))
                                    ).plot_survival_function(ax=ax, color=clr, ci_show=True)
        ax.set_title("%s：外部验证" % args.dataset, fontsize=10)
        ax.set_xlabel("时间（月）")
        ax.set_ylabel("生存率")
        ax.grid(alpha=0.25, linewidth=0.5)
        fig.tight_layout()
        fig.savefig(os.path.join(RES, "fig_external_%s.png" % args.dataset))
        plt.close(fig)
    else:
        print("[info] 未提供临床/生存字段，仅输出 CTS 分组分布")
        out["group_counts"] = sc["CTS_group"].value_counts().to_dict()

    sc.to_csv(os.path.join(PROC, "ext_%s_scores.csv" % args.dataset), encoding="utf-8-sig")
    json.dump(out, open(os.path.join(RES, "05_external_%s.json" % args.dataset), "w"),
              ensure_ascii=False, indent=1, default=str)
    print("[ok] 结果写入 results/05_external_%s.json" % args.dataset)


if __name__ == "__main__":
    main()
