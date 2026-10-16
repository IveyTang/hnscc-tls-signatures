# -*- coding: utf-8 -*-
"""从已下载的原始 STAR-Counts TSV 中补抽 HPV 判别相关基因。

首轮下载只保留了 1291 个免疫微环境相关基因，缺少 CDKN2A(p16) 等
HPV 判别关键基因。本脚本扫描 data/raw/expr/*.tsv，抽取增补基因面板，
生成 data/processed/hnsc_gene_matrix_hpv.csv（避免重复下载 2.3 GB）。

增补面板构成（用于 HPV 转录组判定）：
  A. RB/p16 通路与细胞周期    —— HPV E7 灭活 RB → p16 上调、增殖基因上调
  B. 突变谱相关（HPV+ 特征）  —— PIK3CA/FGFR3 富集，TP53 野生
  C. 鳞状分化与 EGFR 轴       —— HPV- 高表达 CCND1/EGFR/SOX2/AHR
  D. 免疫/干扰素轴            —— 供分类器交叉验证用
"""
import csv
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw")
EXPR = os.path.join(RAW, "expr")
PROC = os.path.join(ROOT, "data", "processed")

PANEL = [
    # A. RB / p16 / 细胞周期 / DNA 复制
    "CDKN2A", "CDKN2B", "CDKN2C", "RB1", "RBL1", "RBL2", "CCND1", "CCND2", "CCND3",
    "CDK4", "CDK6", "CDK2", "E2F1", "E2F2", "E2F3", "E2F4", "MCM2", "MCM3", "MCM4",
    "MCM5", "MCM6", "MCM7", "PCNA", "LIG1", "RPA1", "RPA2", "RPA3", "TYMS", "RFC4",
    "RFC2", "POLE", "PLK1", "TOP2A", "AURKA", "AURKB", "BUB1", "CCNE1", "CDC6",
    "CDC20", "CDKN1A", "CDKN1B", "ORC1", "POLD1", "GINS2", "ASF1B",
    # B. 体细胞突变谱（HPV+ 特征性改变）
    "TP53", "PIK3CA", "PTEN", "NOTCH1", "NOTCH2", "FAT1", "CASP8", "HRAS", "NRAS",
    "NFE2L2", "KEAP1", "FBXW7", "TRAF3", "CYLD", "BRCA2", "ATM", "KMT2D", "NSD1",
    # C. 鳞状分化 / EGFR 轴（HPV- 高表达）
    "SOX2", "TP63", "KRT5", "KRT14", "KRT6A", "KRT16", "KRT17", "AHR", "EGFR",
    "FGFR1", "FGFR3", "LAMA3", "LAMB3", "MMP1", "MMP10", "WNT5A", "ALDH3A1",
    "CYP1B1", "CYP4F11", "G6PD", "TP63", "CRNN",
    # D. 免疫 / 干扰素（用于交叉验证分类器生物学合理性）
    "IFNG", "GZMA", "GZMB", "PRF1", "CXCL9", "CXCL10", "CXCL11", "STAT1", "IRF1",
    "HLA-DRA", "HLA-A", "B2M", "CD8A", "CD8B", "FOXP3", "CD274", "PDCD1", "CTLA4",
    "MZB1", "CXCL13", "LAG3", "TIGIT", "IDO1",
]


def parse(path, keep):
    vals = {}
    header = None
    i_name = i_tpm = None
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            if header is None:
                header = line.rstrip("\n").split("\t")
                try:
                    i_name = header.index("gene_name")
                    i_tpm = header.index("tpm_unstranded")
                except ValueError:
                    return vals
                continue
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) <= i_tpm:
                continue
            g = p[i_name]
            if g in keep:
                try:
                    vals[g] = float(p[i_tpm])
                except ValueError:
                    pass
    return vals


def main():
    keep = set(PANEL)
    files = sorted(f for f in os.listdir(EXPR) if f.endswith(".tsv"))
    print("[info] 原始文件 %d 个；目标基因 %d 个" % (len(files), len(keep)))
    sys.stdout.flush()

    results = {}
    fails = []
    with ThreadPoolExecutor(max_workers=10) as ex:
        futs = {ex.submit(parse, os.path.join(EXPR, f), keep): f for f in files}
        done = 0
        for fu in as_completed(futs):
            f = futs[fu]
            done += 1
            try:
                v = fu.result()
                if v:
                    results[f] = v
            except Exception:
                fails.append(f)
            if done % 150 == 0:
                print("   ... %d/%d" % (done, len(files)))
                sys.stdout.flush()

    genes = sorted({g for v in results.values() for g in v})
    print("[ok] 成功解析 %d 个文件，命中基因 %d 个" % (len(results), len(genes)))
    found = set(genes)
    missing = [g for g in PANEL if g not in found]
    print("     未命中: %s" % (", ".join(missing) if missing else "无"))

    # 列名从 file_id 映射为「样本barcode|样本类型」，与主矩阵 hnsc_gene_matrix.csv 口径一致
    import csv as _csv
    f2s = {}
    with open(os.path.join(RAW, "tcga_hnsc_rnaseq_manifest.tsv"), encoding="utf-8") as fh:
        for row in _csv.DictReader(fh, delimiter="\t"):
            f2s[row["file_id"]] = "%s|%s" % (row["sample_submitter_id"],
                                             row["sample_type"])
    cols = []
    for f in results:
        fid = f[:-4] if f.endswith(".tsv") else f
        cols.append(f2s.get(fid, fid))
    print("     列名映射: %d/%d 命中 manifest" % (sum(1 for c in cols if "|" in c),
                                              len(cols)))

    outp = os.path.join(PROC, "hnsc_gene_matrix_hpv.csv")
    with open(outp, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["gene"] + cols)
        for g in genes:
            w.writerow([g] + [("%.6g" % results[f][g]) if g in results[f] else ""
                              for f in results])
    print("     输出 -> %s" % outp)
    if fails:
        print("     失败 %d 个" % len(fails))


if __name__ == "__main__":
    main()
