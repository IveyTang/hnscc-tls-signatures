# -*- coding: utf-8 -*-
"""
选题四 / Step 1.3 —— 并发下载 TCGA-HNSC STAR-Counts 表达数据

说明：
  * GDC /data/<file_id> 返回的是**纯文本 TSV**（非 gzip），首行为 "# gene-model: GENCODE v36"。
  * 仅保留本研究需要的基因（TLS 签名 + 免疫微环境 + 管家基因），避免落盘 2.3 GB 全量矩阵。

用法：
  python scripts/02_download_expression.py                # 全量
  python scripts/02_download_expression.py --limit 20      # 验证链路
  python scripts/02_download_expression.py --workers 8     # 并发数

输出：
  data/raw/expr/<file_id>.tsv              原始文件（可删除）
  data/processed/hnsc_gene_matrix.csv      基因 × 样本 TPM 子矩阵
  data/processed/expr_download_meta.json   下载元信息
"""
import argparse
import csv
import json
import os
import re
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

# 免疫/基质相关基因的通用命名模式（覆盖 LM22 / xCell 免疫标记基因的大部分）
IMMUNE_PATTERN = re.compile(
    r"^(CD\d|CD1[0-9A-Z]|CD2[0-9A-Z]|CD3[0-9A-Z]|CD4[0-9A-Z]|CD8[0-9A-Z]|"
    r"CCL\d|CXCL\d|CXCR\d|CCR\d|XCL\d|IL\d|IL\d+R|IL\d+[A-Z]|"
    r"HLA-|IGH|IGK|IGL|IGJ|TNF|IFN|IFI|ISG|GZ[ABKM]|PRF1|NKG7|"
    r"KIR\d|TLR\d|OAS\d|MX\d|STAT\d|IRF\d|GBP\d|PSMB\d|TAP\d|"
    r"MS4A\d|PTPRC|LAMP3|AICDA|BCL6|MEF2B|RGS1\d|CR2|LT[AB]\d?|"
    r"FOXP3|CTLA4|LAG3|HAVCR2|TIGIT|PDCD1|VSIR|TOX|TCF7|"
    r"ACTA2|COL\d|FAP|PDGFR|TAGLN|VIM|THY1|PDPN|CLDN5|"
    r"KRT\d|EPCAM|MKI67|TOP2A|MZB1|JCHAIN|XBP1|PRDM1|SDC1)")


def _is_immune(gene):
    return bool(IMMUNE_PATTERN.match(gene))

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw")
PROC = os.path.join(ROOT, "data", "processed")
EXPR = os.path.join(RAW, "expr")
os.makedirs(EXPR, exist_ok=True)
os.makedirs(PROC, exist_ok=True)
API = "https://api.gdc.cancer.gov/data/"


def load_genes():
    cfg = json.load(open(os.path.join(ROOT, "config", "tls_genesets.json"), encoding="utf-8"))
    genes = set()
    for v in cfg.values():
        if isinstance(v, dict) and "genes" in v:
            genes.update(v["genes"])
    genes.update([
        # 免疫检查点
        "CD274", "PDCD1", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "IDO1",
        # B 细胞 / 浆细胞
        "MS4A1", "CD79A", "CD79B", "CD19", "IGHG1", "IGKC", "MZB1", "JCHAIN", "PRDM1", "XBP1",
        # T 细胞
        "CD3D", "CD3E", "CD4", "CD8A", "CD8B", "FOXP3", "IL2RA", "TRAC",
        # 基质 / CAF
        "ACTA2", "FAP", "COL1A1", "COL1A2", "PDGFRB", "VIM", "TAGLN",
        # 上皮 / 对照
        "KRT5", "KRT14", "EPCAM", "KRT13", "KRT4", "MKI67",
        # 管家
        "GAPDH", "ACTB", "B2M", "TAP1", "TBP",
    ])
    # 再以命名模式扩充免疫/基质基因（保证微环境分析可用）
    return sorted(genes), _is_immune


def download_one(file_id, timeout=300):
    path = os.path.join(EXPR, file_id + ".tsv")
    if os.path.exists(path) and os.path.getsize(path) > 10000:
        return path, True
    req = urllib.request.Request(API + file_id)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        data = r.read()
    with open(path, "wb") as f:
        f.write(data)
    return path, False


def parse_star_counts(path, keep, predicate=None):
    """STAR-Counts 纯文本：gene_id / gene_name / gene_type / unstranded / ... / tpm_unstranded / ..."""
    vals = {}
    header = None
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
            if g in keep or (predicate is not None and predicate(g)):
                try:
                    vals[g] = float(p[i_tpm])
                except ValueError:
                    pass
    return vals


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    keep, predicate = load_genes()
    man = list(csv.DictReader(
        open(os.path.join(RAW, "tcga_hnsc_rnaseq_manifest.tsv"), encoding="utf-8"),
        delimiter="\t"))
    if args.limit:
        man = man[:args.limit]

    print("[info] 文件 %d 个；显式基因 %d 个 + 免疫命名模式；并发 %d"
          % (len(man), len(keep), args.workers))
    sys.stdout.flush()

    t0 = time.time()
    results, failed = [], []
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(download_one, r["file_id"]): r for r in man}
        for i, fu in enumerate(as_completed(futs), 1):
            r = futs[fu]
            fid = r["file_id"]
            try:
                p, cached = fu.result()
                sample = "%s|%s" % (r.get("sample_submitter_id", ""), r.get("sample_type", ""))
                results.append((sample, p))
            except Exception as e:
                failed.append((fid, str(e)))
            if i % 50 == 0 or i == len(man):
                el = time.time() - t0
                print("  [%d/%d] %.1f min, est.total %.1f min, 失败 %d"
                      % (i, len(man), el / 60, el / 60 * len(man) / i, len(failed)))
                sys.stdout.flush()

    print("[info] 解析表达值 ...")
    cols, mat = [], {}
    for sample, p in results:
        vals = parse_star_counts(p, keep, predicate)
        if not vals:
            continue
        cols.append(sample)
        for g, v in vals.items():
            mat.setdefault(g, {})[sample] = v

    out = os.path.join(PROC, "hnsc_gene_matrix.csv")
    with open(out, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene"] + cols)
        for g in sorted(mat):
            w.writerow([g] + [("%.4f" % mat[g].get(c, 0.0)) for c in cols])

    meta = {"n_files_ok": len(cols), "n_genes": len(mat),
            "n_failed": len(failed), "elapsed_min": round((time.time() - t0) / 60, 1),
            "genes": sorted(mat.keys()), "failed": failed[:20]}
    json.dump(meta, open(os.path.join(PROC, "expr_download_meta.json"), "w"), indent=1)
    print("\n[ok] 矩阵 %d 基因 × %d 样本 -> %s" % (len(mat), len(cols), out))
    print("[ok] 失败 %d 个" % len(failed))


if __name__ == "__main__":
    main()
