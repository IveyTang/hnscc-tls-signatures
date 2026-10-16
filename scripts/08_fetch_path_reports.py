# -*- coding: utf-8 -*-
"""下载 TCGA-HNSC 全部 Pathology Report PDF 并文本挖掘 p16 / HPV 记录。

目的：BCR XML 中 p16/ISH 标注仅覆盖 115/528 例，其中落在分析集的只有 81 例，
不足以支撑按 HPV 分层的生存分析。TCGA 的病理报告 PDF 常以自由文本记录
p16 免疫组化或 HPV 检测结果，可作为标注的补充来源。

输出：data/raw/path_reports/ 下 PDF + data/raw/tcga_hnsc_path_report_hpv.csv
"""
import csv
import io
import json
import os
import re
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw")
OUT = os.path.join(RAW, "path_reports")
os.makedirs(OUT, exist_ok=True)
GDC = "https://api.gdc.cancer.gov"


def post(path, payload):
    req = urllib.request.Request(
        GDC + path, data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"})
    return json.loads(urllib.request.urlopen(req, timeout=120).read().decode("utf-8"))


def fetch_manifest():
    r = post("/files", {
        "filters": {"op": "and", "content": [
            {"op": "in", "content": {"field": "cases.project.project_id",
                                     "value": ["TCGA-HNSC"]}},
            {"op": "in", "content": {"field": "data_type",
                                     "value": ["Pathology Report"]}}]},
        "fields": "file_id,file_name,cases.submitter_id",
        "size": 2000, "format": "JSON"})
    return [{"file_id": h["file_id"],
             "file_name": h["file_name"],
             "submitter_id": (h.get("cases") or [{}])[0].get("submitter_id")}
            for h in r["data"]["hits"]]


def dl(item):
    p = os.path.join(OUT, item["file_id"] + ".pdf")
    if os.path.exists(p) and os.path.getsize(p) > 500:
        return item, p, None
    err = "no attempt"
    for _ in range(3):
        try:
            d = urllib.request.urlopen(GDC + "/data/" + item["file_id"],
                                       timeout=180).read()
            with open(p, "wb") as fh:
                fh.write(d)
            return item, p, None
        except Exception as e:
            err = str(e)
    return item, None, err


# 关键：只抓含 p16/HPV 语义的句子，避免误判
PAT_P16 = re.compile(r"\bp\s*16\b|\bp16\b|P16\b", re.I)
PAT_HPV = re.compile(r"hpv|human\s+papilloma|papillomavirus", re.I)
PAT_POS = re.compile(r"positive|pos\b|\+", re.I)
PAT_NEG = re.compile(r"negative|neg\b|not\s+detected|non-?reactive|no\s+evidence", re.I)


def classify_text(txt):
    """返回 (hpv_call, evidence) — hpv_call ∈ {Positive, Negative, None}。"""
    hits = []
    for sent in re.split(r"(?<=[.;\n])\s+", txt):
        s = sent.strip()
        if len(s) > 400:
            continue
        if PAT_HPV.search(s) or PAT_P16.search(s):
            hits.append(s)
    if not hits:
        return None, ""
    # 只在明确出现阳性/阴性判词时下结论
    pos = neg = 0
    ev = []
    for s in hits[:25]:
        # 需要同时含检测对象与判词
        if not (PAT_HPV.search(s) or PAT_P16.search(s)):
            continue
        has_pos = bool(re.search(r"\bpositive\b|\bpos\b", s, re.I))
        has_neg = bool(re.search(r"\bnegative\b|\bneg\b|not\s+detected|"
                                 r"non-?reactive|no\s+evidence", s, re.I))
        if has_pos and not has_neg:
            pos += 1
            ev.append(s)
        elif has_neg and not has_pos:
            neg += 1
            ev.append(s)
    if pos and not neg:
        return "Positive", " | ".join(ev[:3])[:600]
    if neg and not pos:
        return "Negative", " | ".join(ev[:3])[:600]
    return None, " | ".join(hits[:3])[:600]


def extract_text(path):
    try:
        from pypdf import PdfReader
        r = PdfReader(path)
        return "\n".join((p.extract_text() or "") for p in r.pages)
    except Exception as e:
        return "__ERR__" + str(e)


def main():
    only_analysis = "--analysis-only" in sys.argv
    man = fetch_manifest()
    print("[info] 病理报告 %d 份" % len(man))
    sys.stdout.flush()

    rows, fails, noext = [], [], 0
    with ThreadPoolExecutor(max_workers=10) as ex:
        futs = {ex.submit(dl, it): it for it in man}
        done = 0
        for fu in as_completed(futs):
            it, p, err = fu.result()
            done += 1
            if not p:
                fails.append(it["file_id"])
                continue
            txt = extract_text(p)
            if txt.startswith("__ERR__"):
                noext += 1
                continue
            call, ev = classify_text(txt)
            rows.append({"submitter_id": it["submitter_id"],
                         "file_id": it["file_id"],
                         "hpv_from_report": call or "",
                         "evidence": ev,
                         "text_len": len(txt)})
            if done % 100 == 0:
                print("   ... %d/%d" % (done, len(man)))
                sys.stdout.flush()

    outp = os.path.join(RAW, "tcga_hnsc_path_report_hpv.csv")
    with open(outp, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=["submitter_id", "file_id",
                                           "hpv_from_report", "text_len",
                                           "evidence"])
        w.writeheader()
        for r in sorted(rows, key=lambda r: r["submitter_id"] or ""):
            w.writerow(r)

    from collections import Counter
    c = Counter(r["hpv_from_report"] or "(未判定)" for r in rows)
    print("\n[ok] 处理 %d 份；PDF 解析失败 %d；下载失败 %d" % (len(rows), noext, len(fails)))
    print("     文本判定分布:", dict(c))
    print("     输出 ->", outp)


if __name__ == "__main__":
    main()
