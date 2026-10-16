# -*- coding: utf-8 -*-
"""下载 TCGA-HNSC 全部 BCR XML 并提取 HPV/p16 状态 + BCR 层面补充字段。

背景：GDC 结构化 API 的 diagnoses.hpv_positive_type 在 TCGA-HNSC 全为空，
但 BCR XML（data_format=bcr xml）中保留了 TSS 上报的原始字段：
  hnsc:hpv_status_by_p16_testing  (CDE 3108263, preferred_name=hpv_status_p16)
  hnsc:hpv_status_by_ish_testing  (CDE 3108261, preferred_name=hpv_status_ish)
本脚本用 ElementTree 解析（避免自闭合标签误匹配），输出单表 CSV。
"""
import csv
import json
import os
import sys
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw")
OUT = os.path.join(RAW, "bcr_xml")
os.makedirs(OUT, exist_ok=True)
GDC = "https://api.gdc.cancer.gov"

# 目标字段：BCR 本地名 -> 输出列名
FIELDS = {
    "hpv_status_by_p16_testing": "hpv_p16",
    "hpv_status_by_ish_testing": "hpv_ish",
    "anatomic_neoplasm_subdivision": "anatomic_subdivision",
    "icd_o_3_site": "icd_o3_site",
    "icd_o_3_histology": "icd_o3_histology",
    "tumor_tissue_site": "tumor_tissue_site",
    "number_pack_years_smoked": "pack_years",
    "tobacco_smoking_history": "tobacco_smoking_history",
    "alcohol_history_documented": "alcohol_history",
    "margin_status": "margin_status",
    "lymphovascular_invasion_present": "lvi",
    "perineural_invasion_present": "pni",
    "lymph_node_examined_count": "ln_examined",
    "number_of_lymphnodes_positive_by_he": "ln_positive",
    "presence_of_pathological_nodal_extracapsular_spread": "ene",
    "neoplasm_histologic_grade": "histologic_grade",
    "person_neoplasm_cancer_status": "bcr_cancer_status",
    "method_of_curative_tx": "curative_tx",
    "history_other_malignancy": "other_malignancy",
}


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
            {"op": "in", "content": {"field": "data_format", "value": ["bcr xml"]}},
            {"op": "in", "content": {"field": "data_type",
                                     "value": ["Clinical Supplement"]}}]},
        "fields": "file_id,file_name,cases.submitter_id",
        "size": 2000, "format": "JSON"})
    return [{"file_id": h["file_id"], "file_name": h["file_name"],
             "submitter_id": (h.get("cases") or [{}])[0].get("submitter_id")}
            for h in r["data"]["hits"]]


def dl(item):
    p = os.path.join(OUT, item["file_id"] + ".xml")
    if os.path.exists(p) and os.path.getsize(p) > 1000:
        return item, p, None
    err = "no attempt"
    for _ in range(3):
        try:
            d = urllib.request.urlopen(GDC + "/data/" + item["file_id"],
                                       timeout=120).read()
            with open(p, "wb") as fh:
                fh.write(d)
            return item, p, None
        except Exception as e:
            err = str(e)
    return item, None, err


def strip_ns(tag):
    return tag.split("}", 1)[1] if "}" in tag else tag


def parse(path, submitter_id):
    """用 ElementTree 解析；同名标签只取第一个非空值。"""
    rec = {"submitter_id": submitter_id}
    try:
        root = ET.parse(path).getroot()
    except Exception as e:
        rec["parse_error"] = str(e)
        return rec
    want = set(FIELDS)
    for el in root.iter():
        local = strip_ns(el.tag)
        if local in want and FIELDS[local] not in rec:
            v = (el.text or "").strip()
            if v:
                rec[FIELDS[local]] = v
        if local == "bcr_patient_barcode" and "barcode" not in rec:
            v = (el.text or "").strip()
            if v:
                rec["barcode"] = v
    for k in FIELDS.values():
        rec.setdefault(k, None)
    return rec


def main():
    man = fetch_manifest()
    print("[info] BCR XML 文件 %d 个" % len(man))
    sys.stdout.flush()

    rows, fails = [], []
    with ThreadPoolExecutor(max_workers=12) as ex:
        futs = {ex.submit(dl, it): it for it in man}
        done = 0
        for fu in as_completed(futs):
            it, p, err = fu.result()
            done += 1
            if p:
                rows.append(parse(p, it["submitter_id"]))
            else:
                fails.append((it["file_id"], err))
            if done % 150 == 0:
                print("   ... %d/%d" % (done, len(man)))
                sys.stdout.flush()

    cols = ["submitter_id", "barcode"] + list(FIELDS.values()) + ["parse_error"]
    outp = os.path.join(RAW, "tcga_hnsc_bcr_hpv.csv")
    with open(outp, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in sorted(rows, key=lambda r: r.get("submitter_id") or ""):
            w.writerow(r)

    print("\n[ok] 解析 %d 例，失败 %d 例 -> %s" % (len(rows), len(fails), outp))
    print("\n--- 关键字段覆盖度 ---")
    for k in ["hpv_p16", "hpv_ish", "anatomic_subdivision", "pack_years",
              "tobacco_smoking_history", "alcohol_history", "margin_status",
              "lvi", "pni", "ln_examined", "ln_positive", "ene",
              "histologic_grade", "bcr_cancer_status", "curative_tx",
              "other_malignancy"]:
        vals = [r.get(k) for r in rows if r.get(k)]
        c = Counter(vals)
        top = ", ".join("%s=%d" % (a, b) for a, b in c.most_common(6))
        print("  %-26s %3d/%d  %s" % (k, len(vals), len(rows), top))

    print("\n--- p16 状态 × 解剖部位（交叉） ---")
    cross = Counter()
    for r in rows:
        p = r.get("hpv_p16")
        if p in ("Positive", "Negative"):
            cross[(p, r.get("anatomic_subdivision"))] += 1
    for (p, site), n in sorted(cross.items(), key=lambda x: -x[1])[:20]:
        print("  %-9s %-26s %d" % (p, site, n))

    if fails:
        print("\n[warn] 失败文件 %d 个" % len(fails))


if __name__ == "__main__":
    main()
