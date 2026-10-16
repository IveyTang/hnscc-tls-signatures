# -*- coding: utf-8 -*-
"""探查 TCGA-HNSC 的 HPV 状态可得性：临床字段覆盖度 + 是否存在病毒类数据。"""
import json
import urllib.request

GDC = "https://api.gdc.cancer.gov"


def post(path, payload):
    req = urllib.request.Request(
        GDC + path, data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"})
    return json.loads(urllib.request.urlopen(req, timeout=90).read().decode("utf-8"))


print("=" * 70)
print("1. TCGA-HNSC 全部 data_category")
r = post("/files", {
    "filters": {"op": "in", "content": {"field": "cases.project.project_id",
                                        "value": ["TCGA-HNSC"]}},
    "facets": "data_category,data_type,experimental_strategy",
    "size": 0})
for k, v in r.get("data", {}).get("aggregations", {}).items():
    print("--", k)
    for b in v["buckets"]:
        print("     %-40s %s" % (b["key"], b["doc_count"]))

print()
print("=" * 70)
print("2. 临床 HPV 相关字段覆盖度")
CASES = ["diagnoses.hpv_positive_type", "diagnoses.primary_diagnosis",
         "diagnoses.tissue_or_organ_of_origin", "diagnoses.ajcc_pathologic_stage",
         "demographic.sex_at_birth"]
cases = post("/cases", {
    "filters": {"op": "in", "content": {"field": "project.project_id",
                                        "value": ["TCGA-HNSC"]}},
    "fields": ",".join(CASES), "format": "JSON", "size": 2000,
    "expand": "diagnoses,demographic"})["data"]["hits"]

n = len(cases)
cnt = {f: 0 for f in CASES}
vals = {}
for h in cases:
    for dx in (h.get("diagnoses") or [{}]):
        for f in CASES:
            if f.startswith("diagnoses."):
                v = dx.get(f.split(".", 1)[1])
                if v not in (None, "", "not reported", "--", "'--", "unknown"):
                    cnt[f] += 1
                    if f == "diagnoses.hpv_positive_type":
                        vals[v] = vals.get(v, 0) + 1
    dm = h.get("demographic") or {}
    for f in CASES:
        if f.startswith("demographic."):
            v = dm.get(f.split(".", 1)[1])
            if v not in (None, "", "not reported", "--", "'--"):
                cnt[f] += 1
print("总病例 %d" % n)
for f in CASES:
    print("  %-42s 非空 %4d (%.1f%%)" % (f, cnt[f], 100.0 * cnt[f] / n))
print("  hpv_positive_type 取值:", vals if vals else "(全空)")
