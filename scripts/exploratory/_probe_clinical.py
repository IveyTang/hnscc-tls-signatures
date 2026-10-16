# -*- coding: utf-8 -*-
"""探查 GDC 的 Clinical Supplement / Pathology Report 是否含 HPV 或 p16 信息。"""
import json
import urllib.request

GDC = "https://api.gdc.cancer.gov"


def post(path, payload):
    req = urllib.request.Request(
        GDC + path, data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"})
    return json.loads(urllib.request.urlopen(req, timeout=90).read().decode("utf-8"))


def get(path, params=None):
    url = GDC + path
    if params:
        url += "?" + urllib.parse.urlencode(params)
    return json.loads(urllib.request.urlopen(url, timeout=90).read().decode("utf-8"))


import urllib.parse  # noqa

print("=" * 70)
print("1. clinical 类别下的文件类型")
r = post("/files", {
    "filters": {"op": "and", "content": [
        {"op": "in", "content": {"field": "cases.project.project_id",
                                 "value": ["TCGA-HNSC"]}},
        {"op": "in", "content": {"field": "data_category", "value": ["clinical"]}}]},
    "facets": "data_type,data_format,access", "size": 0})
for k, v in r["data"]["aggregations"].items():
    print("--", k)
    for b in v["buckets"]:
        print("     %-42s %s" % (b["key"], b["doc_count"]))

print()
print("=" * 70)
print("2. 取一份 Clinical Supplement 看内容")
r = post("/files", {
    "filters": {"op": "and", "content": [
        {"op": "in", "content": {"field": "cases.project.project_id",
                                 "value": ["TCGA-HNSC"]}},
        {"op": "in", "content": {"field": "data_type",
                                 "value": ["Clinical Supplement"]}}]},
    "fields": "file_id,file_name,file_size,access,data_format",
    "size": 3, "format": "JSON"})
for h in r["data"]["hits"]:
    print("  %s | %s | %.1f KB | %s | %s" % (
        h["file_id"], h["file_name"], h["file_size"] / 1024,
        h.get("data_format"), h.get("access")))

print()
print("=" * 70)
print("3. Pathology Report 样例")
r = post("/files", {
    "filters": {"op": "and", "content": [
        {"op": "in", "content": {"field": "cases.project.project_id",
                                 "value": ["TCGA-HNSC"]}},
        {"op": "in", "content": {"field": "data_type",
                                 "value": ["Pathology Report"]}}]},
    "fields": "file_id,file_name,file_size,access,data_format",
    "size": 3, "format": "JSON"})
for h in r["data"]["hits"]:
    print("  %s | %s | %.1f KB | %s | %s" % (
        h["file_id"], h["file_name"], h["file_size"] / 1024,
        h.get("data_format"), h.get("access")))

print()
print("=" * 70)
print("4. 是否有 HPV 相关的实验/注释类数据")
for dt in ["Viral Genotyping", "HPV", "Viral Expression", "Viral Integration"]:
    try:
        r = post("/files", {
            "filters": {"op": "and", "content": [
                {"op": "in", "content": {"field": "cases.project.project_id",
                                         "value": ["TCGA-HNSC"]}},
                {"op": "in", "content": {"field": "data_type", "value": [dt]}}]},
            "size": 0})
        print("  %-24s -> %s" % (dt, r["data"]["pagination"]["total"]))
    except Exception as e:
        print("  %-24s -> ERR %s" % (dt, e))
