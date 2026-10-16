# -*- coding: utf-8 -*-
"""
选题四 / Step 1.1 —— 通过 GDC API 获取 TCGA-HNSC 临床信息与 RNA-seq 文件清单

本机实测：api.gdc.cancer.gov 可达（200，~2.1s）；TCGA-HNSC 共 528 cases、
1132 个 STAR-Counts 转录组文件。

输出：
  data/raw/tcga_hnsc_clinical.json    原始临床 JSON
  data/raw/tcga_hnsc_clinical.tsv     扁平化临床表
  data/raw/tcga_hnsc_rnaseq_manifest.tsv  转录组文件清单（可直接喂给 gdc-client）
"""
import json
import os
import sys
import urllib.request

API = "https://api.gdc.cancer.gov"
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw")
os.makedirs(RAW, exist_ok=True)


def gdc_post(endpoint, payload, timeout=120):
    req = urllib.request.Request(
        API + endpoint,
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode("utf-8"))


CLINICAL_FIELDS = [
    "submitter_id",
    "case_id",
    # 人口学
    # 注意：GDC 已弃用 demographic.gender（返回空），性别字段为 demographic.sex_at_birth
    "demographic.sex_at_birth",
    "demographic.age_at_index",
    "demographic.race",
    "demographic.ethnicity",
    "demographic.vital_status",
    "demographic.days_to_death",
    # 诊断（注意：部分病例存在多个诊断，含第二原发肿瘤，需按
    #  diagnosis_is_primary_disease / 形态学编码筛选，见 01 脚本）
    "diagnoses.age_at_diagnosis",
    "diagnoses.primary_diagnosis",
    "diagnoses.morphology",
    "diagnoses.tumor_grade",
    "diagnoses.tissue_or_organ_of_origin",
    "diagnoses.site_of_resection_or_biopsy",
    "diagnoses.icd_10_code",
    "diagnoses.diagnosis_is_primary_disease",
    "diagnoses.classification_of_tumor",
    "diagnoses.ajcc_pathologic_stage",
    "diagnoses.ajcc_clinical_stage",
    "diagnoses.ajcc_pathologic_t",
    "diagnoses.ajcc_clinical_t",
    "diagnoses.ajcc_pathologic_n",
    "diagnoses.ajcc_clinical_n",
    "diagnoses.ajcc_pathologic_m",
    "diagnoses.ajcc_clinical_m",
    "diagnoses.year_of_diagnosis",
    "diagnoses.days_to_last_follow_up",
    "diagnoses.vital_status",
    "diagnoses.days_to_recurrence",
    "diagnoses.progression_or_recurrence",
    "diagnoses.last_known_disease_status",
    # 暴露
    "exposures.tobacco_smoking_status",
    "exposures.cigarettes_per_day",
    "exposures.pack_years_smoked",
    "exposures.alcohol_history",
    "exposures.alcohol_intensity",
    "exposures.years_smoked",
]


def fetch_clinical():
    def try_fetch(fields):
        payload = {
            "filters": {"op": "in", "content": {
                "field": "project.project_id", "value": ["TCGA-HNSC"]}},
            "fields": ",".join(fields),
            "format": "JSON",
            "size": 2000,
            "expand": "diagnoses,demographic,exposures",
        }
        return gdc_post("/cases", payload)

    try:
        res = try_fetch(CLINICAL_FIELDS)
    except Exception as e:
        print("[warn] 完整字段请求失败（%s），回退到最小字段集" % e)
        res = try_fetch(["submitter_id", "case_id", "demographic.gender",
                         "demographic.age_at_index", "demographic.vital_status",
                         "demographic.days_to_death",
                         "diagnoses.ajcc_pathologic_stage",
                         "diagnoses.days_to_last_follow_up"])

    hits = res["data"]["hits"]
    with open(os.path.join(RAW, "tcga_hnsc_clinical.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=1)

    keys = []
    seen = set()
    for h in hits:
        for k in h:
            if k not in seen:
                seen.add(k)
                keys.append(k)
    keys.sort()
    out = os.path.join(RAW, "tcga_hnsc_clinical.tsv")
    with open(out, "w", encoding="utf-8", newline="") as f:
        f.write("\t".join(keys) + "\n")
        for h in hits:
            f.write("\t".join(str(h.get(k, "")).replace("\t", " ") for k in keys) + "\n")
    print("[ok] 临床样本 %d 例 -> %s" % (len(hits), out))
    return len(hits), keys


def fetch_manifest():
    payload = {
        "filters": {"op": "and", "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-HNSC"]}},
            {"op": "in", "content": {"field": "data_category", "value": ["Transcriptome Profiling"]}},
            {"op": "in", "content": {"field": "analysis.workflow_type", "value": ["STAR - Counts"]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
        ]},
        "fields": "file_id,file_name,md5sum,file_size,data_type,experimental_strategy,"
                  "cases.submitter_id,cases.samples.sample_type,cases.samples.submitter_id",
        "format": "JSON",
        "size": 5000,
    }
    res = gdc_post("/files", payload)
    hits = res["data"]["hits"]
    cols = ["file_id", "file_name", "md5sum", "file_size", "data_type",
            "experimental_strategy", "case_submitter_id", "sample_type", "sample_submitter_id"]
    out = os.path.join(RAW, "tcga_hnsc_rnaseq_manifest.tsv")
    with open(out, "w", encoding="utf-8", newline="") as f:
        f.write("\t".join(cols) + "\n")
        for h in hits:
            cs = (h.get("cases") or [{}])[0]
            sm = ((cs.get("samples") or [{}])[0]) if cs else {}
            f.write("\t".join([
                h.get("file_id", ""), h.get("file_name", ""), h.get("md5sum", ""),
                str(h.get("file_size", "")), h.get("data_type", ""),
                h.get("experimental_strategy", ""),
                cs.get("submitter_id", ""), sm.get("sample_type", ""),
                sm.get("submitter_id", ""),
            ]) + "\n")
    print("[ok] 转录组文件清单 %d 个 -> %s" % (len(hits), out))
    return len(hits)


if __name__ == "__main__":
    print("=== 选题四 Step 1.1：TCGA-HNSC 数据获取 ===")
    n_cases, keys = fetch_clinical()
    n_files = fetch_manifest()
    print("\n临床字段（%d）：%s" % (len(keys), ", ".join(keys)))
    print("\n下一步：将 manifest 交给 gdc-client 批量下载，"
          "或改用 GDC API 分块拉取 expression 文件。")
