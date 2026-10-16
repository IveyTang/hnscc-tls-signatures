# -*- coding: utf-8 -*-
"""
选题四 / Step 1.2 —— 将 GDC 嵌套临床 JSON 扁平化为分析用临床表（含主诊断筛选）

关键修复（v2）：
  1. 【多诊断污染】部分病例含第二原发肿瘤诊断（前列腺、脑、淋巴瘤等）。
     原 v1 直接取 diagnoses[0] 会选错诊断。本版按
     `diagnosis_is_primary_disease` + 鳞癌形态学编码 + 头颈部位 三重条件挑选主诊断。
  2. 【部位字段】GDC 的 `site_of_primary_disease` 为空；正确字段是
     `tissue_or_organ_of_origin`，并进一步映射到解剖亚部位
     （口腔 / 口咽 / 喉 / 下咽 / 其他头颈 / 非头颈）。
  3. 【分期】`ajcc_pathologic_stage` 缺失较多，回退到 `ajcc_clinical_stage`；
     同时补齐 T/N/M。
  4. 【OS】死亡用 days_to_death，存活/失访用 days_to_last_follow_up。

输入：data/raw/tcga_hnsc_clinical.json
输出：data/processed/hnsc_clinical.csv     全部病例
      data/processed/hnsc_analysis_set.csv 分析集（头颈鳞癌 + 可算 OS）
"""
import json
import os
import re
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw")
PROC = os.path.join(ROOT, "data", "processed")
os.makedirs(PROC, exist_ok=True)

# 鳞状细胞癌形态学编码（ICD-O-3）
SQUAMOUS = {"8070/3", "8071/3", "8072/3", "8073/3", "8074/3",
            "8075/3", "8076/3", "8078/3", "8083/3", "8084/3", "8052/3"}

SUBSITE_RULES = [
    ("口咽", ["base of tongue", "tonsil", "oropharyn", "soft palate", "uvula",
              "vallecula", "lingual tonsil", "adenoid"]),
    ("下咽", ["hypopharyn", "pyriform", "postcricoid", "aryepiglottic"]),
    ("喉",   ["larynx", "glottis", "supraglottis", "subglottis", "vocal cord"]),
    ("口腔", ["lip", "tongue", "floor of mouth", "gum", "buccal", "cheek mucosa",
              "hard palate", "mouth", "oral cavity", "retromolar", "palate, nos"]),
    ("其他头颈", ["head, face or neck", "head and neck", "nasal cavity", "sinus",
                  "parotid", "salivary", "pharyn"]),
]
NON_HEADNECK = ["prostate", "brain", "lung", "breast", "kidney", "colon", "rectum",
                "ovary", "uterus", "stomach", "liver", "pancrea", "skin", "bladder",
                "thyroid", "lymph node", "bone marrow", "unknown"]


def map_site(organ):
    if not organ:
        return "未知"
    o = str(organ).lower()
    for name, keys in SUBSITE_RULES:
        for k in keys:
            if k in o:
                return name
    for k in NON_HEADNECK:
        if k in o:
            return "非头颈"
    return "其他头颈"


def norm_stage(path_stage, clin_stage):
    s = path_stage or clin_stage
    if not s or s == "Not Reported":
        return None
    s = str(s).replace("Stage ", "").strip().upper()
    return {"I": "I", "II": "II", "III": "III", "IVA": "IVA",
            "IVB": "IVB", "IVC": "IVC", "IV": "IV", "0": "0", "IS": "0"}.get(s, s)


def norm_tn(x):
    if not x or x in ("Not Reported", "unknown", "Unknown"):
        return None
    return str(x).strip().upper().replace("T", "T", 1) if str(x).startswith(("T", "t")) \
        else str(x).strip().upper()


def pick_primary_diagnosis(dxs):
    """挑选主诊断：优先 主诊断标记 + 鳞癌 + 头颈部位；逐级放宽。"""
    if not dxs:
        return None, "无诊断记录"
    sq_head = [d for d in dxs
               if d.get("diagnosis_is_primary_disease")
               and d.get("morphology") in SQUAMOUS
               and map_site(d.get("tissue_or_organ_of_origin")) not in ("非头颈", "未知")]
    if sq_head:
        return sq_head[0], "primary+鳞癌+头颈"
    sq = [d for d in dxs if d.get("morphology") in SQUAMOUS]
    if sq:
        return sq[0], "鳞癌(未标记主诊断)"
    prim = [d for d in dxs if d.get("diagnosis_is_primary_disease")]
    if prim:
        return prim[0], "主诊断标记(非鳞癌)"
    return dxs[0], "默认首条"


def main():
    d = json.load(open(os.path.join(RAW, "tcga_hnsc_clinical.json"), encoding="utf-8"))
    rows = []
    for h in d["data"]["hits"]:
        demo = h.get("demographic") or {}
        dxs = h.get("diagnoses") or []
        ex = (h.get("exposures") or [{}])[0]
        sel, how = pick_primary_diagnosis(dxs)

        if sel:
            days_death = sel.get("days_to_death") or demo.get("days_to_death")
            days_fu = sel.get("days_to_last_follow_up")
            vital = demo.get("vital_status") or sel.get("vital_status")
            stage = norm_stage(sel.get("ajcc_pathologic_stage"), sel.get("ajcc_clinical_stage"))
            t = sel.get("ajcc_pathologic_t") or sel.get("ajcc_clinical_t")
            n = sel.get("ajcc_pathologic_n") or sel.get("ajcc_clinical_n")
            m = sel.get("ajcc_pathologic_m") or sel.get("ajcc_clinical_m")
            organ = sel.get("tissue_or_organ_of_origin")
            morph = sel.get("morphology")
            pdx = sel.get("primary_diagnosis")
            grade = sel.get("tumor_grade")
            year = sel.get("year_of_diagnosis")
            rec = sel.get("days_to_recurrence")
        else:
            days_death = demo.get("days_to_death")
            days_fu = None
            vital = demo.get("vital_status")
            stage = t = n = m = organ = morph = pdx = grade = year = rec = None

        if vital == "Dead":
            os_time, os_event = days_death, 1
        elif vital == "Alive":
            os_time, os_event = days_fu, 0
        else:
            os_time, os_event = (days_death if days_death is not None else days_fu), None

        n_dx = len(dxs)
        sex = demo.get("sex_at_birth") or demo.get("gender")
        rows.append({
            "case_id": h.get("case_id"),
            "bcr_patient_barcode": h.get("submitter_id"),
            "n_diagnoses": n_dx,
            "multi_primary": n_dx > 1,
            "diagnosis_rule": how,
            "sex_at_birth": sex,
            "gender": sex,
            "age_at_index": demo.get("age_at_index"),
            "race": demo.get("race"),
            "ethnicity": demo.get("ethnicity"),
            "primary_diagnosis": pdx,
            "morphology": morph,
            "is_squamous": morph in SQUAMOUS,
            "tumor_grade": grade,
            "tissue_or_organ_of_origin": organ,
            "subsite": map_site(organ),
            "stage": stage,
            "t_stage": t,
            "n_stage": n,
            "m_stage": m,
            "year_of_diagnosis": year,
            "smoking_status": ex.get("tobacco_smoking_status"),
            "cigarettes_per_day": ex.get("cigarettes_per_day"),
            "pack_years": ex.get("pack_years_smoked"),
            "alcohol_history": ex.get("alcohol_history"),
            "alcohol_intensity": ex.get("alcohol_intensity"),
            "vital_status": vital,
            "days_to_death": days_death,
            "days_to_last_follow_up": days_fu,
            "days_to_recurrence": rec,
            "os_time_days": os_time,
            "os_event": os_event,
        })

    df = pd.DataFrame(rows)
    out = os.path.join(PROC, "hnsc_clinical.csv")
    df.to_csv(out, index=False, encoding="utf-8-sig")

    # ---- 分析集：头颈鳞癌 + 可算 OS ----
    an = df[(df["is_squamous"])
            & (df["subsite"].isin(["口腔", "口咽", "喉", "下咽", "其他头颈"]))
            & (df["os_time_days"].notna()) & (df["os_time_days"] > 0)
            & (df["os_event"].notna())].copy()
    an["os_event"] = an["os_event"].astype(int)
    an_out = os.path.join(PROC, "hnsc_analysis_set.csv")
    an.to_csv(an_out, index=False, encoding="utf-8-sig")

    print("[ok] 全部病例 %d -> %s" % (len(df), out))
    print("[ok] 分析集   %d -> %s" % (len(an), an_out))

    print("\n--- 主诊断筛选规则分布 ---")
    print(df["diagnosis_rule"].value_counts().to_string())
    print("\n--- 多原发病例 ---")
    print("  含多诊断的病例: %d（已按主诊断规则筛选）" % int(df["multi_primary"].sum()))

    print("\n--- 分析集字段完整度 ---")
    for c in ["stage", "t_stage", "n_stage", "m_stage", "smoking_status",
              "age_at_index", "tumor_grade", "subsite"]:
        nn = an[c].notna().sum()
        print("  %-16s %4d/%d (%.1f%%)" % (c, nn, len(an), 100.0 * nn / len(an)))

    print("\n--- 分析集亚部位分布 ---")
    print(an["subsite"].value_counts().to_string())
    print("\n--- 分析集分期分布 ---")
    print(an["stage"].value_counts(dropna=False).to_string())
    print("\n--- 生存 ---")
    print("  可算 OS: %d；死亡事件: %d；中位随访(d): %.0f"
          % (len(an), int(an["os_event"].sum()),
             an.loc[an["os_event"] == 0, "os_time_days"].median()))


if __name__ == "__main__":
    main()
