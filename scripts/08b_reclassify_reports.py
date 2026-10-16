# -*- coding: utf-8 -*-
"""病理报告 HPV 记录：文本缓存 + 改进版分类规则。

背景：首版规则（08）对 530 份报告只判出 40 份，而诊断显示另有 54 份"未判定"
报告中其实含明确结果。漏判原因有三：
  1. 双判词句（"p16 strong marking, and HPV is negative"）被 pos/neg 冲突规则丢弃；
  2. OCR 噪声（~egatlve / HPVis / RELATF.D / hnmunostalns）导致判词匹配失败；
  3. 诊断性表述（"INFILTRATING HPV-RELATED SQUAMOUS CELL CARCINOMA"）未被当作证据。

本脚本：先把全部 PDF 文本缓存为 jsonl（只做一次，约 6 分钟），
再基于缓存用分句 + 就近匹配 + OCR 容错规则重新判定，可反复快速迭代。

输出：
  data/raw/path_reports_text.jsonl          文本缓存
  data/raw/tcga_hnsc_path_report_hpv_v2.csv 改进后的判定
"""
import glob
import json
import os
import re
import sys
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw")
PDFDIR = os.path.join(RAW, "path_reports")
CACHE = os.path.join(RAW, "path_reports_text.jsonl")

# ---------------------------------------------------------------- OCR 容错
# 大写化 + 常见 OCR 讹误归一，便于统一正则
OCR_FIX = [
    (re.compile(r"[~\-_]\s*eg[a-z]{0,3}tive", re.I), "negative"),
    (re.compile(r"\bneg[a-z]{0,2}tlve\b", re.I), "negative"),
    (re.compile(r"\bposltive\b|\bpositlve\b", re.I), "positive"),
    (re.compile(r"\bhnmuno", re.I), "immuno"),
    (re.compile(r"\bslaln", re.I), "stain"),
    (re.compile(r"\bHPVis\b|\bHPV\s*is\b", re.I), "HPV is"),
    (re.compile(r"hpv[\s\-]*relat[a-z.]*", re.I), "hpv-related"),
    (re.compile(r"\bIHPV\b", re.I), "HPV"),
    (re.compile(r"\bP16\b", re.I), "p16"),
    (re.compile(r"\bISHx?\b", re.I), "ISH"),
    (re.compile(r"\bIN\s*SITU\s+HYBRIDIZATION\w*", re.I), "in situ hybridization"),
]


def fix_ocr(s):
    for pat, rep in OCR_FIX:
        s = pat.sub(rep, s)
    return s


# ------------------------------------------------- 非结果句（排除）
# 这些句子提到 HPV/p16 但并未给出检测结果
NONRESULT = re.compile(
    r"testing\s+ordered|test(?:ing)?\s+(?:has\s+been\s+)?request(?:ed)?|"
    r"\bpending\b|will\s+be\s+reported\s+separately|"
    r"ordered\s*:\s*(?:yes|no|per)|"
    r"instructions|"
    r"1\s*x\s*(?:p16|papilloma|hpv)|"
    r"to\s+be\s+performed|not\s+performed|"
    r"please\s+test|"
    r"reagents?\s*\(ASRs?\)",
    re.I)

# ------------------------------------------------- 结果判词（就近匹配）
OBJ = r"(?:p16|hpv(?:\s*[-/ ]?\s*(?:16|18|ish|in\s*situ))?|in\s*situ\s*hybrid\w*|\bish\b|papillomavirus?)"
POSW = r"(?:positive|pos\b|\+|detected|strong\s+(?:dot[- ]like\s+)?nuclear\s+signal|strong\s+marking|diffuse\s+strong|block\s*positive|reactive)"
NEGW = r"(?:negative|neg\b|not\s+detected|non[-\s]?reactive|no\s+evidence\s+of|wild\s*type|absent)"

RE_POS = re.compile(OBJ + r"[^.\n;]{0,60}?" + POSW, re.I)
RE_NEG = re.compile(OBJ + r"[^.\n;]{0,60}?" + NEGW, re.I)
# 反向语序：判词在前（"Negative for HPV"）
RE_NEG_REV = re.compile(r"\bnegative\s+for\s+" + OBJ, re.I)
# 诊断性表述
RE_DIAG_POS = re.compile(r"hpv[-\s]*related|hpv[-\s]*positive|hpv[-\s]*associated", re.I)

# p16 与 ISH 分开记录（p16 是替代标志物，ISH 才是金标准）
RE_P16 = re.compile(r"p16", re.I)
RE_ISH = re.compile(r"\bish\b|in\s*situ\s*hybrid", re.I)


def sentences(txt):
    """按句号/分号/换行切句，保留表格行。"""
    parts = re.split(r"(?<=[.;:])\s+|\n", txt)
    return [p.strip() for p in parts if p.strip()]


def clauses(s):
    """把句子再拆成子句，使 'p16 strong marking, and HPV is negative'
    能分别归到两个检测对象上。"""
    return [c.strip() for c in re.split(r",|;\s*|\band\b|\bwith\b", s, flags=re.I)
            if c.strip()]


def classify(txt):
    """返回 (call, p16_call, ish_call, evidence)。call ∈ {Positive, Negative, None}
    规则：按子句逐个对象判定；ISH 优先（金标准），无 ISH 时回退 p16。"""
    t = fix_ocr(txt)
    p16v = ishv = None
    ev = []
    for s in sentences(t):
        if len(s) > 500:
            continue
        if NONRESULT.search(s) and not (RE_POS.search(s) or RE_NEG.search(s)):
            continue
        if not (RE_P16.search(s) or re.search(r"hpv", s, re.I)):
            continue
        # 诊断性表述（"HPV-related SCC"）直接判阳性
        if RE_DIAG_POS.search(s) and not RE_NEG.search(s):
            ishv = "Positive"
            ev.append(s)
            continue
        for c in clauses(s):
            is_p16 = bool(RE_P16.search(c))
            is_ish = bool(RE_ISH.search(c)) or bool(re.search(r"hpv", c, re.I))
            if not (is_p16 or is_ish):
                # 子句本身无对象，但整句是单对象时可继承
                if len(clauses(s)) == 1:
                    is_p16 = bool(RE_P16.search(s))
                    is_ish = bool(RE_ISH.search(s)) or bool(re.search(r"hpv", s, re.I))
                else:
                    continue
            pos = RE_POS.search(c)
            neg = RE_NEG.search(c) or RE_NEG_REV.search(c)
            if pos and not neg:
                if is_p16:
                    p16v = "Positive"
                if is_ish:
                    ishv = "Positive"
                ev.append(c)
            elif neg and not pos:
                if is_p16:
                    p16v = "Negative"
                if is_ish:
                    ishv = "Negative"
                ev.append(c)
    call = ishv or p16v
    return call, p16v, ishv, " | ".join(ev[:4])[:600]


def dump_texts():
    """把全部 PDF 文本抽到 jsonl 缓存（跳过已缓存项）。"""
    done = set()
    if os.path.exists(CACHE):
        with open(CACHE, encoding="utf-8") as fh:
            for line in fh:
                try:
                    done.add(json.loads(line)["file_id"])
                except Exception:
                    pass
    files = [p for p in sorted(glob.glob(os.path.join(PDFDIR, "*.pdf")))
             if os.path.basename(p)[:-4] not in done]
    print("[info] 待解析 PDF %d 份（已缓存 %d）" % (len(files), len(done)))
    sys.stdout.flush()
    if not files:
        return

    def one(p):
        from pypdf import PdfReader
        fid = os.path.basename(p)[:-4]
        try:
            txt = "\n".join((pg.extract_text() or "") for pg in PdfReader(p).pages)
        except Exception as e:
            txt = ""
            return {"file_id": fid, "text": "", "err": str(e)}
        return {"file_id": fid, "text": txt, "err": ""}

    n = 0
    with open(CACHE, "a", encoding="utf-8") as fh:
        with ThreadPoolExecutor(max_workers=8) as ex:
            for r in ex.map(one, files):
                fh.write(json.dumps(r, ensure_ascii=False) + "\n")
                n += 1
                if n % 100 == 0:
                    print("   ... %d/%d" % (n, len(files)))
                    sys.stdout.flush()
    print("[ok] 文本缓存 -> %s" % CACHE)


def main():
    if "--dump" in sys.argv or not os.path.exists(CACHE):
        dump_texts()

    # 建立 file_id → submitter_id 映射（复用 08 的产物）
    import csv
    f2s = {}
    p = os.path.join(RAW, "tcga_hnsc_path_report_hpv.csv")
    if os.path.exists(p):
        with open(p, encoding="utf-8-sig") as fh:
            for r in csv.DictReader(fh):
                f2s[r["file_id"]] = r["submitter_id"]

    rows = []
    with open(CACHE, encoding="utf-8") as fh:
        for line in fh:
            r = json.loads(line)
            call, p16v, ishv, ev = classify(r["text"])
            rows.append({
                "submitter_id": f2s.get(r["file_id"], ""),
                "file_id": r["file_id"],
                "hpv_from_report": call or "",
                "p16_from_report": p16v or "",
                "ish_from_report": ishv or "",
                "evidence": ev,
                "text_len": len(r["text"]),
            })

    outp = os.path.join(RAW, "tcga_hnsc_path_report_hpv_v2.csv")
    with open(outp, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=["submitter_id", "file_id",
                                           "hpv_from_report", "p16_from_report",
                                           "ish_from_report", "text_len", "evidence"])
        w.writeheader()
        for r in sorted(rows, key=lambda x: (x["submitter_id"] or "", x["file_id"])):
            w.writerow(r)

    from collections import Counter
    c = Counter(r["hpv_from_report"] or "(未判定)" for r in rows)
    print("\n[ok] 重分类 %d 份" % len(rows))
    print("     HPV 判定:", dict(c))
    print("     p16 判定: ", dict(Counter(r["p16_from_report"] or "(无)" for r in rows)))
    print("     ISH 判定: ", dict(Counter(r["ish_from_report"] or "(无)" for r in rows)))
    print("     输出 ->", outp)


if __name__ == "__main__":
    main()
