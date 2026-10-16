# -*- coding: utf-8 -*-
"""诊断：为何部分病理报告提及 p16/HPV 却未能判定。"""
import glob
import os
import re
from concurrent.futures import ThreadPoolExecutor

from pypdf import PdfReader

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "data", "raw", "path_reports")


def get_text(p):
    try:
        return "\n".join((pg.extract_text() or "") for pg in PdfReader(p).pages)
    except Exception:
        return ""


def mentions(t):
    return bool(re.search(r"p\s*16|p16|hpv|papilloma", t, re.I))


def main():
    files = sorted(glob.glob(os.path.join(D, "*.pdf")))
    texts = []
    with ThreadPoolExecutor(max_workers=8) as ex:
        for p, t in zip(files, ex.map(get_text, files)):
            if mentions(t):
                texts.append((p, t))
    print("提及 p16/HPV 的报告 %d 份" % len(texts))

    unresolved = []
    for p, t in texts:
        sents = [s.strip() for s in re.split(r"(?<=[.;\n])\s+", t)
                 if re.search(r"p\s*16|p16|hpv|papilloma", s, re.I) and len(s) < 400]
        pos = [s for s in sents if re.search(r"\bpositive\b|\bpos\b|\+", s, re.I)]
        neg = [s for s in sents if re.search(r"\bnegative\b|\bneg\b|not\s+detected|"
                                             r"non-?reactive|no\s+evidence", s, re.I)]
        if not (pos and not neg) and not (neg and not pos):
            unresolved.append((p, sents[:4]))

    print("未判定 %d 份\n" % len(unresolved))
    for p, sents in unresolved[:30]:
        print("### " + os.path.basename(p))
        for s in sents:
            print("   · " + s[:220].replace("\n", " "))
        print()


if __name__ == "__main__":
    main()
