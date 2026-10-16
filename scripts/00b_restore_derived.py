# -*- coding: utf-8 -*-
"""Restore the pipeline working directory from the committed derived tables.

This lets the statistical steps (03, 04, 06, 11) be reproduced without any
network access. The committed tables under ``data/derived/`` are the exact
outputs of steps 00-02 and 07-10, so downstream steps see the same inputs they
would have seen had the full pipeline been run.

Usage:
    python scripts/00b_restore_derived.py [--force]

Committed tables and the step that normally produces them:

    hnsc_analysis_set.csv       <- step 01
    hnsc_gene_matrix.csv        <- step 02   (1291 genes x 563 samples, 5.1 MB)
    hnsc_gene_matrix_hpv.csv    <- step 09
    hnsc_tls_scores.csv         <- step 03
    hnsc_hpv_status.csv         <- step 10
    hnsc_with_hpv.csv           <- step 10

Steps that still require network access (they read files not committed here):

    00, 01, 02   GDC clinical records and STAR-Counts expression
    05           external validation cohorts (GSE65858, GSE159067)
    07, 08, 08b  BCR clinical XML and diagnostic pathology report PDFs
    09, 10       derive from step 07/08 outputs
"""
import argparse
import os
import shutil
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "data", "derived")
PROC = os.path.join(ROOT, "data", "processed")

# Extra provenance file kept one level up
EXTRA = [("expr_download_meta.json", os.path.join(ROOT, "data", "expr_download_meta.json"))]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--force", action="store_true",
                    help="overwrite files that already exist in data/processed/")
    args = ap.parse_args()

    if not os.path.isdir(SRC):
        sys.exit("ERROR: %s not found. Run this script from a full checkout of the repository." % SRC)

    os.makedirs(PROC, exist_ok=True)

    copied, skipped = [], []
    for name in sorted(os.listdir(SRC)):
        if not name.endswith(".csv"):
            continue
        src = os.path.join(SRC, name)
        dst = os.path.join(PROC, name)
        if os.path.exists(dst) and not args.force:
            skipped.append(name)
            continue
        shutil.copy2(src, dst)
        copied.append((name, os.path.getsize(dst)))

    for name, src in EXTRA:
        if not os.path.exists(src):
            continue
        dst = os.path.join(PROC, name)
        if os.path.exists(dst) and not args.force:
            skipped.append(name)
            continue
        shutil.copy2(src, dst)
        copied.append((name, os.path.getsize(dst)))

    for name, size in copied:
        print("[ok] %-32s %8.1f KB" % (name, size / 1024.0))
    for name in skipped:
        print("[--] %-32s already present (use --force to overwrite)" % name)

    print("\nRestored %d file(s) into %s" % (len(copied), PROC))
    print("\nYou can now run the offline-capable steps:")
    print("    python scripts/03_build_tls_score.py")
    print("    python scripts/04_internal_validation.py")
    print("    python scripts/06_microenvironment_nomogram.py")
    print("    python scripts/11_hpv_stratified.py")
    print("\nNote: step 03 re-reads the expression matrix; step 10 (HPV labelling)")
    print("additionally needs the raw BCR/PDF files from steps 07-08 and cannot be")
    print("run offline. Its output, hnsc_hpv_status.csv, is already restored above.")


if __name__ == "__main__":
    main()
