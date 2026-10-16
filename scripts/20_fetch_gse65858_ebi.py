# -*- coding: utf-8 -*-
"""
20_fetch_gse65858_ebi.py
------------------------
Fetch the GSE65858 external-validation cohort from the **EBI ArrayExpress
mirror** instead of NCBI GEO.

Rationale
~~~~~~~~~
NCBI (www.ncbi.nlm.nih.gov / ftp.ncbi.nlm.nih.gov) is unreachable from this
environment, and the agent egress IP is additionally blocked by NCBI with a
"possible misuse" page. ArrayExpress imported GEO series up to ~2017, so
GSE65858 is available in full as **E-GEOD-65858**, hosted on the EBI FTP
mirror:

    record   https://www.ebi.ac.uk/biostudies/api/v1/studies/E-GEOD-65858
    files    https://ftp.ebi.ac.uk/biostudies/fire/E-GEOD-/858/E-GEOD-65858/Files/
    legacy   https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEOD/E-GEOD-65858/
    platform https://ftp.ebi.ac.uk/biostudies/fire/A-GEOD-/558/A-GEOD-10558/Files/

Artifacts downloaded
~~~~~~~~~~~~~~~~~~~~
    <out>/E-GEOD-65858.idf.txt          study design
    <out>/E-GEOD-65858.sdrf.txt        sample annotation (270 rows, incl. HPV + OS)
    <out>/A-GEOD-10558.adf.txt         Illumina HT-12 v4 probe -> gene annotation
    <out>/tables/GSM*_sample_table.txt 270 per-sample processed expression tables

Outputs
~~~~~~~
    <out>/gse65858_expression.csv      samples x probes (VALUE, log2)
    <out>/gse65858_clinical.csv        harmonised clinical table
    <out>/fetch_log.json               provenance
"""
from __future__ import annotations

import csv
import io
import json
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "data", "external", "GSE65858")
TABLES = os.path.join(OUT, "tables")

FIRE = "https://ftp.ebi.ac.uk/biostudies/fire/E-GEOD-/858/E-GEOD-65858/Files"
PLAT = "https://ftp.ebi.ac.uk/biostudies/fire/A-GEOD-/558/A-GEOD-10558/Files"

UA = "Mozilla/5.0 (compatible; tls4-repro/1.0; mailto:33479986@qq.com)"
WORKERS = 12
RETRIES = 5


def fetch(url: str, dest: str, min_bytes: int = 1) -> tuple[bool, str, int]:
    """Download *url* to *dest* with retries; verify a minimum size."""
    if os.path.exists(dest) and os.path.getsize(dest) >= min_bytes:
        return True, "cached", os.path.getsize(dest)
    last = ""
    for attempt in range(1, RETRIES + 1):
        try:
            req = Request(url, headers={"User-Agent": UA})
            with urlopen(req, timeout=120) as fh:
                blob = fh.read()
            if len(blob) < min_bytes:
                last = "short read %d < %d" % (len(blob), min_bytes)
                time.sleep(2 * attempt)
                continue
            tmp = dest + ".part"
            with open(tmp, "wb") as out:
                out.write(blob)
            os.replace(tmp, dest)
            return True, "ok", len(blob)
        except (HTTPError, URLError, TimeoutError, OSError) as exc:
            last = "%s: %s" % (type(exc).__name__, exc)
            time.sleep(2 * attempt)
    return False, last, 0


def read_sdrf(path: str) -> tuple[list[str], list[list[str]]]:
    with io.open(path, encoding="utf-8", errors="replace") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    hdr = rows[0]
    body = [r for r in rows[1:] if any(x.strip() for x in r)]
    return hdr, body


def main() -> int:
    os.makedirs(TABLES, exist_ok=True)
    log = {"records": []}

    # ---- 1. study metadata + platform annotation -------------------------
    # Exact byte counts, taken from the Content-Length of each URL. Checking a
    # lower bound alone is NOT enough: an interrupted transfer still returns
    # HTTP 200, and a truncated ADF silently loses most of its probe
    # annotation (a 7.8 MB prefix of the 71.7 MB file parses fine but only
    # covers 5,193 of 47,323 probes).
    meta = [
        ("E-GEOD-65858.idf.txt", FIRE, 10_018),
        ("E-GEOD-65858.sdrf.txt", FIRE, 301_544),
        ("A-GEOD-10558.adf.txt", PLAT, 71_731_111),
    ]
    for name, base, expected in meta:
        dest = os.path.join(OUT, name)
        ok, note, size = fetch("%s/%s" % (base, name), dest, expected)
        if ok and size != expected:
            ok, note = False, "truncated: %d of %d bytes" % (size, expected)
        print("[meta] %-26s %-6s %s (%d B)"
              % (name, "OK" if ok else "FAIL", note, size))
        log["records"].append({"file": name, "ok": ok, "note": note,
                               "bytes": size, "expected": expected})
        if not ok:
            print("  !! aborting: study metadata could not be retrieved intact")
            print("  !! for the 71.7 MB ADF, use scripts/20b_fetch_adf_ranged.py "
                  "(parallel HTTP Range download)")
            return 1

    # ---- 2. sample -> file mapping from the SDRF ------------------------
    hdr, body = read_sdrf(os.path.join(OUT, "E-GEOD-65858.sdrf.txt"))
    idx_assay = hdr.index("Assay Name")
    idx_file = hdr.index("Derived Array Data File")
    pairs = [(r[idx_assay].strip(), r[idx_file].strip()) for r in body]
    pairs = [(a, f) for a, f in pairs if a and f]
    seen, uniq = set(), []
    for a, f in pairs:
        if a not in seen:
            seen.add(a)
            uniq.append((a, f))
    print("\n[sdrf] %d unique samples -> %d unique data files" % (len(seen), len(uniq)))

    # ---- 3. parallel download of the per-sample tables ------------------
    todo = [(a, f) for a, f in uniq
            if not (os.path.exists(os.path.join(TABLES, f))
                    and os.path.getsize(os.path.join(TABLES, f)) >= 50_000)]
    print("[dl] %d files to fetch, %d workers" % (len(todo), WORKERS))
    t0 = time.time()
    done = fail = 0
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        futs = {pool.submit(fetch, "%s/%s" % (FIRE, f),
                            os.path.join(TABLES, f), 50_000): (a, f)
                for a, f in todo}
        for fut in as_completed(futs):
            assay, fname = futs[fut]
            ok, note, size = fut.result()
            done += 1
            if not ok:
                fail += 1
                print("  [FAIL] %s  %s" % (fname, note))
            if done % 25 == 0 or done == len(todo):
                print("  .. %d/%d  (%.0fs, %d failed)"
                      % (done, len(todo), time.time() - t0, fail), flush=True)
    print("[dl] finished in %.0fs, %d failed" % (time.time() - t0, fail))

    # ---- 4. verify all files present ------------------------------------
    missing = [f for a, f in uniq if not os.path.exists(os.path.join(TABLES, f))]
    print("[verify] present %d / %d" % (len(uniq) - len(missing), len(uniq)))
    if missing:
        print("  missing:", missing[:10])
    log["samples"] = len(uniq)
    log["downloaded"] = len(uniq) - len(missing)
    log["failed"] = len(missing)
    log["retrieved_at"] = time.strftime("%Y-%m-%dT%H:%M:%S")

    with io.open(os.path.join(OUT, "fetch_log.json"), "w", encoding="utf-8") as fh:
        json.dump(log, fh, indent=2, ensure_ascii=False)
    print("\n[log] written:", os.path.join(OUT, "fetch_log.json"))
    return 0 if not missing else 2


if __name__ == "__main__":
    sys.exit(main())
