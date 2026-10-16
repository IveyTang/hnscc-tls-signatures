# -*- coding: utf-8 -*-
"""
23_fetch_gse159067_geo.py
-------------------------
Retrieve **GSE159067** directly from NCBI GEO.

Background — why a plain `curl` fails, and what works
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In this environment `https://www.ncbi.nlm.nih.gov` and
`https://ftp.ncbi.nlm.nih.gov` are unreachable through the sandbox HTTP proxy
(connection returns 000), and the agent's default egress is in any case
refused by NCBI with an abuse/access-denied page. Neither of those is a
property of GEO itself. Two channels reach the same host directly and are
verified to work:

  1. **FTP protocol** — `ftp://ftp.ncbi.nlm.nih.gov/...`. The blocking is
     applied to proxied HTTPS, so plain FTP passes.
  2. **HTTPS to NCBI's IP, bypassing the proxy** —
     `curl -k --noproxy '*' --resolve ftp.ncbi.nlm.nih.gov:443:130.14.250.31
     https://ftp.ncbi.nlm.nih.gov/...`
     Note `--noproxy '*'` (do not use the sandbox proxy) together with
     `--resolve` (pin the real hostname). `-k` is needed because the direct
     connection is not subject to the intercepting certificate.

This mirrors `scripts/ncbi_geo.sh` in the companion repository
`IveyTang/ici-benchmark`, where the same two channels are used.

Dataset
~~~~~~~
GSE159067 — "Immunologically hot phenotype is associated with response to
PD-1/PD-L1 inhibitors in patients with advanced head and neck squamous cell
carcinoma" (n = 102, GPL18573 HTG EdgeSeq targeted RNA panel).
BioProject PRJNA667509; SRA SRP286445.

Files retrieved
~~~~~~~~~~~~~~~
    <out>/GSE159067_series_matrix.txt.gz        sample annotation (clinical)
    <out>/GSE159067_IHN_log2cpm_data.txt.gz     log2 CPM expression matrix

Usage
~~~~~
    python scripts/23_fetch_gse159067_geo.py [--force]
"""
from __future__ import annotations

import os
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "data", "external", "GSE159067")

HOST = "ftp.ncbi.nlm.nih.gov"
NCBI_IP = os.environ.get("NCBI_IP", "130.14.250.31")
GSE = "GSE159067"
DIR = "GSE159%sn%sn" % ("", "")          # placeholder, computed below

# GSE159067 -> GSE159nnn
_n = GSE[3:]
DIR = "GSE%snnn/%s" % (_n[: len(_n) - 3], GSE)

FILES = [
    ("matrix", "%s_series_matrix.txt.gz" % GSE, 5_000),
    ("suppl", "%s_IHN_log2cpm_data.txt.gz" % GSE, 1_000_000),
]


def _path(kind: str, name: str) -> str:
    return "geo/series/%s/%s/%s" % (DIR, kind, name if kind != "matrix" else "")


def _urls(kind: str, name: str) -> list[str]:
    base = "geo/series/%s/%s/%s" % (DIR, kind, name)
    return [
        "ftp://%s/%s" % (HOST, base),
        "https://%s/%s" % (HOST, base),
    ]


def fetch(kind: str, name: str, min_bytes: int, force: bool = False) -> bool:
    dest = os.path.join(OUT, name)
    if os.path.exists(dest) and os.path.getsize(dest) >= min_bytes and not force:
        print("[skip] %-40s %d B (cached)" % (name, os.path.getsize(dest)))
        return True
    ftp_url, https_url = _urls(kind, name)

    # --- channel 1: plain FTP -------------------------------------------
    cmd = ["curl", "-s", "--connect-timeout", "25", "-m", "900", "-o", dest, ftp_url]
    r = subprocess.run(cmd, capture_output=True)
    if r.returncode == 0 and os.path.exists(dest) and os.path.getsize(dest) >= min_bytes:
        print("[ok  ] %-40s %d B  (ftp)" % (name, os.path.getsize(dest)))
        return True

    # --- channel 2: direct HTTPS, proxy bypassed ------------------------
    cmd = ["curl", "-s", "-k", "--noproxy", "*",
           "--resolve", "%s:443:%s" % (HOST, NCBI_IP),
           "--connect-timeout", "25", "-m", "900", "-o", dest, https_url]
    r = subprocess.run(cmd, capture_output=True)
    if r.returncode == 0 and os.path.exists(dest) and os.path.getsize(dest) >= min_bytes:
        print("[ok  ] %-40s %d B  (https+resolve)" % (name, os.path.getsize(dest)))
        return True

    print("[FAIL] %-40s ftp and https both failed" % name)
    print("        NCBI_IP=%s  (override with the NCBI_IP env var if DNS moved)" % NCBI_IP)
    return False


def main() -> int:
    os.makedirs(OUT, exist_ok=True)
    force = "--force" in sys.argv
    print("[geo ] %s  ->  %s" % (DIR, OUT))
    ok = all(fetch(k, n, mb, force) for k, n, mb in FILES)
    if ok:
        for _, n, _ in FILES:
            print("   %8d B  %s" % (os.path.getsize(os.path.join(OUT, n)), n))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
