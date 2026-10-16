# -*- coding: utf-8 -*-
"""
20b_fetch_adf_ranged.py
-----------------------
Download the full Illumina HT-12 v4 array-design annotation from the EBI FTP
mirror using parallel HTTP Range requests.

Why this separate step
~~~~~~~~~~~~~~~~~~~~~~
`A-GEOD-10558.adf.txt` is 71,731,111 bytes. A single stream through the local
proxy sustains only ~175 KB/s, and a plain download is prone to silent
truncation (the earlier attempt stopped at 10.8% while still returning HTTP
200). Ranged requests let us fetch it quickly and *verify the length*, which
is the property we actually need: the file must contain all 47,323 probes so
that every gene in the TLS gene sets has a probe.

Layout of the file
~~~~~~~~~~~~~~~~~~
Lines 0-10   header block (key <tab> value)
Line  11     blank
Line  12     column header row (30 fields, includes Comment[Symbol])
Line  13+    one row per probe: Reporter Name, ..., Comment[Symbol], ...
"""
from __future__ import annotations

import io
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from urllib.request import Request, urlopen

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "data", "external", "GSE65858")
URL = ("https://ftp.ebi.ac.uk/biostudies/fire/A-GEOD-/558/"
       "A-GEOD-10558/Files/A-GEOD-10558.adf.txt")
DEST = os.path.join(OUT, "A-GEOD-10558.adf.txt")
EXPECTED = 71_731_111          # Content-Length observed via HEAD
CHUNKS = 12
UA = "Mozilla/5.0 (compatible; tls4-repro/1.0; mailto:33479986@qq.com)"


def get_range(start: int, end: int, attempt_max: int = 5) -> bytes:
    last = ""
    for attempt in range(1, attempt_max + 1):
        try:
            req = Request(URL, headers={"User-Agent": UA,
                                        "Range": "bytes=%d-%d" % (start, end)})
            with urlopen(req, timeout=180) as fh:
                blob = fh.read()
            if blob:
                return blob
            last = "empty body"
        except Exception as exc:                       # noqa: BLE001
            last = "%s: %s" % (type(exc).__name__, exc)
        time.sleep(1.5 * attempt)
    raise RuntimeError("range %d-%d failed: %s" % (start, end, last))


def main() -> int:
    os.makedirs(OUT, exist_ok=True)

    if os.path.exists(DEST) and os.path.getsize(DEST) == EXPECTED:
        print("[adf] already complete (%d bytes) - nothing to do" % EXPECTED)
        return 0

    size = EXPECTED
    step = size // CHUNKS
    spans = []
    for i in range(CHUNKS):
        s = i * step
        e = (size - 1) if i == CHUNKS - 1 else (s + step - 1)
        spans.append((i, s, e))

    print("[adf] %s" % URL)
    print("[adf] %d bytes in %d chunks (%d workers)"
          % (size, CHUNKS, CHUNKS))
    t0 = time.time()
    parts: dict[int, bytes] = {}
    with ThreadPoolExecutor(max_workers=CHUNKS) as pool:
        futs = {pool.submit(get_range, s, e): i for i, s, e in spans}
        for fut in as_completed(futs):
            i = futs[fut]
            blob = fut.result()
            parts[i] = blob
            print("    chunk %2d  %9d bytes  (%.0fs)"
                  % (i, len(blob), time.time() - t0), flush=True)

    blob = b"".join(parts[i] for i in sorted(parts))
    print("[adf] assembled %d bytes" % len(blob))
    if len(blob) != size:
        print("  !! size mismatch: got %d, expected %d" % (len(blob), size))
        return 2

    tmp = DEST + ".part"
    with open(tmp, "wb") as fh:
        fh.write(blob)
    os.replace(tmp, DEST)
    print("[adf] written -> %s  (%.0fs total)"
          % (os.path.basename(DEST), time.time() - t0))

    # --- structural sanity check ----------------------------------------
    text = blob.decode("utf-8", errors="replace")
    lines = text.split("\n")
    hdr_row = None
    for i, ln in enumerate(lines[:40]):
        if ln.startswith("Reporter Name\t"):
            hdr_row = i
            break
    n_rows = len([l for l in lines if l.startswith("ILMN_")])
    print("[adf] header row at line %s ; ILMN_ rows = %d"
          % (hdr_row, n_rows))
    return 0 if n_rows > 40_000 else 3


if __name__ == "__main__":
    sys.exit(main())
