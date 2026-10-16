#!/usr/bin/env python3
"""Publish the archived release of this repository.

Two things have to happen for the manuscript's "Availability of data and
materials" statement to be complete:

    1. A GitHub release has to exist, tagged with the version in CITATION.cff.
       If the Zenodo-GitHub integration is switched on, Zenodo archives that
       release automatically and mints a DOI.
    2. Either that DOI, or one minted through the Zenodo API directly, has to be
       written back into README.md, CITATION.cff and the manuscript.

This script covers both routes. It stores no credentials: the GitHub token is
read from the Git Credential Manager (the same store `git push` uses), and the
Zenodo token from the ZENODO_TOKEN environment variable.

Usage
-----
    # 1. GitHub release only (route A, after enabling Zenodo integration)
    python zenodo/create_release.py --github-release

    # 2. Zenodo deposit via API instead (route B; needs ZENODO_TOKEN)
    set ZENODO_TOKEN=...            # Windows
    export ZENODO_TOKEN=...         # bash
    python zenodo/create_release.py --zenodo-deposit

    # 3. Both
    python zenodo/create_release.py --github-release --zenodo-deposit

    # Dry run: print what would happen, change nothing
    python zenodo/create_release.py --github-release --dry-run

Environment
-----------
    GITHUB_REPO     overrides the detected "owner/name" (default: parsed from
                    the `origin` remote)
    ZENODO_TOKEN    Zenodo personal access token with deposit:write and
                    deposit:actions
    ZENODO_SANDBOX  set to 1 to target sandbox.zenodo.org instead
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import urllib.error
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
UA = "hnscc-tls-signatures-release-tool"
TIMEOUT = 120


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def log(msg: str) -> None:
    print(msg, flush=True)


def run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


def api(url, token=None, method="GET", payload=None, raw=None,
        content_type="application/json", ok=(200, 201)):
    data = None
    headers = {"User-Agent": UA, "Accept": "application/vnd.github+json"}
    if token:
        headers["Authorization"] = "Bearer " + token
    if raw is not None:
        data = raw
        headers["Content-Type"] = content_type
    elif payload is not None:
        data = json.dumps(payload).encode("utf-8")
        headers["Content-Type"] = "application/json"
    req = urllib.request.Request(url, data=data, headers=headers, method=method)
    try:
        with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
            body = r.read()
            if not body:
                return r.status, {}
            try:
                return r.status, json.loads(body.decode("utf-8"))
            except json.JSONDecodeError:
                return r.status, {"_raw": body[:400].decode("utf-8", "replace")}
    except urllib.error.HTTPError as e:
        detail = e.read().decode("utf-8", "replace")[:600]
        raise SystemExit("HTTP %s from %s\n%s" % (e.code, url, detail))


def find_gcm():
    """Locate the Git Credential Manager shipped with the local Git."""
    for name in ("git-credential-manager", "git-credential-manager.exe"):
        p = shutil.which(name)
        if p:
            return p
    probe = run(["git", "--exec-path"])
    if probe.returncode == 0:
        for cand in ("git-credential-manager", "git-credential-manager.exe"):
            p = os.path.join(os.path.dirname(probe.stdout.strip()), cand)
            if os.path.exists(p):
                return p
    return None


def github_token():
    """Read the github.com credential back out of the system credential store."""
    env = os.environ.get("GITHUB_TOKEN") or os.environ.get("GH_TOKEN")
    if env:
        return env.strip()
    gcm = find_gcm()
    if not gcm:
        raise SystemExit(
            "No GitHub token available. Either set GITHUB_TOKEN, or run\n"
            "  git push\n"
            "once so that Git Credential Manager caches a credential.")
    p = run([gcm, "get"],
            input="protocol=https\nhost=github.com\n\n")
    token = ""
    for line in p.stdout.splitlines():
        if line.startswith("password="):
            token = line.split("=", 1)[1].strip()
    if not token:
        raise SystemExit(
            "Git Credential Manager did not return a token. Run `git push` in "
            "this repository and complete the browser sign-in, then retry.")
    return token


def remote_slug():
    if os.environ.get("GITHUB_REPO"):
        return os.environ["GITHUB_REPO"]
    p = run(["git", "remote", "get-url", "origin"], cwd=ROOT)
    url = p.stdout.strip()
    m = re.search(r"github\.com[:/]+([^/]+)/([^/]+?)(?:\.git)?$", url)
    if not m:
        raise SystemExit("Could not parse owner/name from origin remote: %r" % url)
    return "%s/%s" % (m.group(1), m.group(2))


def read_cff_version():
    path = os.path.join(ROOT, "CITATION.cff")
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            m = re.match(r'version:\s*"?([^"\s]+)"?', line)
            if m:
                return m.group(1)
    return "1.0.0"


def read_zenodo_metadata():
    path = os.path.join(ROOT, ".zenodo.json")
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def make_archive(version):
    """Zip the working tree the same way `git archive` would, minus .git."""
    out = os.path.join(ROOT, "..",
                       "%s-v%s.zip" % (os.path.basename(ROOT), version))
    out = os.path.abspath(out)
    skip_dirs = {".git", "__pycache__", ".venv", "venv", ".idea", ".vscode"}
    skip_dirs.add(os.path.join("data", "processed"))
    import zipfile
    n = 0
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        for base, dirs, files in os.walk(ROOT):
            rel = os.path.relpath(base, ROOT)
            dirs[:] = [d for d in dirs if d not in skip_dirs]
            if rel in skip_dirs or "__pycache__" in base:
                continue
            for f in sorted(files):
                p = os.path.join(base, f)
                if "__pycache__" in p:
                    continue
                z.write(p, os.path.join(os.path.basename(ROOT),
                                        os.path.relpath(p, ROOT)))
                n += 1
    log("  archive: %s (%d files, %.1f MB)"
        % (out, n, os.path.getsize(out) / 1024 / 1024))
    return out


# --------------------------------------------------------------------------
# route A - GitHub release
# --------------------------------------------------------------------------

def github_release(slug, token, version, tag, dry_run):
    log("[A] GitHub release %s on %s" % (tag, slug))
    if dry_run:
        log("  dry run - would create tag %s and release %s" % (tag, version))
        return None

    _, ref = api("https://api.github.com/repos/%s/git/ref/heads/main" % slug,
                 token)
    sha = ref["object"]["sha"]
    log("  main is at %s" % sha[:12])

    status, _ = api("https://api.github.com/repos/%s/git/refs" % slug, token)
    existing = {t["ref"] for t in _} if isinstance(_, list) else set()
    want = "refs/tags/" + tag
    if want in existing:
        log("  tag %s already exists - leaving it alone" % tag)
    else:
        api("https://api.github.com/repos/%s/git/refs" % slug, token,
            method="POST",
            payload={"ref": want, "sha": sha})
        log("  created tag %s -> %s" % (tag, sha[:12]))

    body = (
        "First archived release accompanying the manuscript\n"
        "*Tertiary lymphoid structure gene signatures are not interchangeable in "
        "head and neck squamous cell carcinoma: a retrospective multi-signature "
        "transcriptomic analysis* (Journal of Translational Medicine, submitted).\n\n"
        "Contents: the full analysis pipeline, derived patient-level tables, and "
        "all result tables and figures for 438 TCGA-HNSC patients.\n\n"
        "Code is MIT licensed; derived data are CC BY 4.0. Primary data are "
        "publicly available from the Genomic Data Commons (TCGA-HNSC)."
    )
    try:
        status, rel = api("https://api.github.com/repos/%s/releases" % slug,
                          token, method="POST",
                          payload={"tag_name": tag, "name": "v" + version,
                                   "body": body, "draft": False,
                                   "prerelease": False})
    except SystemExit as e:
        if "HTTP 422" in str(e):
            log("  a release for %s already exists - skipping" % tag)
            return None
        raise
    log("  release published: %s" % rel.get("html_url"))
    return rel.get("html_url")


# --------------------------------------------------------------------------
# route B - Zenodo deposit
# --------------------------------------------------------------------------

def zenodo_deposit(version, token, do_publish, dry_run):
    base = ("https://sandbox.zenodo.org/api"
            if os.environ.get("ZENODO_SANDBOX") == "1"
            else "https://zenodo.org/api")
    log("[B] Zenodo deposit (%s)" % base)
    meta = read_zenodo_metadata()
    meta.pop("communities", None)

    archive = make_archive(version)
    if dry_run:
        log("  dry run - would upload %s" % archive)
        return None

    _, dep = api(base + "/deposit/depositions", token, method="POST",
                 payload={"metadata": meta})
    dep_id = dep["id"]
    log("  created deposition %s" % dep_id)

    bucket = dep["links"]["bucket"]
    fname = os.path.basename(archive)
    with open(archive, "rb") as fh:
        api("%s/%s" % (bucket, fname), token, method="PUT", raw=fh.read(),
            content_type="application/octet-stream", ok=(200, 201))
    log("  uploaded %s" % fname)

    _, done = api("%s/deposit/depositions/%s" % (base, dep_id), token,
                  method="PUT", payload={"metadata": meta})
    log("  metadata set")

    if do_publish:
        _, pub = api("%s/deposit/depositions/%s/actions/publish"
                     % (base, dep_id), token, method="POST", ok=(202,))
        doi = pub.get("doi") or pub.get("metadata", {}).get("doi")
        log("  published - DOI %s" % doi)
        log("  record   %s" % pub.get("links", {}).get("record_html", ""))
        return doi
    log("  left as a draft: review at %s"
        % dep["links"].get("html", ""))
    return None


# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--github-release", action="store_true",
                    help="create tag + GitHub release (route A)")
    ap.add_argument("--zenodo-deposit", action="store_true",
                    help="create a Zenodo deposit via the API (route B)")
    ap.add_argument("--no-publish", action="store_true",
                    help="with --zenodo-deposit, leave the deposit as a draft")
    ap.add_argument("--tag", default=None, help="override the tag name")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the plan without changing anything")
    args = ap.parse_args()

    if not (args.github_release or args.zenodo_deposit):
        ap.error("choose --github-release, --zenodo-deposit, or both")

    version = read_cff_version()
    tag = args.tag or ("v" + version)
    log("repository : %s" % ROOT)
    log("version    : %s   tag: %s" % (version, tag))

    if args.github_release:
        slug = remote_slug()
        tok = github_token()
        log("github     : %s (token %d chars, from %s)"
            % (slug, len(tok),
               "environment" if os.environ.get("GITHUB_TOKEN") else "credential store"))
        github_release(slug, tok, version, tag, args.dry_run)

    if args.zenodo_deposit:
        tok = os.environ.get("ZENODO_TOKEN")
        if not tok:
            raise SystemExit(
                "ZENODO_TOKEN is not set. Create a personal access token at\n"
                "  https://zenodo.org/account/settings/applications/tokens/new/\n"
                "with the scopes deposit:write and deposit:actions, then\n"
                "  set ZENODO_TOKEN=...      (Windows)\n"
                "  export ZENODO_TOKEN=...   (bash)")
        zenodo_deposit(version, tok, not args.no_publish, args.dry_run)

    log("done")


if __name__ == "__main__":
    sys.exit(main())
