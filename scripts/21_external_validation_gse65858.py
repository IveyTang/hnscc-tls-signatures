# -*- coding: utf-8 -*-
"""
21_external_validation_gse65858.py
----------------------------------
External validation of the consensus TLS score (CTS) in **GSE65858**
(Wichmann/Rosolowski et al.), retrieved from the EBI ArrayExpress mirror as
`E-GEOD-65858` because NCBI GEO is unreachable from this environment.

Design
~~~~~~
The scoring method is copied *exactly* from the discovery pipeline
(`03_build_tls_score.py`) so that the validation is like-for-like and
platform-agnostic:

    per-gene z-score across samples  ->  gene-set mean (z_mean)
    CTS = mean percentile rank of the five TLS gene sets' z_mean

Nothing is re-tuned on the validation cohort: gene sets are unchanged, the
score definition is unchanged, and both the *discovery-derived* cut-off (0.249)
and a *cohort-internal* optimal cut-off are reported.

External cohort
~~~~~~~~~~~~~~~
    270 HNSCC tumours, Illumina HumanHT-12 v4 (GPL10558)
    94 deaths (34.8%)
    HPV16 RNA status available for 250/270 (DNA+RNA+ 35, DNA+RNA- 19, DNA- 196)
    UICC stage, T/N category, tumour site, age, sex, and the published
    expression cluster (incl. the immune-response "Atypical IR" cluster)

Inputs
~~~~~~
    data/external/GSE65858/E-GEOD-65858.sdrf.txt
    data/external/GSE65858/A-GEOD-10558.adf.txt
    data/external/GSE65858/tables/GSM*_sample_table.txt
    config/tls_genesets.json

Outputs
~~~~~~~
    results/21_ext_cts_scores.csv
    results/21_ext_convergence.csv
    results/21_ext_programme_correlation.csv
    results/21_ext_survival.csv
    results/21_ext_summary.json
    results/fig_ext_km_cts.png
    results/fig_ext_convergence.png
    results/fig_ext_cts_by_hpv.png
"""
from __future__ import annotations

import csv
import glob
import io
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
EXT = os.path.join(ROOT, "data", "external", "GSE65858")
TABLES = os.path.join(EXT, "tables")
RES = os.path.join(ROOT, "results")
os.makedirs(RES, exist_ok=True)

ADF = os.path.join(EXT, "A-GEOD-10558.adf.txt")
SDRF = os.path.join(EXT, "E-GEOD-65858.sdrf.txt")
GENESETS = os.path.join(ROOT, "config", "tls_genesets.json")

CONSENSUS_SETS = ["G1_Coppola_12chemokine", "G2_11chemokine_HCC",
                  "G3_HNSCC_13chemokine", "G4_TLS_maturation_GC",
                  "G6_HEV_lymphoid_stroma"]
CHEMOKINE_SETS = ["G1_Coppola_12chemokine", "G2_11chemokine_HCC",
                  "G3_HNSCC_13chemokine"]

# Microenvironmental programmes, identical to 06_microenvironment_nomogram.py
MARKERS = {
    "B cells":        ["MS4A1", "CD19", "CD79A", "CD79B", "TCL1A", "BANK1"],
    "GC B cells":     ["AICDA", "BCL6", "MEF2B", "RGS13", "MKI67"],
    "Plasma cells":   ["MZB1", "JCHAIN", "IGHG1", "IGKC", "XBP1", "PRDM1", "SDC1"],
    "Tfh":            ["CXCL13", "IL21", "BCL6", "ICOS", "PDCD1", "SH2D1A"],
    "CD8 T cells":    ["CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "NKG7"],
    "Treg":           ["FOXP3", "IL2RA", "CTLA4", "IKZF2"],
    "CAF":            ["ACTA2", "FAP", "COL1A1", "COL1A2", "PDGFRB", "TAGLN", "THY1"],
    "Endothelium/HEV": ["PECAM1", "VWF", "CLDN5", "PDPN", "CD34"],
    "Macrophages":    ["CD68", "CD163", "CSF1R", "MRC1"],
    "IFN-gamma":      ["IFNG", "STAT1", "IDO1", "CXCL9", "CXCL10", "GBP1", "IRF1"],
    "Checkpoint":     ["CD274", "PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT"],
    "Epithelial":     ["EPCAM", "KRT5", "KRT14", "KRT13", "KRT4"],
}

# Discovery-cohort optimal cut-off, carried over unchanged
DISCOVERY_CUT = 0.2489


# ---------------------------------------------------------------- helpers
def load_gene_sets():
    d = json.load(io.open(GENESETS, encoding="utf-8"))
    out = {}
    for k, v in d.items():
        if isinstance(v, dict) and isinstance(v.get("genes"), list):
            out[k] = [g.strip().upper() for g in v["genes"]]
    return out


def load_probe_map():
    """Reporter Name (ILMN_*) -> HGNC symbol, from the ArrayExpress ADF."""
    with io.open(ADF, encoding="utf-8", errors="replace", newline="") as fh:
        lines = fh.read().split("\n")
    hdr_i = next(i for i, ln in enumerate(lines[:60])
                 if ln.startswith("Reporter Name\t"))
    cols = lines[hdr_i].split("\t")
    i_rep = cols.index("Reporter Name")
    i_sym = None
    for cand in ("Comment[Symbol]", "Comment[ILMN_Gene]"):
        if cand in cols:
            i_sym = cols.index(cand)
            break
    if i_sym is None:
        raise RuntimeError("no gene-symbol column in ADF header")
    m, n = {}, 0
    for ln in lines[hdr_i + 1:]:
        if not ln.startswith("ILMN_"):
            continue
        f = ln.split("\t")
        if len(f) <= max(i_rep, i_sym):
            continue
        n += 1
        probe, sym = f[i_rep].strip(), f[i_sym].strip().upper()
        if probe and sym and sym not in ("", "NA"):
            m[probe] = sym
    return m, n, cols[hdr_i]


def load_expression(probe_map):
    files = sorted(glob.glob(os.path.join(TABLES, "GSM*_sample_table.txt")))
    if not files:
        raise RuntimeError("no sample tables found in %s" % TABLES)
    frames, kept = [], []
    for fp in files:
        gsm = os.path.basename(fp).split("_")[0]
        try:
            df = pd.read_csv(fp, sep="\t", index_col=0,
                             usecols=[0, 1], encoding="utf-8")
        except Exception:
            continue
        if df.shape[0] < 1000:
            continue
        s = pd.to_numeric(df.iloc[:, 0], errors="coerce")
        s.name = gsm
        frames.append(s)
        kept.append(gsm)
    mat = pd.concat(frames, axis=1)
    mat.index = [str(i).strip() for i in mat.index]
    mat = mat[~mat.index.duplicated(keep="first")]
    return mat, kept


def load_clinical():
    with io.open(SDRF, encoding="utf-8", errors="replace") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    hdr = rows[0]
    body = [r for r in rows[1:] if any(x.strip() for x in r)]

    def g(name, r):
        i = hdr.index(name) if name in hdr else None
        return r[i].strip() if (i is not None and i < len(r)) else ""

    recs = []
    for r in body:
        stage = g("Characteristics [uicc_stage]", r)
        stage_ord = {"I": 1, "II": 2, "III": 3, "IVA": 4, "IVB": 5, "IVC": 6}
        try:
            os_days = float(g("Characteristics [os]", r))
        except ValueError:
            os_days = np.nan
        ev = g("Characteristics [os_event]", r).upper()
        try:
            age = float(g("Characteristics [age]", r))
        except ValueError:
            age = np.nan
        recs.append({
            "gsm": g("Assay Name", r),
            "patient_id": g("Characteristics [ID]", r),
            "data_file": g("Derived Array Data File", r),
            "hpv_dna": g("Characteristics [hpv_dna]", r),
            "hpv_rna": g("Characteristics [hpv16_dna_rna]", r),
            "os_days": os_days,
            "os_event": (1 if ev == "TRUE" else 0) if ev else np.nan,
            "stage": stage,
            "stage_ord": stage_ord.get(stage, np.nan),
            "t_category": g("Characteristics [t_category]", r),
            "n_category": g("Characteristics [n_category]", r),
            "tumor_site": g("Characteristics [tumor_site]", r),
            "sex": g("Characteristics [sex]", r),
            "age": age,
            "cluster": g("Characteristics [consensus_cluster]", r),
            "tp53": g("Characteristics [tp53_mutation]", r),
            "treatment": g("Characteristics [treatment]", r),
        })
    df = pd.DataFrame(recs).set_index("gsm")
    return df


def zscore(mat):
    mu = mat.mean(axis=1)
    sd = mat.std(axis=1).replace(0, np.nan)
    return mat.sub(mu, axis=0).div(sd, axis=0)


def zmean(z, genes):
    present = [g for g in genes if g in z.index]
    if not present:
        return None, 0
    return z.loc[present].mean(axis=0), len(present)


def best_cutpoint(score, time, event, n_cand=60, min_grp=15):
    """Vectorised log-rank scan over 10-90% quantiles of the score."""
    qs = np.linspace(0.10, 0.90, n_cand)
    cands = np.unique(np.quantile(score, qs))
    best = (None, -np.inf)
    for c in cands:
        hi = score > c
        if hi.sum() < min_grp or (~hi).sum() < min_grp:
            continue
        t, e = time[hi], event[hi]
        u, f = time[~hi], event[~hi]
        chi2, _ = _logrank(t, e, u, f)
        if chi2 > best[1]:
            best = (float(c), float(chi2))
    return best


def _logrank(t1, e1, t2, e2):
    times = np.unique(np.concatenate([t1, t2]))
    O1 = E1 = V = 0.0
    for tt in times:
        n1 = (t1 >= tt).sum()
        n2 = (t2 >= tt).sum()
        n = n1 + n2
        d1 = ((t1 == tt) & (e1 == 1)).sum()
        d2 = ((t2 == tt) & (e2 == 1)).sum()
        d = d1 + d2
        if n <= 1 or d == 0:
            continue
        O1 += d1
        E1 += d * n1 / n
        V += d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1)
    if V <= 0:
        return 0.0, 1.0
    chi2 = (O1 - E1) ** 2 / V
    p = float(stats.chi2.sf(chi2, 1))
    return float(chi2), p


def main():
    print("=" * 74)
    print("GSE65858 external validation (via EBI ArrayExpress mirror)")
    print("=" * 74)

    gs = load_gene_sets()
    print("[genesets] %d sets: %s" % (len(gs), ", ".join(gs)))

    probe_map, n_probe, hdrline = load_probe_map()
    print("[adf] %d probes carry a gene symbol" % n_probe)

    raw, gsms = load_expression(probe_map)
    print("[expr] %d probes x %d samples" % raw.shape)

    # ---- collapse probes -> genes -------------------------------------
    sym = raw.index.map(lambda p: probe_map.get(p))
    keep = ~pd.isna(sym)
    sub = raw.loc[keep].copy()
    sub.index = sym[keep]
    gene_mat = sub.groupby(level=0).mean()
    print("[expr] %d genes after collapsing" % gene_mat.shape[0])

    clin = load_clinical()
    common = [s for s in gene_mat.columns if s in clin.index]
    gene_mat = gene_mat[common]
    clin = clin.loc[common]
    print("[merge] %d samples with expression + clinical" % len(common))

    z = zscore(gene_mat)
    print("[score] z-matrix %d x %d" % z.shape)

    # ---- gene set scores ----------------------------------------------
    scores, cover = {}, {}
    for name, genes in gs.items():
        s, n = zmean(z, genes)
        if s is not None:
            scores[name] = s
            cover[name] = "%d/%d" % (n, len(genes))
    sc = pd.DataFrame(scores)
    print("[score] coverage: %s" % json.dumps(cover))

    # ---- consensus TLS score (identical definition) --------------------
    have = [c for c in CONSENSUS_SETS if c in sc.columns]
    ranks = sc[have].rank(pct=True)
    cts = ranks.mean(axis=1)
    sc["CTS"] = cts

    # ---- programme scores ---------------------------------------------
    prog = {}
    for name, genes in MARKERS.items():
        s, n = zmean(z, genes)
        if s is not None and n >= 2:
            prog[name] = s
    prog_df = pd.DataFrame(prog)
    print("[score] %d microenvironmental programmes" % prog_df.shape[1])

    out = clin.join(sc)
    out = out.join(prog_df.add_prefix("prog_"))
    out.to_csv(os.path.join(RES, "21_ext_cts_scores.csv"),
               encoding="utf-8-sig")

    summary = {
        "cohort": "GSE65858 (E-GEOD-65858, EBI ArrayExpress mirror)",
        "n_samples": int(len(common)),
        "n_deaths": int(np.nansum(out["os_event"])),
        "geneset_coverage": cover,
        "n_genes_scored": int(gene_mat.shape[0]),
        "discovery_cutoff": DISCOVERY_CUT,
        "os_unit": "days (as submitted to ArrayExpress)",
    }

    # ================= A. signature convergence =========================
    conv = []
    chem = [c for c in CHEMOKINE_SETS if c in sc.columns]
    pairs = [(a, b) for i, a in enumerate(chem) for b in chem[i + 1:]]
    pairs += [(c, "G7_Tcell_inflamed_cytolytic")
              for c in chem + ["G4_TLS_maturation_GC", "G6_HEV_lymphoid_stroma"]
              if c in sc.columns and "G7_Tcell_inflamed_cytolytic" in sc.columns]
    for a, b in pairs:
        rho = stats.spearmanr(sc[a], sc[b]).statistic
        conv.append({"set_a": a, "set_b": b, "spearman_rho": round(float(rho), 4)})
    conv_df = pd.DataFrame(conv)
    conv_df.to_csv(os.path.join(RES, "21_ext_convergence.csv"),
                   index=False, encoding="utf-8-sig")
    cc = conv_df[conv_df["set_a"].isin(chem) & conv_df["set_b"].isin(chem)]
    summary["convergence"] = {
        "chemokine_pairs_min_rho": float(cc["spearman_rho"].min()),
        "chemokine_pairs_max_rho": float(cc["spearman_rho"].max()),
        "chemokine_vs_Tcell_rho": [
            {"set": r.set_a if r.set_a != "G7_Tcell_inflamed_cytolytic" else r.set_b,
             "rho": r.spearman_rho}
            for r in conv_df.itertuples()
            if "G7_Tcell_inflamed_cytolytic" in (r.set_a, r.set_b)],
    }

    # ================= B. programme correlations ========================
    prows = []
    for name in prog_df.columns:
        rho = stats.spearmanr(cts, prog_df[name]).statistic
        prows.append({"programme": name, "spearman_rho_with_CTS":
                      round(float(rho), 4)})
    prog_corr = pd.DataFrame(prows).sort_values(
        "spearman_rho_with_CTS", ascending=False)
    prog_corr.to_csv(os.path.join(RES, "21_ext_programme_correlation.csv"),
                     index=False, encoding="utf-8-sig")
    summary["programme_correlation"] = prog_corr.to_dict("records")

    # ================= C. survival ======================================
    d = out.dropna(subset=["os_days", "os_event", "CTS"]).copy()
    d = d[d["os_days"] > 0]
    time = d["os_days"].to_numpy(float)
    event = d["os_event"].to_numpy(int)
    score = d["CTS"].to_numpy(float)

    surv_rows = []
    # C1 discovery-derived cut-off
    hi = score > DISCOVERY_CUT
    chi2, p = _logrank(time[hi], event[hi], time[~hi], event[~hi])
    med_hi = _median_surv(time[hi], event[hi])
    med_lo = _median_surv(time[~hi], event[~hi])
    surv_rows.append({"analysis": "discovery cut-off %.4f" % DISCOVERY_CUT,
                      "n_high": int(hi.sum()), "n_low": int((~hi).sum()),
                      "chi2": round(chi2, 3), "p": round(p, 5),
                      "median_os_high_days": med_hi, "median_os_low_days": med_lo})
    print("[surv] discovery cut-off: chi2=%.2f p=%.4g  (hi n=%d, lo n=%d)"
          % (chi2, p, hi.sum(), (~hi).sum()))

    # C2 cohort-internal optimal cut-off
    cut, chi2b = best_cutpoint(score, time, event)
    if cut is not None:
        hb = score > cut
        _, pb = _logrank(time[hb], event[hb], time[~hb], event[~hb])
        surv_rows.append({"analysis": "cohort-optimal cut-off %.4f" % cut,
                          "n_high": int(hb.sum()), "n_low": int((~hb).sum()),
                          "chi2": round(chi2b, 3), "p": round(pb, 5),
                          "median_os_high_days": _median_surv(time[hb], event[hb]),
                          "median_os_low_days": _median_surv(time[~hb], event[~hb])})
        print("[surv] cohort-optimal cut=%.4f chi2=%.2f p=%.4g" % (cut, chi2b, pb))
        summary["ext_optimal_cutoff"] = float(cut)

    # C3 median cut-off
    med = float(np.median(score))
    hm = score > med
    chi2m, pm = _logrank(time[hm], event[hm], time[~hm], event[~hm])
    surv_rows.append({"analysis": "median cut-off %.4f" % med,
                      "n_high": int(hm.sum()), "n_low": int((~hm).sum()),
                      "chi2": round(chi2m, 3), "p": round(pm, 5),
                      "median_os_high_days": _median_surv(time[hm], event[hm]),
                      "median_os_low_days": _median_surv(time[~hm], event[~hm])})
    print("[surv] median cut-off: chi2=%.2f p=%.4g" % (chi2m, pm))

    # C4 continuous, and adjusted for stage/age
    hr, p_uni = _cox1(score, time, event)
    surv_rows.append({"analysis": "continuous CTS (univariable Cox)",
                      "hr": hr, "p": round(p_uni, 5)})
    print("[surv] univariable Cox HR=%.3f p=%.4g" % (hr, p_uni))

    cov = d[["stage_ord", "age"]].to_numpy(float)
    mask = ~np.isnan(cov).any(axis=1)
    hr_m, p_m = _cox_multi(score[mask], cov[mask], time[mask], event[mask])
    surv_rows.append({"analysis": "continuous CTS (adjusted stage + age)",
                      "hr": hr_m, "p": round(p_m, 5)})
    print("[surv] adjusted Cox HR=%.3f p=%.4g" % (hr_m, p_m))

    # C5 positive controls and individual signatures ----------------------
    # The pipeline itself must be able to detect *some* prognostic signal in
    # the validation cohort, otherwise a null CTS result is uninterpretable.
    # G7 (generic T cell-inflamed) is the pre-specified positive control.
    per_set = []
    for name in list(sc.columns):
        if name == "CTS":
            continue
        v = sc[name].reindex(d.index).to_numpy(float)
        if np.isnan(v).any():
            continue
        cut_i, chi2_i = best_cutpoint(v, time, event)
        hr_i, p_hr = _cox1(v, time, event)
        row = {"score": name, "hr_per_sd": round(hr_i, 4),
               "p_cox": round(p_hr, 5),
               "optimal_cut": round(cut_i, 4) if cut_i is not None else None,
               "chi2_at_cut": round(chi2_i, 3)}
        if cut_i is not None:
            hb = v > cut_i
            _, pb = _logrank(time[hb], event[hb], time[~hb], event[~hb])
            row["p_at_cut"] = round(pb, 5)
        per_set.append(row)
    per_set_df = pd.DataFrame(per_set).sort_values("p_cox")
    per_set_df.to_csv(os.path.join(RES, "21_ext_per_signature_survival.csv"),
                      index=False, encoding="utf-8-sig")
    summary["per_signature_survival"] = per_set
    print("\n[surv] per-signature cox / optimal cut-off (external cohort):")
    for r in per_set:
        print("    %-32s HR=%6.3f p=%.4g   cut=%s p_at_cut=%s"
              % (r["score"], r["hr_per_sd"], r["p_cox"], r["optimal_cut"],
                 r.get("p_at_cut", "-")))

    surv_df = pd.DataFrame(surv_rows)
    surv_df.to_csv(os.path.join(RES, "21_ext_survival.csv"),
                   index=False, encoding="utf-8-sig")
    summary["survival"] = surv_rows

    # ================= D. HPV association ===============================
    hpv = out.dropna(subset=["hpv_rna"])
    grp = hpv["hpv_rna"].replace({"NA": np.nan}).dropna()
    groups = {}
    for lab in ["DNA-", "DNA+RNA-", "DNA+RNA+"]:
        v = hpv.loc[grp.index][hpv.loc[grp.index, "hpv_rna"] == lab]["CTS"]
        if len(v):
            groups[lab] = v
    hpv_res = {"n_annotated": int(len(grp)),
               "groups": {k: {"n": int(len(v)), "median_cts": round(float(v.median()), 4)}
                          for k, v in groups.items()}}
    if "DNA-" in groups and "DNA+RNA+" in groups:
        u, pu = stats.mannwhitneyu(groups["DNA+RNA+"], groups["DNA-"],
                                   alternative="two-sided")
        hpv_res["mannwhitney_DNAposRNApos_vs_DNAneg"] = {
            "U": float(u), "p": float(pu)}
        print("[hpv] DNA+RNA+ median CTS=%.4f vs DNA- %.4f  p=%.4g"
              % (groups["DNA+RNA+"].median(), groups["DNA-"].median(), pu))
    # binary: any HPV DNA (consistent with the discovery-cohort definition)
    hpv["hpv_bin"] = np.where(hpv["hpv_dna"].isin(["HPV16", "Other HPV"]),
                              "positive",
                              np.where(hpv["hpv_dna"] == "Negative", "negative", ""))
    b = hpv[hpv["hpv_bin"] != ""]
    gpos = b[b["hpv_bin"] == "positive"]["CTS"]
    gneg = b[b["hpv_bin"] == "negative"]["CTS"]
    if len(gpos) and len(gneg):
        u2, p2 = stats.mannwhitneyu(gpos, gneg, alternative="two-sided")
        hpv_res["dna_binary"] = {"n_pos": int(len(gpos)), "n_neg": int(len(gneg)),
                                 "median_pos": round(float(gpos.median()), 4),
                                 "median_neg": round(float(gneg.median()), 4),
                                 "p": float(p2)}
        print("[hpv] DNA-positive median CTS=%.4f (n=%d) vs negative %.4f (n=%d) p=%.4g"
              % (gpos.median(), len(gpos), gneg.median(), len(gneg), p2))
    # HPV-stratified survival
    strat = []
    for lab in ["DNA-", "DNA+RNA-", "DNA+RNA+"]:
        sd = out[out["hpv_rna"] == lab]
        if len(sd) < 20:
            continue
        st = sd["os_days"].to_numpy(float)
        se = sd["os_event"].to_numpy(int)
        ss = sd["CTS"].to_numpy(float)
        h = ss > DISCOVERY_CUT
        if h.sum() >= 5 and (~h).sum() >= 5:
            chi2s, ps = _logrank(st[h], se[h], st[~h], se[~h])
            strat.append({"stratum": lab, "n": int(len(sd)),
                          "n_cts_high": int(h.sum()), "n_cts_low": int((~h).sum()),
                          "chi2": round(chi2s, 3), "p": round(ps, 5)})
            print("[hpv-strat] %-10s n=%3d  high=%3d low=%3d  p=%.4g"
                  % (lab, len(sd), h.sum(), (~h).sum(), ps))
        else:
            strat.append({"stratum": lab, "n": int(len(sd)),
                          "n_cts_high": int(h.sum()), "n_cts_low": int((~h).sum()),
                          "note": "too few CTS-low for a stratified comparison"})
            print("[hpv-strat] %-10s n=%3d  high=%3d low=%3d  -> not evaluable"
                  % (lab, len(sd), h.sum(), (~h).sum()))
    hpv_res["stratified"] = strat
    summary["hpv"] = hpv_res

    # ================= E. published immune cluster ======================
    cl = out.dropna(subset=["cluster"])
    byc = cl.groupby("cluster")["CTS"].agg(["count", "median"]).round(4)
    summary["cluster_cts"] = byc.reset_index().to_dict("records")
    print("[cluster] CTS by published expression cluster:")
    for r in summary["cluster_cts"]:
        print("    %-20s n=%3d  median CTS=%.4f"
              % (r["cluster"], r["count"], r["median"]))
    # is the immune-response cluster enriched for high CTS?
    grps = [g["CTS"].to_numpy(float) for _, g in cl.groupby("cluster")]
    if len(grps) > 2:
        kw = stats.kruskal(*grps)
        top = byc["median"].idxmax()
        summary["cluster_kruskal"] = {"H": float(kw.statistic),
                                      "p": float(kw.pvalue),
                                      "highest_median_cluster": top}
        print("[cluster] Kruskal-Wallis H=%.2f p=%.4g ; highest median = %s"
              % (kw.statistic, kw.pvalue, top))

    # ================= F. figures =======================================
    step_grid = np.unique(np.concatenate([time[event == 1], [time.max()]]))

    def km(t, e):
        ts = np.unique(t[e == 1])
        surv, at = 1.0, []
        for tt in ts:
            n = (t >= tt).sum()
            dd = ((t == tt) & (e == 1)).sum()
            if n > 0:
                surv *= (1 - dd / n)
            at.append((tt, surv))
        return np.array([0.0] + [a for a, _ in at]), np.array([1.0] + [b for _, b in at])

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.rcParams["font.sans-serif"] = ["DejaVu Sans"]
        plt.rcParams["axes.unicode_minus"] = False

        # Fig: KM for CTS (null) vs the positive control G7, side by side
        ctrl = "G7_Tcell_inflamed_cytolytic"
        ctrl_row = next((r for r in per_set if r["score"] == ctrl), None)
        panels = [("Consensus TLS score (CTS)", score,
                   DISCOVERY_CUT, "%.3f" % DISCOVERY_CUT, p)]
        if ctrl_row and ctrl_row.get("optimal_cut") is not None and ctrl in sc.columns:
            cv = sc[ctrl].reindex(d.index).to_numpy(float)
            panels.append(("T cell-inflamed signature (positive control)", cv,
                           ctrl_row["optimal_cut"],
                           "%.3f" % ctrl_row["optimal_cut"],
                           ctrl_row.get("p_at_cut", np.nan)))

        fig, axes = plt.subplots(1, len(panels), figsize=(5.4 * len(panels), 4.9),
                                 dpi=300, squeeze=False)
        for ax, (title, vec, cutv, cutlab, pv) in zip(axes[0], panels):
            m = vec > cutv
            for mask, lab, col in [(m, "high (n=%d)" % m.sum(), "#A32D2D"),
                                   (~m, "low (n=%d)" % (~m).sum(), "#0F6E56")]:
                xs, ys = km(time[mask], event[mask])
                ax.step(xs, ys, where="post", color=col, lw=1.8, label=lab)
            ax.set_xlabel("Overall survival (days)")
            ax.set_ylabel("Survival probability")
            ax.set_title("%s\ncut-off %s" % (title, cutlab), fontsize=9)
            ax.text(0.98, 0.95, "log-rank p = %.3f" % pv, transform=ax.transAxes,
                    ha="right", va="top", fontsize=9)
            ax.legend(frameon=False, fontsize=8, loc="lower left")
            ax.set_ylim(0, 1.02)
        fig.suptitle("GSE65858 external validation (n = %d, %d deaths)"
                     % (len(d), int(event.sum())), fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        fig.savefig(os.path.join(RES, "fig_ext_km_cts.png"))
        plt.close(fig)

        # Fig: signature convergence heatmap
        sets_all = [c for c in list(sc.columns) if c != "CTS"]
        M = np.ones((len(sets_all), len(sets_all)))
        for i, a in enumerate(sets_all):
            for j, b in enumerate(sets_all):
                M[i, j] = stats.spearmanr(sc[a], sc[b]).statistic
        fig, ax = plt.subplots(figsize=(7.6, 6.4), dpi=300)
        im = ax.imshow(M, cmap="RdYlBu_r", vmin=0.5, vmax=1.0)
        lab = [s.split("_", 1)[1] if "_" in s else s for s in sets_all]
        ax.set_xticks(range(len(lab)))
        ax.set_yticks(range(len(lab)))
        ax.set_xticklabels(lab, rotation=40, ha="right", fontsize=7)
        ax.set_yticklabels(lab, fontsize=7)
        for i in range(len(lab)):
            for j in range(len(lab)):
                ax.text(j, i, "%.2f" % M[i, j], ha="center", va="center",
                        fontsize=6.2,
                        color="white" if M[i, j] > 0.88 else "black")
        fig.colorbar(im, ax=ax, shrink=0.8, label="Spearman rho")
        ax.set_title("GSE65858: gene set score convergence", fontsize=10)
        fig.tight_layout()
        fig.savefig(os.path.join(RES, "fig_ext_convergence.png"))
        plt.close(fig)

        # Fig: CTS by HPV group
        order = [k for k in ["DNA-", "DNA+RNA-", "DNA+RNA+"] if k in groups]
        if order:
            fig, ax = plt.subplots(figsize=(5.4, 4.6), dpi=300)
            data = [groups[k].values for k in order]
            bp = ax.boxplot(data, patch_artist=True, showfliers=False,
                            widths=0.55)
            ax.set_xticks(range(1, len(order) + 1))
            ax.set_xticklabels(order)
            for patch, col in zip(bp["boxes"], ["#4A6FA5", "#C58C3A", "#A32D2D"]):
                patch.set_facecolor(col)
                patch.set_alpha(0.75)
            for i, k in enumerate(order):
                y = groups[k].values
                x = np.random.default_rng(0).normal(i + 1, 0.06, len(y))
                ax.scatter(x, y, s=5, color="#333333", alpha=0.4, zorder=3)
            ax.axhline(DISCOVERY_CUT, color="grey", ls="--", lw=1)
            ax.text(0.02, DISCOVERY_CUT, " discovery cut-off",
                    transform=ax.get_yaxis_transform(), fontsize=7, va="bottom",
                    color="grey")
            ax.set_ylabel("Consensus TLS score")
            ax.set_xlabel("HPV16 status")
            pp = hpv_res.get("mannwhitney_DNAposRNApos_vs_DNAneg", {}).get("p")
            if pp is not None:
                ax.set_title("GSE65858: CTS by HPV16 RNA status\n"
                             "DNA+RNA+ vs DNA-  p = %.4g" % pp, fontsize=10)
            fig.tight_layout()
            fig.savefig(os.path.join(RES, "fig_ext_cts_by_hpv.png"))
            plt.close(fig)
        figs_ok = True
    except Exception as exc:                                  # noqa: BLE001
        print("[fig] skipped: %s: %s" % (type(exc).__name__, exc))
        figs_ok = False
    summary["figures_written"] = figs_ok

    with io.open(os.path.join(RES, "21_ext_summary.json"), "w",
                 encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
    print("\n[ok] wrote 21_ext_* into %s" % RES)
    return 0


def _median_surv(t, e):
    ts = np.unique(t[e == 1])
    surv = 1.0
    for tt in ts:
        n = (t >= tt).sum()
        dd = ((t == tt) & (e == 1)).sum()
        if n > 0:
            surv *= (1 - dd / n)
        if surv <= 0.5:
            return float(tt)
    return None


def _cox1(score, time, event):
    """Single-covariate Cox via Newton-Raphson (no lifelines dependency)."""
    return _cox_multi(score.reshape(-1, 1), None, time, event, single=True)


def _cox_multi(score, cov, time, event, single=False):
    X = score.reshape(-1, 1) if single else np.column_stack([score, cov])
    X = X - X.mean(axis=0)
    sd = X.std(axis=0)
    sd[sd == 0] = 1.0
    X = X / sd
    n, k = X.shape
    beta = np.zeros(k)
    for _ in range(60):
        eta = X @ beta
        eta -= eta.max()
        w = np.exp(eta)
        order = np.argsort(-time)
        Xo, to, eo, wo = X[order], time[order], event[order], w[order]
        # risk-set cumulative sums
        W = np.cumsum(wo)
        WX = np.cumsum(Xo * wo[:, None], axis=0)
        WXX = np.cumsum(Xo[:, None, :] * Xo[:, :, None] * wo[:, None, None], axis=0)
        grad = np.zeros(k)
        hess = np.zeros((k, k))
        for i in range(n):
            if eo[i] != 1:
                continue
            s0 = W[i]
            s1 = WX[i] / s0
            s2 = WXX[i] / s0 - np.outer(s1, s1)
            grad += Xo[i] - s1
            hess += s2
        try:
            delta = np.linalg.solve(hess, grad)
        except np.linalg.LinAlgError:
            break
        beta = beta + delta
        if np.max(np.abs(delta)) < 1e-8:
            break
    # covariance = inverse of observed information
    eta = X @ beta
    eta -= eta.max()
    w = np.exp(eta)
    order = np.argsort(-time)
    Xo, to, eo, wo = X[order], time[order], event[order], w[order]
    W = np.cumsum(wo)
    WX = np.cumsum(Xo * wo[:, None], axis=0)
    WXX = np.cumsum(Xo[:, None, :] * Xo[:, :, None] * wo[:, None, None], axis=0)
    hess = np.zeros((k, k))
    for i in range(n):
        if eo[i] != 1:
            continue
        s0 = W[i]
        s1 = WX[i] / s0
        s2 = WXX[i] / s0 - np.outer(s1, s1)
        hess += s2
    try:
        covb = np.linalg.inv(hess)
    except np.linalg.LinAlgError:
        return float(np.nan), float(np.nan)
    se = float(np.sqrt(covb[0, 0]))
    zstat = beta[0] / se if se > 0 else np.nan
    pval = float(2 * stats.norm.sf(abs(zstat)))
    return float(np.exp(beta[0])), pval


if __name__ == "__main__":
    sys.exit(main())
