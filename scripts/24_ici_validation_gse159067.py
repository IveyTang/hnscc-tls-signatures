# -*- coding: utf-8 -*-
"""
24_ici_validation_gse159067.py
------------------------------
Immune-checkpoint-inhibitor (ICI) validation of the consensus TLS score (CTS)
in **GSE159067**, using the identical gene sets and scoring algorithm applied
in the discovery (TCGA-HNSC) and prognostic-validation (GSE65858) cohorts,
with no re-tuning.

Cohort
~~~~~~
102 patients with advanced HNSCC treated with PD-1/PD-L1 inhibitors.
  RECIST best response : PD 64, SD 27, PR 6, CR 5  (objective response 11/102)
  PFS / OS             : available with events (96 / 92)
  HPV status           : 0 = negative (49), 1 = positive (11), 2 = unknown (42)
  Published HOT score  : a continuous immune "hot phenotype" score with a
                         HOT/COLD call (47 / 55) from the original study.

Why this cohort matters for this manuscript
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The paper argues that chemokine-based TLS signatures largely measure general
immune inflammation. GSE159067 provides an **independent, externally derived
immune-phenotype label** (the published HOT score) plus an immunotherapy
response endpoint, so it tests two things at once:
  (i)  does CTS track an independently published immune-inflammation score?
  (ii) does CTS predict benefit from checkpoint blockade?

Platform caveat
~~~~~~~~~~~~~~~
GPL18573 is the HTG EdgeSeq Oncology Biomarker Panel (2,567 genes), so
coverage is partial. Of the 81 genes across the seven sets, 70 are present
(86%), including 12/12 and 11/11 for the Coppola and HCC chemokine signatures
and 12/13 for the HNSCC 13-chemokine signature. The three chemokine
signatures that are the subject of this paper are therefore essentially fully
represented. The plasma-cell programme is the one weak spot (2 of 7 markers).

Outputs
~~~~~~~
    results/24_ici_scores.csv
    results/24_ici_summary.json
    results/24_ici_response.csv
    results/24_ici_survival.csv
    results/fig_ici_validation.png
    results/table7_ici_validation.csv
"""
from __future__ import annotations

import gzip
import io
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
EXT = os.path.join(ROOT, "data", "external", "GSE159067")
RES = os.path.join(ROOT, "results")
os.makedirs(RES, exist_ok=True)

MATRIX = os.path.join(EXT, "GSE159067_IHN_log2cpm_data.txt.gz")
SERIES = os.path.join(EXT, "GSE159067_series_matrix.txt.gz")
GENESETS = os.path.join(ROOT, "config", "tls_genesets.json")

CONSENSUS_SETS = ["G1_Coppola_12chemokine", "G2_11chemokine_HCC",
                  "G3_HNSCC_13chemokine", "G4_TLS_maturation_GC",
                  "G6_HEV_lymphoid_stroma"]
CHEMOKINE_SETS = ["G1_Coppola_12chemokine", "G2_11chemokine_HCC",
                  "G3_HNSCC_13chemokine"]
CTRL = "G7_Tcell_inflamed_cytolytic"
DISCOVERY_CUT = 0.2489

MARKERS = {
    "B cells":         ["MS4A1", "CD19", "CD79A", "CD79B", "TCL1A", "BANK1"],
    "GC B cells":      ["AICDA", "BCL6", "MEF2B", "RGS13", "MKI67"],
    "Plasma cells":    ["MZB1", "JCHAIN", "IGHG1", "IGKC", "XBP1", "PRDM1", "SDC1"],
    "Tfh":             ["CXCL13", "IL21", "BCL6", "ICOS", "PDCD1", "SH2D1A"],
    "CD8 T cells":     ["CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "NKG7"],
    "Treg":            ["FOXP3", "IL2RA", "CTLA4", "IKZF2"],
    "CAF":             ["ACTA2", "FAP", "COL1A1", "COL1A2", "PDGFRB", "TAGLN", "THY1"],
    "Endothelium/HEV": ["PECAM1", "VWF", "CLDN5", "PDPN", "CD34"],
    "Macrophages":     ["CD68", "CD163", "CSF1R", "MRC1"],
    "IFN-gamma":       ["IFNG", "STAT1", "IDO1", "CXCL9", "CXCL10", "GBP1", "IRF1"],
    "Checkpoint":      ["CD274", "PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT"],
    "Epithelial":      ["EPCAM", "KRT5", "KRT14", "KRT13", "KRT4"],
}


def load_gene_sets():
    d = json.load(io.open(GENESETS, encoding="utf-8"))
    return {k: [g.strip().upper() for g in v["genes"]]
            for k, v in d.items()
            if isinstance(v, dict) and isinstance(v.get("genes"), list)}


def load_expression():
    with gzip.open(MATRIX, "rt", encoding="utf-8", errors="replace") as fh:
        df = pd.read_csv(fh, sep="\t", index_col=0)
    df.index = [str(i).strip().upper() for i in df.index]
    df = df[~df.index.duplicated(keep="first")]
    # drop control probes reported alongside the panel
    df = df[~df.index.str.startswith(("NEG_CTRL", "POS_CTRL", "HK_"))]
    return df


def load_clinical():
    """Parse the GEO series matrix.

    All `!Sample_*` rows share one column order, but the rows are not emitted
    in a fixed sequence (characteristics come before description), so the
    table has to be built **positionally** rather than line by line.
    """
    blocks: dict[str, list[list[str]]] = {}
    with gzip.open(SERIES, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith("!Sample_"):
                continue
            p = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")]
            blocks.setdefault(p[0], []).append(p[1:])

    key_head = ("!Sample_description"
                if "!Sample_description" in blocks else "!Sample_geo_accession")
    keys = blocks[key_head][0]
    n = len(keys)
    rows = [{"sample": k} for k in keys]

    for head, store in (("!Sample_geo_accession", "gsm"),
                        ("!Sample_title", "title"),
                        ("!Sample_source_name_ch1", "source")):
        if head in blocks:
            for i, v in enumerate(blocks[head][0][:n]):
                rows[i][store] = v

    for vals in blocks.get("!Sample_characteristics_ch1", []):
        for i, v in enumerate(vals[:n]):
            if ":" in v:
                f, val = v.split(":", 1)
                rows[i][f.strip()] = val.strip()

    out = []
    for d in rows:
        out.append({
            "sample": d["sample"],
            "gsm": d.get("gsm"),
            "site": d.get("site"),
            "age": _f(d.get("age")),
            "sex": d.get("Sex"),
            "tobacco": _f(d.get("tobacco status")),
            "alcohol": _f(d.get("alcohol status")),
            "hpv": _f(d.get("hpv status (0 or 1 or 2 if unknown)")),
            "line": _f(d.get("immunotherapy line")),
            "response": d.get("best response on immunotherapy (recist)"),
            "os_m": _f(d.get("os (month)")),
            "os_e": _f(d.get("os (event)")),
            "pfs_m": _f(d.get("pfs (month)")),
            "pfs_e": _f(d.get("pfs (event)")),
            "hot_phenotype": _f(d.get("hot phenotype")),
            "hot_score": d.get("hot score"),
        })
    return pd.DataFrame(out).set_index("sample")


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


def zscore(mat):
    mu = mat.mean(axis=1)
    sd = mat.std(axis=1).replace(0, np.nan)
    return mat.sub(mu, axis=0).div(sd, axis=0)


def zmean(z, genes):
    present = [g for g in genes if g in z.index]
    if not present:
        return None, 0
    return z.loc[present].mean(axis=0), len(present)


def logit(y, x):
    """Simple logistic regression with one standardised covariate."""
    x = (x - x.mean()) / x.std()
    X = np.column_stack([np.ones_like(x), x])
    beta = np.zeros(2)
    for _ in range(80):
        p = 1 / (1 + np.exp(-(X @ beta)))
        g = X.T @ (y - p)
        W = p * (1 - p) + 1e-9
        H = X.T @ (X * W[:, None])
        try:
            step = np.linalg.solve(H, g)
        except np.linalg.LinAlgError:
            break
        beta += step
        if np.max(np.abs(step)) < 1e-9:
            break
    p = 1 / (1 + np.exp(-(X @ beta)))
    W = p * (1 - p) + 1e-9
    H = X.T @ (X * W[:, None])
    try:
        se = np.sqrt(np.linalg.inv(H).diagonal())
    except np.linalg.LinAlgError:
        se = np.array([np.nan, np.nan])
    z = beta[1] / se[1] if se[1] > 0 else np.nan
    return float(np.exp(beta[1])), float(2 * stats.norm.sf(abs(z)))


def logrank(t1, e1, t2, e2):
    times = np.unique(np.concatenate([t1, t2]))
    O = E = V = 0.0
    for tt in times:
        n1, n2 = (t1 >= tt).sum(), (t2 >= tt).sum()
        n = n1 + n2
        d1 = ((t1 == tt) & (e1 == 1)).sum()
        d2 = ((t2 == tt) & (e2 == 1)).sum()
        dd = d1 + d2
        if n <= 1 or dd == 0:
            continue
        O += d1
        E += dd * n1 / n
        V += dd * (n1 / n) * (1 - n1 / n) * (n - dd) / (n - 1)
    if V <= 0:
        return 0.0, 1.0
    chi2 = (O - E) ** 2 / V
    return float(chi2), float(stats.chi2.sf(chi2, 1))


def cox1(time, event, x):
    x = (x - x.mean()) / x.std()
    X = np.column_stack([x])
    beta = np.zeros(1)
    for _ in range(60):
        eta = X @ beta
        eta -= eta.max()
        w = np.exp(eta)
        order = np.argsort(-time)
        Xo, to, eo, wo = X[order], time[order], event[order], w[order]
        W = np.cumsum(wo)
        WX = np.cumsum(Xo * wo[:, None], axis=0)
        WXX = np.cumsum(Xo[:, None, :] * Xo[:, :, None] * wo[:, None, None], axis=0)
        g = np.zeros(1)
        H = np.zeros((1, 1))
        for i in range(len(to)):
            if eo[i] != 1:
                continue
            s1 = WX[i] / W[i]
            s2 = WXX[i] / W[i] - np.outer(s1, s1)
            g += Xo[i] - s1
            H += s2
        try:
            d = np.linalg.solve(H, g)
        except np.linalg.LinAlgError:
            return np.nan, np.nan
        beta += d
        if np.max(np.abs(d)) < 1e-9:
            break
    try:
        se = float(np.sqrt(np.linalg.inv(H)[0, 0]))
    except np.linalg.LinAlgError:
        return float(np.exp(beta[0])), np.nan
    z = beta[0] / se if se > 0 else np.nan
    return float(np.exp(beta[0])), float(2 * stats.norm.sf(abs(z)))


def auc(score, label):
    score = np.asarray(score, float)
    label = np.asarray(label, int)
    r = stats.rankdata(score)
    n1 = label.sum()
    n0 = len(label) - n1
    if n1 == 0 or n0 == 0:
        return np.nan
    return float((r[label == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def main():
    print("=" * 76)
    print("GSE159067 - ICI validation (direct GEO retrieval)")
    print("=" * 76)

    gs = load_gene_sets()
    expr = load_expression()
    clin = load_clinical()
    print("[expr] %d genes x %d samples" % expr.shape)
    print("[clin] %d samples" % len(clin))

    common = [s for s in expr.columns if s in clin.index]
    expr = expr[common]
    clin = clin.loc[common]
    print("[merge] %d samples with expression + annotation" % len(common))

    z = zscore(expr)
    scores, cover = {}, {}
    for name, genes in gs.items():
        s, n = zmean(z, genes)
        if s is not None:
            scores[name] = s
            cover[name] = "%d/%d" % (n, len(genes))
    sc = pd.DataFrame(scores)
    print("[score] coverage: %s" % json.dumps(cover))

    have = [c for c in CONSENSUS_SETS if c in sc.columns]
    sc["CTS"] = sc[have].rank(pct=True).mean(axis=1)
    prog = {}
    for name, genes in MARKERS.items():
        s, n = zmean(z, genes)
        if s is not None and n >= 2:
            prog[name] = s
    prog_df = pd.DataFrame(prog)

    out = clin.join(sc).join(prog_df.add_prefix("prog_"))
    out.to_csv(os.path.join(RES, "24_ici_scores.csv"), encoding="utf-8-sig")

    summary = {
        "cohort": "GSE159067 (advanced HNSCC on PD-1/PD-L1 inhibitors)",
        "n": int(len(out)),
        "geneset_coverage": cover,
        "n_panel_genes": int(expr.shape[0]),
        "discovery_cutoff": DISCOVERY_CUT,
    }

    # ---------- A. convergence (replication of the structural claim) ----
    chem = [c for c in CHEMOKINE_SETS if c in sc.columns]
    pairs = [(a, b) for i, a in enumerate(chem) for b in chem[i + 1:]]
    conv = [{"set_a": a, "set_b": b,
             "spearman_rho": round(float(stats.spearmanr(sc[a], sc[b]).statistic), 4)}
            for a, b in pairs]
    if CTRL in sc.columns:
        for c in chem:
            conv.append({"set_a": c, "set_b": CTRL,
                         "spearman_rho": round(float(stats.spearmanr(sc[c], sc[CTRL]).statistic), 4)})
    pd.DataFrame(conv).to_csv(os.path.join(RES, "24_ici_convergence.csv"),
                              index=False, encoding="utf-8-sig")
    cc = [r["spearman_rho"] for r in conv if CTRL not in (r["set_a"], r["set_b"])]
    ctrl_rho = [r["spearman_rho"] for r in conv if CTRL in (r["set_a"], r["set_b"])]
    summary["convergence"] = {
        "chemokine_pairs_min_rho": min(cc) if cc else None,
        "chemokine_pairs_max_rho": max(cc) if cc else None,
        "chemokine_vs_Tcell_rho": ctrl_rho,
    }

    # ---------- B. programme correlations -------------------------------
    prows = [{"programme": n,
              "spearman_rho_with_CTS": round(float(stats.spearmanr(out["CTS"], out["prog_" + n]).statistic), 4)}
             for n in prog_df.columns]
    pc = pd.DataFrame(prows).sort_values("spearman_rho_with_CTS", ascending=False)
    pc.to_csv(os.path.join(RES, "24_ici_programme_correlation.csv"),
              index=False, encoding="utf-8-sig")
    summary["programme_correlation"] = pc.to_dict("records")

    # ---------- C. association with the published HOT score -------------
    hp = out.dropna(subset=["hot_phenotype"])
    rho_hot = float(stats.spearmanr(hp["CTS"], hp["hot_phenotype"]).statistic)
    hot = out.copy()
    hot["hot_bin"] = (hot["hot_score"] == "HOT").astype(int)
    a_hot = auc(hot["CTS"], hot["hot_bin"])
    a_ctrl = auc(hot[CTRL], hot["hot_bin"]) if CTRL in hot.columns else np.nan
    summary["hot_phenotype"] = {
        "n_with_score": int(len(hp)),
        "spearman_rho_CTS_vs_hot_phenotype": round(rho_hot, 4),
        "auc_CTS_for_HOT_call": round(float(a_hot), 4),
        "auc_Tcell_inflamed_for_HOT_call": round(float(a_ctrl), 4) if a_ctrl == a_ctrl else None,
    }
    print("[hot] CTS vs published hot phenotype rho=%.3f ; AUC(HOT)=%.3f (T-cell control %.3f)"
          % (rho_hot, a_hot, a_ctrl))

    # ---------- D. response to checkpoint blockade ----------------------
    d = out.dropna(subset=["response"]).copy()
    d["orr"] = d["response"].isin(["CR", "PR"]).astype(int)
    d["dcr"] = d["response"].isin(["CR", "PR", "SD"]).astype(int)
    print("[resp] n=%d  ORR=%d (%.1f%%)  DCR=%d (%.1f%%)"
          % (len(d), d["orr"].sum(), 100 * d["orr"].mean(),
             d["dcr"].sum(), 100 * d["dcr"].mean()))

    rrows = []
    for label, vec in [("CTS", d["CTS"])] + \
                      [(c, d[c]) for c in list(sc.columns) if c != "CTS"]:
        for ep in ("orr", "dcr"):
            orr, p = logit(d[ep].to_numpy(float), vec.to_numpy(float))
            rrows.append({"score": label, "endpoint": ep.upper(),
                          "odds_ratio_per_sd": round(orr, 4),
                          "p": round(p, 5),
                          "auc": round(auc(vec, d[ep]), 4)})
    rd = pd.DataFrame(rrows).sort_values(["endpoint", "p"])
    rd.to_csv(os.path.join(RES, "24_ici_response.csv"), index=False, encoding="utf-8-sig")
    summary["response"] = {
        "n": int(len(d)), "orr": int(d["orr"].sum()), "dcr": int(d["dcr"].sum()),
        "details": rd[rd["score"].isin(["CTS", CTRL])].to_dict("records"),
        "all": rrows,
    }
    for r in rd[rd["score"].isin(["CTS", CTRL])].to_dict("records"):
        print("[resp] %-28s %-4s OR=%.3f p=%.4g AUC=%.3f"
              % (r["score"], r["endpoint"], r["odds_ratio_per_sd"], r["p"], r["auc"]))

    # ---------- E. survival on immunotherapy ----------------------------
    srows = []
    for ep, tc, ec in [("PFS", "pfs_m", "pfs_e"), ("OS", "os_m", "os_e")]:
        sd = out.dropna(subset=[tc, ec, "CTS"])
        sd = sd[sd[tc] > 0]
        hr, p = cox1(sd[tc].to_numpy(float), sd[ec].to_numpy(int), sd["CTS"].to_numpy(float))
        srows.append({"endpoint": ep, "score": "CTS", "n": int(len(sd)),
                      "events": int(sd[ec].sum()), "hr_per_sd": round(hr, 4), "p": round(p, 5)})
        if CTRL in sd.columns:
            hr2, p2 = cox1(sd[tc].to_numpy(float), sd[ec].to_numpy(int), sd[CTRL].to_numpy(float))
            srows.append({"endpoint": ep, "score": CTRL, "n": int(len(sd)),
                          "events": int(sd[ec].sum()), "hr_per_sd": round(hr2, 4), "p": round(p2, 5)})
        hi = sd["CTS"] > DISCOVERY_CUT
        if hi.sum() >= 5 and (~hi).sum() >= 5:
            chi2, pl = logrank(sd.loc[hi, tc].to_numpy(float), sd.loc[hi, ec].to_numpy(int),
                               sd.loc[~hi, tc].to_numpy(float), sd.loc[~hi, ec].to_numpy(int))
            srows.append({"endpoint": ep, "score": "CTS (discovery cut-off)",
                          "n": int(len(sd)), "events": int(sd[ec].sum()),
                          "chi2": round(chi2, 3), "p": round(pl, 5)})
    sdf = pd.DataFrame(srows)
    sdf.to_csv(os.path.join(RES, "24_ici_survival.csv"), index=False, encoding="utf-8-sig")
    summary["survival"] = srows
    print()
    for r in srows:
        print("[surv] %-4s %-28s HR=%s p=%s"
              % (r["endpoint"], r["score"], r.get("hr_per_sd", "-"), r["p"]))

    # ---------- F. HPV association --------------------------------------
    hv = out.dropna(subset=["hpv"])
    g0 = hv[hv["hpv"] == 0]["CTS"]
    g1 = hv[hv["hpv"] == 1]["CTS"]
    if len(g0) and len(g1):
        u, p = stats.mannwhitneyu(g1, g0, alternative="two-sided")
        summary["hpv"] = {"n_neg": int(len(g0)), "n_pos": int(len(g1)),
                          "median_cts_neg": round(float(g0.median()), 4),
                          "median_cts_pos": round(float(g1.median()), 4),
                          "p": float(p)}
        print("[hpv] median CTS pos=%.4f (n=%d) vs neg=%.4f (n=%d) p=%.4g"
              % (g1.median(), len(g1), g0.median(), len(g0), p))

    # ---------- G. figure ------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.rcParams["font.sans-serif"] = ["DejaVu Sans"]
        plt.rcParams["axes.unicode_minus"] = False
        fig, axes = plt.subplots(1, 3, figsize=(11.6, 4.3), dpi=300)

        ax = axes[0]
        ax.scatter(hp["CTS"], hp["hot_phenotype"], s=14, alpha=0.7, color="#1f3b73")
        ax.set_xlabel("Consensus TLS score", fontsize=8)
        ax.set_ylabel("Published HOT phenotype score", fontsize=8)
        ax.set_title("A  CTS vs published hot phenotype\nSpearman $\\rho$ = %.3f" % rho_hot, fontsize=9, loc="left")
        ax.tick_params(labelsize=7)

        ax = axes[1]
        order = ["CR", "PR", "SD", "PD"]
        data = [d.loc[d["response"] == o, "CTS"].dropna().values for o in order if (d["response"] == o).sum()]
        labs = [o for o in order if (d["response"] == o).sum()]
        bp = ax.boxplot(data, patch_artist=True, showfliers=False, widths=0.55)
        ax.set_xticks(range(1, len(labs) + 1))
        ax.set_xticklabels(labs)
        for b in bp["boxes"]:
            b.set_facecolor("#4A6FA5"); b.set_alpha(0.75)
        ax.set_xlabel("RECIST best response", fontsize=8)
        ax.set_ylabel("Consensus TLS score", fontsize=8)
        rr = next(r for r in rrows if r["score"] == "CTS" and r["endpoint"] == "ORR")
        ax.set_title("B  CTS by response\nORR OR = %.2f per SD, p = %.3f" % (rr["odds_ratio_per_sd"], rr["p"]),
                     fontsize=9, loc="left")
        ax.tick_params(labelsize=7)

        ax = axes[2]
        sets_all = [c for c in sc.columns if c != "CTS"]
        M = np.array([[stats.spearmanr(sc[a], sc[b]).statistic for b in sets_all] for a in sets_all])
        im = ax.imshow(M, cmap="RdYlBu_r", vmin=0.3, vmax=1.0)
        lab = [s.split("_", 1)[1] for s in sets_all]
        ax.set_xticks(range(len(lab))); ax.set_yticks(range(len(lab)))
        ax.set_xticklabels(lab, rotation=45, ha="right", fontsize=6)
        ax.set_yticklabels(lab, fontsize=6)
        for i in range(len(lab)):
            for j in range(len(lab)):
                ax.text(j, i, "%.2f" % M[i, j], ha="center", va="center", fontsize=5.2,
                        color="white" if M[i, j] > 0.88 else "black")
        ax.set_title("C  Score convergence in the ICI cohort", fontsize=9, loc="left")
        fig.colorbar(im, ax=ax, shrink=0.8)

        fig.suptitle("GSE159067: ICI-treated HNSCC (n = %d) - external validation" % len(out), fontsize=10)
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        fig.savefig(os.path.join(RES, "fig_ici_validation.png"))
        plt.close(fig)
        summary["figure"] = True
    except Exception as exc:                                   # noqa: BLE001
        print("[fig] skipped: %s: %s" % (type(exc).__name__, exc))
        summary["figure"] = False

    # ---------- H. summary table ----------------------------------------
    rows = [
        ["Cohort", "GSE159067 - advanced HNSCC on PD-1/PD-L1 inhibitors"],
        ["Patients", str(len(out))],
        ["Panel", "HTG EdgeSeq Oncology Biomarker Panel (2,567 genes)"],
        ["Gene set coverage", "; ".join("%s %s" % (k.split("_")[0], v) for k, v in cover.items())],
        ["Chemokine pairs, rho (min-max)",
         "%.3f-%.3f" % (min(cc), max(cc)) if cc else "n/a"],
        ["Chemokine vs T cell-inflamed, rho",
         "%.3f-%.3f" % (min(ctrl_rho), max(ctrl_rho)) if ctrl_rho else "n/a"],
        ["CTS vs published hot phenotype, rho", "%.3f" % rho_hot],
        ["CTS AUC for published HOT call", "%.3f" % a_hot],
        ["Objective response rate", "%d/%d (%.1f%%)" % (d["orr"].sum(), len(d), 100 * d["orr"].mean())],
        ["Disease control rate", "%d/%d (%.1f%%)" % (d["dcr"].sum(), len(d), 100 * d["dcr"].mean())],
    ]
    for r in rrows:
        if r["score"] == "CTS":
            rows.append(["CTS, %s (OR per SD)" % r["endpoint"], "%.2f (p = %.3f)" % (r["odds_ratio_per_sd"], r["p"])])
    for r in srows:
        if r["score"] == "CTS":
            rows.append(["CTS, %s (HR per SD)" % r["endpoint"], "%.3f (p = %.3f)" % (r["hr_per_sd"], r["p"])])
    if "hpv" in summary:
        h = summary["hpv"]
        rows.append(["Median CTS, HPV+ vs HPV-", "%.3f vs %.3f (p = %.4g)" % (h["median_cts_pos"], h["median_cts_neg"], h["p"])])
    t7 = pd.DataFrame(rows, columns=["Item", "Value"])
    t7.to_csv(os.path.join(RES, "table7_ici_validation.csv"), index=False, encoding="utf-8-sig")
    for r in rows:
        print("   %-44s %s" % (r[0], r[1]))

    with io.open(os.path.join(RES, "24_ici_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
    print("\n[ok] wrote 24_ici_* and table7 into %s" % RES)
    return 0


if __name__ == "__main__":
    sys.exit(main())
