"""
03b_cutpoint_stability.py
-------------------------
Quantify how stable the "optimal" CTS cut-off actually is.

Background
----------
`03_build_tls_score.py` selects the cut-off that maximises the log-rank
statistic (the "optimal cut-off", CTS = 0.2489) and reports a 200-resample
bootstrap of that cut-off. The bootstrap median (0.437) against the point
estimate (0.249) already hinted at instability, but a median alone does not
say whether the *prognostic claim* is stable. This script answers that
directly: inside each bootstrap resample it re-runs the **entire** selection
procedure and records

  1. the distribution of the selected cut-off (median, IQR, range);
  2. the proportion of resamples whose selected cut-off lies within +/-0.05
     of the discovery value;
  3. **the proportion of resamples in which the selected cut-off yields a
     significant log-rank test (p < 0.05)** -- how often a "significant"
     survival difference would have been found at all;
  4. the proportion that are significant *and* land in the same terzile as
     the discovery value.

Statistic 3 is the interpretable one: a survival association detectable in
only a minority of resamples should not be presented as a validated
threshold. Degenerate resamples are counted and reported, not silently
dropped.

The cut-off search reproduces `03_build_tls_score.py` exactly (same candidate
grid: 60 values evenly spaced between the 10th and 90th percentile,
`min_grp = 15`), so the point estimate agrees with the discovery pipeline.
The log-rank statistic uses the same vectorised risk-set formula; the whole
candidate grid is evaluated as one matrix operation, which makes 1000
resamples run in seconds rather than hours.

Input   data/processed/hnsc_with_hpv.csv   (CTS, os_time_days, os_event)
Output  results/03b_cutpoint_stability.json
        results/03b_cutpoint_stability.csv
"""

import json
import os

import numpy as np
import pandas as pd
from scipy.stats import chi2

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data", "processed")
RES = os.path.join(ROOT, "results")

N_BOOT = 1000
SEED = 20261016
N_CAND = 60          # same as 03_build_tls_score.py
MIN_GRP = 15         # same as 03_build_tls_score.py
P_LO, P_HI = 10.0, 90.0
WITHIN = 0.05


def cut_scan(score, time, event):
    """Vectorised log-rank scan over the full candidate grid.

    Returns (cut, chi2, p) for the most significant candidate, or
    (None, None, None) when no candidate satisfies the minimum group size.
    Identical candidate grid and risk-set formula to 03_build_tls_score.py.
    """
    lo, hi = np.percentile(score, P_LO), np.percentile(score, P_HI)
    cands = np.unique(np.round(np.linspace(lo, hi, N_CAND), 6))
    order = np.argsort(time)
    t_s, e_s, s_s = time[order], event[order], score[order]
    n = len(t_s)

    at_risk = np.arange(n, 0, -1).astype(float)        # risk set size
    d_ind = (e_s == 1).astype(float)                   # event indicator
    ev = e_s == 1

    G = s_s[None, :] > cands[:, None]                  # (C, n) group-1 mask
    n1 = G.sum(axis=1)
    ok = (n1 >= MIN_GRP) & ((n - n1) >= MIN_GRP)

    at1 = np.cumsum(G[:, ::-1], axis=1)[:, ::-1].astype(float)   # group-1 at risk

    with np.errstate(divide="ignore", invalid="ignore"):
        exp = d_ind[None, :] * at1 / at_risk[None, :]
    O1 = (G & ev[None, :]).sum(axis=1).astype(float)
    E1 = np.nansum(exp[:, ev], axis=1)

    keep = at_risk > 1
    frac = at1[:, keep] / at_risk[keep]
    v = np.nansum(
        frac * (1 - frac) * d_ind[keep][None, :]
        * (at_risk[keep] - d_ind[keep])[None, :]
        / (at_risk[keep] - 1)[None, :], axis=1)

    stat = np.full(len(cands), -np.inf)
    good = ok & (v > 1e-9) & np.isfinite(v)
    stat[good] = (O1[good] - E1[good]) ** 2 / v[good]

    k = int(np.argmax(stat))
    if not np.isfinite(stat[k]):
        return None, None, None
    return float(cands[k]), float(stat[k]), float(chi2.sf(stat[k], 1))


def logrank_at(score, time, event, c):
    """Log-rank chi2 / p for a FIXED threshold c (no scanning).

    This is the honest reproducibility measure for a reported threshold:
    the selection step is not repeated, so the test cannot benefit from
    maximisation.
    """
    g = score > c
    n1, n0 = int(g.sum()), int((~g).sum())
    if n1 == 0 or n0 == 0:
        return None, None
    order = np.argsort(time)
    t_s, e_s, g_s = time[order], event[order], g[order]
    n = len(t_s)
    at_risk = np.arange(n, 0, -1).astype(float)
    d_ind = (e_s == 1).astype(float)
    ev = e_s == 1
    at1 = np.cumsum(g_s[::-1])[::-1].astype(float)
    if at1[0] == 0 or at1[0] == n:
        return None, None
    with np.errstate(divide="ignore", invalid="ignore"):
        exp = d_ind * at1 / at_risk
    O1 = float((g_s & ev).sum())
    E1 = float(np.nansum(exp[ev]))
    keep = at_risk > 1
    frac = at1[keep] / at_risk[keep]
    v = float(np.nansum(frac * (1 - frac) * d_ind[keep]
                        * (at_risk[keep] - d_ind[keep])
                        / (at_risk[keep] - 1)))
    if v <= 1e-9 or not np.isfinite(v):
        return None, None
    stat = (O1 - E1) ** 2 / v
    return float(stat), float(chi2.sf(stat, 1))


def main():
    df = pd.read_csv(os.path.join(DATA, "hnsc_with_hpv.csv"))
    need = ["CTS", "os_time_days", "os_event"]
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise SystemExit("missing columns: %s" % missing)
    df = df.dropna(subset=need).reset_index(drop=True)

    score = df["CTS"].values.astype(float)
    time = df["os_time_days"].values.astype(float)
    event = df["os_event"].values.astype(int)
    n = len(df)

    pt_cut, pt_chi2, pt_p = cut_scan(score, time, event)
    pt_median = float(np.median(score))
    print("[data] n = %d, deaths = %d" % (n, int(event.sum())))
    print("[point] optimal cut-off = %.4f  chi2 = %.4f  p = %.4g"
          % (pt_cut, pt_chi2, pt_p))

    # verify the vectorised scan against the discovery pipeline
    DISCO_CUT, DISCO_CHI2 = 0.248926, 5.721466392604056
    if abs(pt_cut - DISCO_CUT) > 1e-6 or abs(pt_chi2 - DISCO_CHI2) > 1e-6:
        print("[warn] point estimate differs from the discovery pipeline: "
              "%.6f/%.6f vs %.6f/%.6f"
              % (pt_cut, pt_chi2, DISCO_CUT, DISCO_CHI2))
    else:
        print("[ok] reproduces the discovery pipeline exactly "
              "(cut %.6f, chi2 %.6f)" % (pt_cut, pt_chi2))
    p33, p67 = np.percentile(score, [100 / 3, 200 / 3])
    pt_terz = int(np.searchsorted([p33, p67], pt_cut))

    rng = np.random.default_rng(SEED)
    rows = []
    degenerate = 0
    fixed_p_opt = []      # p at the discovery threshold, applied unchanged
    fixed_p_med = []      # p at the median threshold, applied unchanged
    for b in range(N_BOOT):
        idx = rng.integers(0, n, n)
        s_b, t_b, e_b = score[idx], time[idx], event[idx]
        c, chi2v, p = cut_scan(s_b, t_b, e_b)
        if c is None:
            degenerate += 1
            continue
        tz = int(np.searchsorted([p33, p67], c))
        _, p_opt = logrank_at(s_b, t_b, e_b, pt_cut)
        _, p_med = logrank_at(s_b, t_b, e_b, pt_median)
        if p_opt is not None:
            fixed_p_opt.append(p_opt)
        if p_med is not None:
            fixed_p_med.append(p_med)
        rows.append({"bootstrap": b, "cut": round(c, 6),
                     "chi2": round(chi2v, 6), "p": p, "terzile": tz,
                     "same_terzile": int(tz == pt_terz),
                     "sig": int(p < 0.05),
                     "within_005": int(abs(c - pt_cut) <= 0.05),
                     "p_at_fixed_cut": p_opt if p_opt is not None else np.nan,
                     "p_at_median_cut": p_med if p_med is not None else np.nan})

    bdf = pd.DataFrame(rows)
    if bdf.empty:
        raise SystemExit("no valid bootstrap resamples")
    bdf.to_csv(os.path.join(RES, "03b_cutpoint_stability.csv"),
               index=False, encoding="utf-8-sig")

    fpo = np.array(fixed_p_opt)
    fpm = np.array(fixed_p_med)

    q25, q50, q75 = np.percentile(bdf["cut"], [25, 50, 75])
    summary = {
        "n": n,
        "n_deaths": int(event.sum()),
        "n_bootstrap_requested": N_BOOT,
        "n_bootstrap_valid": int(len(bdf)),
        "n_bootstrap_degenerate": int(degenerate),
        "seed": SEED,
        "n_candidates": N_CAND,
        "candidate_percentile_range": [P_LO, P_HI],
        "min_per_group": MIN_GRP,
        "within_tolerance": WITHIN,
        "point_cut": round(float(pt_cut), 6),
        "point_chi2": round(float(pt_chi2), 6),
        "point_p": float(pt_p),
        "point_terzile": pt_terz,
        "terzile_boundaries": [round(float(p33), 6), round(float(p67), 6)],
        "cut_median": round(float(q50), 6),
        "cut_iqr": [round(float(q25), 6), round(float(q75), 6)],
        "cut_range": [round(float(bdf["cut"].min()), 6),
                      round(float(bdf["cut"].max()), 6)],
        "frac_within_0.05_of_point": round(float(bdf["within_005"].mean()), 4),
        "frac_same_terzile": round(float(bdf["same_terzile"].mean()), 4),
        "frac_significant_p05": round(float(bdf["sig"].mean()), 4),
        "frac_sig_and_same_terzile": round(
            float((bdf["sig"] & (bdf["terzile"] == pt_terz)).mean()), 4),
        "p_median": round(float(bdf["p"].median()), 4),
        "p_q25": round(float(np.percentile(bdf["p"], 25)), 4),
        "p_q75": round(float(np.percentile(bdf["p"], 75)), 4),
        # --- honest reproducibility of the REPORTED threshold ---
        "median_cut": round(pt_median, 6),
        "fixed_cut_used": round(float(pt_cut), 6),
        "fixed_cut_sig_rate_p05": round(float((fpo < 0.05).mean()), 4),
        "fixed_cut_p_median": round(float(np.median(fpo)), 4),
        "fixed_cut_p_iqr": [round(float(np.percentile(fpo, 25)), 4),
                            round(float(np.percentile(fpo, 75)), 4)],
        "median_cut_sig_rate_p05": round(float((fpm < 0.05).mean()), 4),
        "median_cut_p_median": round(float(np.median(fpm)), 4),
    }
    with open(os.path.join(RES, "03b_cutpoint_stability.json"), "w",
              encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=2)

    print()
    print("--- cut-off stability over %d resamples (%d degenerate) ---"
          % (len(bdf), degenerate))
    print("  median selected cut-off : %.4f   (IQR %.4f-%.4f, range %.4f-%.4f)"
          % (q50, q25, q75, bdf["cut"].min(), bdf["cut"].max()))
    print("  within +/-%.2f of point  : %.1f%%" % (WITHIN, 100 * bdf["within_005"].mean()))
    print("  same terzile as point   : %.1f%%" % (100 * bdf["same_terzile"].mean()))
    print("  median log-rank p       : %.4f" % bdf["p"].median())
    print("  log-rank p < 0.05       : %.1f%%   <-- the interpretable number"
          % (100 * bdf["sig"].mean()))
    print("  significant AND same terzile : %.1f%%"
          % (100 * (bdf["sig"] & (bdf["terzile"] == pt_terz)).mean()))
    print()
    print("--- reproducibility of the REPORTED threshold (no scanning) ---")
    print("  discovery cut %.4f applied unchanged : significant in %.1f%% "
          "of resamples (median p = %.4f, IQR %.4f-%.4f)"
          % (pt_cut, 100 * (fpo < 0.05).mean(), np.median(fpo),
             np.percentile(fpo, 25), np.percentile(fpo, 75)))
    print("  median cut    %.4f applied unchanged : significant in %.1f%% "
          "(median p = %.4f)"
          % (pt_median, 100 * (fpm < 0.05).mean(), np.median(fpm)))


if __name__ == "__main__":
    main()
