"""
25_hpv_classifier_external.py
-----------------------------
Cross-platform external validation of the transcriptomic HPV classifier.

Why this matters
----------------
The discovery paper annotates 84 of 438 TCGA-HNSC patients for HPV directly
(ISH / p16 / pathology text) and derives a transcriptomic classifier for the
remainder, validated by nested leave-one-out cross-validation within those
same 84 patients (LOO-AUC 0.990). That is an *internal* validation on the same
platform (RNA-seq, log2 TPM+1). Reviewers will reasonably ask whether the
classifier is a TCGA artefact.

GSE65858 provides an independent answer: 270 HNSCC patients profiled on the
Illumina HumanHT-12 v4 array with HPV16 status resolved at *both* the DNA and
the viral E6*I RNA level (DNA+RNA+ 35, DNA+RNA- 19, DNA- 196, NA 20).

What is transplanted, and what is not
-------------------------------------
The eight genes selected by the nested LOO procedure and their fitted
logistic-regression coefficients are used **unchanged**. The training
`StandardScaler` mean and the model intercept were not recorded in the
discovery summary, so the external application standardises each gene within
GSE65858 and evaluates the linear predictor

    L = sum_g w_g * z_g

Discrimination (ROC-AUC) is invariant to additive constants, so L supports an
exact AUC comparison without the intercept. Sensitivity and specificity
*are* threshold dependent; those are reported at Youden thresholds re-derived
within GSE65858 and labelled as such, never as a transplanted operating point.

The comparison is deliberately conservative: RNA-seq log2 TPM versus array
log2 intensity are different measurement scales, and the classifier was
selected on TCGA (predominantly non-oropharyngeal, largely HPV-negative),
whereas GSE65858 is enriched for oropharyngeal disease.

Input   results/10_hpv_summary.json          (genes + coefficients)
        data/external/GSE65858/**            (array + SDRF metadata)
Output  results/25_hpv_classifier_external.json
        results/25_hpv_classifier_external.csv
        results/fig_hpv_classifier_external.png
"""

import glob
import importlib.util
import json
import os
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "results")
EXT = os.path.join(ROOT, "data", "external", "GSE65858")

N_BOOT = 2000
SEED = 20261016


def _load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_hpv_labels():
    """HPV16 DNA/RNA status per sample, from the ArrayExpress SDRF."""
    import csv
    sdrf = os.path.join(EXT, "E-GEOD-65858.sdrf.txt")
    with open(sdrf, encoding="utf-8") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    hdr = rows[0]
    i_assay = hdr.index("Assay Name")
    i_rna = hdr.index("Characteristics [hpv16_dna_rna]")
    rec = {}
    for r in rows[1:]:
        if len(r) <= max(i_assay, i_rna):
            continue
        gsm = r[i_assay].split("_")[0].strip()
        rec[gsm] = r[i_rna].strip()
    return pd.Series(rec, name="hpv16_dna_rna")


def auc_ci(y, s, n_boot=N_BOOT, seed=SEED):
    """Bootstrap 95% CI for the AUC."""
    rng = np.random.default_rng(seed)
    y = np.asarray(y)
    s = np.asarray(s, dtype=float)
    n = len(y)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        if len(np.unique(y[idx])) < 2:
            continue
        vals.append(roc_auc_score(y[idx], s[idx]))
    vals = np.array(vals)
    return (round(float(np.percentile(vals, 2.5)), 4),
            round(float(np.percentile(vals, 97.5)), 4),
            len(vals))


def sens_spec(y, s, thr):
    pred = (s >= thr).astype(int)
    y = np.asarray(y)
    tp = int(((pred == 1) & (y == 1)).sum())
    fn = int(((pred == 0) & (y == 1)).sum())
    tn = int(((pred == 0) & (y == 0)).sum())
    fp = int(((pred == 1) & (y == 0)).sum())
    return (tp / (tp + fn) if tp + fn else np.nan,
            tn / (tn + fp) if tn + fp else np.nan)


def main():
    sys.path.insert(0, HERE)
    print("=" * 74)
    print("Cross-platform validation of the transcriptomic HPV classifier")
    print("=" * 74)

    # ---- model to transplant -----------------------------------------
    m10 = json.load(open(os.path.join(RES, "10_hpv_summary.json"),
                         encoding="utf-8"))
    genes = list(m10["top_genes"])
    coefs = dict(m10["coefficients"])
    print("[model] %d genes: %s" % (len(genes), ", ".join(genes)))
    print("[model] discovery nested-LOO AUC = %.3f (within TCGA, same platform)"
          % m10["loo_auc"])
    print("[model] coefficients: %s"
          % ", ".join("%s=%+.3f" % (g, coefs[g]) for g in genes))

    # ---- external expression ----------------------------------------
    ext21 = _load_module(os.path.join(HERE, "21_external_validation_gse65858.py"),
                         "ext21")
    probe_map, n_probe, _ = ext21.load_probe_map()
    raw, gsms = ext21.load_expression(probe_map)
    sym = raw.index.map(lambda p: probe_map.get(p))
    keep = ~pd.isna(sym)
    sub = raw.loc[keep].copy()
    sub.index = sym[keep]
    gene_mat = sub.groupby(level=0).mean()
    print("[expr] GSE65858: %d genes x %d samples" % gene_mat.shape)

    present = [g for g in genes if g in gene_mat.index]
    absent = [g for g in genes if g not in gene_mat.index]
    print("[expr] classifier gene coverage: %d/%d%s"
          % (len(present), len(genes),
             "" if not absent else "  (absent: %s)" % ", ".join(absent)))
    if absent:
        print("[stop] cannot transplant a model with missing features")
        return

    # ---- labels ------------------------------------------------------
    hpv = load_hpv_labels()
    common = [s for s in gene_mat.columns if s in hpv.index]
    gene_mat = gene_mat[common]
    hpv = hpv.loc[common]
    print("[labels] %d samples with expression + HPV status" % len(common))
    print("[labels] %s" % hpv.value_counts().to_dict())

    # ---- standardise within the external cohort, apply weights --------
    z = gene_mat.loc[genes]
    mu = z.mean(axis=1)
    sd = z.std(axis=1).replace(0, np.nan)
    z = z.sub(mu, axis=0).div(sd, axis=0)
    w = np.array([coefs[g] for g in genes], dtype=float)
    L = (z.T.values * w).sum(axis=1)          # intercept-free linear predictor
    Lmean = z.T.values.mean(axis=1)           # unweighted panel mean
    cdkn2a = (gene_mat.loc["CDKN2A"].T.values if "CDKN2A" in gene_mat.index
              else None)

    # ---- comparisons -------------------------------------------------
    comparisons = {
        "primary_dnaRnaPos_vs_dnaNeg": ({"DNA+RNA+"}, {"DNA-"}),
        "secondary_anyDnaPos_vs_dnaNeg": ({"DNA+RNA+", "DNA+RNA-"}, {"DNA-"}),
    }

    site = None
    try:
        clin = ext21.load_clinical()
        site = clin.reindex(common)["tumor_site"]
    except Exception as exc:                                  # pragma: no cover
        print("[warn] tumor site not available: %s" % exc)

    rows = []

    def evaluate(name, pos, neg, mask_extra=None):
        sel = hpv.isin(pos | neg)
        if mask_extra is not None:
            sel = sel & mask_extra
        y = hpv[sel].isin(pos).astype(int).values
        if y.sum() == 0 or (y == 0).sum() == 0:
            print("  [skip] %s: degenerate" % name)
            return
        outs = {"linear_predictor_8genes": L[sel.values],
                "unweighted_mean_8genes": Lmean[sel.values]}
        if cdkn2a is not None:
            outs["CDKN2A_alone"] = cdkn2a[sel.values]
        print("  --- %s  (n=%d: %d pos / %d neg) ---"
              % (name, len(y), int(y.sum()), int((y == 0).sum())))
        for label, s in outs.items():
            a = float(roc_auc_score(y, s))
            lo, hi, nb = auc_ci(y, s)
            fpr, tpr, thr = roc_curve(y, s)
            j = np.argmax(tpr - fpr)
            se, sp = sens_spec(y, s, thr[j])
            print("      %-26s AUC %.3f [%.3f-%.3f]  Youden-cut sens %.3f spec %.3f"
                  % (label, a, lo, hi, se, sp))
            rows.append({"comparison": name, "n": int(len(y)),
                         "n_pos": int(y.sum()), "n_neg": int((y == 0).sum()),
                         "score": label, "auc": round(a, 4),
                         "auc_lo": lo, "auc_hi": hi,
                         "sens_at_youden": round(float(se), 4),
                         "spec_at_youden": round(float(sp), 4),
                         "n_boot_valid": nb})

    print()
    print("--- discrimination of the transplanted model ---")
    for name, (pos, neg) in comparisons.items():
        evaluate(name, pos, neg)
    # oropharyngeal subset (the anatomical site where HPV status matters most)
    if site is not None:
        oro = site.astype(str).str.strip().str.lower().eq("oropharynx")
        print("  [site] oropharynx n = %d" % int(oro.sum()))
        pos, neg = comparisons["primary_dnaRnaPos_vs_dnaNeg"]
        evaluate("primary_oropharynx_only", pos, neg, mask_extra=oro)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(RES, "25_hpv_classifier_external.csv"),
              index=False, encoding="utf-8-sig")

    # ---- summary -----------------------------------------------------
    def _auc(comp, score="linear_predictor_8genes"):
        m = df[(df["comparison"] == comp) & (df["score"] == score)]
        return None if m.empty else m.iloc[0].to_dict()

    prim = _auc("primary_dnaRnaPos_vs_dnaNeg")
    sec = _auc("secondary_anyDnaPos_vs_dnaNeg")
    oro = _auc("primary_oropharynx_only")
    base = _auc("primary_dnaRnaPos_vs_dnaNeg", "CDKN2A_alone")
    unw = _auc("primary_dnaRnaPos_vs_dnaNeg", "unweighted_mean_8genes")

    summary = {
        "cohort": "GSE65858 (Illumina HumanHT-12 v4)",
        "n_with_expression_and_hpv": int(len(hpv)),
        "label_distribution": {str(k): int(v)
                               for k, v in hpv.value_counts().items()},
        "model_genes": genes,
        "coefficients": coefs,
        "discovery_loo_auc": round(float(m10["loo_auc"]), 4),
        "standardisation": "per-gene z-score within GSE65858 samples",
        "intercept_note": ("the discovery intercept was not recorded; ROC-AUC is "
                           "invariant to additive constants, so discrimination is "
                           "comparable, whereas sensitivity/specificity at a fixed "
                           "operating point are not and are re-derived here"),
        "primary_dnaRnaPos_vs_dnaNeg": prim,
        "secondary_anyDnaPos_vs_dnaNeg": sec,
        "primary_oropharynx_only": oro,
        "baseline_CDKN2A_alone": base,
        "baseline_unweighted_mean": unw,
    }
    with open(os.path.join(RES, "25_hpv_classifier_external.json"), "w",
              encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=2)

    # ---- figure ------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        pos, neg = comparisons["primary_dnaRnaPos_vs_dnaNeg"]
        sel = hpv.isin(pos | neg)
        y = hpv[sel].isin(pos).astype(int).values
        fig, ax = plt.subplots(figsize=(5.2, 5.0), dpi=300)
        for label, s, col in (("8-gene classifier", L[sel.values], "#A32D2D"),
                              ("unweighted panel mean", Lmean[sel.values], "#4A6FA5"),
                              ("CDKN2A alone", (cdkn2a[sel.values]
                                                if cdkn2a is not None else None),
                               "#0F6E56")):
            if s is None:
                continue
            fpr, tpr, _ = roc_curve(y, s)
            a = roc_auc_score(y, s)
            ax.plot(fpr, tpr, lw=1.8, color=col,
                    label="%s (AUC %.3f)" % (label, a))
        ax.plot([0, 1], [0, 1], ls=":", lw=1, color="#888888")
        ax.set_xlabel("1 - specificity")
        ax.set_ylabel("Sensitivity")
        ax.set_title("HPV classifier transplanted to GSE65858\n"
                     "DNA+RNA+ vs DNA- (n = %d)" % len(y), fontsize=9.5)
        ax.legend(frameon=False, fontsize=8, loc="lower right")
        ax.set_xlim(-0.01, 1.01)
        ax.set_ylim(-0.01, 1.01)
        fig.tight_layout()
        fig.savefig(os.path.join(RES, "fig_hpv_classifier_external.png"))
        plt.close(fig)
        print("\n[fig] fig_hpv_classifier_external.png")
    except Exception as exc:                                  # pragma: no cover
        print("\n[warn] figure failed: %s" % exc)

    print()
    print("--- headline numbers ---")
    if prim:
        print("  discovery (same-platform LOO)      AUC %.3f" % m10["loo_auc"])
        print("  external, DNA+RNA+ vs DNA-         AUC %.3f [%.3f-%.3f]"
              % (prim["auc"], prim["auc_lo"], prim["auc_hi"]))
    if sec:
        print("  external, any HPV16 DNA+ vs DNA-   AUC %.3f [%.3f-%.3f]"
              % (sec["auc"], sec["auc_lo"], sec["auc_hi"]))
    if oro:
        print("  external, oropharynx only          AUC %.3f [%.3f-%.3f]"
              % (oro["auc"], oro["auc_lo"], oro["auc_hi"]))
    if base:
        print("  external, CDKN2A alone             AUC %.3f" % base["auc"])
    if unw:
        print("  external, unweighted panel mean    AUC %.3f" % unw["auc"])


if __name__ == "__main__":
    main()
