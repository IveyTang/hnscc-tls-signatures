"""
26_supp_figure_hpv_roc.py
-------------------------
Supplementary Figure S1: three-panel ROC figure for the cross-platform
external validation of the transcriptomic HPV classifier.

Panel A  DNA+RNA+ vs DNA-          (primary target, n = 231)
Panel B  any HPV16 DNA+ vs DNA-    (secondary target, n = 250)
Panel C  oropharynx only, DNA+RNA+ vs DNA-   (n = 75)

The three curves in each panel are the transplanted 8-gene classifier
(unchanged coefficients), the unweighted mean of the same eight genes, and
CDKN2A alone as a single-gene benchmark.

Input   results/25_hpv_classifier_external.json / .csv   (numbers, for the legend)
        results/10_hpv_summary.json                      (genes + coefficients)
        data/external/GSE65858/**                        (array + SDRF)
Output  results/fig_supp_hpv_roc.png
"""
import importlib.util
import json
import os

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "results")


def _load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    j25 = json.load(open(os.path.join(RES, "25_hpv_classifier_external.json"),
                         encoding="utf-8"))
    m10 = json.load(open(os.path.join(RES, "10_hpv_summary.json"),
                         encoding="utf-8"))
    genes = list(m10["top_genes"])
    coefs = dict(m10["coefficients"])

    ext21 = _load_module(os.path.join(HERE, "21_external_validation_gse65858.py"),
                         "ext21")
    mod25 = _load_module(os.path.join(HERE, "25_hpv_classifier_external.py"),
                         "mod25")

    probe_map, _, _ = ext21.load_probe_map()
    raw, _ = ext21.load_expression(probe_map)
    sym = raw.index.map(lambda p: probe_map.get(p))
    keep = ~pd.isna(sym)
    sub = raw.loc[keep].copy()
    sub.index = sym[keep]
    gene_mat = sub.groupby(level=0).mean()

    hpv = mod25.load_hpv_labels()
    common = [s for s in gene_mat.columns if s in hpv.index]
    gene_mat = gene_mat[common]
    hpv = hpv.loc[common]

    z = gene_mat.loc[genes]
    z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1).replace(0, np.nan), axis=0)
    w = np.array([coefs[g] for g in genes], dtype=float)
    L = pd.Series((z.T.values * w).sum(axis=1), index=common)
    Lmean = pd.Series(z.T.values.mean(axis=1), index=common)
    cdkn2a = pd.Series(gene_mat.loc["CDKN2A"], index=common)

    clin = ext21.load_clinical()
    site = clin.reindex(common)["tumor_site"].astype(str).str.strip().str.lower()
    oro = site.eq("oropharynx")

    panels = [
        ("A", "DNA+RNA+ versus DNA-",
         hpv.isin({"DNA+RNA+", "DNA-"}), hpv.isin({"DNA+RNA+"}), oro | True),
        ("B", "any HPV16 DNA+ versus DNA-",
         hpv.isin({"DNA+RNA+", "DNA+RNA-", "DNA-"}), hpv.isin({"DNA+RNA+", "DNA+RNA-"}), oro | True),
        ("C", "Oropharynx only, DNA+RNA+ versus DNA-",
         hpv.isin({"DNA+RNA+", "DNA-"}), hpv.isin({"DNA+RNA+"}), oro),
    ]

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"8-gene classifier": "#A32D2D",
              "unweighted panel mean": "#4A6FA5",
              "CDKN2A alone": "#0F6E56"}
    series = {"8-gene classifier": L,
              "unweighted panel mean": Lmean,
              "CDKN2A alone": cdkn2a}

    fig, axes = plt.subplots(1, 3, figsize=(10.6, 3.7), dpi=300)
    for ax, (tag, title, sel, pos, mask) in zip(axes, panels):
        sel = sel & mask
        y = pos[sel].astype(int).values
        n = len(y)
        for label, s in series.items():
            v = s[sel].values
            fpr, tpr, _ = roc_curve(y, v)
            a = roc_auc_score(y, v)
            ax.plot(fpr, tpr, lw=1.8, color=colors[label],
                    label="%s (%.3f)" % (label, a))
        ax.plot([0, 1], [0, 1], ls=":", lw=1, color="#888888")
        ax.set_xlim(-0.01, 1.01)
        ax.set_ylim(-0.01, 1.01)
        ax.set_xlabel("1 - specificity", fontsize=9)
        ax.set_ylabel("Sensitivity", fontsize=9)
        ax.set_title("%s. %s\n(n = %d; %d positive)" % (tag, title, n, int(y.sum())),
                     fontsize=8.6)
        ax.tick_params(labelsize=8)
        ax.legend(frameon=False, fontsize=7.2, loc="lower right")
    fig.tight_layout()
    out = os.path.join(RES, "fig_supp_hpv_roc.png")
    fig.savefig(out)
    plt.close(fig)
    print("[fig] %s" % out)

    for ax_tag, (tag, title, sel, pos, mask) in zip("ABC", panels):
        sel = sel & mask
        y = pos[sel].astype(int).values
        line = ["  %s  %s  n=%d (%d pos)" % (tag, title, len(y), int(y.sum()))]
        for label, s in series.items():
            line.append("%s=%.3f" % (label, roc_auc_score(y, s[sel].values)))
        print("  ".join(line))


if __name__ == "__main__":
    main()
