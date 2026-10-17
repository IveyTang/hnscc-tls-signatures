# Reproducibility

## Environment

| Item | Value |
|---|---|
| Operating system | Windows 11 (also tested on Linux) |
| Python | 3.13.12 (requires ≥ 3.10) |
| Package versions | see `requirements-pinned.txt` |
| Random seed | `20261016`, fixed for all bootstrap resampling |
| Network requirement | Steps 00–02, 07, and 08 only |

Install:

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements-pinned.txt     # exact versions
# or: pip install -r requirements.txt       # minimum versions
```

## Reproduction paths

### Path A — full reproduction from the Genomic Data Commons

Requires network access and about 2.5 GB of free disk space (2.3 GB of which is raw STAR-Counts files that can be deleted afterwards).

```bash
python scripts/00_fetch_gdc_manifest.py       # ~2 min   clinical records + file manifest
python scripts/01_build_clinical_table.py     # ~1 min   cohort construction
python scripts/02_download_expression.py --workers 10   # ~8 min   expression download
python scripts/03_build_tls_score.py          # ~1 min   scoring + convergence
python scripts/04_internal_validation.py      # ~2 min   survival + discrimination
python scripts/06_microenvironment_nomogram.py # ~2 min  microenvironment + nomogram
python scripts/07_fetch_bcr_hpv.py            # ~2 min   BCR XML HPV fields
python scripts/08_fetch_path_reports.py       # ~10 min  pathology report PDFs
python scripts/08b_reclassify_reports.py      # ~2 min   clause-level p16 reclassification
python scripts/09_extract_hpv_panel.py        # ~1 min   107-gene panel
python scripts/10_hpv_status.py               # ~2 min   HPV labels + classifier
python scripts/11_hpv_stratified.py           # ~1 min   HPV-stratified analysis
```

Total runtime roughly **35 minutes**, dominated by the expression download and the pathology report retrieval.

### Path B — reproduce the statistics only (no network)

Every table and figure in the manuscript can be regenerated from the committed derived tables with no network access. A helper script restores them into the pipeline's working directory:

```bash
python scripts/00b_restore_derived.py       # copies data/derived/*.csv -> data/processed/

python scripts/03_build_tls_score.py            # scoring + convergence
python scripts/04_internal_validation.py        # survival + discrimination
python scripts/06_microenvironment_nomogram.py  # microenvironment + nomogram
python scripts/11_hpv_stratified.py             # HPV-stratified analysis
```

Total runtime under **10 minutes**.

Step **10** (`10_hpv_status.py`) cannot run offline: it reads the raw BCR clinical XML and pathology-report CSV produced by steps 07–08, which are not committed because they are retrieved from the GDC at runtime. Its output, `hnsc_hpv_status.csv`, is committed and restored by `00b_restore_derived.py`, so step 11 — the HPV-stratified analysis reported in the manuscript — runs without it.

Steps **03** and **06** need `hnsc_gene_matrix.csv` (5.1 MB, 1291 genes × 563 samples), which **is** committed under `data/derived/`. If you would rather not download the repository at that size, delete that one file and skip steps 03 and 06; steps 04 and 11 will still run against the committed `hnsc_tls_scores.csv` and `hnsc_with_hpv.csv`.

### Path C — reproduce the cut-off stability analysis and the HPV cross-platform validation (no network)

Both of these reuse derived tables that are committed, so neither needs network access:

```bash
python scripts/03b_cutpoint_stability.py        # ~3 min   1,000 bootstrap resamples of the cut-off search
python scripts/25_hpv_classifier_external.py    # ~1 min   transplant the 8-gene HPV classifier into GSE65858
python scripts/26_supp_figure_hpv_roc.py        # ~1 min   Supplementary Figure S1 (three-panel ROC)
```

`03b` reproduces the discovery-cohort point estimate exactly (cut-off 0.248926, χ² 5.721466) before resampling. `25` and `26` need the GSE65858 array and SDRF under `data/external/GSE65858/`; see the retrieval commands in `scripts/README.md`.

## Determinism

| Source of variation | Handling |
|---|---|
| Bootstrap resampling | `numpy.random.default_rng(20261016)`, fixed (1,000 resamples for the cut-off; 2,000 for AUC confidence intervals) |
| Classifier cross-validation | Deterministic leave-one-out; feature selection repeated within each fold |
| Logistic regression solver | `lbfgs`, default tolerance; convergence checked |
| Optimal cut-off search | Exhaustive deterministic scan over percentile grid |
| Figure rendering | Matplotlib version affects anti-aliasing and font metrics only; all numerical results are stable |

Re-running the full pipeline on the same input reproduces every value in `results/` exactly, except for PNG byte-level differences arising from font rendering.

## Verification

Cross-checks performed between implementations:

- The vectorised log-rank statistic used for cut-off selection matches `lifelines.statistics.logrank_test` exactly at the selected threshold (χ² = 5.72 by both methods).
- Time-dependent AUC uses a cumulative/dynamic definition; scores are negated to define a risk score, with 0.5 as the null value.
- The multivariable Cox model is fitted with a ridge penalty of 0.1 and zero-variance covariates are removed before fitting. Unpenalised refits produced equivalent conclusions.
- The HPV classifier was trained with feature selection inside each cross-validation fold specifically to test for and exclude information leakage.

## Known limitations on reproducibility

1. **External validation requires an initial retrieval step.** Both external cohorts are committed as derived tables, so every reported statistic reproduces offline; re-running the retrieval itself (`scripts/20_fetch_gse65858_ebi.py`, `scripts/23_fetch_gse159067_geo.py`) needs network access. Note that the earlier release of this repository stated that NCBI/GEO was unreachable from the analysis environment; that claim was wrong (all requests were being routed through an HTTP proxy that returned HTTP 000 for NCBI hosts) and has been retracted. GEO is reachable over plain FTP, or over HTTPS to NCBI's IP with the proxy bypassed — the pattern used by `23_fetch_gse159067_geo.py`, and the same one used by the team's `IveyTang/ici-benchmark` repository.
2. **Pathology report parsing depends on OCR of PDF text.** `scripts/08b_reclassify_reports.py` applies clause-level splitting and correction of common optical-character-recognition errors. A different PDF extraction library may yield slightly different text and therefore a different annotation for the small number of patients (3) whose label derives from this source.
3. **Annotation coverage.** Direct HPV labels exist for only 84 of 438 patients. The classifier is validated against those same labels, so its performance estimate is not independent of the annotation set; performance outside the oropharynx is less certain.
