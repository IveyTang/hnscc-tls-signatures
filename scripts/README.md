# Pipeline scripts

Scripts are numbered in execution order. Each reads from `data/` and writes to `results/`; paths are resolved relative to the repository root, so scripts can be run from anywhere.

## Acquisition (network required)

| Script | Purpose | Runtime |
|---|---|---|
| `00_fetch_gdc_manifest.py` | Retrieve case-level clinical records and the file manifest for TCGA-HNSC from the GDC API | ~2 min |
| `01_build_clinical_table.py` | Build the analysis cohort: three-way primary-tumour filter, stage fallback, OS computation | ~1 min |
| `02_download_expression.py` | Download STAR-Counts files and extract `tpm_unstranded` into a gene × sample matrix | ~8 min, 2.3 GB |
| `07_fetch_bcr_hpv.py` | Parse HPV fields (`hpv_status_by_ish_testing`, `hpv_status_by_p16_testing`) from BCR clinical XML | ~2 min |
| `08_fetch_path_reports.py` | Retrieve diagnostic pathology report PDFs | ~10 min |
| `08b_reclassify_reports.py` | Clause-level p16/HPV reclassification with optical-character-recognition correction | ~2 min |

## Analysis (offline once the above have run)

| Script | Purpose | Key outputs |
|---|---|---|
| `03_build_tls_score.py` | Score seven gene sets under three algorithms; build the consensus TLS score; signature convergence matrix | `03_score_summary.json`, `signature_convergence.csv`, `hnsc_tls_scores.csv` |
| `04_internal_validation.py` | Kaplan–Meier, log-rank, univariable and multivariable Cox, time-dependent AUC, C-index | `04_*.csv/json`, `fig_km_cts.png`, `fig_forest_multivariate.png`, `fig_roc_timedependent.png`, `fig_signature_convergence.png` |
| `05_external_validation.py` | External validation on an external cohort using the same platform-agnostic rank scoring, without re-tuning | *blocked — see below* |
| `06_microenvironment_nomogram.py` | Marker-based microenvironmental programmes, nomogram, calibration, decision curve analysis | `06_*.csv/json`, `fig_nomogram.png`, `fig_calibration.png`, `fig_dca.png`, `fig_microenv_correlation.png` |
| `09_extract_hpv_panel.py` | Extract the 107-gene HPV host-expression panel | `hnsc_gene_matrix_hpv.csv` |
| `10_hpv_status.py` | Three-source HPV labelling (ISH > p16 > pathology report) and nested leave-one-out transcriptomic classifier | `10_hpv_summary.json`, `hnsc_hpv_status.csv`, `hnsc_with_hpv.csv` |
| `11_hpv_stratified.py` | HPV-stratified Kaplan–Meier and Cox analysis, confounding assessment, and interaction testing | `11_*.csv/json`, `fig_km_cts_by_hpv.png` |

## Utility

| Script | Purpose |
|---|---|
| `00b_restore_derived.py` | Copy the committed derived tables into `data/processed/` so the offline steps can run without network access |

## Exploratory

`exploratory/` contains ad-hoc probes used while establishing the data schema:

- `_probe_clinical.py` — inspect available GDC clinical fields and their completeness
- `_probe_hpv.py` — locate HPV-related fields in the BCR XML
- `_diag_reports.py` — diagnose pathology-report text extraction

These are retained because they document how the GDC schema pitfalls described in `../docs/data_dictionary.md` were discovered. They are **not** part of the analysis pipeline and do not need to be run.

## Blocked step

`05_external_validation.py` is complete and ready to run, but the external cohorts could not be retrieved in the analysis environment because NCBI/GEO was unreachable (HTTP 000). To use it:

1. Download **GSE65858** (HPV-annotated HNSCC surgical cohort with survival) and **GSE159067** (ICI-treated cohort).
2. Place them in `data/raw/geo/<dataset>/` with a sample annotation file.
3. Run:
   ```bash
   python scripts/05_external_validation.py --dataset GSE65858 --time-col OS_time --event-col OS_event
   ```

The script applies the same platform-agnostic rank-based scoring and the TCGA-derived cut-off, without re-tuning. The manuscript therefore describes the findings as discovery-only.
