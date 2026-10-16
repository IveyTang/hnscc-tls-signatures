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
| `05_external_validation.py` | Generic external-validation harness for a local cohort with a sample annotation file | — |
| `06_microenvironment_nomogram.py` | Marker-based microenvironmental programmes, nomogram, calibration, decision curve analysis | `06_*.csv/json`, `fig_nomogram.png`, `fig_calibration.png`, `fig_dca.png`, `fig_microenv_correlation.png` |
| `09_extract_hpv_panel.py` | Extract the 107-gene HPV host-expression panel | `hnsc_gene_matrix_hpv.csv` |
| `10_hpv_status.py` | Three-source HPV labelling (ISH > p16 > pathology report) and nested leave-one-out transcriptomic classifier | `10_hpv_summary.json`, `hnsc_hpv_status.csv`, `hnsc_with_hpv.csv` |
| `11_hpv_stratified.py` | HPV-stratified Kaplan–Meier and Cox analysis, confounding assessment, and interaction testing | `11_*.csv/json`, `fig_km_cts_by_hpv.png` |
| `20_fetch_gse65858_ebi.py` | Fetch the external cohort (GSE65858) from the **EMBL-EBI ArrayExpress mirror** (`E-GEOD-65858`) plus the platform annotation, bypassing NCBI | `data/external/GSE65858/**`, `fetch_log.json` |
| `20b_fetch_adf_ranged.py` | Download the 71.7 MB Illumina HT-12 v4 array design with parallel HTTP Range requests and verify its exact length | `A-GEOD-10558.adf.txt` |
| `21_external_validation_gse65858.py` | External validation: identical gene sets and scoring, no re-tuning; convergence, programme correlation, survival, HPV association, published-cluster check | `21_ext_*.csv/json`, `fig_ext_*.png` |
| `22_ext_figure_panel.py` | Compose the four external-validation panels into the manuscript figure and build the summary table | `fig_ext_validation.png`, `table6_external_validation.csv` |

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

## External cohort: why ArrayExpress instead of GEO

NCBI (`www.ncbi.nlm.nih.gov`, `ftp.ncbi.nlm.nih.gov`) is unreachable from the analysis
environment used here, and the agent egress IP is additionally blocked by NCBI with an
access-denied page. ArrayExpress imported GEO series up to about 2017, so **GSE65858 is
available in full as `E-GEOD-65858`** on the EMBL-EBI mirror, together with the matching
array design `A-GEOD-10558`:

```
record     https://www.ebi.ac.uk/biostudies/api/v1/studies/E-GEOD-65858
files      https://ftp.ebi.ac.uk/biostudies/fire/E-GEOD-/858/E-GEOD-65858/Files/
legacy     https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEOD/E-GEOD-65858/
platform   https://ftp.ebi.ac.uk/biostudies/fire/A-GEOD-/558/A-GEOD-10558/Files/
```

Run the retrieval and the validation in this order:

```bash
python scripts/20_fetch_gse65858_ebi.py       # 270 per-sample tables + SDRF + IDF (~4 min)
python scripts/20b_fetch_adf_ranged.py        # 71.7 MB array design, range-parallel (~1 min)
python scripts/21_external_validation_gse65858.py
python scripts/22_ext_figure_panel.py
```

Two practical notes that cost time to discover:

- A truncated transfer still returns HTTP 200. The first attempt at the array design
  stopped at 10.8% (7.8 MB of 71.7 MB) and parsed cleanly while silently covering only
  5,193 of 47,323 probes. Verify exact byte lengths, not minimum sizes; `20b` does this
  and also demonstrates the parallel-range pattern, which is ~6× faster than a single
  stream through the local proxy.
- `GSE159067`, the immunotherapy-treated cohort used in the original analysis plan, was
  submitted after ArrayExpress stopped importing GEO and has **no** ArrayExpress mirror.
  It is therefore not part of this repository's external validation, and no
  immunotherapy-response claim rests on it.
