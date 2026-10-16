# Data directory

Everything in `derived/` was computed by the scripts in `../scripts/` from publicly available TCGA-HNSC data obtained through the Genomic Data Commons. No controlled-access or identifiable data are included.

## Files

| File | Rows | Size | Description |
|---|---|---|---|
| `derived/hnsc_analysis_set.csv` | 438 | 119 KB | Final analysis cohort after the three-way primary-tumour filter (one row per patient) |
| `derived/hnsc_gene_matrix.csv` | 1291 | 5.1 MB | Expression matrix for all analysed genes (log2 TPM + 1), samples in columns |
| `derived/hnsc_tls_scores.csv` | 438 | 246 KB | Analysis cohort plus the seven gene set scores, their first principal components, the consensus TLS score, and clinical covariates |
| `derived/hnsc_hpv_status.csv` | 438 | 52 KB | Analysis cohort plus HPV label, label source, classifier probability, and predicted status |
| `derived/hnsc_with_hpv.csv` | 438 | 298 KB | Complete merged table used by the survival analyses: gene set scores + CTS + full clinical annotation + HPV status |
| `derived/hnsc_gene_matrix_hpv.csv` | 107 | 467 KB | Wide matrix of the 107-gene HPV host-expression panel (log2 TPM + 1); columns are TCGA aliquots |
| `expr_download_meta.json` | — | 18 KB | Provenance record for the expression download: file UUIDs, sample barcodes, and sample types |

## Recommended entry point

For reanalysing the survival results, start from **`derived/hnsc_with_hpv.csv`** — it contains every variable used in the manuscript's Tables 1–5 and Figures 1–9.

```python
import pandas as pd
df = pd.read_csv("data/derived/hnsc_with_hpv.csv")
print(df.shape)          # (438, 69)
print(df["CTS"].describe())
```

## Notes on the columns

- **`CTS`** — the consensus TLS score: the mean percentile rank of the five TLS-related gene sets scored by the per-gene z-score algorithm. Scale 0–1. Higher values indicate a more immune-inflamed phenotype.
- **`CTS_group_opt`** — dichotomisation at the optimal cut-off (0.249), determined by maximising the log-rank statistic across candidate thresholds between the 10th and 90th percentiles.
- **`CTS_group_median`** — dichotomisation at the median (sensitivity analysis).
- **`os_time_days` / `os_event`** — overall survival in days from diagnosis and the event indicator (1 = death, 0 = censored).
- **`hpv_final`** — the HPV status used in the manuscript: direct histological/immunohistochemical annotation where available, otherwise the transcriptomic classifier.
- **`hpv_label`** — direct annotation only (`Positive` / `Negative` for 84 patients, empty otherwise); `hpv_label_source` records whether it came from ISH, p16, or pathology-report text.
- **`hpv_pred` / `hpv_prob`** — classifier output for all 438 patients.
- **`pack_years_x` / `pack_years_y`** — pack-years from two different GDC exposure encodings; `_x` is the value used in the analysis.
- **`multi_primary` / `diagnosis_rule`** — provenance for the primary-tumour filter described in the manuscript Methods.

See `../docs/data_dictionary.md` for a column-by-column reference including all clinical fields.

## Regenerating the raw data

The raw STAR-Counts expression files (about 2.3 GB, 563 files) are **not** committed. To rebuild them:

```bash
python scripts/00_fetch_gdc_manifest.py
python scripts/01_build_clinical_table.py
python scripts/02_download_expression.py --workers 10
```

This takes roughly 10 minutes on a typical connection. Because the files are cached, re-running the pipeline is idempotent.

## Licence

Derived data in this directory are released under **CC BY 4.0** (see `../LICENSE-DATA`). Upstream TCGA data are subject to the NIH Genomic Data Sharing Policy; users should also cite the TCGA-HNSC marker paper (Lawrence et al., *Nature* 2015;517:576–582).
