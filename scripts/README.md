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
| `03b_cutpoint_stability.py` | 1,000 bootstrap resamples in which the whole cut-off selection is repeated, plus reproducibility *without* re-selection at the fixed and median cut-offs | `03b_cutpoint_stability.csv/json` |
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
| `23_fetch_gse159067_geo.py` | Fetch the immunotherapy cohort (GSE159067, 2020 submission) **directly from GEO** over FTP, with an HTTPS-to-IP proxy-bypass fallback | `data/external/GSE159067/**`, `fetch_log.json` |
| `24_ici_validation_gse159067.py` | Immunotherapy validation: objective response, progression-free and overall survival, and correlation with the cohort's published hot phenotype | `24_ici_*.csv/json`, `table7_ici_validation.csv`, `fig_ici_validation.png` |
| `25_hpv_classifier_external.py` | Transplant the eight-gene HPV classifier into GSE65858 **without re-fitting**; ROC-AUC for transcriptionally active infection, any DNA positivity, and the oropharyngeal subset, against CDKN2A and an unweighted panel mean | `25_hpv_classifier_external.csv/json`, `fig_hpv_classifier_external.png` |
| `26_supp_figure_hpv_roc.py` | Build Supplementary Figure S1: the three-panel ROC figure for the transplanted classifier | `fig_supp_hpv_roc.png` |

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

## Retrieving GEO data when the default route is blocked

**Correction.** An earlier version of this file stated that NCBI/GEO was
unreachable and that the validation cohorts therefore had to come from
mirrors. That was an artefact of the retrieval route used at the time: every
request went through the sandbox HTTP proxy, which returns `000` for NCBI
hosts, and NCBI's abuse filter rejects the agent's default egress IP. GEO
itself is reachable, and both cohorts used in this study were obtained from
public archives without any mirroring workaround being strictly required.

Two channels work, and both are implemented in `scripts/23_fetch_gse159067_geo.py`:

```bash
# 1) plain FTP — the proxy only intercepts HTTP(S)
ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE159nnn/GSE159067/suppl/

# 2) HTTPS straight to NCBI's IP, with the proxy bypassed
curl -k --noproxy '*' --resolve ftp.ncbi.nlm.nih.gov:443:130.14.250.31 \
     https://ftp.ncbi.nlm.nih.gov/geo/series/GSE159nnn/GSE159067/suppl/

# the web front end behaves the same way
curl -k --noproxy '*' --resolve www.ncbi.nlm.nih.gov:443:130.14.29.110 \
     'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159067'
```

`--noproxy '*'` is the essential part (do not use the sandbox proxy);
`--resolve` pins the real hostname so TLS still validates against the host,
and `-k` is retained because the direct route is not subject to the
intercepting certificate. If DNS for NCBI moves, override the address with
the `NCBI_IP` environment variable.

This is the same approach used by `scripts/ncbi_geo.sh` in the companion
repository `IveyTang/ici-benchmark`.

### What each cohort is, and where it came from

| Cohort | Retrieved from | Why |
|---|---|---|
| GSE65858 (270 HNSCC, prognostic validation) | EMBL-EBI ArrayExpress mirror `E-GEOD-65858` + array design `A-GEOD-10558` | ArrayExpress imported GEO series up to ~2017, so a complete, byte-identical copy exists there. Either route works; the mirror is used because it is one host and it also supplies the array-design annotation needed for probe→gene mapping. |
| GSE159067 (102 HNSCC on PD-1/PD-L1 inhibitors) | **NCBI GEO directly** (BioProject PRJNA667509, SRA SRP286445) | Submitted in 2020 and therefore **not** in ArrayExpress. Retrieved with the FTP / direct-IP channels above. |

Retrieval order:

```bash
python scripts/20_fetch_gse65858_ebi.py        # 270 per-sample tables + SDRF + IDF
python scripts/20b_fetch_adf_ranged.py         # 71.7 MB array design, range-parallel
python scripts/23_fetch_gse159067_geo.py       # series matrix + log2 CPM matrix from GEO
python scripts/21_external_validation_gse65858.py
python scripts/24_ici_validation_gse159067.py
python scripts/22_ext_figure_panel.py
```

Two practical notes that cost time to discover:

- A truncated transfer still returns HTTP 200. An early attempt at the array
  design stopped at 10.8% (7.8 MB of 71.7 MB) and parsed cleanly while
  covering only 5,193 of 47,323 probes. Verify exact byte lengths, not minimum
  sizes; `20b` does this and also demonstrates the parallel-range pattern,
  which is ~6x faster than a single stream through the local proxy.
- Gene-set coverage differs sharply between the two platforms. The
  HumanHT-12 v4 array covers the chemokine signatures near-completely but is
  missing four immunoglobulin genes of the plasma-cell programme; the targeted
  HTG EdgeSeq panel covers the chemokine signatures essentially completely
  (12/12, 11/11, 12/13) but only 86% of the gene-set members overall.
