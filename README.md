# hnscc-tls-signatures

**Tertiary lymphoid structure gene signatures are not interchangeable in head and neck squamous cell carcinoma: a retrospective multi-signature transcriptomic analysis**

This repository contains the complete, reproducible analysis pipeline, derived data, and result tables and figures for the manuscript above, submitted to *Journal of Translational Medicine*.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23415040.svg)](https://doi.org/10.5281/zenodo.23415040)

| | |
|---|---|
| **Repository** | https://github.com/IveyTang/hnscc-tls-signatures |
| **Archived release** | [10.5281/zenodo.23415040](https://doi.org/10.5281/zenodo.23415040) — concept DOI, always resolves to the latest version<br>v1.0.0: [10.5281/zenodo.23415041](https://doi.org/10.5281/zenodo.23415041) |
| **Version** | 1.0.0 (see [`CHANGELOG.md`](CHANGELOG.md)) |
| **Licence** | Code MIT · derived data CC BY 4.0 |

---

## What this study shows

Chemokine-based transcriptional "TLS signatures" are widely used as scalable surrogates for tertiary lymphoid structures (TLS) in bulk tumour expression data. Using 438 patients from TCGA-HNSC and seven gene sets scored under three platform-agnostic algorithms, and then validating the result in an independent 270-patient cohort (GSE65858, retrieved from the EMBL-EBI ArrayExpress mirror as E-GEOD-65858), we show that:

1. **The three chemokine-based TLS signatures are nearly redundant** with one another (Spearman ρ ≥ 0.987).
2. **They correlate almost as strongly with a generic T cell-inflamed signature** (ρ = 0.89) as with each other, indicating that they predominantly capture general immune inflammation rather than TLS neogenesis.
3. **A consensus TLS score (CTS) tracks checkpoint (ρ = 0.898) and regulatory T cell (ρ = 0.889) programmes far more strongly** than the structural features of TLS — germinal-centre B cells (ρ = 0.462), plasma cells (ρ = 0.496), or high endothelial venules (ρ = 0.507).
4. **Prognostically, CTS is cut-off dependent and not independent of stage and age** (multivariable HR 0.757, 95% CI 0.463–1.238, p = 0.268), and its discrimination (C-index 0.536) does not exceed that of a generic T cell-inflamed signature (0.546).
5. **HPV status does not confound the score's behaviour.** CTS is markedly higher in HPV-positive tumours (median 0.672 vs 0.472, p = 2.9 × 10⁻⁵), but adjusting for HPV does not attenuate the CTS hazard ratio (0.660 → 0.739) and the interaction test is null (p = 0.881). Notably, only 5 of 68 HPV-positive patients fall below the optimal cut-off, so a threshold-based chemokine score cannot discriminate within HPV-positive disease.
6. **The structural findings replicate externally while the prognostic claim does not.** In GSE65858 the chemokine signatures remained nearly redundant (ρ = 0.974–0.990) and again tracked checkpoint expression (ρ = 0.827) far more than germinal-centre B cells (ρ = 0.146) or plasma cells (ρ = 0.077); the immune-response expression cluster published with that cohort carried the highest CTS (median 0.742; Kruskal–Wallis p = 3 × 10⁻¹⁵); and the HPV association reproduced almost exactly (median 0.670 vs 0.445, p = 0.0012). Prognosis, however, did **not** replicate (HR 1.002, p = 0.985) — even though the individual gene sets composing the score, and the generic T cell-inflamed positive control, each retained cut-off-dependent associations in the same cohort (p = 0.008, 0.046 and 0.010 respectively).

7. **In an immunotherapy cohort the score tracks a published immune phenotype and predicts benefit, but no better than generic inflammation.** In 102 patients with advanced HNSCC treated with PD-1/PD-L1 inhibitors (GSE159067) the chemokine signatures were again nearly interchangeable (ρ = 0.990–0.993) and correlated as strongly with a generic T cell-inflamed signature (ρ = 0.880–0.883) as in the discovery cohort. The consensus score correlated strongly with the immune "hot phenotype" score published with that cohort (ρ = 0.737; AUC 0.874 for its HOT/COLD call) and predicted objective response (OR 2.52 per SD, p = 0.017) and longer progression-free survival (HR 0.735, p = 0.002). But the generic T cell-inflamed signature matched it everywhere and beat it on objective response (OR 2.71, AUC 0.767) and overall survival (HR 0.686).

Together, these findings provide quantitative support for the long-recognised lack of standardisation in TLS scoring and indicate that chemokine signatures should not be treated as interchangeable surrogates of TLS.

---

## Repository layout

```
config/
  tls_genesets.json          Gene set definitions (7 sets) and scoring strategy
scripts/
  00_fetch_gdc_manifest.py   Retrieve clinical records and file manifest from GDC
  01_build_clinical_table.py Build cohort with the three-way primary-tumour filter
  02_download_expression.py  Download STAR-Counts files and extract TPM values
  03_build_tls_score.py      Score seven gene sets by three algorithms; build CTS
  03b_cutpoint_stability.py  Bootstrap the cut-off selection (1,000 resamples)
  04_internal_validation.py  Kaplan-Meier, Cox, time-dependent AUC, C-index, convergence
  05_external_validation.py  Generic external-validation harness (superseded by 21 / 24)
  06_microenvironment_nomogram.py  Marker programmes, nomogram, calibration, DCA
  07_fetch_bcr_hpv.py        Parse HPV fields from BCR clinical XML
  08_fetch_path_reports.py   Retrieve diagnostic pathology report PDFs
  08b_reclassify_reports.py  Clause-level p16/HPV reclassification with OCR correction
  09_extract_hpv_panel.py    Extract the 107-gene HPV host-expression panel
  10_hpv_status.py           Three-source HPV labelling + nested LOO classifier
  11_hpv_stratified.py       HPV-stratified survival analysis and interaction testing
  20_fetch_gse65858_ebi.py   Fetch GSE65858 from the EBI ArrayExpress mirror
  20b_fetch_adf_ranged.py    Parallel-range download of the HT-12 v4 array design
  21_external_validation_gse65858.py  Prognostic + phenotypic validation, no re-tuning
  22_ext_figure_panel.py     Compose the external-validation figure and summary table
  23_fetch_gse159067_geo.py  Fetch the immunotherapy cohort directly from GEO
  24_ici_validation_gse159067.py  Response, survival and phenotype validation
  25_hpv_classifier_external.py   Transplant the HPV classifier into GSE65858
  26_supp_figure_hpv_roc.py  Build Supplementary Figure S1 (three-panel ROC)
  exploratory/               Ad-hoc GDC field-discovery probes (not part of the pipeline)
data/
  expr_download_meta.json    Provenance for the expression download
  derived/                   Derived, patient-level analysis tables (see data/README.md)
results/                     All outputs, flat, exactly as the scripts write them
  *.csv                        Result tables
  *.json                     Machine-readable result summaries
  fig_*.png                  Figures
  *.log                      Pipeline run logs
docs/
  analysis_plan.md           Pre-specified analysis plan
  data_dictionary.md         Column-by-column data dictionary
  reproducibility.md         How to reproduce from scratch, including runtime
  supplementary.md           REMARK checklist, bootstrap cut-off table, Figure S1 legend
zenodo/
  metadata.json              Deposit metadata for Zenodo
  DEPOSIT_GUIDE.md           Step-by-step archiving instructions
```

> The scripts write directly into `results/` with no subdirectories. The committed layout mirrors that exactly, so re-running any step overwrites its own output in place rather than scattering files — `git status` stays clean when the pipeline is re-run on the committed data.

---

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

The pipeline is numbered and should be run in order. Steps 00–02 acquire data from the Genomic Data Commons and require network access plus roughly 2.5 GB of free disk space; steps 03–11 operate on the local derived tables.

```bash
python scripts/00_fetch_gdc_manifest.py
python scripts/01_build_clinical_table.py
python scripts/02_download_expression.py --workers 10   # ~8 min, ~2.3 GB
python scripts/03_build_tls_score.py
python scripts/04_internal_validation.py
python scripts/06_microenvironment_nomogram.py
python scripts/07_fetch_bcr_hpv.py
python scripts/08_fetch_path_reports.py
python scripts/08b_reclassify_reports.py
python scripts/09_extract_hpv_panel.py
python scripts/10_hpv_status.py
python scripts/11_hpv_stratified.py
```

If you only want to reproduce the statistical results, restore the committed derived tables and skip the network-dependent steps:

```bash
python scripts/00b_restore_derived.py           # data/derived/ -> data/processed/
python scripts/03_build_tls_score.py
python scripts/04_internal_validation.py
python scripts/06_microenvironment_nomogram.py
python scripts/11_hpv_stratified.py
```

See `docs/reproducibility.md` for runtime, random seeds, and software versions.

---

## Data availability

All primary data are publicly available from the **Genomic Data Commons** (https://portal.gdc.cancer.gov/) under project **TCGA-HNSC**. No controlled-access or identifiable data were used, and no additional ethics approval was required.

Derived patient-level tables generated by this pipeline are included under `data/derived/` and are released under CC BY 4.0. The raw STAR-Counts expression matrix (2.3 GB) is not committed but can be regenerated in about eight minutes with `scripts/02_download_expression.py`.

**External validation data.** Both validation cohorts are public. **GSE65858** (270 HNSCC) is available from NCBI GEO and, identically, from the **EMBL-EBI ArrayExpress** mirror as **E-GEOD-65858** with the array design **A-GEOD-10558**; **GSE159067** (102 HNSCC on PD-1/PD-L1 inhibitors) is available from NCBI GEO directly (BioProject PRJNA667509) and has no ArrayExpress mirror, having been submitted after ArrayExpress stopped importing GEO. An earlier version of this README described NCBI/GEO as unreachable; that was an artefact of routing all requests through the sandbox HTTP proxy, and `scripts/README.md` documents the two channels (plain FTP, and HTTPS with the proxy bypassed) that work. See `scripts/README.md` for the retrieval order and the two pitfalls (silent truncation of very large files, and the fact that `GSE159067` has no ArrayExpress mirror because it was submitted after ArrayExpress stopped importing GEO).

The harmonised clinical table and the per-patient validation scores are committed under `data/derived/gse65858_clinical_and_cts.csv`; the 270 per-sample expression tables (~225 MB) are not committed but are regenerated by `scripts/20_fetch_gse65858_ebi.py` in about four minutes.

---

## Data handling notes

Two aspects of the GDC schema required specific handling and are implemented in `scripts/01_build_clinical_table.py`. They are worth flagging because they can silently corrupt a TCGA-HNSC cohort:

- **Multiple primary diagnoses.** 172 of 528 participants (32.6%) carry more than one recorded diagnosis, including second primaries of the prostate, brain, and haematopoietic system. Selecting `diagnoses[0]` contaminates the cohort. We retain only the diagnosis flagged `diagnosis_is_primary_disease` that also has a squamous morphology code (8070/3, 8071/3, 8072/3, 8073/3, 8074/3, 8075/3, 8076/3, 8078/3, 8083/3, 8084/3, 8052/3) and a head-and-neck site of origin.
- **Deprecated sex field.** The legacy field `demographic.gender` returns no values; `demographic.sex_at_birth` is complete and is used instead. A zero-variance covariate from the deprecated field will cause `lifelines` to fail with a `delta contains nan` error.

A third note concerns HPV annotation: direct labels exist for only 84 of 438 patients (in situ hybridisation 60, p16 immunohistochemistry 21, pathology report text 3). The remaining patients are classified by a 107-gene transcriptomic classifier trained under nested leave-one-out cross-validation (LOO-AUC 0.990, concordant with direct annotation in 82 of 84 cases). Predicted status is used for stratification and adjustment only, never as an outcome.

---

## Citation

If you use this code or data, please cite the manuscript and the archived release:

```
Tang X, Zhang Y, He Y, et al. Tertiary lymphoid structure gene signatures are not
interchangeable in head and neck squamous cell carcinoma: a retrospective
multi-signature transcriptomic analysis. Journal of Translational Medicine. (in submission)
```

The repository is citable at https://github.com/IveyTang/hnscc-tls-signatures, and release 1.0.0 is archived at https://doi.org/10.5281/zenodo.23415040 (concept DOI). See `CITATION.cff` for machine-readable metadata.

```
Tang X, Zhang Y, He Y. hnscc-tls-signatures: analysis code and derived data for
"Tertiary lymphoid structure gene signatures are not interchangeable in head and
neck squamous cell carcinoma" (v1.0.0). Zenodo. 2026.
https://doi.org/10.5281/zenodo.23415040
```

To archive a revision, bump `version` in `.zenodo.json` and `CITATION.cff`, commit, then run `python zenodo/create_release.py --github-release` with a new tag — Zenodo files it under the same concept DOI.

---

## Licence

- **Code** (`scripts/`, `config/`): MIT — see `LICENSE`.
- **Derived data and results** (`data/derived/`, `results/`): CC BY 4.0 — see `LICENSE-DATA`.

---

## Contact

Xi Tang, MD — Department of Head and Neck Surgery, Chongqing University Cancer Hospital, Chongqing, China
Email: 33479986@qq.com · ORCID: [0000-0003-2740-1368](https://orcid.org/0000-0003-2740-1368)
