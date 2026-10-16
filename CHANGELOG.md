# Changelog

All notable changes to this repository are documented here. This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] — 2026-10-16

Initial release accompanying the manuscript *"Tertiary lymphoid structure gene signatures are not interchangeable in head and neck squamous cell carcinoma: a retrospective multi-signature transcriptomic analysis"* (submitted to *Journal of Translational Medicine*).

**Archived release.** Concept DOI [`10.5281/zenodo.23415040`](https://doi.org/10.5281/zenodo.23415040) · version DOI for 1.0.0 [`10.5281/zenodo.23415041`](https://doi.org/10.5281/zenodo.23415041) · record <https://zenodo.org/records/23415041>. Archived file: `IveyTang/hnscc-tls-signatures-v1.0.0.zip` (3,251,526 bytes, md5 `4bf07ad6f2d04be5b2d43e9d5429c976`).

### Pipeline

- **Cohort construction** — GDC retrieval, three-way primary-tumour filter (primary-disease flag + squamous morphology code + head-and-neck site), and stage fallback, yielding 438 of 528 TCGA-HNSC patients. Handles two schema pitfalls that silently corrupt the cohort: multiple primary diagnoses and the deprecated `demographic.gender` field.
- **Expression** — STAR-Counts `tpm_unstranded` extraction with concurrency and local caching, producing a 1291-gene × 563-sample matrix.
- **Scoring** — seven gene sets scored under three platform-agnostic algorithms (z-score, ssGSEA, PCA1), combined into a consensus TLS score from the five TLS-related sets.
- **Internal validation** — Kaplan–Meier and log-rank analysis, univariable and multivariable Cox regression, time-dependent AUC, concordance index, and signature convergence analysis.
- **Microenvironment and clinical utility** — twelve marker-based programmes, nomogram, calibration, and decision curve analysis.
- **HPV status determination** — hierarchical annotation from in situ hybridisation, p16 immunohistochemistry, and clause-level pathology-report parsing with OCR correction (84 patients), plus a 107-gene nested leave-one-out transcriptomic classifier for the remainder (leave-one-out AUC 0.990).
- **HPV-stratified analysis** — stratified Kaplan–Meier and Cox analysis, confounding assessment, and formal interaction testing.

### Results included

- All result tables (`results/*.csv`), figures (`results/fig_*.png`), machine-readable summaries (`results/*.json`), and pipeline run logs (`results/*.log`).
- Derived patient-level analysis tables (`data/derived/`).

The `results/` directory is committed flat, matching the scripts' own output paths, so re-running the pipeline overwrites in place and leaves `git status` clean.

### Known issues

- **External validation not included.** The external cohorts (GSE65858, GSE159067) could not be retrieved because NCBI/GEO was unreachable from the analysis environment (HTTP 000). `scripts/05_external_validation.py` is complete and ready to run once the data are supplied. Findings are reported as discovery-only.
- **`10_hpv_status.py` cannot run offline.** It depends on raw BCR clinical XML and pathology-report files retrieved from the GDC at runtime. Its output is committed, so downstream steps run without it.
- **HPV classifier scope.** Validated only against the 84 directly annotated patients; performance outside the oropharynx is less certain. Intended for stratification and adjustment, never as a diagnostic test or outcome variable.

### Documentation

- `docs/analysis_plan.md` — pre-specified analysis plan with post-hoc decisions flagged.
- `docs/data_dictionary.md` — column-by-column reference for all derived tables.
- `docs/reproducibility.md` — environment, runtimes, determinism, and verification steps.
- `zenodo/DEPOSIT_GUIDE.md` — archival instructions for all three deposit routes.
