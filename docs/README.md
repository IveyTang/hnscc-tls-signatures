# Repository structure

This file is a map of the repository. For the study summary and quick start, see the top-level [README](../README.md).

```
hnscc-tls-signatures/
|
|-- README.md                  Study summary, quick start, data availability, citation
|-- CHANGELOG.md               Release history
|-- CITATION.cff               Machine-readable citation metadata
|-- LICENSE                    MIT  - applies to code
|-- LICENSE-DATA               CC BY 4.0 - applies to derived data and results
|-- .zenodo.json               Zenodo metadata for the GitHub integration
|-- requirements.txt           Minimum-version dependencies
|-- requirements-pinned.txt    Exact versions used for the published analysis
|-- .gitignore
|
|-- config/
|   `-- tls_genesets.json      The seven gene sets, their membership, and the scoring strategy
|
|-- scripts/
|   |-- 00_fetch_gdc_manifest.py          GDC clinical records + file manifest
|   |-- 00b_restore_derived.py            Restore committed tables for offline reproduction
|   |-- 01_build_clinical_table.py        Cohort construction (three-way primary filter)
|   |-- 02_download_expression.py         STAR-Counts download and TPM extraction
|   |-- 03_build_tls_score.py             Gene set scoring + consensus score
|   |-- 04_internal_validation.py         Survival, discrimination, convergence
|   |-- 05_external_validation.py         External validation (blocked - see README)
|   |-- 06_microenvironment_nomogram.py   Programmes, nomogram, calibration, DCA
|   |-- 07_fetch_bcr_hpv.py               HPV fields from BCR clinical XML
|   |-- 08_fetch_path_reports.py          Diagnostic pathology report retrieval
|   |-- 08b_reclassify_reports.py         Clause-level p16 reclassification
|   |-- 09_extract_hpv_panel.py           107-gene HPV host-expression panel
|   |-- 10_hpv_status.py                  HPV labelling + nested LOO classifier
|   |-- 11_hpv_stratified.py              HPV-stratified survival and interaction
|   |-- README.md                         Per-script reference
|   `-- exploratory/                      Ad-hoc GDC schema probes (not part of the pipeline)
|
|-- data/
|   |-- README.md              Data directory guide
|   |-- expr_download_meta.json  Expression download provenance
|   `-- derived/               Derived patient-level tables (committed, CC BY 4.0)
|       |-- hnsc_analysis_set.csv
|       |-- hnsc_gene_matrix.csv
|       |-- hnsc_gene_matrix_hpv.csv
|       |-- hnsc_hpv_status.csv
|       |-- hnsc_tls_scores.csv
|       `-- hnsc_with_hpv.csv
|
|-- results/                   All outputs, flat, exactly as the scripts write them
|   |-- 03_score_summary.json, 04_*.csv, 06_*.csv, 10_*.json, 11_*.csv, ...
|   |-- signature_convergence.csv, table1_baseline.txt
|   |-- fig_*.png             All figures
|   `-- *.log                 Pipeline run logs
|
|-- docs/
|   |-- analysis_plan.md       Pre-specified plan; post-hoc decisions flagged
|   |-- data_dictionary.md     Column-by-column reference
|   `-- reproducibility.md     Environment, runtimes, determinism, verification
|
`-- zenodo/
    |-- metadata.json          Deposit-ready metadata (web form or REST API)
    `-- DEPOSIT_GUIDE.md       Archival instructions for all three routes
```

## Where to look for what

| I want to… | Start here |
|---|---|
| Understand the study in two minutes | [`../README.md`](../README.md) |
| Reproduce the figures and tables | [`../docs/reproducibility.md`](../docs/reproducibility.md) |
| Know what a column means | [`../docs/data_dictionary.md`](../docs/data_dictionary.md) |
| Understand why an analysis choice was made | [`analysis_plan.md`](analysis_plan.md) |
| Reanalyse the data myself | `../data/derived/hnsc_with_hpv.csv` |
| Archive a release on Zenodo | [`../zenodo/DEPOSIT_GUIDE.md`](../zenodo/DEPOSIT_GUIDE.md) |
| Know which script does what | [`../scripts/README.md`](../scripts/README.md) |

## Conventions

- **Paths** are resolved relative to the repository root inside every script, so scripts run correctly from any working directory.
- **Numbering** reflects execution order. `00b` was inserted after `00` rather than renumbering the pipeline, so that the numbering in the manuscript's Methods continues to match.
- **Determinism** — the random seed is fixed at `20261016` for all resampling.
- **Outputs** are written to `results/` and are safe to regenerate; committed copies exist so the results can be verified without network access.
