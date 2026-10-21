# Changelog

All notable changes to this repository are documented here. This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.4.0]

Analysis release accompanying the pre-submission audit of the manuscript. Adds one
multivariable model used throughout, a clinical-only comparator for the decision
curve analysis, bootstrap intervals for discrimination, and English versions of
Figures 1-8.

### Added

- `scripts/30_unified_models.py` — refits the HPV-adjusted estimate on the same
  covariates as the primary multivariable model (age, sex, stage dummy-coded,
  subsite, smoking status; ridge 0.1), adds a clinical-only model to the decision
  curve analysis, and bootstraps the C-index and 36-month AUC differences
  (1,000 resamples). Outputs `results/30_*`.
- `scripts/31_figures_english.py` — redraws Figures 1-8 with English labels at
  300 dpi. The previously deposited figures carried Chinese axis labels and
  in-figure titles, which are not usable in an English-language submission.
  Outputs to `results/fig_en/`.

### Changed

- **HPV-adjusted hazard ratio.** The manuscript now reports one multivariable
  model. Adding HPV status to that model moves the consensus score HR from 0.757
  to 0.810 (95% CI 0.495-1.325); HPV positivity is independently protective
  (HR 0.511, 95% CI 0.335-0.779). The earlier 0.660/0.739 pair came from a
  parsimonious model (age, sex, stage IV) reported under the same label.
- **Decision curve analysis.** A clinical-only curve has been added. Net benefit
  for the combined and clinical-only models is essentially identical (0.225 vs
  0.221, thresholds 0.10-0.65; the combined model is higher at 47% of
  thresholds), which supports the paper's conclusion rather than the earlier
  reading of the figure.
- **Discrimination.** Reported with bootstrap intervals (dC-index -0.010, 95% CI
  -0.034 to 0.014; dAUC at 36 months -0.013, 95% CI -0.050 to 0.018) and
  described as indistinguishable rather than inferior.
- **GSE65858 probe accounting.** Corrected to the measured values: 47,231 probes
  in the array design, 43,961 with a gene symbol, 31,330 present in the sample
  tables, 22,000 genes after collapsing.
- `README.md` finding 5 now reports the unified-model estimate.
- `docs/reproducibility.md` records the model change and the DCA/bootstrap
  additions as cross-checks.

## [1.3.1]

Metadata release accompanying a change of target journal. No analysis, result
table or figure changes; the archived 1.3.0 files are unchanged.

### Changed

- **Manuscript title.** The manuscript is retitled "Chemokine-based tertiary
  lymphoid structure signatures in head and neck cancer are redundant with
  generic inflammation and do not replicate prognostically". The previous title
  claimed only that the signatures are not interchangeable with one another;
  the study also shows that they are redundant with a generic inflammation
  signature and that the prognostic association does not replicate in an
  external cohort. The new title states both. Updated in `README.md`,
  `CITATION.cff`, `.zenodo.json`, `zenodo/metadata.json` and
  `docs/supplementary.md`. The title recorded for the 1.3.0 deposition on Zenodo
  is historical and is **not** altered.
- **Target journal.** The manuscript is submitted to *Frontiers in Immunology*
  (section Cancer Immunity and Immunotherapy) rather than to *Journal of
  Translational Medicine*, following three desk rejections of related work on
  novelty/scope grounds. No content change follows from this.

### Added

- `docs/analysis_plan.md` now carries an explicit record of one deviation from
  the frozen plan: the cut-off stability analysis was pre-specified with 200
  bootstrap resamples and was run with 1,000. The pre-specified text is left
  verbatim and the deviation is noted beneath it. The manuscript discloses the
  same deviation in the Methods.

## [1.3.0]

Adds the cut-off stability analysis, the cross-platform validation of the
transcriptomic HPV classifier, and the supplementary reporting materials.

### Added

- **Cut-off stability** (`scripts/03b_cutpoint_stability.py`): 1,000 bootstrap
  resamples in which the *entire* cut-off selection procedure is repeated
  (60 candidate thresholds, minimum 15 patients per group). Reproducibility is
  also assessed *without* re-selection, by applying the discovery cut-off
  unchanged to each resample, because a scan over many thresholds will nearly
  always locate some nominally significant cut-off. Outputs
  `results/03b_cutpoint_stability.csv` (per-resample) and `.json` (summary).
  The vectorised log-rank implementation reproduces the discovery point
  estimate exactly (cut-off 0.248926, χ² 5.721466).
- **Cross-platform validation of the HPV classifier**
  (`scripts/25_hpv_classifier_external.py`): the eight genes and coefficients
  selected under nested leave-one-out in TCGA-HNSC are transplanted *without
  re-fitting* into GSE65858, where HPV16 status was resolved independently at
  both the DNA and the viral E6\*I RNA level. Because ROC-AUC is invariant to
  the intercept (which was not recorded in the discovery analysis), the
  comparison is exact; sensitivity and specificity are re-derived within the
  external cohort at Youden thresholds and labelled as such. Outputs
  `results/25_hpv_classifier_external.csv` / `.json` and
  `results/fig_hpv_classifier_external.png`.
- **Supplementary Figure S1** (`scripts/26_supp_figure_hpv_roc.py`): the
  three-panel ROC figure (primary DNA+RNA+ versus DNA−; any HPV16 DNA+
  versus DNA−; oropharynx only), each with the transplanted classifier, the
  unweighted panel mean and CDKN2A alone. Output
  `results/fig_supp_hpv_roc.png`.
- **`docs/supplementary.md`**: the manuscript's supplementary material — the
  20-item REMARK reporting checklist mapped item by item to the manuscript,
  the bootstrap cut-off distribution, and the Supplementary Figure S1 legend.

### Changed

- `docs/reproducibility.md`: corrected the statement that NCBI/GEO was
  unreachable. That claim was an artefact of routing all requests through an
  HTTP proxy that returns HTTP 000 for NCBI hosts; GEO is reachable over plain
  FTP or over HTTPS with the proxy bypassed. Added a Path C covering the
  cut-off stability and HPV cross-platform analyses.

## [1.2.0]

Adds the immunotherapy (ICI) validation cohort and corrects the earlier claim
that NCBI/GEO was unreachable.

### Added

- **GSE159067** (102 advanced HNSCC patients on PD-1/PD-L1 inhibitors),
  retrieved **directly from GEO**:
  - `scripts/23_fetch_gse159067_geo.py` — retrieval over plain FTP and, as a
    fallback, HTTPS to NCBI's IP with the proxy bypassed.
  - `scripts/24_ici_validation_gse159067.py` — response, survival and
    phenotype validation with the same gene sets and scoring, no re-tuning.
- Derived data `data/derived/gse159067_clinical_and_cts.csv`,
  `results/24_ici_*`, `results/table7_ici_validation.csv`,
  `results/fig_ici_validation.png`.

### Findings

- Chemokine signatures were again nearly interchangeable (ρ = 0.990–0.993) and
  correlated as strongly with a generic T cell-inflamed signature
  (ρ = 0.880–0.883) as in the discovery cohort.
- CTS tracked the immune "hot phenotype" published with this cohort
  (ρ = 0.737; AUC 0.874 for the HOT/COLD call), and predicted objective
  response (OR 2.52 per SD, p = 0.017) and longer PFS (HR 0.735, p = 0.002).
- The generic T cell-inflamed control matched it on every endpoint and was
  better for objective response (OR 2.71; AUC 0.767) and overall survival
  (HR 0.686).

### Corrected

- **The claim that NCBI/GEO was unreachable was wrong.** All earlier requests
  went through the sandbox HTTP proxy, which returns `000` for NCBI hosts and
  whose egress IP NCBI's abuse filter rejects. Plain FTP, or HTTPS with
  `--noproxy '*'` and `--resolve` to NCBI's IP, both work. `scripts/README.md`
  now documents this instead of the mirror workaround.
- References 36–49 had been mis-numbered as 1–14 in an earlier edit; the whole
  list is now renumbered sequentially and a numbering invariant was added to
  the compliance scan.

---

## [1.1.0]

Adds the external validation cohort and removes the last blocked analysis step.

### Added

- **External validation in GSE65858** (270 HNSCC patients, 94 deaths, Illumina HumanHT-12 v4),
  retrieved from the **EMBL-EBI ArrayExpress mirror** as `E-GEOD-65858` because NCBI was
  unreachable from the analysis environment.
  - `scripts/20_fetch_gse65858_ebi.py` — parallel retrieval of 270 per-sample tables, the
    SDRF and the IDF, with exact byte-length verification.
  - `scripts/20b_fetch_adf_ranged.py` — 71.7 MB array design via parallel HTTP Range
    requests, verifying the full 47,323-probe coverage.
  - `scripts/21_external_validation_gse65858.py` — identical gene sets and scoring, no
    re-tuning: convergence, microenvironmental programme correlation, survival,
    HPV association, and a check against the expression clusters published with the cohort.
  - `scripts/22_ext_figure_panel.py` — four-panel validation figure and summary table.
- Derived validation data: `data/derived/gse65858_clinical_and_cts.csv`,
  `data/gse65858_sdrf.txt`, `results/21_ext_*`, `results/table6_external_validation.csv`,
  `results/fig_ext_*.png`.

### Findings

- **Replicated:** chemokine signature redundancy (ρ = 0.974–0.990); the score's coupling to
  checkpoint expression (ρ = 0.827) rather than to germinal-centre B cells (ρ = 0.146) or
  plasma cells (ρ = 0.077); and the HPV association (median 0.670 in HPV16 DNA+RNA+ vs 0.445
  in DNA−, p = 0.0012, against 0.672 vs 0.472 in the discovery cohort).
- **Did not replicate:** the prognostic association (HR 1.002 per SD, p = 0.985; all cut-offs
  null), although the individual gene sets and the positive control retained cut-off-dependent
  associations in the same cohort (p = 0.008, 0.046 and 0.010).
- The immune-response expression cluster published with GSE65858 carried the highest
  consensus score (median 0.742; Kruskal–Wallis p = 3 × 10⁻¹⁵).

### Notes

- `GSE159067` (immunotherapy-treated cohort in the original analysis plan) has **no**
  ArrayExpress mirror, having been submitted after ArrayExpress stopped importing GEO. It is
  not part of this validation and no immunotherapy-response claim rests on it.
- A truncated HTTP transfer still returns status 200: an initial attempt at the array design
  silently covered only 5,193 of 47,323 probes. All retrieval steps now verify exact sizes.

---

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
