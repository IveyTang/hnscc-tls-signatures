# Analysis plan

The analysis plan was specified before the survival analyses were run. Points at which a decision was made in response to the data are flagged explicitly as **post hoc**.

## Objectives

1. Quantify the concordance among published TLS-related gene signatures in HNSCC.
2. Determine whether a consensus TLS score (CTS) is prognostic and independent of established clinicopathological factors.
3. Characterise which microenvironmental programmes the signatures actually track, using a generic T cell-inflamed signature as a benchmark.

## Cohort

All 528 TCGA-HNSC participants with expanded clinical entities (diagnoses, demographic, exposures) were retrieved from the GDC API. Eligibility required:

1. A diagnosis flagged `diagnosis_is_primary_disease`.
2. A squamous morphology code: 8070/3, 8071/3, 8072/3, 8073/3, 8074/3, 8075/3, 8076/3, 8078/3, 8083/3, 8084/3, or 8052/3.
3. A head-and-neck site of origin (`tissue_or_organ_of_origin`), mapped to oral cavity, oropharynx, larynx, or hypopharynx.
4. Computable overall survival (OS).
5. A primary tumour RNA sequencing profile.

> **Rationale for criterion 1.** 172 of 528 participants (32.6%) carry more than one recorded diagnosis, including second primaries of the prostate, brain, and haematopoietic system. Selecting `diagnoses[0]` would contaminate the cohort with non-HNSCC primaries.

OS was defined as the interval from diagnosis to death from any cause, censored at last follow-up.

## Gene sets

Seven sets were assembled (`config/tls_genesets.json`):

**Entered into the consensus (five TLS-related sets)**

| ID | Set | Genes |
|---|---|---|
| G1 | Coppola 12-chemokine | 12 |
| G2 | 11-chemokine (as applied to hepatocellular carcinoma) | 11 |
| G3 | Published HNSCC 13-chemokine TLS signature | 13 |
| G4 | TLS maturation / germinal centre | 14 |
| G6 | High endothelial venule / lymphoid stroma | 9 |

**Excluded from the consensus by pre-specified decision (two attribution/benchmark sets)**

| ID | Set | Genes | Role |
|---|---|---|---|
| G5 | T helper 1 / interferon-γ effector programme | 12 | Attribution |
| G7 | Generic T cell-inflamed / cytolytic signature | 10 | Positive control |

> **Why exclude G5 and G7.** Inclusion would have conflated general immune inflammation with TLS-specific biology — which is precisely the question under test. G7 in particular serves as the benchmark against which the TLS-specific content of G1–G4 and G6 is judged.

## Scoring

Three algorithms were applied to each gene set:

1. **z-score** — per-gene z-score across samples, then average.
2. **ssGSEA** — single-sample gene set enrichment analysis.
3. **PCA1** — first principal component of the z-scored gene matrix, sign-aligned to the z-score mean.

All three are designed to be platform-agnostic, so that the scores transfer to microarray-based external cohorts without re-fitting. Gene-level TPM values were obtained from GDC STAR-Counts files (`tpm_unstranded`) and transformed as log2(TPM + 1); only primary tumour samples were retained, and samples from the same participant were averaged.

**The consensus TLS score (CTS)** is the mean percentile rank of the five TLS-related sets under the z-score algorithm, which is invariant to platform-specific scaling.

## Cut-off determination

The optimal cut-off maximises the log-rank statistic over candidate thresholds between the 10th and 90th percentiles, subject to a minimum of 15 patients per group. The vectorised implementation was validated against `lifelines` and matched exactly (χ² = 5.72 by both methods).

Stability was assessed by 200 bootstrap resamples with a fixed seed (20261016). The median cut-off and its bootstrap distribution are reported as sensitivity analyses. **The cut-off is not proposed as a validated clinical threshold.**

> **Deviation, recorded after the fact.** This plan specified 200 resamples. The implemented analysis (`scripts/03b_cutpoint_stability.py`) uses **1,000**, so that the 10th and 90th percentiles of the selected threshold could be estimated with useful precision. No other element of the procedure was changed, and the text above is left exactly as pre-specified; this note is the record of the deviation. The vectorised log-rank implementation reproduces the discovery-cohort point estimate exactly (cut-off 0.248926, χ² 5.721466). The manuscript states the deviation and its reason in the Methods.

## Statistical analysis

- Kaplan–Meier estimation with the log-rank test.
- Univariable and multivariable Cox proportional hazards models (`lifelines`), with a ridge penalty of 0.1 and pre-fitting removal of zero-variance covariates. The multivariable model included CTS, age, sex, stage (dummy-coded), subsite, and smoking status.
- Proportional hazards inspected by Schoenfeld residuals.
- Time-dependent AUC at 12, 36, and 60 months using a cumulative/dynamic definition, and the concordance index. CTS is protective in direction, so scores were negated to define a risk score for AUC estimation, with 0.5 as the null value.
- Microenvironmental programmes approximated by marker-based mean z-scores for 12 programmes (B cells, germinal-centre B cells, plasma cells, T follicular helper cells, CD8 T cells, regulatory T cells, cancer-associated fibroblasts, endothelium/high endothelial venules, macrophages, interferon-γ, immune checkpoints, and epithelium as a negative control); associations quantified by Spearman correlation.
- Nomogram with points scaled so that the largest absolute coefficient equals 100; calibration by comparing Kaplan–Meier estimates with mean predicted survival within quartiles of predicted risk; decision curve analysis over threshold probabilities 0.05–0.80 for 36-month survival.
- All tests two-sided, α = 0.05. The signature concordance analysis is descriptive and no multiplicity adjustment was applied; prognostic analyses are reported with effect sizes and confidence intervals and are exploratory.

## HPV analysis

### Why it was added

HPV status is a strong prognostic factor in HNSCC and is known to be associated with TLS enrichment. Any survival advantage attributed to a TLS signature could therefore be confounded by HPV status unless that status is modelled explicitly. This concern is raised in the published literature and was the reason the analysis was undertaken.

### Status determination

HPV status was resolved from three complementary sources, applied hierarchically (ISH > p16 > pathology report):

1. `hpv_status_by_ish_testing` — in situ hybridisation; the most specific assay.
2. `hpv_status_by_p16_testing` — p16 immunohistochemistry; a reliable surrogate in the oropharynx but less so elsewhere.
3. Diagnostic pathology report PDFs — parsed at clause level with optical-character-recognition correction to limit false negatives.

This yielded direct labels for **84 of 438 patients** (23 positive; 60 by ISH, 21 by p16, 3 by report).

For the remaining patients, a transcriptomic classifier was trained on a 107-gene host-expression panel covering the HPV-driven axis (CDKN2A, RB1 and E2F pathway members, cell-cycle and DNA-replication genes), squamous differentiation, and immune programmes. Feature selection and model fitting both occur inside each leave-one-out fold, so no information leaks from the held-out sample.

**Pre-specified sanity checks** — all three had to pass before the classifier could be used:

| Check | Expected | Observed |
|---|---|---|
| Enrichment in the oropharynx | Predicted HPV+ concentrated in the oropharynx | 37/56 (66.1%) vs 31/382 (8.1%); Fisher p = 2.5 × 10⁻²¹ |
| Longer OS in predicted HPV+ | Positive prognostic association | log-rank p = 0.0031; HR 0.506 |
| Higher CTS in predicted HPV+ | More immune-inflamed | median 0.651 vs 0.473; Mann–Whitney p = 1.1 × 10⁻⁴ |

Reported performance: leave-one-out AUC **0.990** (sensitivity 0.957, specificity 0.951) against the 84 annotated patients, versus 0.934 for single-gene CDKN2A; concordant with direct annotation in 82 of 84 cases.

**Intended use.** Predicted status is used for stratification and covariate adjustment only, and never as an outcome.

### Confounding and interaction

Three analyses were specified:

1. Compare the CTS distribution between HPV strata (Mann–Whitney test).
2. Repeat Kaplan–Meier analysis and Cox regression within HPV-positive / HPV-negative strata and within oropharyngeal / non-oropharyngeal strata.
3. Fit the multivariable model with and without HPV status to observe the change in the CTS hazard ratio, and test a first-order interaction term (CTS × HPV) by likelihood-ratio test — in the full cohort using predicted status and in the directly annotated subset.

> **Post hoc note.** The directly annotated subset (n = 84, 33 deaths) was added as a sensitivity analysis after the classifier was built, because predicted status is model-derived. The full-cohort result is primary; the annotated-subset result is reported as a check.

## Software

Python 3.13 with pandas, NumPy, SciPy, statsmodels, scikit-learn, `lifelines`, and `gseapy`. Exact versions in `../requirements-pinned.txt`.
