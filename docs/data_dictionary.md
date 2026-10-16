# Data dictionary

Column-by-column reference for the derived tables in `data/derived/`. All tables are keyed on `barcode` (the TCGA participant barcode) unless noted.

---

## 1. Gene set scores

Present in `hnsc_tls_scores.csv` and `hnsc_with_hpv.csv`.

| Column | Type | Description |
|---|---|---|
| `G1_Coppola_12chemokine` | float | Coppola 12-chemokine signature, z-score algorithm |
| `G2_11chemokine_HCC` | float | 11-chemokine signature as applied to hepatocellular carcinoma, z-score algorithm |
| `G3_HNSCC_13chemokine` | float | Published HNSCC 13-chemokine TLS signature, z-score algorithm |
| `G4_TLS_maturation_GC` | float | TLS maturation / germinal-centre set, z-score algorithm |
| `G5_Th1_IFNgamma` | float | T helper 1 / interferon-γ effector programme (attribution only; excluded from the consensus) |
| `G6_HEV_lymphoid_stroma` | float | High endothelial venule / lymphoid stroma set, z-score algorithm |
| `G7_Tcell_inflamed_cytolytic` | float | Generic T cell-inflamed / cytolytic signature (positive control; excluded from the consensus) |
| `*__pca1` | float | First principal component of the same gene set, sign-aligned to the z-score mean |

Gene set membership is defined in `config/tls_genesets.json`.

## 2. Consensus TLS score and groupings

| Column | Type | Description |
|---|---|---|
| `CTS` | float | **Consensus TLS score.** Mean percentile rank of the five TLS-related gene sets (G1, G2, G3, G4, G6) under the z-score algorithm. Range 0–1; higher = more immune-inflamed |
| `CTS_group_opt` | str | `High` / `Low` at the optimal cut-off (CTS ≥ 0.249) |
| `CTS_group_median` | str | `High` / `Low` at the median cut-off (sensitivity analysis) |

The optimal cut-off was selected by maximising the log-rank statistic over candidate thresholds between the 10th and 90th percentiles, subject to at least 15 patients per group. It is unstable under bootstrap resampling (median 0.437, IQR 0.245–0.697) and is reported for descriptive purposes only.

## 3. Cohort construction provenance

| Column | Type | Description |
|---|---|---|
| `barcode` | str | TCGA participant barcode (12 characters) |
| `case_id` | str | GDC case UUID |
| `n_diagnoses` | int | Number of diagnoses recorded for the participant in the GDC |
| `multi_primary` | bool | True if `n_diagnoses` > 1 |
| `diagnosis_rule` | str | Filter that selected the retained diagnosis; `primary+鳞癌+头颈` = primary-disease flag + squamous morphology + head-and-neck site |
| `is_squamous` | bool | Whether the retained diagnosis carries a squamous morphology code |
| `morphology` | str | ICD-O-3 morphology code (e.g. `8070/3`, `8071/3`) |
| `primary_diagnosis` | str | Free-text primary diagnosis |
| `tissue_or_organ_of_origin` | str | GDC site-of-origin field (used because `site_of_primary_disease` is empty) |
| `subsite` | str | Harmonised subsite: 口腔 (oral cavity), 口咽 (oropharynx), 喉 (larynx), 下咽 (hypopharynx) |

## 4. Demographics and exposures

| Column | Type | Description |
|---|---|---|
| `sex_at_birth` | str | `male` / `female`. **Note:** the legacy field `gender` is deprecated and returns no values |
| `gender` | str | Deprecated field, retained for transparency; expected to be empty |
| `age_at_index` | float | Age in years at diagnosis |
| `race` | str | Self-reported race |
| `ethnicity` | str | Self-reported ethnicity |
| `year_of_diagnosis` | int | Calendar year of diagnosis |
| `smoking_status` | str | GDC smoking status category |
| `cigarettes_per_day` | float | Cigarettes per day |
| `pack_years_x` | float | Pack-years (primary encoding; used in the analysis) |
| `pack_years_y` | float | Pack-years (alternate encoding, present in the HPV tables only) |
| `tobacco_smoking_history` | float | Coded tobacco smoking history (present in HPV tables only) |
| `alcohol_history` | str | Alcohol use history |
| `alcohol_intensity` | str | Alcohol use intensity |

## 5. Tumour and staging

| Column | Type | Description |
|---|---|---|
| `tumor_grade` | str | Pathological grade (G1–G4, GX) |
| `histologic_grade` | str | Histological grade as recorded in the BCR XML (HPV tables only) |
| `stage` | str | AJCC stage. Pathological stage where available (70.6% of patients), otherwise clinical stage |
| `t_stage`, `n_stage`, `m_stage` | str | AJCC T, N, and M classification |
| `margin_status` | str | Surgical margin status (HPV tables only) |
| `lvi` | str | Lymphovascular invasion (HPV tables only) |
| `pni` | str | Perineural invasion (HPV tables only) |
| `ln_examined`, `ln_positive` | float | Lymph nodes examined / positive (HPV tables only) |
| `ene` | str | Extranodal extension (HPV tables only) |
| `anatomic_subdivision` | str | BCR anatomic subdivision (HPV tables only) |
| `bcr_cancer_status` | str | BCR cancer status at last follow-up (HPV tables only) |

## 6. Survival endpoints

| Column | Type | Description |
|---|---|---|
| `vital_status` | str | `Alive` / `Dead` |
| `days_to_death` | float | Days from diagnosis to death (empty if alive) |
| `days_to_last_follow_up` | float | Days from diagnosis to last follow-up |
| `days_to_recurrence` | float | Days from diagnosis to recurrence, where recorded |
| `os_time_days` | float | **Overall survival time in days** (death or last follow-up) |
| `os_event` | int | **Event indicator** — 1 = death, 0 = censored |

## 7. HPV status

Present in `hnsc_hpv_status.csv` and `hnsc_with_hpv.csv`.

| Column | Type | Description |
|---|---|---|
| `hpv_ish` | str | HPV status by in situ hybridisation, from BCR XML (`hnsc:hpv_status_by_ish_testing`) |
| `hpv_p16` | str | HPV status by p16 immunohistochemistry, from BCR XML (`hnsc:hpv_status_by_p16_testing`) |
| `hpv_from_report` | str | HPV/p16 call parsed from the diagnostic pathology report PDF text |
| `hpv_label` | str | Direct annotation, resolved hierarchically (ISH > p16 > report). Populated for 84 patients (23 positive, 61 negative) |
| `hpv_label_source` | str | `ISH` (60), `p16` (21), or `pathology report` (3) |
| `hpv_prob` | float | Classifier probability of HPV positivity (all 438 patients) |
| `hpv_pred` | str | Classifier call, `Positive` / `Negative` at a threshold of 0.392 (68 positive, 370 negative) |
| `hpv_final` | str | **Status used in the manuscript** — direct annotation where available, otherwise the classifier |
| `hpv_source` | str | Provenance of `hpv_final`: `ISH`, `p16`, `pathology report`, or `classifier` |

### The HPV classifier

Trained on a 107-gene host-expression panel covering the HPV-driven axis (CDKN2A, RB1 and E2F pathway members, cell-cycle and DNA-replication genes), squamous differentiation, and immune programmes. Feature selection and fitting both occur **inside** each leave-one-out fold, so no information leaks from the held-out sample.

- Leave-one-out AUC **0.990** (sensitivity 0.957, specificity 0.951) against the 84 directly annotated patients
- Single-gene CDKN2A reference baseline AUC 0.934
- Concordant with direct annotation in 82 of 84 cases (97.6%)

**Pre-specified sanity checks** (all passed): predicted HPV-positive tumours must be enriched in the oropharynx (37/56, 66.1% vs 31/382, 8.1%; Fisher p = 2.5 × 10⁻²¹), have longer OS (log-rank p = 0.0031), and have higher CTS (median 0.651 vs 0.473; Mann–Whitney p = 1.1 × 10⁻⁴).

**Intended use:** stratification and covariate adjustment only. The classifier output is **not** a validated diagnostic test and must not be used as an outcome variable or in a clinical setting.

The 8 selected genes and their logistic-regression coefficients are recorded in `results/10_hpv_summary.json`.

---

## Derived expression matrix

`hnsc_gene_matrix_hpv.csv` is a transposed, wide matrix: one row per gene (107 rows), one column per TCGA aliquot, with the sample type appended to the column name (e.g. `TCGA-CV-6961-01A|Primary Tumor`). Values are log2(TPM + 1).
