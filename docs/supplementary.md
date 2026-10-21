# Supplementary material

**Manuscript title:** Tertiary lymphoid structure gene signatures are not interchangeable in head and neck squamous cell carcinoma: a retrospective multi-signature transcriptomic analysis

This file contains one reporting checklist table, one additional results table, and one supplementary figure.

---

## Supplementary Table S1. REMARK reporting checklist

Checklist items are reproduced from the REporting recommendations for tumour MARKer prognostic studies (REMARK) [51]. Section headings and item numbering follow the original publication.

| Item | Recommendation | Where addressed |
|---|---|---|
| **Introduction** | | |
| 1 | State the marker examined, the study objectives, and any prespecified hypotheses. | Introduction, final paragraph (four stated objectives: inter-signature concordance; prognostic independence of the consensus TLS score; attribution of the signals to specific microenvironmental programmes; and reproducibility of the cut-off under resampling). Prespecified hypotheses are also flagged in Methods (gene sets; cut-off determination and stability) and in the abstract. |
| **Materials and Methods** | | |
| **Patients** | | |
| 2 | Describe the characteristics of the study patients, including their source and inclusion and exclusion criteria. | Methods, "Study design and data source" and "Cohort construction" (source: TCGA-HNSC via the GDC API; inclusion: primary head-and-neck squamous diagnosis, computable overall survival, primary tumour RNA sequencing profile; exclusions: non-primary or non-squamous diagnoses, second primaries used as the index lesion). Results, "Cohort characteristics"; Table 1. |
| 3 | Describe treatments received and how chosen. | Not available. Treatment information is incomplete in the TCGA clinical file and was not modelled; this is stated in Methods ("Statistical analysis") and listed as a limitation in Discussion. |
| **Specimen characteristics** | | |
| 4 | Describe type of biological material used, and methods of preservation and storage. | Methods, "Expression processing and scoring" (primary tumour tissue, one sample per participant, averaged where multiple aliquots existed). Material type and preservation follow the contributing TCGA tissue source sites; no additional specimen handling was performed by the present authors. |
| **Assay methods** | | |
| 5 | Specify the assay method used and provide a detailed protocol, including reagents, quality control, reproducibility, quantitation, scoring and reporting protocols. Specify blinding. | Methods, "Expression processing and scoring" (GDC STAR-Counts, `tpm_unstranded`, log2(TPM + 1)) and "External validation cohorts" (GSE65858: Illumina HumanHT-12 v4, RSN normalisation with ComBat batch correction by the original submitters; GSE159067: HTG EdgeSeq Oncology Biomarker Panel). All scoring rules are defined analytically in Methods ("Gene sets"; "Expression processing and scoring"), and the analysis code is publicly archived. No laboratory assays were performed, so blinding is not applicable; the transcriptomic HPV classifier was fitted and evaluated without reference to outcome, and was used for stratification and adjustment only. |
| **Study design** | | |
| 6 | State the method of case selection, whether prospective or retrospective, stratification or matching, the time period from which cases were taken, the end of follow-up, and the median follow-up. | Methods, "Study design and data source" (retrospective secondary analysis of publicly available de-identified data; no matching or stratification at selection) and "Cohort construction". Case accrual and follow-up windows are those of the originating studies. Median follow-up among censored patients was 29.4 months (Results; Table 1). |
| 7 | Precisely define all clinical end points examined. | Methods, "Cohort construction" (overall survival from diagnosis to death from any cause, censoring at last follow-up) and "External validation cohorts" (GSE159067: RECIST best response, objective response, progression-free survival, overall survival). |
| 8 | List all candidate variables initially examined or considered for inclusion in models. | Methods, "Statistical analysis" (consensus TLS score, age, sex, stage, subsite, smoking status, HPV status); Table 4; Table 5. |
| 9 | Give rationale for sample size; if designed to detect a specified effect size, give the target power and effect size. | No formal power calculation was performed. The analysis uses all 438 eligible TCGA-HNSC patients and both external cohorts in full; this is a retrospective analysis of existing data sets and the sample size was fixed by the source cohorts. The limited power of subgroup and response analyses is acknowledged in Discussion. |
| **Statistical analysis methods** | | |
| 10 | Specify all statistical methods, including variable selection, model-building, verification of assumptions, and handling of missing data. | Methods, "Statistical analysis" (Kaplan–Meier with log-rank; univariable and multivariable Cox models with a ridge penalty of 0.1 and pre-fitting removal of zero-variance covariates; Schoenfeld residuals for the proportional-hazards assumption; time-dependent AUC; concordance index; logistic regression for response) and "Cut-off determination and stability". Missing data: pathological stage was incomplete and clinical stage was used as a fallback; HPV status was completed by a validated transcriptomic classifier, with the residual uncertainty quantified and discussed. |
| 11 | Clarify how marker values were handled; if relevant, describe methods used for cutpoint determination. | Methods, "Expression processing and scoring" (three algorithms; platform-agnostic z-score consensus) and "Cut-off determination and stability" (log-rank scan over thresholds between the 10th and 90th percentiles with a minimum of 15 patients per group; median cut-off as sensitivity analysis; 1,000 bootstrap resamples with a fixed seed; reproducibility assessed both with and without re-selection). Results, "The cut-off that generates the survival difference is not reproducible" and "Prognostic performance of the consensus TLS score"; Supplementary Table S2. |
| **Results** | | |
| **Data** | | |
| 12 | Describe the flow of patients through the study, including numbers at each stage and reasons for dropout; report numbers of patients and events overall and per subgroup. | Methods, "Cohort construction" (528 participants → 448 with a primary head-and-neck squamous diagnosis, complete stage and computable survival → 438 with primary tumour expression). Numbers and events are reported for the whole cohort and for each stratum in Table 1, Table 5 and Table 6, and in the external validation paragraphs of Methods. |
| 13 | Report distributions of basic demographic characteristics, standard prognostic variables, and the tumour marker, including numbers of missing values. | Table 1 (age, sex, subsite, stage, T and N classification, grade, smoking, deaths, follow-up, HPV status, with missing categories shown explicitly); Results, "Cohort characteristics" (consensus TLS score range, median and interquartile range). |
| **Analysis and presentation** | | |
| 14 | Show the relation of the marker to standard prognostic variables. | Results, "Prognostic performance of the consensus TLS score" (between-group balance tests across sex, N classification, stage, grade, smoking, age, T classification and subsite); Table 4; Table 5. |
| 15 | Present univariate analyses showing the relation between the marker and outcome, with the estimated effect; preferably provide similar analyses for all other variables; a Kaplan–Meier plot is recommended for time-to-event outcomes. | Table 4 (univariable hazard ratios for the consensus score, age, sex and stage) and Figure 1 (Kaplan–Meier estimates by consensus score group). Concordance indices and time-dependent AUC for the marker and for all comparator signatures are given in Results and Figure 3. |
| 16 | For key multivariable analyses, report estimated effects with confidence intervals for the marker and, at least for the final model, all other variables in the model. | Table 4 (multivariable model including the consensus score, age, sex, stage, subsite and smoking status) and Figure 2 (forest plot with 95% confidence intervals). |
| 17 | Provide estimated effects with confidence intervals from an analysis in which the marker and standard prognostic variables are included, regardless of significance. | Table 4 and Figure 2 report the full multivariable model irrespective of statistical significance; the confidence interval for the consensus score spans unity and is reported as such. |
| 18 | If done, report results of further investigations such as checking assumptions, sensitivity analyses, and internal validation. | Methods and Results (proportional-hazards inspection; median versus optimal cut-off; 1,000-resample bootstrap of the cut-off selection and of the fixed cut-off; HPV-stratified and interaction analyses; external validation in two independent cohorts; cross-platform validation of the HPV classifier in Supplementary Figure S1). |
| **Discussion** | | |
| 19 | Interpret the results in the context of the prespecified hypotheses and other relevant studies; include a discussion of limitations of the study. | Discussion, "Comparison with previous work" and "Strengths and limitations". |
| 20 | Discuss implications for future research and clinical value. | Discussion, final paragraph; Conclusion. |

*REMARK, REporting recommendations for tumour MARKer prognostic studies; GDC, Genomic Data Commons; TCGA-HNSC, The Cancer Genome Atlas Head and Neck Squamous Cell Carcinoma; TPM, transcripts per million; RSN, robust spline normalisation; AUC, area under the receiver operating characteristic curve; RECIST, Response Evaluation Criteria in Solid Tumours.*

---

## Supplementary Table S2. Bootstrap distribution of the selected cut-off

In each of 1,000 bootstrap resamples of the discovery cohort (438 patients, 218 deaths; seed 20261016) the entire cut-off selection procedure was repeated across 60 candidate thresholds between the 10th and 90th percentiles of the consensus TLS score, subject to a minimum of 15 patients per group. The table reports the quantiles of the selected threshold, of the maximised log-rank statistic, and of the log-rank p value obtained when the discovery-derived cut-off (0.249) was instead applied unchanged to the resample.

| Quantile | Selected cut-off | Log-rank χ² at the selected cut-off | p at the fixed cut-off (0.249) |
|---|---|---|---|
| Minimum | 0.089 | 1.03 | < 0.001 |
| 10th | 0.209 | 3.46 | 0.0002 |
| 25th | 0.249 | 5.28 | 0.0019 |
| Median | 0.403 | 8.55 | 0.0168 |
| 75th | 0.692 | 12.18 | 0.1106 |
| 90th | 0.727 | 16.07 | 0.3271 |
| Maximum | 0.860 | 34.63 | 0.995 |

Summary statistics: mean 0.444, standard deviation 0.222. The distribution of the selected threshold is bimodal rather than centred on the discovery value: 45.7% of resamples selected a threshold below 0.30 and 40.8% selected one above 0.50, while only 37.6% fell within 0.05 of the discovery value (0.249) and 48.1% fell within the same terzile. A nominally significant log-rank result was obtained by the selection procedure in 87.2% of resamples, which reflects the behaviour of a scan over 60 candidate thresholds and is the reason reproducibility was also assessed without re-selection: holding the cut-off fixed at 0.249, the survival difference remained significant in 63.8% of resamples (median p = 0.0168, interquartile range 0.0019–0.1106); holding it at the sample median, in 14.5%.

*The fixed cut-off is the cut-off that maximised the log-rank statistic in the original cohort (0.249); the sample median cut-off was 0.500. All quantiles are empirical quantiles of the 1,000 resamples.*

---

## Supplementary Figure S1 legend

**Supplementary Figure S1. Cross-platform external validation of the transcriptomic HPV classifier in GSE65858.** The eight genes selected by the nested leave-one-out procedure in TCGA-HNSC and their fitted logistic-regression coefficients were transplanted without any re-fitting into GSE65858 (Illumina HumanHT-12 v4; n = 270 with expression and HPV16 status). Each gene was standardised within the external cohort and the intercept-free linear predictor evaluated. Because the area under the receiver operating characteristic curve is invariant to additive constants, discrimination is directly comparable across the two platforms, whereas sensitivity and specificity at a fixed operating point are not and were re-derived within GSE65858 at Youden thresholds. **(A)** Transcriptionally active HPV16 infection (DNA+RNA+, n = 35) versus HPV16-negative tumours (DNA−, n = 196): area under the curve 0.921 (95% confidence interval 0.872–0.966) for the transplanted classifier, 0.923 for the unweighted mean of the same eight genes, and 0.894 for CDKN2A alone. **(B)** Any HPV16 DNA positivity (DNA+RNA+ and DNA+RNA−, n = 54) versus DNA− (n = 196): discrimination falls to 0.781 (95% confidence interval 0.700–0.856) for the transplanted classifier, and to 0.716 for CDKN2A alone, indicating that the classifier reports on E7-mediated cell-cycle entry (transcriptional activity) rather than on the mere presence of viral DNA. **(C)** Oropharyngeal tumours only (n = 75; 28 positive), the anatomic site in which HPV status is clinically decisive and where the classifier was not trained: area under the curve 0.914 for the transplanted classifier, 0.919 for the unweighted panel mean, and 0.884 for CDKN2A alone. Diagonal dotted lines represent no discrimination. HPV, human papillomavirus; AUC, area under the receiver operating characteristic curve.

---

*Supplementary material status: Supplementary Table S1 (REMARK checklist, 20 items), Supplementary Table S2 (bootstrap cut-off distribution, 1,000 resamples) and Supplementary Figure S1 (three-panel ROC, cross-platform HPV classifier validation) were generated from the same analysis pipeline as the main manuscript and are reproducible from the archived scripts. Supplementary Figure S1 is supplied as a separate file (`FigureS1_fig_supp_hpv_roc.tiff`, 300 dpi).*
