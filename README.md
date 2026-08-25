# Causal Survival Forest Analysis for Individual Treatment Effect Estimation

## Overview

This repository contains the R code used to implement and evaluate a causal survival forest (CSF) model for estimating individual treatment effects (ITEs) in a time-to-event setting.

The analysis framework was designed to characterize treatment effect heterogeneity and identify patient subgroups with potentially differential treatment effects. The repository includes code for multiple imputation, CSF model construction, ITE estimation and visualization, model calibration and discrimination assessment, effect modification analyses, subgroup analyses, sensitivity analyses, and SHAP-based model interpretation.

## Analysis Workflow

The overall analytical workflow includes:

1. Multiple imputation of missing baseline covariates.
2. Construction of the causal survival forest model.
3. Estimation of individual treatment effects (ITEs).
4. Evaluation of model discrimination using the targeting operator characteristic (TOC) curve.
5. Evaluation of model calibration.
6. Assessment of the deviation of conditional average treatment effects (CATEs) from the average treatment effect (ATE).
7. Exploration of associations between baseline characteristics and estimated ITEs.
8. Identification of treatment-effect heterogeneity using effect modification and subgroup analyses.
9. Two-step subgroup discovery to identify clinically relevant patient subgroups.
10. Sensitivity analyses to assess the robustness of the findings.
11. SHAP-based interpretation of the CSF model.

## Code Description

| File                                 | Description                                                                                                                                          |
| ------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------- |
| `code_imputation.R`                  | Performs multiple imputation of missing data.                                                                                                        |
| `code_csf_model_construction.R`      | Fits the causal survival forest (CSF) model and estimates individual treatment effects (ITEs).                                                       |
| `code_ite_plots.R`                   | Visualizes the distribution of estimated ITEs.                                                                                                       |
| `code_categorical_distributions.R`   | Examines the distribution of ITEs across categorical baseline variables.                                                                             |
| `code_pdp_plots.R`                   | Evaluates and visualizes the relationships between continuous baseline variables and estimated ITEs using partial dependence plots.                  |
| `code_calibration_plot.R`            | Generates calibration plots to assess the calibration of the CSF model.                                                                              |
| `code_toc_qini.R`                    | Generates targeting operator characteristic (TOC) curves to evaluate the discriminative ability of the CSF model for treatment-effect heterogeneity. |
| `code_cate_minus_ate_bootstrap_ci.R` | Estimates confidence intervals for the deviation of CATEs from the overall ATE using bootstrap procedures.                                           |
| `code_effect_modifier_plot.R`        | Performs effect modification analyses and visualizes treatment-effect heterogeneity.                                                                 |
| `code_clinical_rule_analysis.R`      | Performs two-step subgroup discovery to identify clinically relevant subgroups with differential treatment effects.                                  |
| `code_subgroup_analysis.R`           | Performs prespecified subgroup analyses to evaluate heterogeneity of treatment effects across clinically relevant subgroups.                         |
| `code_tertile_asmd.R`                | Calculates absolute standardized mean differences (ASMDs) before and after weighting across ITE tertiles to assess covariate balance.                |
| `code_sensitivity_analyses.R`        | Performs sensitivity analyses to assess the robustness of the estimated treatment effects.                                                           |
| `code_shap_analysis.R`               | Performs SHAP analysis to characterize global and local contributions of baseline variables to the estimated treatment effects.                      |

## Software

The analyses were performed using R.

The specific R packages and package versions required to reproduce the analyses are specified within the corresponding R scripts.

## Data Availability

The patient-level data used in this study are not included in this repository because of data access, privacy, and/or institutional restrictions.

The analysis code is provided to facilitate transparency and reproducibility of the reported analyses.

Researchers interested in reproducing the analyses should obtain access to the underlying data through the appropriate data providers or according to the applicable institutional data-access procedures.

## Reproducibility

The R scripts are organized according to the major analytical components of the study. Because the underlying patient-level datasets are not distributed with this repository, the scripts require access to the corresponding analysis dataset and appropriate preprocessing before execution.

The results reported in the manuscript were generated using the analytical procedures implemented in these scripts.

## Citation

If you use this code or analytical framework, please cite the corresponding study.

## Contact

For questions regarding the analysis code or study methodology, please refer to the corresponding author of the associated publication.
