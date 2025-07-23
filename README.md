# The Effect of Smokers Transitioning to E-cigarettes on Physical and Mental Health: An Emulated Trial using Longitudinal Data
Project Overview
This repository contains the R scripts and analysis for a dissertation project that emulates a target trial using longitudinal data from the UK Household Longitudinal Study (UKHLS). The primary objective is to estimate the causal effect on physical and mental health when smokers transition to e-cigarette use, dual-use, or complete cessation, compared to continuing to smoke.

The analysis uses several statistical methods to control for confounding and selection bias, including:

Propensity Score Matching (PSM) to create comparable groups.

Doubly Robust Regression on the matched cohort for the main analysis.

Fixed-Effects Models and Subgroup Analyses as sensitivity checks.

The primary outcomes are the SF-12 Physical Component Summary (PCS), the SF-12 Mental Component Summary (MCS), and the EQ-5D health utility index.


Scripts
The analysis is performed by running the following R scripts in sequential order:

01_data_preparation.R: Loads the raw Stata data, performs all cleaning and variable recoding, and saves a final complete-case dataset.

02_descriptive_unmatched_analysis.R: Creates the "paired" cohort for identifying smoking transitions and runs descriptive and regression analyses on the full unmatched data.

03_propensity_score_matching.R: Performs 1:3 propensity score matching, assesses covariate balance with Love plots, and saves the final matched datasets.

04_matched_regression_analysis.R: Runs non-doubly robust and doubly robust regression models on the matched data and performs detailed model diagnostics.

05_sensitivity_subgroup_analysis.R: Conducts sensitivity analysis using Fixed-Effects models and performs subgroup analyses by age, education, and income.

06_final_plots_tables.R: Generates comprehensive summary tables and forest plots comparing all analytical approaches.


