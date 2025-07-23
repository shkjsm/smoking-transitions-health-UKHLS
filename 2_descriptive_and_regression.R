################################################################################
# Dissertation Project: The Effect of Smokers transitioning to E-cigarette in
#                       Physical and Mental health : An Emulated Trial using Longitudinal Data
#
# SCRIPT 2: DESCRIPTIVE AND UNMATCHED ANALYSIS
#
# What this script does:
# 1. Loads the clean, complete case data.
# 2. Creates a "paired" dataset to identify smoking transitions.
# 3. Generates a table of baseline characteristics.
# 4. Runs unadjusted and fully adjusted OLS regression models.
# 5. Saves the final table as a PNG image and docx.
#
################################################################################


#===============================================================================
# 1. SETUP
#===============================================================================
library(tidyverse)
library(tableone) # For descriptive statistics table
library(estimatr) # For robust standard errors
library(broom)    # For tidying model output
library(gt)       # For publication-quality tables

# --- Define paths ---
INPUT_DIR <- "./output_data/"
OUTPUT_DIR <- "./output_results/"
dir.create(OUTPUT_DIR, showWarnings = FALSE)


#===============================================================================
# 2. LOAD DATA AND CREATE PAIRED COHORT
#===============================================================================

# --- Load the clean dataset from Script 1 ---
ukhls_complete <- readRDS(paste0(INPUT_DIR, "ukhls_complete_cases.rds"))

# --- Function to create paired data between two waves ---
create_paired_data <- function(data, t0_wave, t1_wave) {
  # Baseline data (t0): Smokers not using e-cigarettes
  df_t0 <- data %>%
    filter(wave == t0_wave, age_dv >= 16, smoker == "Smoker", ecigs_use == "Non-user") %>%
    rename_with(~ paste0(., "_t0"), .cols = -pidp)
  
  # Follow-up data (t1)
  df_t1 <- data %>%
    filter(wave == t1_wave) %>%
    rename_with(~ paste0(., "_t1"), .cols = -pidp)
  
  # Join to find individuals present in both waves
  paired_df <- inner_join(df_t0, df_t1, by = "pidp")
  
  return(paired_df)
}

# --- Create the full paired dataset across all waves (7-14) ---
wave_pairs <- 7:13
ukhls_paired <- map_dfr(wave_pairs, ~ create_paired_data(ukhls_complete, .x, .x + 1))

# --- Define transition groups and map EQ-5D scores ---
# Function to map SF-12 to EQ-5D
map_eq5d <- function(pcs, mcs) {
  pmin(-1.6984 + (pcs * 0.07927) + (mcs * 0.02859) + (pcs * mcs * -0.000126) +
         (pcs^2 * -0.00141) + (mcs^2 * -0.00014) + (pcs^3 * 0.0000107), 1)
}

ukhls_paired <- ukhls_paired %>%
  mutate(
    smoking_group = case_when(
      smoker_t0 == "Smoker" & smoker_t1 == "Smoker" & ecigs_use_t1 == "Non-user"     ~ "Continued Smoker",
      smoker_t0 == "Smoker" & smoker_t1 == "Non-smoker" & ecigs_use_t1 == "Current user" ~ "Switcher",
      smoker_t0 == "Smoker" & smoker_t1 == "Smoker" & ecigs_use_t1 == "Current user"     ~ "Dual User",
      smoker_t0 == "Smoker" & smoker_t1 == "Non-smoker" & ecigs_use_t1 == "Non-user"     ~ "Quitter"
    ),
    wave_pair = paste0(wave_t0, "-", wave_t1),
    # Map EQ-5D scores for baseline and follow-up
    eq5d_t0 = map_eq5d(sf12pcs_dv_t0, sf12mcs_dv_t0),
    eq5d_t1 = map_eq5d(sf12pcs_dv_t1, sf12mcs_dv_t1)
  ) %>%
  filter(!is.na(smoking_group)) %>%
  mutate(smoking_group = factor(smoking_group, levels = c("Continued Smoker", "Switcher", "Dual User", "Quitter")))


#===============================================================================
# 3. DESCRIPTIVE ANALYSIS (TABLE)
#===============================================================================

# --- Define variable lists for Table 1 ---
cat_vars <- c("sex_dv_t0", "ethn_dv_t0", "hiqual_dv_t0", "jbstat_t0", "health_t0", "hl2gp_t0", "gor_dv_t0")
cont_vars <- c("age_dv_t0", "log_real_income_t0", "sf12pcs_dv_t0", "sf12mcs_dv_t0", "nkids015_t0")

# --- Create Table 1 ---
table1_obj <- CreateTableOne(
  vars = c(cont_vars, cat_vars),
  strata = "smoking_group",
  data = ukhls_paired,
  test = TRUE
)

# Convert to a data frame for use with gt
table1_df <- as.data.frame(print(table1_obj, printToggle = FALSE, smd = TRUE)) %>%
  rownames_to_column(var = "Characteristic")


#===============================================================================
# 4. UNMATCHED REGRESSION ANALYSIS
#===============================================================================

# --- Define model formulas ---
formula_pcs <- as.formula(sf12pcs_dv_t1 ~ smoking_group + sf12pcs_dv_t0 + sf12mcs_dv_t0 + age_dv_t0 + age_squared_t0 +
                                 sex_dv_t0 + ethn_dv_t0 + hiqual_dv_t0 + health_t0 + jbstat_t0 + log_real_income_t0 +
                                 hl2gp_t0 + gor_dv_t0 + nkids015_t0 + factor(wave_t0))

formula_mcs <- as.formula(sf12mcs_dv_t1 ~ smoking_group + sf12pcs_dv_t0 + sf12mcs_dv_t0 + age_dv_t0 + age_squared_t0 +
                                 sex_dv_t0 + ethn_dv_t0 + hiqual_dv_t0 + health_t0 + jbstat_t0 + log_real_income_t0 +
                                 hl2gp_t0 + gor_dv_t0 + nkids015_t0 + factor(wave_t0))

formula_eq5d <- as.formula(eq5d_t1 ~ smoking_group + eq5d_t0 + age_dv_t0 + age_squared_t0 + sex_dv_t0 + ethn_dv_t0 +
                                  hiqual_dv_t0 + health_t0 + jbstat_t0 + log_real_income_t0 + hl2gp_t0 + gor_dv_t0 +
                                  nkids015_t0 + factor(wave_t0))

# --- Run Unadjusted Models ---
lm_unadj_pcs <- tidy(lm_robust(sf12pcs_dv_t1 ~ smoking_group, data = ukhls_paired, clusters = pidp))
lm_unadj_mcs <- tidy(lm_robust(sf12mcs_dv_t1 ~ smoking_group, data = ukhls_paired, clusters = pidp))
lm_unadj_eq5d <- tidy(lm_robust(eq5d_t1 ~ smoking_group, data = ukhls_paired, clusters = pidp))

# --- Run Fully Adjusted Models ---
lm_adj_pcs <- tidy(lm_robust(formula_pcs, data = ukhls_paired, clusters = pidp))
lm_adj_mcs <- tidy(lm_robust(formula_mcs, data = ukhls_paired, clusters = pidp))
lm_adj_eq5d <- tidy(lm_robust(formula_eq5d, data = ukhls_paired, clusters = pidp))


#===============================================================================
# 5. GENERATE AND SAVE FINAL TABLES
#===============================================================================

# --- Save Table 1: Baseline Characteristics ---
table1_gt <- gt(table1_df) %>%
  tab_header(
    title = md("**Table 1: Baseline Characteristics of the Unmatched Cohort**"),
    subtitle = "Stratified by Smoking Transition Group"
  ) %>%
  tab_options(table.width = pct(100), container.overflow.x = TRUE)

gtsave(table1_gt, filename = paste0(OUTPUT_DIR, "Table1_Baseline_Characteristics.png"))

# --- Helper function for formatting regression result  ---
format_results <- function(estimate, conf.low, conf.high, p.value) {
  stars <- case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE ~ ""
  )
  p_formatted <- case_when(
    p.value < 0.001 ~ "*p*<0.001", TRUE ~ paste0("*p*=", sprintf("%.3f", p.value))
  )
  paste0(
    sprintf("%.3f", estimate), stars, " [",
    sprintf("%.3f", conf.low), ", ",
    sprintf("%.3f", conf.high), "], ",
    p_formatted
  )
}

# --- Prepare data for the regression results table ---
unmatched_results_data <- bind_rows(
  lm_unadj_pcs %>% mutate(Model = "Unadjusted", Outcome = "Physical Health (SF-12 PCS)"),
  lm_adj_pcs %>% mutate(Model = "Fully Adjusted", Outcome = "Physical Health (SF-12 PCS)"),
  lm_unadj_mcs %>% mutate(Model = "Unadjusted", Outcome = "Mental Health (SF-12 MCS)"),
  lm_adj_mcs %>% mutate(Model = "Fully Adjusted", Outcome = "Mental Health (SF-12 MCS)"),
  lm_unadj_eq5d %>% mutate(Model = "Unadjusted", Outcome = "Health-Related Quality of Life (EQ-5D-3L)"),
  lm_adj_eq5d %>% mutate(Model = "Fully Adjusted", Outcome = "Health-Related Quality of Life (EQ-5D-3L)")
) %>%
  filter(str_detect(term, "smoking_group")) %>%
  mutate(
    Comparison = str_remove(term, "smoking_group"),
    Result = format_results(estimate, conf.low, conf.high, p.value)
  ) %>%
  select(Outcome, Comparison, Model, Result) %>%
  pivot_wider(names_from = Model, values_from = Result)

# --- Create and save the regression results table ---
unmatched_results_gt <- gt(unmatched_results_data, groupname_col = "Outcome", rowname_col = "Comparison") %>%
  tab_header(
    title = md("**Associations Between Smoking Transitions and Health Outcomes (Unmatched Cohort)**"),
    subtitle = "Results from Models on the Unmatched Cohort"
  ) %>%
  cols_label(
    `Unadjusted` = md("**Unadjusted Model**<br>Est. [95% CI], *p*-value"),
    `Fully Adjusted` = md("**Fully Adjusted Model**<br>Est. [95% CI], *p*-value")
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups())%>%
  fmt_markdown(columns = everything()) %>%
  tab_options(table.width = pct(100)) %>%
  tab_footnote(footnote = md("Reference group is 'Continued Smoker'. <br>**p*<0.05, ***p*<0.01, ****p*<0.001"))

print(unmatched_results_gt)
gtsave(unmatched_results_gt, filename = paste0(OUTPUT_DIR, "Table2_Unmatched_Results.png"))
