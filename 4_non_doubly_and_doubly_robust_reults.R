################################################################################
# Dissertation Project: The Effect of Smokers transitioning to E-cigarette in
#                       Physical and Mental health : An Emulated Trial using Longitudinal Data
#
# SCRIPT 4: REGRESSION ANALYSIS ON MATCHED DATA
#
# What this script does:
# 1. Loads the three 1:3 matched datasets.
# 2. Maps EQ-5D scores for the matched data.
# 3. Runs non-doubly robust and doubly robust models.
# 4. Performs and saves detailed model diagnostics for the main models.
# 5. Generates and saves a final summary table as a PNG image.
#
################################################################################


#===============================================================================
# 1. SETUP
#===============================================================================
library(tidyverse)
library(estimatr) # For robust standard errors
library(broom)    # For tidying model output
library(gt)       # For publication-quality tables
library(car)      # For VIF (diagnostics)
library(lmtest)   # For Breusch-Pagan test (diagnostics)

# --- Define paths ---
MATCH_DATA_DIR <- "./output_data/matched/"
OUTPUT_DIR <- "./output_results/"
DIAGNOSTICS_DIR <- "./output_results/diagnostics/"

dir.create(OUTPUT_DIR, showWarnings = FALSE)
dir.create(DIAGNOSTICS_DIR, showWarnings = FALSE)


#===============================================================================
# 2. LOAD AND PREPARE MATCHED DATA
#===============================================================================

# --- Load the 1:3 matched datasets from Script 3 ---
matched_data_switcher <- readRDS(paste0(MATCH_DATA_DIR, "matched_data_switcher_1to3.rds"))
matched_data_dual     <- readRDS(paste0(MATCH_DATA_DIR, "matched_data_dual_1to3.rds"))
matched_data_quitter  <- readRDS(paste0(MATCH_DATA_DIR, "matched_data_quitter_1to3.rds"))

# --- Function to map SF-12 to EQ-5D ---
map_eq5d <- function(pcs, mcs) {
  pmin(-1.6984 + (pcs * 0.07927) + (mcs * 0.02859) + (pcs * mcs * -0.000126) +
         (pcs^2 * -0.00141) + (mcs^2 * -0.00014) + (pcs^3 * 0.0000107), 1)
}

# --- Apply function to all three datasets ---
matched_data_switcher <- matched_data_switcher %>% mutate(eq5d_t0 = map_eq5d(sf12pcs_dv_t0, sf12mcs_dv_t0), eq5d_t1 = map_eq5d(sf12pcs_dv_t1, sf12mcs_dv_t1))
matched_data_dual     <- matched_data_dual %>% mutate(eq5d_t0 = map_eq5d(sf12pcs_dv_t0, sf12mcs_dv_t0), eq5d_t1 = map_eq5d(sf12pcs_dv_t1, sf12mcs_dv_t1))
matched_data_quitter  <- matched_data_quitter %>% mutate(eq5d_t0 = map_eq5d(sf12pcs_dv_t0, sf12mcs_dv_t0), eq5d_t1 = map_eq5d(sf12pcs_dv_t1, sf12mcs_dv_t1))


#===============================================================================
# 3. REGRESSION MODEL ESTIMATION
#===============================================================================

# --- Define Model Equations ---
pcs_formula <- as.formula(sf12pcs_dv_t1 ~ smoking_group + sf12pcs_dv_t0 + sf12mcs_dv_t0 + age_dv_t0 + age_squared_t0 + sex_dv_t0 + ethn_dv_t0 + hiqual_dv_t0 + health_t0 + jbstat_t0 + log_real_income_t0 + hl2gp_t0 + gor_dv_t0 + nkids015_t0 + factor(wave_t0))
mcs_formula <- as.formula(sf12mcs_dv_t1 ~ smoking_group + sf12pcs_dv_t0 + sf12mcs_dv_t0 + age_dv_t0 + age_squared_t0 + sex_dv_t0 + ethn_dv_t0 + hiqual_dv_t0 + health_t0 + jbstat_t0 + log_real_income_t0 + hl2gp_t0 + gor_dv_t0 + nkids015_t0 + factor(wave_t0))
eq5d_formula <- as.formula(eq5d_t1 ~ smoking_group + eq5d_t0 + age_dv_t0 + age_squared_t0 + sex_dv_t0 + ethn_dv_t0 + hiqual_dv_t0 + health_t0 + jbstat_t0 + log_real_income_t0 + hl2gp_t0 + gor_dv_t0 + nkids015_t0 + factor(wave_t0))

# --- Non-Doubly Robust Models (Unadjusted on Matched Data) ---
ndr_switcher_pcs  <- tidy(lm_robust(sf12pcs_dv_t1 ~ smoking_group, data = matched_data_switcher, clusters = pidp))
ndr_switcher_mcs  <- tidy(lm_robust(sf12mcs_dv_t1 ~ smoking_group, data = matched_data_switcher, clusters = pidp))
ndr_switcher_eq5d <- tidy(lm_robust(eq5d_t1 ~ smoking_group,         data = matched_data_switcher, clusters = pidp))

ndr_dual_pcs      <- tidy(lm_robust(sf12pcs_dv_t1 ~ smoking_group, data = matched_data_dual, clusters = pidp))
ndr_dual_mcs      <- tidy(lm_robust(sf12mcs_dv_t1 ~ smoking_group, data = matched_data_dual, clusters = pidp))
ndr_dual_eq5d     <- tidy(lm_robust(eq5d_t1 ~ smoking_group,         data = matched_data_dual, clusters = pidp))

ndr_quitter_pcs   <- tidy(lm_robust(sf12pcs_dv_t1 ~ smoking_group, data = matched_data_quitter, clusters = pidp))
ndr_quitter_mcs   <- tidy(lm_robust(sf12mcs_dv_t1 ~ smoking_group, data = matched_data_quitter, clusters = pidp))
ndr_quitter_eq5d  <- tidy(lm_robust(eq5d_t1 ~ smoking_group,         data = matched_data_quitter, clusters = pidp))


# --- Doubly Robust Models (Adjusted on Matched Data) ---
dr_switcher_pcs  <- tidy(lm_robust(pcs_formula, data = matched_data_switcher, clusters = pidp))
dr_switcher_mcs  <- tidy(lm_robust(mcs_formula, data = matched_data_switcher, clusters = pidp))
dr_switcher_eq5d <- tidy(lm_robust(eq5d_formula, data = matched_data_switcher, clusters = pidp))

dr_dual_pcs      <- tidy(lm_robust(pcs_formula, data = matched_data_dual, clusters = pidp))
dr_dual_mcs      <- tidy(lm_robust(mcs_formula, data = matched_data_dual, clusters = pidp))
dr_dual_eq5d     <- tidy(lm_robust(eq5d_formula, data = matched_data_dual, clusters = pidp))

dr_quitter_pcs   <- tidy(lm_robust(pcs_formula, data = matched_data_quitter, clusters = pidp))
dr_quitter_mcs   <- tidy(lm_robust(mcs_formula, data = matched_data_quitter, clusters = pidp))
dr_quitter_eq5d  <- tidy(lm_robust(eq5d_formula, data = matched_data_quitter, clusters = pidp))


#===============================================================================
# 4. MODEL DIAGNOSTICS FOR DOUBLY ROBUST MODELS
#===============================================================================

# --- Re-run DR models with standard lm() to access diagnostic functions ---
  DR_Switcher_PCS = lm(pcs_formula, data = matched_data_switcher)
  DR_Switcher_MCS = lm(mcs_formula, data = matched_data_switcher)
  DR_Switcher_eq5d = lm(eq5d_formula, data = matched_data_switcher)
  DR_Dual_PCS = lm(pcs_formula, data = matched_data_dual)
  DR_Dual_MCS = lm(mcs_formula, data = matched_data_dual)
  DR_Dual_eq5d = lm(eq5d_formula, data = matched_data_dual)
  DR_Quitter_PCS = lm(pcs_formula, data = matched_data_quitter)
  DR_Quitter_MCS = lm(mcs_formula, data = matched_data_quitter)
  DR_Quitter_eq5d = lm(eq5d_formula, data = matched_data_quitter)


# --- Function to perform and save diagnostics for a given model ---
# --- Plots for SF-12 PCS Models ---
png("Diagnostic_Plots_PCS_vs_Switcher.png", width = 1000, height = 1000, res = 120)
par(mfrow = c(2, 2))
plot(DR_Switcher_PCS, main = "Diagnostics: PCS vs. Switcher")
dev.off()

png("Diagnostic_Plots_PCS_vs_DualUser.png", width = 1000, height = 1000, res = 120)
par(mfrow = c(2, 2))
plot(DR_Dual_PCS, main = "Diagnostics: PCS vs. Dual User")
dev.off()

png("Diagnostic_Plots_PCS_vs_Quitter.png", width = 1000, height = 1000, res = 120)
par(mfrow = c(2, 2))
plot(DR_Quitter_PCS, main = "Diagnostics: PCS vs. Quitter")
dev.off()


# --- Plots for SF-12 MCS Models ---
png("Diagnostic_Plots_MCS_vs_Switcher.png", width = 1000, height = 1000, res = 120)
par(mfrow = c(2, 2))
plot(DR_Switcher_MCS, main = "Diagnostics: MCS vs. Switcher")
dev.off()

png("Diagnostic_Plots_MCS_vs_DualUser.png", width = 1000, height = 1000, res = 120)
par(mfrow = c(2, 2))
plot(DR_Dual_MCS, main = "Diagnostics: MCS vs. Dual User")
dev.off()

png("Diagnostic_Plots_MCS_vs_Quitter.png", width = 1000, height = 1000, res = 120)
par(mfrow = c(2, 2))
plot(DR_Quitter_MCS, main = "Diagnostics: MCS vs. Quitter")
dev.off()


# --- Plots for EQ-5D Models ---
png("Diagnostic_Plots_EQ5D_vs_Switcher.png", width = 1000, height = 1000, res = 120)
par(mfrow = c(2, 2))
plot(DR_Switcher_eq5d, main = "Diagnostics: EQ-5D vs. Switcher")
dev.off()

png("Diagnostic_Plots_EQ5D_vs_DualUser.png", width = 1000, height = 1000, res = 120)
par(mfrow = c(2, 2))
plot(DR_Dual_eq5d, main = "Diagnostics: EQ-5D vs. Dual User")
dev.off()

png("Diagnostic_Plots_EQ5D_vs_Quitter.png", width = 1000, height = 1000, res = 120)
par(mfrow = c(2, 2))
plot(DR_Quitter_eq5d, main = "Diagnostics: EQ-5D vs. Quitter")
dev.off()

#===============================================================================
# 5. GENERATE AND SAVE FINAL RESULTS TABLE
#===============================================================================

# --- Helper function for formatting result strings ---
format_results <- function(estimate, conf.low, conf.high, p.value) {
  stars <- case_when(p.value < 0.001 ~ "***", p.value < 0.01  ~ "**", p.value < 0.05  ~ "*", TRUE ~ "")
  p_formatted <- case_when(p.value < 0.001 ~ "*p*<0.001", TRUE ~ paste0("*p*=", sprintf("%.3f", p.value)))
  paste0(sprintf("%.3f", estimate), stars, " [", sprintf("%.3f", conf.low), ", ", sprintf("%.3f", conf.high), "], ", p_formatted)
}

# --- Combine all results into a single dataframe ---
matched_results_data <- bind_rows(
  "Physical Health (SF-12 PCS)" = bind_rows(
    "Unadjusted" = bind_rows("Switcher" = ndr_switcher_pcs, "Dual User" = ndr_dual_pcs, "Quitter" = ndr_quitter_pcs, .id="Comparison"),
    "Fully Adjusted (Doubly Robust)" = bind_rows("Switcher" = dr_switcher_pcs, "Dual User" = dr_dual_pcs, "Quitter" = dr_quitter_pcs, .id="Comparison"), .id="Model"),
  "Mental Health (SF-12 MCS)" = bind_rows(
    "Unadjusted" = bind_rows("Switcher" = ndr_switcher_mcs, "Dual User" = ndr_dual_mcs, "Quitter" = ndr_quitter_mcs, .id="Comparison"),
    "Fully Adjusted (Doubly Robust)" = bind_rows("Switcher" = dr_switcher_mcs, "Dual User" = dr_dual_mcs, "Quitter" = dr_quitter_mcs, .id="Comparison"), .id="Model"),
  "Health-Related Quality of Life (EQ-5D)" = bind_rows(
    "Unadjusted" = bind_rows("Switcher" = ndr_switcher_eq5d, "Dual User" = ndr_dual_eq5d, "Quitter" = ndr_quitter_eq5d, .id="Comparison"),
    "Fully Adjusted (Doubly Robust)" = bind_rows("Switcher" = dr_switcher_eq5d, "Dual User" = dr_dual_eq5d, "Quitter" = dr_quitter_eq5d, .id="Comparison"), .id="Model"),
  .id = "Outcome"
) %>%
  filter(str_detect(term, "smoking_group")) %>%
  mutate(Result = format_results(estimate, conf.low, conf.high, p.value)) %>%
  select(Outcome, Comparison, Model, Result) %>%
  pivot_wider(names_from = Model, values_from = Result)

# --- Create and save the final table ---
matched_results_gt <- gt(matched_results_data, groupname_col = "Outcome", rowname_col = "Comparison") %>%
  tab_header(
    title = md("**Associations from Propensity Score Matched Cohort (1:3 Matching)**"),
    subtitle = "Results from Models on the Matched Cohort"
  ) %>%
  cols_label(
    `Unadjusted` = md("**Unadjusted Model**<br>Est. [95% CI], *p*-value"),
    `Fully Adjusted (Doubly Robust)` = md("**Doubly Robust Model**<br>Est. [95% CI], *p*-value")
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups())%>%
  fmt_markdown(columns = everything()) %>%
  tab_options(table.width = pct(100)) %>%
  tab_footnote(footnote = md("*Note: Reference group is 'Continued Smoker'*.<br>**p*<0.05, ***p*<0.01, ****p*<0.001"))

print(matched_results_gt)
gtsave(matched_results_gt, filename = paste0(OUTPUT_DIR, "Table3_Matched_Results.png"))