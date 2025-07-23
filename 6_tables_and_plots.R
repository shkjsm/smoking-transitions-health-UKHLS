################################################################################
# Dissertation Project: The Effect of Smokers transitioning to E-cigarette in
#                       Physical and Mental health : An Emulated Trial using Longitudinal Data
#
# SCRIPT 6: COMPREHENSIVE TABLE AND FOREST PLOT GENERATION
#
# What this script does:
# 1. Loads all necessary clean and matched data.
# 2. Re-runs all primary models (unmatched, matched, FE) to have results
#    available in the environment.
# 3. Creates and saves a comprehensive summary table comparing all models.
# 4. Creates and saves comprehensive forest plots visualizing all models.
#
################################################################################


#===============================================================================
# 1. SETUP
#===============================================================================
library(tidyverse)
library(estimatr)
library(broom)
library(plm)
library(lmtest)
library(gt)

# --- Define Paths ---
DATA_DIR <- "./output_data/"
MATCH_DATA_DIR <- "./output_data/matched/"
OUTPUT_DIR <- "./output_results/"
dir.create(OUTPUT_DIR, showWarnings = FALSE)


#===============================================================================
# 2. LOAD DATA & RE-RUN MODELS
# This section re-runs the main models to make their results available for
# plotting and creating the final comprehensive table.
#===============================================================================

# --- Load Data ---
ukhls_complete <- readRDS(paste0(DATA_DIR, "ukhls_complete_cases.rds"))
matched_switcher <- readRDS(paste0(MATCH_DATA_DIR, "matched_data_switcher_1to3.rds"))
matched_dual     <- readRDS(paste0(MATCH_DATA_DIR, "matched_data_dual_1to3.rds"))
matched_quitter  <- readRDS(paste0(MATCH_DATA_DIR, "matched_data_quitter_1to3.rds"))

# --- Helper functions to create paired data and map EQ5D ---
create_paired_data <- function(data, t0_wave, t1_wave) {
  df_t0 <- data %>% filter(wave == t0_wave, smoker == "Smoker", ecigs_use == "Non-user") %>% rename_with(~ paste0(., "_t0"), -pidp)
  df_t1 <- data %>% filter(wave == t1_wave) %>% rename_with(~ paste0(., "_t1"), -pidp)
  inner_join(df_t0, df_t1, by = "pidp") %>%
    mutate(smoking_group = case_when(
      smoker_t1 == "Smoker" & ecigs_use_t1 == "Non-user" ~ "Continued Smoker",
      smoker_t1 == "Non-smoker" & ecigs_use_t1 == "Current user" ~ "Switcher",
      smoker_t1 == "Smoker" & ecigs_use_t1 == "Current user" ~ "Dual User",
      smoker_t1 == "Non-smoker" & ecigs_use_t1 == "Non-user" ~ "Quitter")) %>%
    filter(!is.na(smoking_group)) %>%
    mutate(smoking_group = factor(smoking_group, levels=c("Continued Smoker", "Switcher", "Dual User", "Quitter")))
}
map_eq5d <- function(pcs, mcs) { pmin(-1.6984 + (pcs*0.07927) + (mcs*0.02859) + (pcs*mcs*-0.000126) + (pcs^2*-0.00141) + (mcs^2*-0.00014) + (pcs^3*0.0000107), 1) }

# --- Create unmatched paired data ---
unmatched_paired <- map_dfr(7:13, ~create_paired_data(ukhls_complete, .x, .x+1)) %>%
  mutate(eq5d_t0 = map_eq5d(sf12pcs_dv_t0, sf12mcs_dv_t0), eq5d_t1 = map_eq5d(sf12pcs_dv_t1, sf12mcs_dv_t1))

# --- Add EQ5D scores to matched data ---
all_matched <- list(Switcher=matched_switcher, `Dual User`=matched_dual, Quitter=matched_quitter) %>%
  map(~ .x %>% mutate(eq5d_t0 = map_eq5d(sf12pcs_dv_t0, sf12mcs_dv_t0), eq5d_t1 = map_eq5d(sf12pcs_dv_t1, sf12mcs_dv_t1)))

# --- Define all model formulas ---
f_adj_pcs <- as.formula(sf12pcs_dv_t1 ~ smoking_group + sf12pcs_dv_t0 + sf12mcs_dv_t0 + age_dv_t0 + age_squared_t0 + sex_dv_t0 + ethn_dv_t0 + hiqual_dv_t0 + health_t0 + jbstat_t0 + log_real_income_t0 + hl2gp_t0 + gor_dv_t0 + nkids015_t0 + factor(wave_t0))
f_adj_mcs <- as.formula(sf12mcs_dv_t1 ~ smoking_group + sf12pcs_dv_t0 + sf12mcs_dv_t0 + age_dv_t0 + age_squared_t0 + sex_dv_t0 + ethn_dv_t0 + hiqual_dv_t0 + health_t0 + jbstat_t0 + log_real_income_t0 + hl2gp_t0 + gor_dv_t0 + nkids015_t0 + factor(wave_t0))
f_adj_eq5d <- as.formula(eq5d_t1 ~ smoking_group + eq5d_t0 + age_dv_t0 + age_squared_t0 + sex_dv_t0 + ethn_dv_t0 + hiqual_dv_t0 + health_t0 + jbstat_t0 + log_real_income_t0 + hl2gp_t0 + gor_dv_t0 + nkids015_t0 + factor(wave_t0))
f_fe_pcs <- as.formula(sf12pcs_dv_t1 ~ smoking_group + sf12pcs_dv_t0 + age_dv_t0 + health_t0 + jbstat_t0 + log_real_income_t0 + hl2gp_t0 + nkids015_t0)
f_fe_mcs <- as.formula(sf12mcs_dv_t1 ~ smoking_group + sf12mcs_dv_t0 + age_dv_t0 + health_t0 + jbstat_t0 + log_real_income_t0 + hl2gp_t0 + nkids015_t0)
f_fe_eq5d <- as.formula(eq5d_t1 ~ smoking_group + eq5d_t0 + age_dv_t0 + health_t0 + jbstat_t0 + log_real_income_t0 + hl2gp_t0 + nkids015_t0)

# --- Function to run a set of models and return a tidy data frame ---
run_all_models_for_outcome <- function(outcome_name, unmatched_data, matched_data_list, f_unadj, f_adj, f_fe) {
  # Unmatched Models
  unmatched_unadj <- tidy(lm_robust(f_unadj, data = unmatched_data, clusters = pidp))
  unmatched_adj <- tidy(lm_robust(f_adj, data = unmatched_data, clusters = pidp))
  
  # Matched Models
  matched_results <- imap_dfr(matched_data_list, function(data, name) {
    ndr_model <- tidy(lm_robust(f_unadj, data = data, clusters = pidp))
    dr_model <- tidy(lm_robust(f_adj, data = data, clusters = pidp))
    pdata <- pdata.frame(data, index = c("pidp", "wave_t1"))
    fe_model <- coeftest(plm(f_fe, data=pdata, model="within"), vcov=vcovHC(plm(f_fe, data=pdata, model="within"), type="HC1", cluster="group")) %>% tidy(conf.int=TRUE)
    bind_rows(
      "Unadjusted (Matched)" = ndr_model,
      "Fully Adjusted (DR)" = dr_model,
      "Fixed Effects" = fe_model,
      .id = "model_type"
    ) %>% mutate(comparison = name)
  })
   
  bind_rows(
    unmatched_unadj %>% mutate(model_type = "Unadjusted (Unmatched)", comparison = "All"),
    unmatched_adj %>% mutate(model_type = "Fully Adjusted (Unmatched)", comparison = "All"),
    matched_results
  ) %>% mutate(outcome = outcome_name)
}

# --- Execute for all three outcomes ---
pcs_results <- run_all_models_for_outcome("Physical Health (SF-12 PCS)", unmatched_paired, all_matched, as.formula(sf12pcs_dv_t1~smoking_group), f_adj_pcs, f_fe_pcs)
mcs_results <- run_all_models_for_outcome("Mental Health (SF-12 MCS)", unmatched_paired, all_matched, as.formula(sf12mcs_dv_t1~smoking_group), f_adj_mcs, f_fe_mcs)
eq5d_results <- run_all_models_for_outcome("Health-Related Quality of Life (EQ-5D)", unmatched_paired, all_matched, as.formula(eq5d_t1~smoking_group), f_adj_eq5d, f_fe_eq5d)

all_results <- bind_rows(pcs_results, mcs_results, eq5d_results)

#===============================================================================
# 3. GENERATE COMPREHENSIVE TABLE
#===============================================================================

# --- Helper function for formatting ---
format_results <- function(estimate, conf.low, conf.high, p.value) {
  stars <- case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**", p.value < 0.05 ~ "*", TRUE ~ "")
  p_fmt <- case_when(p.value < 0.001 ~ "*p*<0.001", TRUE ~ paste0("*p*=", sprintf("%.3f", p.value)))
  paste0(sprintf("%.3f", estimate), stars, " [", sprintf("%.3f", conf.low), ", ", sprintf("%.3f", conf.high), "], ", p_fmt)
}

# --- Prepare data for the comprehensive table ---
comp_table_data <- all_results %>%
  filter(str_detect(term, "smoking_group")) %>%
  mutate(
    term = str_remove(term, "smoking_group"),
    Result = format_results(estimate, conf.low, conf.high, p.value),
    # For models run on specific comparisons, use that. For unmatched, use the term.
    Comparison = if_else(comparison != "All", comparison, term)
  ) %>%
  select(outcome, Comparison, model_type, Result) %>%
  pivot_wider(names_from = model_type, values_from = Result) %>%
  select(outcome, Comparison, `Unadjusted (Unmatched)`, `Fully Adjusted (Unmatched)`, `Unadjusted (Matched)`, `Fully Adjusted (DR)`, `Fixed Effects`)

# --- Create and save the gt table ---
comprehensive_gt <- gt(comp_table_data, groupname_col = "outcome", rowname_col = "Comparison") %>%
  tab_header(title = md("**Comprehensive Comparison of Model Estimates**"), subtitle = "Reference group is 'Continued Smoker'") %>%
  tab_spanner(label = md("**Unmatched Models**<br>Est. [95% CI], *p*-value"), columns = c("Unadjusted (Unmatched)", "Fully Adjusted (Unmatched)")) %>%
  tab_spanner(label = md("**Matched Models (1:3 Matching)**<br>Est. [95% CI], *p*-value"), columns = c("Unadjusted (Matched)", "Fully Adjusted (DR)", "Fixed Effects")) %>%
  cols_label(
    `Unadjusted (Unmatched)` = md("**Unadjusted**"), `Fully Adjusted (Unmatched)` = md("**Fully Adjusted**"),
    `Unadjusted (Matched)` = md("**Unadjusted**"), `Fully Adjusted (DR)` = md("**Doubly Robust**"), `Fixed Effects` = md("**Fixed Effects**")
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups())%>%
  fmt_markdown(columns = everything()) %>%
  tab_options(table.width = pct(100)) %>%
  tab_footnote(footnote = md("*Note: Reference group is 'Continued Smoker'*.<br>**p*<0.05, ***p*<0.01, ****p*<0.001"))

print(comprehensive_gt)
gtsave(comprehensive_gt, filename = paste0(OUTPUT_DIR, "Table6_Comprehensive_Results.png"))


#===============================================================================
# 4. GENERATE FOREST PLOTS
#===============================================================================

# --- Prepare data for plotting ---
plot_data <- all_results %>%
  filter(str_detect(term, "smoking_group")) %>%
  mutate(
    term = str_remove(term, "smoking_group"),
    Comparison = if_else(comparison != "All", paste("vs.", comparison), paste("vs.", term)),
    Comparison = factor(Comparison, levels = c("vs. Switcher", "vs. Dual User", "vs. Quitter")),
    model_type = factor(model_type, levels = c(
      "Fixed Effects", "Fully Adjusted (DR)", "Unadjusted (Matched)",
      "Fully Adjusted (Unmatched)", "Unadjusted (Unmatched)"
    ))
  )

# --- Function to create a forest plot ---
create_forest_plot <- function(data, plot_title, x_label) {
  ggplot(data, aes(x = estimate, y = model_type)) +
    facet_wrap(~ Comparison, scales = "free_x") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.25, linewidth = 0.8) +
    geom_point(size = 3) +
    labs(title = plot_title, x = x_label, y = "") +
    scale_color_viridis_d(guide = "none") + # No legend, colors distinguish points
    theme_minimal(base_size = 14) +
    theme(panel.spacing = unit(1.5, "lines"), strip.text = element_text(face = "bold"))
}

# --- Generate and save each plot ---
plot_pcs <- create_forest_plot(filter(plot_data, outcome == "Physical Health (SF-12 PCS)"), "Estimated Effect on Physical Health (SF-12 PCS)", "Change in PCS Score (95% CI)")
plot_mcs <- create_forest_plot(filter(plot_data, outcome == "Mental Health (SF-12 MCS)"), "Estimated Effect on Mental Health (SF-12 MCS)", "Change in MCS Score (95% CI)")
plot_eq5d <- create_forest_plot(filter(plot_data, outcome == "Health-Related Quality of Life (EQ-5D)"), "Estimated Effect on Quality of Life (EQ-5D)", "Change in EQ-5D Score (95% CI)")

print(plot_pcs)
print(plot_mcs)
print(plot_eq5d)

ggsave(paste0(OUTPUT_DIR, "Figure_ForestPlot_PCS.png"), plot = plot_pcs, width = 12, height = 5, dpi = 300, bg = "white")
ggsave(paste0(OUTPUT_DIR, "Figure_ForestPlot_MCS.png"), plot = plot_mcs, width = 12, height = 5, dpi = 300, bg = "white")
ggsave(paste0(OUTPUT_DIR, "Figure_ForestPlot_EQ5D.png"), plot = plot_eq5d, width = 12, height = 5, dpi = 300, bg = "white")
