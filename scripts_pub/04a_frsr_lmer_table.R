# ============================================================================ # 
# Linear Mixed-Effects Models: Execution Script
# Run interaction analyses with different covariate configurations
# ============================================================================ # 

# Load functions (assumes lmer_interaction_analysis_functions.R is sourced)
source("scripts_pub/function/lmer_Table_func.R")

# This script documents the statistical analysis workflow used in the manuscript.
# Individual-level clinical and biomarker data are not publicly available due to ethics and privacy restrictions.
# The script assumes that the following preprocessed objects have been loaded:
#   visitdf: a longitudinal dataframe with one row per visit, including ALSFRS-R scores & time from baseline
#   csfdf: a baseline CSF database with one row per participant.
#   etc.
#
# Main analysis:
#   Outcome: ALSFRS-R total score
#   Time variable: months from first visit
#   Random effects: patient-specific intercept and slope for time
#   Main test: time-by-variable interaction

# Configuration ------------------------------------------------------------

# Define predictor sets
core_covariates_lmer <- c("Age_init", "VC_Percent", "deltaFS", "REEC_definite", "OnsetSite")

lmer_predictors <- c("NfL_z", "GFAP_z", "pTau181_z", "pTau217_z","Ab3840_z","Ab4240_z", "Ab_status_bl",
                     "Cr_z", "CK_z", "Alb_z", "CrCys_ratio_z")

lmer_predictors_core <- c(
  "Age_init", "Sex", "VC_Percent", "deltaFS", "REEC_definite", "OnsetSite",
  "genetic_label","BMI_init","ever_Edaravone","ever_Riluzole",
  "APOE_e4_carrier","APOE_e2_carrier"
)

# Specify which covariates should interact with time
interaction_covariates_config <- c("Age_init", "deltaFS", "VC_Percent", "REEC_definite", "OnsetSite") 


# Data Preparation ---------------------------------------------------------

df_vis <- build_df_vis_baseline(
  visitdf %>% filter(!visit_post_event),
  csfdf,
  NotValidSample = NotValidSample,
  keep_multi_visit_only = TRUE,
  min_baseline_ALSFRS = 20
) %>%
  mutate(
    Ab_status_bl = factor(Ab_status_bl, levels = c("negative", "positive")),
    REEC_binary = as.numeric(REEC_definite == "definite")
    )

# Check that baseline variables are invariant within each participant after merging
baseline_check <- df_vis %>%
  group_by(SSBS_ID) %>%
  summarise(
    n_NfL_z = n_distinct(NfL_z, na.rm = TRUE),
    n_GFAP_z = n_distinct(GFAP_z, na.rm = TRUE),
    n_Age_init = n_distinct(Age_init, na.rm = TRUE),
    .groups = "drop"
  )

stopifnot(all(baseline_check$n_NfL_z <= 1, na.rm = TRUE))
stopifnot(all(baseline_check$n_GFAP_z <= 1, na.rm = TRUE))
stopifnot(all(baseline_check$n_Age_init <= 1, na.rm = TRUE))

## sanity check ----------------------------------------------------------

z_check <- df_vis %>%
  distinct(SSBS_ID, NfL_z, GFAP_z, pTau181_z, pTau217_z, Ab4240_z, Ab3840_z) %>%
  summarise(
    across(
      c(NfL_z, GFAP_z, pTau181_z, pTau217_z, Ab4240_z, Ab3840_z),
      list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE))
    )
  )

z_check

z_check_demo <- df_vis %>%
  distinct(SSBS_ID, Age_init, VC_Percent, BMI_init, deltaFS) %>%
  summarise(
    across(
      c(Age_init, VC_Percent, BMI_init, deltaFS),
      list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE))
    )
  )

z_check_demo

# Run Analyses -------------------------------------------------------------

# Biomarker predictors
all_lmer_results <- purrr::map_dfr(lmer_predictors, function(pred) {
  bind_rows(
    fit_unadjusted_lmer(pred),
    fit_age_lmer(pred, interaction_covariates = "Age_init"),
    fit_age_nfl_lmer(pred, interaction_covariates = "Age_init"),
    fit_core_lmer(pred, core_covariates = core_covariates_lmer, 
                  interaction_covariates = interaction_covariates_config)
  )
})

# Core variables (demographic, clinical)
all_corevar_lmer_results <- purrr::map_dfr(lmer_predictors_core, function(pred) {
  bind_rows(
    fit_unadjusted_lmer(pred),
    fit_age_lmer(pred, interaction_covariates = "Age_init"),
    fit_age_nfl_lmer(pred, interaction_covariates = "Age_init"),
    fit_core_lmer(pred, core_covariates = core_covariates_lmer,
                  interaction_covariates = interaction_covariates_config)
  )
})

# Combine all results
all_lmer_results <- bind_rows(
  all_corevar_lmer_results,
  all_lmer_results
)

# Format Results -----------------------------------------------------------

display_long_lmer <- all_lmer_results %>%
  mutate(
    beta_disp = fmt_ci2(beta, beta_lcl, beta_ucl),
    p_disp    = fmt_p2(p_value),
    aic_disp  = fmt_2d(AIC_full),
    daic_disp = fmt_2d(delta_AIC),
    bic_disp  = fmt_2d(BIC_full),
    r2_disp = case_when(
      is.na(R2m_full) ~ NA_character_,
      R2m_full < 0    ~ "<0",
      TRUE            ~ fmt_3d(R2m_full)
    ),
    dr2_disp = case_when(
      is.na(delta_R2m) ~ NA_character_,
      delta_R2m < 0    ~ "<0",
      delta_R2m <= 0.001 ~ "~0",
      TRUE             ~ fmt_3d(delta_R2m)
    ),
    model_key = case_when(
      base_model == "Null"    ~ "Unadj",
      base_model == "Age"     ~ "Age",
      base_model == "Age+NfL" ~ "Age+NfL",
      base_model == "Core"    ~ "Core",
      TRUE ~ base_model
    )
  ) %>%
  mutate(
    dr2_disp = if_else(model_key == "Unadj", NA_character_, dr2_disp),
    daic_disp = if_else(model_key == "Unadj", NA_character_, daic_disp)
  )

# Pivot to Wide Format -----------------------------------------------------

long_to_wide_lmer <- display_long_lmer %>%
  select(biomarker_raw, biomarker_label, model_key,
         beta_disp, p_disp, aic_disp, daic_disp, bic_disp, r2_disp, dr2_disp) %>%
  pivot_longer(
    cols = c(beta_disp, p_disp, aic_disp, daic_disp, bic_disp, r2_disp, dr2_disp),
    names_to = "metric", values_to = "value"
  ) %>%
  mutate(metric = case_when(
    metric == "beta_disp" ~ "β (95% CI)",
    metric == "p_disp"    ~ "p-value",
    metric == "aic_disp"  ~ "AIC",
    metric == "daic_disp" ~ "ΔAIC",
    metric == "bic_disp"  ~ "BIC",
    metric == "r2_disp"   ~ "Marginal R²",
    metric == "dr2_disp"  ~ "ΔMarginal R²",
    TRUE ~ metric
  )) %>%
  filter(!(model_key == "Unadj" & metric %in% c("ΔMarginal R²", "ΔAIC"))) %>%
  mutate(final_col = paste(model_key, metric, sep = " | "))

master_wide_lmer <- long_to_wide_lmer %>%
  select(biomarker_label, final_col, value) %>%
  distinct() %>%
  pivot_wider(names_from = final_col, values_from = value) %>%
  rename(Biomarker = biomarker_label)

# Add N columns for all adjustment levels
for (blk in c("Unadj", "Age", "Age+NfL", "Core")) {
  base_tag <- case_when(
    blk == "Unadj"   ~ "Null",
    blk == "Age"     ~ "Age",
    blk == "Age+NfL" ~ "Age+NfL",
    blk == "Core"    ~ "Core"
  )
  
  n_cols_block <- all_lmer_results %>%
    filter(base_model == base_tag) %>%
    distinct(biomarker_label, n_subjects, n_obs) %>%
    mutate(
      !!paste0(blk, " | N (subj)") := as.character(n_subjects),
      !!paste0(blk, " | N (obs)")  := as.character(n_obs)
    ) %>%
    select(biomarker_label, starts_with(blk))
  
  master_wide_lmer <- master_wide_lmer %>%
    left_join(n_cols_block, by = c("Biomarker" = "biomarker_label"))
}

# Generate Tables ----------------------------------------------------------

# Main analysis table (Unadj + Age+NfL + Core)
master_wide_lmer_main <- select_lmer_columns(
  master_wide_lmer,
  blocks = c("Unadj", "Age+NfL", "Core"),
  metrics = c("Beta", "p", "DeltaAIC"),
  show_n_cols = "subj"
)

ft_lmer_main <- build_lmer_flextable(
  master_wide_lmer_main, 
  include_p = "all", 
  metrics_to_show = c("DeltaAIC", "DeltaR2"),
  core_covariates = core_covariates_lmer,
  label_map = T1label  # Define T1label as needed
)

ft_lmer_main


# 🔍Supple T8 ---------------------------------------------------------------

SuppleT8_pub <- ft_lmer_main %>% 
  flextable::width(j = 1:ncol_keys(.), width = 0.55) %>% 
  flextable::width(j = c(1,2,5,9), width = 1.5) %>%
  flextable::padding(padding.left = 0, padding.right = 0, part = "all") %>% 
  set_caption(caption = flextable::as_paragraph(
    flextable::as_chunk("Supplementary Table 8. Sensitivity analyses of baseline predictors for ALSFRS-R decline",
                        props = officer::fp_text(font.size   = 12, bold = TRUE, 
                                                 font.family = "Times New Roman"))),
    word_stylename = "Table Caption", 
    fp_p = officer::fp_par(text.align = "left", padding= 5)) %>%
  flextable::font(fontname = "Times New Roman", part = "all") %>%
  flextable::fontsize(size=9, part="all")   

SuppleT8_pub

# SuppleT8_pub %>%
#   flextable::save_as_docx(path = "output/Tables/SuppleT8_lmer.docx",
#                            pr_section = officer::prop_section(
#                              page_size = officer::page_size(orient = "landscape",
#                                                             width = 11, height = 8.5),
#                              page_margins = officer::page_mar(bottom = 0.5, top = 0.5,
#                                                               left = 0.5, right = 0.5)
#                            ))


# QC - Check for singular fits in the models ---------------------------------
singularity_check <- all_lmer_results %>%
  filter(singular_full == TRUE) %>%
  select(biomarker_label, model_label, n_subjects, n_obs)

if (nrow(singularity_check) > 0) {
  message("Warning: Some models have singular fit. Consider simplifying random effects.")
  print(singularity_check)
}

singularity_check 
