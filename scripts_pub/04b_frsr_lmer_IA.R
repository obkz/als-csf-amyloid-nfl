# ========================================================================== #
# NfL Interaction Models: Execution Script
# Test whether NfL modifies biomarker-decline associations
# ========================================================================== #

# Load functions
source("scripts_pub/function/lmer_Table_IA_func.R")

# Configuration ------------------------------------------------------------

# Define covariate sets
core_covariates_nfl <- c("Age_init", "VC_Percent", "deltaFS", "REEC_definite", "OnsetSite")

# Specify which covariates should interact with time
# Example configurations:
# Age + NfL model specification
age_nfl_main_effects <- character(0)  # No additional main effects beyond NfL
age_nfl_time_interactions <- c("Age_init")

# Core adjusted model specification  
core_main_effects <- c("Age_init", "VC_Percent", "deltaFS", "REEC_definite", "OnsetSite")

# Subset that interacts with time
core_time_interactions <- c("Age_init", "deltaFS", "VC_Percent", "OnsetSite", "REEC_definite") 

# Predictor lists (exclude NfL itself from testing)
nfl_interaction_predictors <- c(
  "GFAP_z", "pTau181_z", "pTau217_z",
  "Ab3840_z", "Ab4240_z", 
  "Ab_status_bl",
  "APOE_e4_carrier")

nfl_interaction_core_vars <- "Age_init"


# Data Preparation ---------------------------------------------------------

# Assumes df_vis is already prepared with NfL_z and other predictors


# Run Analyses -------------------------------------------------------------

all_nfl_interaction_results <- purrr::map_dfr(nfl_interaction_predictors, function(pred) {
  bind_rows(
    fit_age_nfl_interaction(pred),
    fit_core_nfl_interaction(pred)
  )
})

# Age_init as NfL modifier: use sensitivity wrappers that bypass the guard clause
all_corevar_nfl_results <- purrr::map_dfr(nfl_interaction_core_vars, function(pred) {
  dplyr::bind_rows(
    fit_age_nfl_interaction_sens(pred),
    fit_core_nfl_interaction_sens(pred)
  )
})

 
# Combine all results
all_nfl_interaction_results <- bind_rows(
  all_corevar_nfl_results,
  all_nfl_interaction_results
)

# Format Results -----------------------------------------------------------

display_long_nfl <- all_nfl_interaction_results %>%
  mutate(
    main_beta_disp = fmt_ci2(main_beta, main_beta_lcl, main_beta_ucl),
    main_p_disp    = fmt_p2(main_p_value),
    int_beta_disp  = fmt_ci2(int_beta, int_beta_lcl, int_beta_ucl),
    int_p_disp     = fmt_p2(int_p_value),
    nfl_beta_disp  = fmt_ci2(nfl_beta, nfl_beta_lcl, nfl_beta_ucl),
    nfl_p_disp     = fmt_p2(nfl_p_value),
    aic_disp       = fmt_2d(AIC_full),
    daic_disp      = fmt_2d(delta_AIC),
    bic_disp       = fmt_2d(BIC_full),
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
      base_model == "Age+NfL"  ~ "Age+NfL",
      base_model == "NfL+Core" ~ "NfL+Core",
      TRUE ~ base_model
    )
  )

# Pivot to Wide Format -----------------------------------------------------

long_to_wide_nfl <- display_long_nfl %>%
  select(biomarker_raw, biomarker_label, model_key,
         main_beta_disp, main_p_disp, 
         int_beta_disp, int_p_disp,
         nfl_beta_disp, nfl_p_disp,
         aic_disp, daic_disp, bic_disp, r2_disp, dr2_disp) %>%
  pivot_longer(
    cols = c(main_beta_disp, main_p_disp, int_beta_disp, int_p_disp,
             nfl_beta_disp, nfl_p_disp, aic_disp, daic_disp, bic_disp, 
             r2_disp, dr2_disp),
    names_to = "metric", values_to = "value"
  ) %>%
  mutate(metric = case_when(
    metric == "main_beta_disp" ~ "Main β (95% CI)",
    metric == "main_p_disp"    ~ "Main p-value",
    metric == "int_beta_disp"  ~ "Interaction β (95% CI)",
    metric == "int_p_disp"     ~ "Interaction p-value",
    metric == "nfl_beta_disp"  ~ "NfL β (95% CI)",
    metric == "nfl_p_disp"     ~ "NfL p-value",
    metric == "aic_disp"       ~ "AIC",
    metric == "daic_disp"      ~ "ΔAIC",
    metric == "bic_disp"       ~ "BIC",
    metric == "r2_disp"        ~ "Marginal R²",
    metric == "dr2_disp"       ~ "ΔMarginal R²",
    TRUE ~ metric
  )) %>%
  mutate(final_col = paste(model_key, metric, sep = " | "))

master_wide_nfl <- long_to_wide_nfl %>%
  select(biomarker_label, final_col, value) %>%
  distinct() %>%
  pivot_wider(names_from = final_col, values_from = value) %>%
  rename(Biomarker = biomarker_label)

# Add N columns for all blocks
for (blk in c("Age+NfL", "NfL+Core")) {
  n_cols_block <- all_nfl_interaction_results %>%
    filter(base_model == blk) %>%
    distinct(biomarker_label, n_subjects, n_obs) %>%
    mutate(
      !!paste0(blk, " | N (subj)") := as.character(n_subjects),
      !!paste0(blk, " | N (obs)")  := as.character(n_obs)
    ) %>%
    select(biomarker_label, starts_with(blk))
  
  master_wide_nfl <- master_wide_nfl %>%
    left_join(n_cols_block, by = c("Biomarker" = "biomarker_label"))
}

# Generate Tables ----------------------------------------------------------

master_wide_nfl_core <- select_nfl_interaction_columns(
  master_wide_nfl,
  blocks = c("NfL+Core"),
  metrics = c("NflBeta", "InteractionBeta", "InteractionP","DeltaAIC"),
  show_n_cols = "subj"
)

ft_nfl_core <- build_nfl_interaction_flextable(
  master_wide_nfl_core, 
  include_main_p = "none",
  include_int_p = "all",
  include_nfl = "all",
  metrics_to_show = c("DeltaAIC"),
  core_covariates = core_covariates_nfl,
  label_map = T1label
)

ft_nfl_core


# 📊 Table2 ----------------------------------------------------------------

Table2_pub <- ft_nfl_core %>% 
  set_caption(caption = flextable::as_paragraph(
    flextable::as_chunk("Table 2: NfL × Biomarker Interaction Analyses on Functional Decline",
                        props = officer::fp_text(font.size   = 12, bold = TRUE, 
                                                 font.family = "Times New Roman"))),
    word_stylename = "Table Caption", 
    fp_p = officer::fp_par(text.align = "left", padding= 5)) %>%
  flextable::font(fontname = "Times New Roman", part = "all") %>%
  flextable::fontsize(size=9, part="all") 

Table2_pub 

# Table2_pub %>%
#   flextable::save_as_docx(path = "output/Tables/Table2_r.docx",
#                            pr_section = officer::prop_section(
#                              page_size = officer::page_size(orient = "landscape",
#                                                             width = 11, height = 8.5),
#                              page_margins = officer::page_mar(bottom = 0.5, top = 0.5,
#                                                               left = 0.5, right = 0.5)
#                            ))

# Additional Checks --------------------------------------------------------

# Optional: Check for convergence issues
singularity_check_nfl <- all_nfl_interaction_results %>%
  filter(singular_full == TRUE) %>%
  select(biomarker_label, model_label, n_subjects, n_obs)

if (nrow(singularity_check_nfl) > 0) {
  message("Warning: Some NfL interaction models have singular fit.")
  print(singularity_check_nfl)
}

# Optional: Identify significant 3-way interactions
significant_interactions <- all_nfl_interaction_results %>%
  filter(int_p_value < 0.05) %>%
  arrange(int_p_value) %>%
  select(biomarker_label, model_label, int_beta, int_p_value, delta_AIC)

if (nrow(significant_interactions) > 0) {
  message("\nSignificant 3-way interactions (Time*NfL*Biomarker, p < 0.05):")
  print(significant_interactions)
}


# NULL CHECK: Identify any models where interaction results are missing (NA)
all_nfl_interaction_results %>%
  filter(is.na(int_beta) | is.na(int_p_value))


all_nfl_interaction_results %>%
  summarise(
    missing_int_beta = sum(is.na(int_beta)),
    missing_int_p    = sum(is.na(int_p_value)),
    missing_nfl_beta = sum(is.na(nfl_beta)),
    singular_n       = sum(singular_full)
  )
