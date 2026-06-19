
# Exploratory analysis for revision ---------------------------------------

# Age x NfL x Time and APOE4 x NfL x Time sensitivity analyses ------------
# Model 1 (NfL only):  Time * NfL + core covariates
# Model 2 (3-way):     Time * NfL * modifier + core covariates

ctrl <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))

run_modifier_comparison <- function(modifier_col, modifier_label) {
  
  af <- df_vis %>%
    select(ALSFRSR_Total, SSBS_ID, MonthsFromFirstVisit,
           NfL_z, all_of(modifier_col),
           Age_init, deltaFS, VC_Percent, OnsetSite, REEC_definite) %>%
    drop_na() %>%
    group_by(SSBS_ID) %>%
    filter(n_distinct(MonthsFromFirstVisit) >= 2) %>%
    ungroup()
  
  # Base: Time * NfL + core covariates
  # Full: Time * NfL * modifier + core covariates
  fit_nfl_interaction(
    biomarker_col          = modifier_col,
    base_covariates        = character(0),
    interaction_covariates = core_time_interactions,
    base_label             = "NfL+Core",
    model_label            = "NfL+Core adjusted",
    data                   = af
  ) %>%
    mutate(biomarker_label = modifier_label)
}


# Sanity check: Run models to make sure it works in the same way as the main analysis. 

# modifier_results <- bind_rows(
#   run_modifier_comparison("Age_init",        "Age at baseline"),
#   run_modifier_comparison("APOE_e4_carrier", "APOE ε4 status"),
#   run_modifier_comparison("Ab3840_z",        "Aβ38/40 ratio"),
#   run_modifier_comparison("Ab4240_z",        "Aβ42/40 ratio"),
#   run_modifier_comparison("Ab_status_bl",    "Aβ status (positive vs negative)")
# )


# For Revision
modifier_results <- bind_rows(
  run_modifier_comparison("Age_init",        "Age at baseline"),
  run_modifier_comparison("APOE_e4_carrier", "APOE ε4 status")
)

modifier_results %>%
  select(biomarker_label, nfl_beta, nfl_beta_lcl, nfl_beta_ucl, nfl_p_value,
         int_beta, int_beta_lcl, int_beta_ucl, int_p_value,
         delta_AIC, n_subjects, n_obs, singular_full)


# Age & ApoE ---------------------------------------------------------------


display_long_modifier <- modifier_results %>%
  mutate(
    int_beta_disp = fmt_ci2(int_beta, int_beta_lcl, int_beta_ucl),
    int_p_disp    = fmt_p2(int_p_value),
    nfl_beta_disp = fmt_ci2(nfl_beta, nfl_beta_lcl, nfl_beta_ucl),
    nfl_p_disp    = fmt_p2(nfl_p_value),
    daic_disp     = fmt_2d(delta_AIC),
    model_key     = "NfL+Core"
  )

long_to_wide_modifier <- display_long_modifier %>%
  select(biomarker_label, model_key,
         int_beta_disp, int_p_disp,
         nfl_beta_disp, nfl_p_disp,
         daic_disp) %>%
  pivot_longer(
    cols = c(int_beta_disp, int_p_disp,
             nfl_beta_disp, nfl_p_disp,
             daic_disp),
    names_to = "metric", values_to = "value"
  ) %>%
  mutate(metric = case_when(
    metric == "int_beta_disp" ~ "Interaction β (95% CI)",
    metric == "int_p_disp"    ~ "Interaction p-value",
    metric == "nfl_beta_disp" ~ "NfL β (95% CI)",
    metric == "nfl_p_disp"    ~ "NfL p-value",
    metric == "daic_disp"     ~ "ΔAIC",
    TRUE ~ metric
  )) %>%
  mutate(final_col = paste(model_key, metric, sep = " | "))

master_wide_modifier <- long_to_wide_modifier %>%
  select(biomarker_label, final_col, value) %>%
  distinct() %>%
  pivot_wider(names_from = final_col, values_from = value) %>%
  rename(Biomarker = biomarker_label)

# Add N column
master_wide_modifier <- master_wide_modifier %>%
  left_join(
    modifier_results %>%
      mutate(`NfL+Core | N (subj)` = as.character(n_subjects)) %>%
      select(biomarker_label, `NfL+Core | N (subj)`) %>%
      rename(Biomarker = biomarker_label),
    by = "Biomarker"
  )

# Table using existing functions
master_wide_modifier_sel <- select_nfl_interaction_columns(
  master_wide_modifier,
  blocks      = c("NfL+Core"),
  metrics     = c("NflBeta", "InteractionBeta", "InteractionP", "DeltaAIC"),
  show_n_cols = "subj"
)

ft_modifier <- build_nfl_interaction_flextable(
  master_wide_modifier_sel,
  include_main_p  = "none",
  include_int_p   = "all",
  include_nfl     = "all",
  metrics_to_show = c("DeltaAIC"),
  core_covariates = core_covariates_nfl) %>%
  flextable::compose(
    i = 2, j = 1,
    value = flextable::as_paragraph("Modifier"),
    part = "header"
  ) %>%
  flextable::font(fontname = "Times New Roman", part = "all") %>%
  flextable::fontsize(size = 9, part = "all") %>% 
  # title caption
  set_caption(caption = flextable::as_paragraph(
    flextable::as_chunk("Table R1 for Reviewer: Sensitivity analysis for Time × NfL × Modifier interactions.",
                        props = officer::fp_text(font.size   = 12, bold = TRUE, 
                                                 font.family = "Times New Roman"))),
    word_stylename = "Table Caption", 
    fp_p = officer::fp_par(text.align = "left", padding= 5))

ft_modifier


# Time*NfL*age ------------------------------------------------------------

sens_predictors_all <- c("Ab3840_z", "Ab4240_z", "Ab_status_bl")

# Extract both 3-way interaction terms from REML fit
sens_dual3way_results_full <- purrr::map_dfr(sens_predictors_all, function(bio) {
  
  af <- df_vis %>%
    select(ALSFRSR_Total, SSBS_ID, MonthsFromFirstVisit,
           NfL_z, Age_init,
           deltaFS, VC_Percent, OnsetSite, REEC_definite,
           all_of(bio)) %>%
    drop_na() %>%
    group_by(SSBS_ID) %>%
    filter(n_distinct(MonthsFromFirstVisit) >= 2) %>%
    ungroup()
  
  ctrl <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
  
  base_formula <- as.formula(paste(
    "ALSFRSR_Total ~",
    "MonthsFromFirstVisit * NfL_z * Age_init +",
    "MonthsFromFirstVisit * deltaFS +",
    "MonthsFromFirstVisit * VC_Percent +",
    "MonthsFromFirstVisit * OnsetSite +",
    "MonthsFromFirstVisit * REEC_definite +",
    "(MonthsFromFirstVisit | SSBS_ID)"
  ))
  
  full_formula <- as.formula(paste(
    "ALSFRSR_Total ~",
    "MonthsFromFirstVisit * NfL_z *", bio, "+",
    "MonthsFromFirstVisit * NfL_z * Age_init +",
    "MonthsFromFirstVisit * deltaFS +",
    "MonthsFromFirstVisit * VC_Percent +",
    "MonthsFromFirstVisit * OnsetSite +",
    "MonthsFromFirstVisit * REEC_definite +",
    "(MonthsFromFirstVisit | SSBS_ID)"
  ))
  
  fit_base_ML   <- lmerTest::lmer(base_formula, data = af, REML = FALSE, control = ctrl)
  fit_full_ML   <- lmerTest::lmer(full_formula, data = af, REML = FALSE, control = ctrl)
  fit_full_REML <- lmerTest::lmer(full_formula, data = af, REML = TRUE,  control = ctrl)
  
  is_factor <- is.factor(af[[bio]])
  if (is_factor) {
    lvs           <- levels(af[[bio]])
    bio_int_term  <- paste0("MonthsFromFirstVisit:NfL_z:", bio, lvs[2])
    display_label <- sprintf("%s (%s vs %s)",
                             stringr::str_remove(bio, "_bl$|_z$"),
                             lvs[2], lvs[1])
  } else {
    bio_int_term  <- paste0("MonthsFromFirstVisit:NfL_z:", bio)
    display_label <- stringr::str_remove(bio, "_bl$|_z$")
  }
  age_int_term <- "MonthsFromFirstVisit:NfL_z:Age_init"
  
  coef_tbl <- broom.mixed::tidy(fit_full_REML, effects = "fixed", conf.int = TRUE)
  
  bio_coef <- coef_tbl %>% filter(term == bio_int_term)
  age_coef <- coef_tbl %>% filter(term == age_int_term)
  
  tibble::tibble(
    biomarker_label  = display_label,
    # Time x NfL x Biomarker (primary)
    bio_beta         = bio_coef$estimate,
    bio_beta_lcl     = bio_coef$conf.low,
    bio_beta_ucl     = bio_coef$conf.high,
    bio_p_value      = bio_coef$p.value,
    # Time x NfL x Age (shown for transparency)
    age_beta         = age_coef$estimate,
    age_beta_lcl     = age_coef$conf.low,
    age_beta_ucl     = age_coef$conf.high,
    age_p_value      = age_coef$p.value,
    delta_AIC        = AIC(fit_full_ML) - AIC(fit_base_ML),
    n_subjects       = n_distinct(af$SSBS_ID),
    n_obs            = nrow(af),
    singular_full    = lme4::isSingular(fit_full_REML)
  )
})

# Build flextable
sens_table_data <- sens_dual3way_results_full %>%
  filter(biomarker_label %in% c("Ab4240", "Ab3840",
                                "Ab_status (positive vs negative)")) %>%
  mutate(
    Biomarker     = case_when(
      biomarker_label == "Ab4240"                           ~ "Aβ42/40 ratio",
      biomarker_label == "Ab3840"                           ~ "Aβ38/40 ratio",
      biomarker_label == "Ab_status (positive vs negative)" ~ "Aβ status (positive vs negative)",
      TRUE ~ biomarker_label
    ),
    bio_beta_disp = fmt_ci2(bio_beta, bio_beta_lcl, bio_beta_ucl),
    bio_p_disp    = fmt_p2(bio_p_value),
    age_beta_disp = fmt_ci2(age_beta, age_beta_lcl, age_beta_ucl),
    age_p_disp    = fmt_p2(age_p_value),
    daic_disp     = fmt_2d(delta_AIC),
    n_subj_disp   = as.character(n_subjects)
  ) %>%
  select(Biomarker, 
         bio_beta_disp, bio_p_disp, 
         age_beta_disp, age_p_disp,
         daic_disp, n_subj_disp)

sens_table_data



## Flextable ---------------------------------------------------------------

brd <- officer::fp_border(color = "black", width = 1)

ft_sens <- flextable::flextable(sens_table_data) %>%
  flextable::set_header_labels(
    Biomarker     = "Biomarker",
    bio_beta_disp = "β (95% CI)",
    bio_p_disp    = "p-value",
    age_beta_disp = "β (95% CI)",
    age_p_disp    = "p-value",
    daic_disp     = "ΔAIC",
    n_subj_disp   = "N (subjects)"
  ) %>%
  flextable::add_header_row(
    values     = c("", 
                   "Time × NfL × Biomarker", "",
                   "Time × NfL × Age", "",
                   "", ""),
    colwidths  = c(1, 1, 1, 1, 1, 1, 1)
  ) %>%
  flextable::merge_h(part = "header") %>% 
  flextable::align(align = "center", part = "header") %>%
  flextable::align(align = "center", j = 2:7) %>%
  flextable::bold(part = "header") %>%
  flextable::bold(j = 1) %>%
  flextable::border_remove() %>%
  flextable::hline_top(border = brd, part = "header") %>%
  flextable::hline(i = 1, border = brd, part = "header") %>%
  flextable::hline(i = 2, border = brd, part = "header") %>%
  flextable::hline_bottom(border = brd, part = "body") %>%
  flextable::autofit() %>%
  flextable::add_footer_lines(values = c(
    paste(
      "Sensitivity analysis: both Time × NfL × Biomarker and Time × NfL × Age were",
      "included simultaneously in the same model.",
      "The base model includes Time × NfL × Age three-way interactions with core clinical features all interacted with time.",
      "Full model adds Time × NfL × Biomarker to the base model.",
      "ΔAIC = AIC(full) − AIC(base); values < −2 indicate improved model fit.",
      "Continuous variables standardised (per 1-SD). N = 218 (complete cases, ≥2 visits)."
    ),
    "Abbreviations: Aβ, amyloid-β; AIC, Akaike information criterion; ALSFRS-R, Amyotrophic Lateral Sclerosis Functional Rating Scale–Revised; β, regression coefficient; CI, confidence interval; ΔFS, disease progression rate from symptom onset to baseline; NfL, neurofilament light chain; REEC, Revised El Escorial Criteria; %VC, percent predicted vital capacity"
  )) %>%
  flextable::hline_top(border = brd, part = "footer")

ft_sens_pub <- ft_sens %>%
  flextable::font(fontname = "Times New Roman", part = "all") %>%
  flextable::fontsize(size = 9, part = "all") %>% 
  # title caption
  set_caption(caption = flextable::as_paragraph(
    flextable::as_chunk("Table R2 for Reviewer: Sensitivity analysis including both Time × NfL × Biomarker and Time × NfL × Age",
                        props = officer::fp_text(font.size   = 12, bold = TRUE, 
                                                 font.family = "Times New Roman"))),
    word_stylename = "Table Caption", 
    fp_p = officer::fp_par(text.align = "left", padding= 5)) 


ft_sens_pub


# Time*NfL*APOE4 -------------------------------------------------------------------

# Sensitivity analysis for reviewer: Time x NfL x APOE4 as confounder
# Same structure as Age sensitivity but with APOE_e4_carrier instead of Age_init

sens_dual3way_apoe_results <- purrr::map_dfr(
  c("Ab3840_z", "Ab4240_z", "Ab_status_bl"),
  function(bio) {
    
    af <- df_vis %>%
      select(ALSFRSR_Total, SSBS_ID, MonthsFromFirstVisit,
             NfL_z, APOE_e4_carrier, Age_init,
             deltaFS, VC_Percent, OnsetSite, REEC_definite,
             all_of(bio)) %>%
      drop_na() %>%
      group_by(SSBS_ID) %>%
      filter(n_distinct(MonthsFromFirstVisit) >= 2) %>%
      ungroup()
    
    ctrl <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
    
    # Base model: Time*NfL*APOE4 + core covariates (Age as 2-way only)
    base_formula <- as.formula(paste(
      "ALSFRSR_Total ~",
      "MonthsFromFirstVisit * NfL_z * APOE_e4_carrier +",
      "MonthsFromFirstVisit * Age_init +",
      "MonthsFromFirstVisit * deltaFS +",
      "MonthsFromFirstVisit * VC_Percent +",
      "MonthsFromFirstVisit * OnsetSite +",
      "MonthsFromFirstVisit * REEC_definite +",
      "(MonthsFromFirstVisit | SSBS_ID)"
    ))
    
    # Full model: adds Time*NfL*Biomarker on top
    full_formula <- as.formula(paste(
      "ALSFRSR_Total ~",
      "MonthsFromFirstVisit * NfL_z *", bio, "+",
      "MonthsFromFirstVisit * NfL_z * APOE_e4_carrier +",
      "MonthsFromFirstVisit * Age_init +",
      "MonthsFromFirstVisit * deltaFS +",
      "MonthsFromFirstVisit * VC_Percent +",
      "MonthsFromFirstVisit * OnsetSite +",
      "MonthsFromFirstVisit * REEC_definite +",
      "(MonthsFromFirstVisit | SSBS_ID)"
    ))
    
    fit_base_ML   <- lmerTest::lmer(base_formula, data = af, REML = FALSE, control = ctrl)
    fit_full_ML   <- lmerTest::lmer(full_formula, data = af, REML = FALSE, control = ctrl)
    fit_full_REML <- lmerTest::lmer(full_formula, data = af, REML = TRUE,  control = ctrl)
    
    is_factor <- is.factor(af[[bio]])
    if (is_factor) {
      lvs          <- levels(af[[bio]])
      bio_int_term <- paste0("MonthsFromFirstVisit:NfL_z:", bio, lvs[2])
      display_label <- sprintf("%s (%s vs %s)",
                               stringr::str_remove(bio, "_bl$|_z$"),
                               lvs[2], lvs[1])
    } else {
      bio_int_term  <- paste0("MonthsFromFirstVisit:NfL_z:", bio)
      display_label <- stringr::str_remove(bio, "_bl$|_z$")
    }
    apoe_int_term <- "MonthsFromFirstVisit:NfL_z:APOE_e4_carrier"
    
    coef_tbl <- broom.mixed::tidy(fit_full_REML, effects = "fixed", conf.int = TRUE)
    
    bio_coef  <- coef_tbl %>% filter(term == bio_int_term)
    apoe_coef <- coef_tbl %>% filter(term == apoe_int_term)
    
    tibble::tibble(
      biomarker_label = display_label,
      bio_beta        = bio_coef$estimate,
      bio_beta_lcl    = bio_coef$conf.low,
      bio_beta_ucl    = bio_coef$conf.high,
      bio_p_value     = bio_coef$p.value,
      apoe_beta       = apoe_coef$estimate,
      apoe_beta_lcl   = apoe_coef$conf.low,
      apoe_beta_ucl   = apoe_coef$conf.high,
      apoe_p_value    = apoe_coef$p.value,
      delta_AIC       = AIC(fit_full_ML) - AIC(fit_base_ML),
      n_subjects      = n_distinct(af$SSBS_ID),
      n_obs           = nrow(af),
      singular_full   = lme4::isSingular(fit_full_REML)
    )
  }
)

# Quick check
sens_dual3way_apoe_results %>%
  select(biomarker_label, bio_beta, bio_beta_lcl, bio_beta_ucl,
         bio_p_value, apoe_beta, apoe_p_value, delta_AIC,
         n_subjects, n_obs, singular_full)

## Flextable ------------------------------------------------------------------------

sens_table_data_apoe <- sens_dual3way_apoe_results %>%
  filter(biomarker_label %in% c("Ab4240", "Ab3840",
                                "Ab_status (positive vs negative)")) %>%
  mutate(
    Biomarker      = case_when(
      biomarker_label == "Ab4240"                           ~ "Aβ42/40 ratio",
      biomarker_label == "Ab3840"                           ~ "Aβ38/40 ratio",
      biomarker_label == "Ab_status (positive vs negative)" ~ "Aβ status (positive vs negative)",
      TRUE ~ biomarker_label
    ),
    bio_beta_disp  = fmt_ci2(bio_beta,  bio_beta_lcl,  bio_beta_ucl),
    bio_p_disp     = fmt_p2(bio_p_value),
    apoe_beta_disp = fmt_ci2(apoe_beta, apoe_beta_lcl, apoe_beta_ucl),
    apoe_p_disp    = fmt_p2(apoe_p_value),
    daic_disp      = fmt_2d(delta_AIC),
    n_subj_disp    = as.character(n_subjects)
  ) %>%
  select(Biomarker,
         bio_beta_disp, bio_p_disp,
         apoe_beta_disp, apoe_p_disp,
         daic_disp, n_subj_disp)

brd <- officer::fp_border(color = "black", width = 1)

ft_sens_apoe <- flextable::flextable(sens_table_data_apoe) %>%
  flextable::set_header_labels(
    Biomarker      = "Biomarker",
    bio_beta_disp  = "β (95% CI)",
    bio_p_disp     = "p-value",
    apoe_beta_disp = "β (95% CI)",
    apoe_p_disp    = "p-value",
    daic_disp      = "ΔAIC",
    n_subj_disp    = "N (subjects)"
  ) %>%
  flextable::add_header_row(
    values    = c("",
                  "Time × NfL × Biomarker", "",
                  "Time × NfL × APOE4", "",
                  "", ""),
    colwidths = c(1, 1, 1, 1, 1, 1, 1)
  ) %>%
  flextable::merge_h(part = "header") %>%
  flextable::align(align = "center", part = "header") %>%
  flextable::align(align = "center", j = 2:7) %>%
  flextable::bold(part = "header") %>%
  flextable::bold(j = 1) %>%
  flextable::border_remove() %>%
  flextable::hline_top(border = brd, part = "header") %>%
  flextable::hline(i = 1, border = brd, part = "header") %>%
  flextable::hline(i = 2, border = brd, part = "header") %>%
  flextable::hline_bottom(border = brd, part = "body") %>%
  flextable::autofit() %>%
  flextable::add_footer_lines(values = c(
    paste(
      "Sensitivity analysis: both Time × NfL × Biomarker and Time × NfL × APOE4 were",
      "included simultaneously in the same model.",
      "Base model: Time × NfL × APOE4 + Time × Age + Time × deltaFS + Time × VC% +",
      "Time × onset site + Time × diagnostic certainty + (Time | patient).",
      "Full model adds Time × NfL × Biomarker to the base model.",
      "ΔAIC = AIC(full) − AIC(base); values < −2 indicate improved model fit.",
      "Continuous variables standardised (per 1-SD).",
      "N = 172 (complete cases with APOE4 data available and ≥2 visits)."
    ),
    "Aβ, amyloid-β; AIC, Akaike information criterion; APOE, apolipoprotein E; β, regression coefficient; CI, confidence interval; ΔAIC, change in AIC; NfL, neurofilament light."
  )) %>%
  flextable::hline_top(border = brd, part = "footer")

ft_sens_apoe_pub <- ft_sens_apoe %>%
  flextable::font(fontname = "Times New Roman", part = "all") %>%
  flextable::fontsize(size = 9, part = "all") %>% 
  # title caption
  set_caption(caption = flextable::as_paragraph(
    flextable::as_chunk("Table R3 for Reviewer: Sensitivity analysis including both Time × NfL × Biomarker and Time × NfL × APOE",
                        props = officer::fp_text(font.size   = 12, bold = TRUE, 
                                                 font.family = "Times New Roman"))),
    word_stylename = "Table Caption", 
    fp_p = officer::fp_par(text.align = "left", padding= 5)) %>%
  flextable::font(fontname = "Times New Roman", part = "all") %>%
  flextable::fontsize(size=9, part="all")   


ft_sens_apoe_pub



# Sanity check ------------------------------------------------------------

# n=218
sens_dual3way_results_full %>%
  select(biomarker_label, n_subjects, n_obs, singular_full)

# n=172
sens_dual3way_apoe_results %>%
  select(biomarker_label, n_subjects, n_obs, singular_full)


sens_dual3way_results_full %>%
  summarise(
    missing_bio_beta = sum(is.na(bio_beta)),
    missing_bio_p    = sum(is.na(bio_p_value)),
    missing_age_beta = sum(is.na(age_beta)),
    singular_n       = sum(singular_full)
  )

sens_dual3way_apoe_results %>%
  summarise(
    missing_bio_beta  = sum(is.na(bio_beta)),
    missing_bio_p     = sum(is.na(bio_p_value)),
    missing_apoe_beta = sum(is.na(apoe_beta)),
    singular_n        = sum(singular_full)
  )


# Export ------------------------------------------------------------------

library(readxl)
library(officer)

ft_modifier
ft_sens_pub
ft_sens_apoe_pub 

# # Create Word document
# RSuppleT_doc <- read_docx()
# 
# # Add first table
# RSuppleT_doc <- body_add_flextable(RSuppleT_doc, value = ft_modifier)
# RSuppleT_doc <- body_end_section_landscape(RSuppleT_doc)
# cat("Successfully added: ft_modifier (landscape)\n")
# 
# # Add remaining tables
# RSuppleT_doc <- add_flex_with_orientation(RSuppleT_doc, ft_sens_pub, "landscape", "SuppleT2_pub")
# RSuppleT_doc <- add_flex_with_orientation(RSuppleT_doc, ft_sens_apoe_pub, "landscape", "SuppleT3_pub")
# 
# # Save
# print(RSuppleT_doc, target = "For_Reviewer_Tables_1.docx")
# cat("\n=== Document saved successfully ===\n")
# 
