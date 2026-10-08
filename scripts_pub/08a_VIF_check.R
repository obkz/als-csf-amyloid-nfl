# Multicollinearity Diagnostics: Role-Based Term Extraction ------------------
# Self-contained diagnostic script for VIF reporting 
# Requires the following objects to already exist in the session (built by
# the main LME execution scripts):
#   df_vis                      - analysis dataframe
#   lmer_predictors              - biomarker list for Supple T8 models
#   nfl_interaction_predictors   - biomarker list for Table 2 (3-way) models
#   core_covariates_lmer         - Core covariate set for Supple T8 models
#   interaction_covariates_config - covariates interacting with time (Supple T8)
#   core_main_effects             - Core covariate set for Table 2 models
#   core_time_interactions        - covariates interacting with time (Table 2)

library(performance)

time_var <- "MonthsFromFirstVisit"

# --- Model-only fitters (mirror fit_core_lmer() / fit_core_nfl_interaction()) --

fit_core_lmer_model_only <- function(predictor_col,
                                     core_covariates = core_covariates_lmer,
                                     interaction_covariates = interaction_covariates_config,
                                     data = df_vis) {
  
  re_term <- paste0("(", time_var, " | SSBS_ID)")
  main_effects <- setdiff(core_covariates, interaction_covariates)
  
  rhs_parts <- c(paste(time_var, "*", predictor_col))
  if (length(interaction_covariates) > 0) {
    rhs_parts <- c(rhs_parts, paste(time_var, "*", interaction_covariates))
  }
  if (length(main_effects) > 0) {
    rhs_parts <- c(rhs_parts, main_effects)
  }
  rhs <- paste(rhs_parts, collapse = " + ")
  full_formula <- as.formula(paste("ALSFRSR_Total ~", rhs, "+", re_term))
  
  model_vars <- c("ALSFRSR_Total", "SSBS_ID", time_var, predictor_col,
                  core_covariates, interaction_covariates)
  
  analysis_frame <- data %>%
    dplyr::select(dplyr::all_of(model_vars)) %>%
    tidyr::drop_na() %>%
    dplyr::group_by(SSBS_ID) %>%
    dplyr::filter(dplyr::n_distinct(.data[[time_var]]) >= 2) %>%
    dplyr::ungroup()
  
  ctrl <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
  lmerTest::lmer(full_formula, data = analysis_frame, REML = TRUE, control = ctrl)
}

fit_core_nfl_interaction_model_only <- function(biomarker_col,
                                                base_covariates = core_main_effects,
                                                interaction_covariates = core_time_interactions,
                                                data = df_vis) {
  
  nfl_var <- "NfL_z"
  re_term <- paste0("(", time_var, " | SSBS_ID)")
  
  rhs_parts <- c(paste(time_var, "*", nfl_var, "*", biomarker_col))
  if (length(interaction_covariates) > 0) {
    rhs_parts <- c(rhs_parts, paste(time_var, "*", interaction_covariates))
  }
  if (length(base_covariates) > 0) {
    rhs_parts <- c(rhs_parts, base_covariates)
  }
  rhs <- paste(rhs_parts, collapse = " + ")
  full_formula <- as.formula(paste("ALSFRSR_Total ~", rhs, "+", re_term))
  
  model_vars <- c("ALSFRSR_Total", "SSBS_ID", time_var, nfl_var, biomarker_col,
                  base_covariates, interaction_covariates)
  
  analysis_frame <- data %>%
    dplyr::select(dplyr::all_of(model_vars)) %>%
    tidyr::drop_na() %>%
    dplyr::group_by(SSBS_ID) %>%
    dplyr::filter(dplyr::n_distinct(.data[[time_var]]) >= 2) %>%
    dplyr::ungroup()
  
  ctrl <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
  lmerTest::lmer(full_formula, data = analysis_frame, REML = TRUE, control = ctrl)
}

# --- Table A: Supple T8 companion (single-biomarker Core models, 2-way) ----

extract_vif_lmer <- function(biomarker_col) {
  fit <- tryCatch(fit_core_lmer_model_only(biomarker_col), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  
  vif_tbl <- tryCatch(performance::check_collinearity(fit), error = function(e) NULL)
  if (is.null(vif_tbl)) return(NULL)
  
  main_term <- biomarker_col
  time_term <- paste0(time_var, ":", biomarker_col)
  predictor_terms <- c(main_term, time_term)
  
  role_row <- tibble::tibble(
    biomarker = biomarker_col,
    "Biomarker (main effect)" = vif_tbl$VIF[vif_tbl$Term == main_term][1],
    "Time × Biomarker"        = vif_tbl$VIF[vif_tbl$Term == time_term][1]
  )
  
  covariate_vifs <- vif_tbl$VIF[!vif_tbl$Term %in% predictor_terms]
  role_row$covariate_VIF_range <- sprintf("%.2f–%.2f", min(covariate_vifs), max(covariate_vifs))
  
  role_row
}

SuppleT8_vif_table <- purrr::map_dfr(lmer_predictors, extract_vif_lmer)

SuppleT8_vif_table




# Global max VIF check across all Supple T8 fully adjusted (Core) models ----

# Predictors that actually have a Core-adjusted model in Supple T8
supple8_core_predictors <- all_lmer_results %>%
  dplyr::filter(base_model == "Core") %>%
  dplyr::distinct(biomarker_raw) %>%
  dplyr::pull(biomarker_raw)

supple8_core_predictors
length(supple8_core_predictors)

supple8_all_vifs <- purrr::map_dfr(supple8_core_predictors, function(pred) {
  
  fit <- tryCatch(
    fit_core_lmer_model_only(pred),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(tibble::tibble(
      predictor = pred,
      term = NA_character_,
      vif = NA_real_,
      status = "fit failed"
    ))
  }
  
  vif_tbl <- tryCatch(
    performance::check_collinearity(fit),
    error = function(e) NULL
  )
  
  if (is.null(vif_tbl)) {
    return(tibble::tibble(
      predictor = pred,
      term = NA_character_,
      vif = NA_real_,
      status = "VIF failed"
    ))
  }
  
  tibble::tibble(
    predictor = pred,
    term = vif_tbl$Term,
    vif = vif_tbl$VIF,
    status = "OK"
  )
})

supple8_all_vifs %>%
  dplyr::filter(status != "OK")

supple8_all_vifs %>%
  dplyr::filter(status == "OK") %>%
  dplyr::distinct(predictor) %>%
  dplyr::summarise(
    n_models_checked = dplyr::n(),
    max_vif = max(supple8_all_vifs$vif, na.rm = TRUE)
  )

# Top rows by VIF, across all terms and all Supple T8 models
supple8_all_vifs %>% dplyr::filter(status == "OK") %>% dplyr::slice_max(vif, n = 10)



# --- Table B: Table 2 companion (Time*NfL*Biomarker three-way models) ------

extract_vif_nfl <- function(biomarker_col) {
  fit <- tryCatch(fit_core_nfl_interaction_model_only(biomarker_col), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  
  vif_tbl <- tryCatch(performance::check_collinearity(fit), error = function(e) NULL)
  if (is.null(vif_tbl)) return(NULL)
  
  lookup <- function(term) {
    val <- vif_tbl$VIF[vif_tbl$Term == term]
    if (length(val) == 0) NA_real_ else val[1]
  }
  
  predictor_terms <- c(
    "NfL_z", paste0(time_var, ":NfL_z"),
    biomarker_col, paste0(time_var, ":", biomarker_col),
    paste0("NfL_z:", biomarker_col), paste0(time_var, ":NfL_z:", biomarker_col)
  )
  
  role_row <- tibble::tibble(
    biomarker = biomarker_col,
    "NfL (main effect)"      = lookup(predictor_terms[1]),
    "Time × NfL"             = lookup(predictor_terms[2]),
    "Biomarker (main effect)" = lookup(predictor_terms[3]),
    "Time × Biomarker"        = lookup(predictor_terms[4]),
    "NfL × Biomarker"         = lookup(predictor_terms[5]),
    "Time × NfL × Biomarker"  = lookup(predictor_terms[6])
  )
  
  covariate_vifs <- vif_tbl$VIF[!vif_tbl$Term %in% predictor_terms]
  role_row$covariate_VIF_range <- sprintf("%.2f–%.2f", min(covariate_vifs), max(covariate_vifs))
  
  role_row
}


Table2_vif_table <- purrr::map_dfr(nfl_interaction_predictors, extract_vif_nfl)

Table2_vif_table



# ============================================================ #
# QC ----------------------------------------------------------
# ============================================================ #

qc_core_model_only <- purrr::map_dfr(lmer_predictors, function(pred) {
  
  fit_qc <- tryCatch(fit_core_lmer_model_only(pred), error = function(e) NULL)
  if (is.null(fit_qc)) {
    return(tibble::tibble(predictor = pred, status = "fit_qc failed"))
  }
  
  # Build the interaction term name, matching fit_lmer_interaction()'s
  # factor-vs-continuous branching
  analysis_frame_qc <- model.frame(fit_qc)
  is_factor_pred <- is.factor(analysis_frame_qc[[pred]])
  
  if (is_factor_pred) {
    lvs <- levels(analysis_frame_qc[[pred]])
    if (length(lvs) != 2) {
      return(tibble::tibble(predictor = pred, status = "non-binary factor, skipped"))
    }
    interaction_term <- paste0(time_var, ":", pred, lvs[2])
  } else {
    interaction_term <- paste0(time_var, ":", pred)
  }
  
  coef_qc <- broom.mixed::tidy(fit_qc, effects = "fixed", conf.int = TRUE) %>%
    dplyr::filter(term == interaction_term)
  
  ref_row <- all_lmer_results %>%
    dplyr::filter(biomarker_raw == pred, base_model == "Core")
  
  if (nrow(ref_row) == 0 || nrow(coef_qc) == 0) {
    return(tibble::tibble(predictor = pred, status = "term or reference row missing"))
  }
  
  tibble::tibble(
    predictor = pred,
    beta_match = isTRUE(all.equal(coef_qc$estimate, ref_row$beta, tolerance = 1e-8)),
    p_match    = isTRUE(all.equal(coef_qc$p.value, ref_row$p_value, tolerance = 1e-8)),
    nobs_match = nobs(fit_qc) == ref_row$n_obs,
    nsubj_match = dplyr::n_distinct(model.frame(fit_qc)$SSBS_ID) == ref_row$n_subjects,
    status = "compared"
  )
})

qc_nfl_model_only <- purrr::map_dfr(nfl_interaction_predictors, function(pred) {
  
  fit_qc <- tryCatch(fit_core_nfl_interaction_model_only(pred), error = function(e) NULL)
  if (is.null(fit_qc)) {
    return(tibble::tibble(predictor = pred, status = "fit_qc failed"))
  }
  
  # Build the 3-way interaction term name, matching fit_nfl_interaction()'s
  # factor-vs-continuous branching
  analysis_frame_qc <- model.frame(fit_qc)
  is_factor_bio <- is.factor(analysis_frame_qc[[pred]])
  
  if (is_factor_bio) {
    lvs <- levels(analysis_frame_qc[[pred]])
    if (length(lvs) != 2) {
      return(tibble::tibble(predictor = pred, status = "non-binary factor, skipped"))
    }
    interaction_term <- paste0(time_var, ":NfL_z:", pred, lvs[2])
  } else {
    interaction_term <- paste0(time_var, ":NfL_z:", pred)
  }
  
  coef_qc <- broom.mixed::tidy(fit_qc, effects = "fixed", conf.int = TRUE) %>%
    dplyr::filter(term == interaction_term)
  
  ref_row <- all_nfl_interaction_results %>%
    dplyr::filter(biomarker_raw == pred, base_model == "NfL+Core")
  
  if (nrow(ref_row) == 0 || nrow(coef_qc) == 0) {
    return(tibble::tibble(predictor = pred, status = "term or reference row missing"))
  }
  
  tibble::tibble(
    predictor = pred,
    beta_match = isTRUE(all.equal(coef_qc$estimate, ref_row$int_beta, tolerance = 1e-8)),
    p_match    = isTRUE(all.equal(coef_qc$p.value, ref_row$int_p_value, tolerance = 1e-8)),
    nobs_match = nobs(fit_qc) == ref_row$n_obs,
    nsubj_match = dplyr::n_distinct(model.frame(fit_qc)$SSBS_ID) == ref_row$n_subjects,
    status = "compared"
  )
})

# ok if empty, otherwise investigate any mismatches
qc_core_model_only %>%
  dplyr::filter(
    status != "compared" |
      !dplyr::coalesce(beta_match, FALSE) |
      !dplyr::coalesce(p_match, FALSE) |
      !dplyr::coalesce(nobs_match, FALSE) |
      !dplyr::coalesce(nsubj_match, FALSE)
  )

qc_nfl_model_only %>%
  dplyr::filter(
    status != "compared" |
      !dplyr::coalesce(beta_match, FALSE) |
      !dplyr::coalesce(p_match, FALSE) |
      !dplyr::coalesce(nobs_match, FALSE) |
      !dplyr::coalesce(nsubj_match, FALSE)
  )



fit_test_for_qc <- fit_core_nfl_interaction_model_only("Ab_status_bl")
performance::check_collinearity(fit_test_for_qc)




# Reviewer-Facing VIF Table for Response Letter (Reviewer 5, comment #4) ----
# Formats Table2_vif_table (Time x NfL x Biomarker three-way interaction
# models) as a grouped-header flextable for the response letter docx.
# Confirmed maximum VIF across all fixed-effect terms in all models: 1.92
# (MonthsFromFirstVisit main effect, APOE_e4_carrier model).

library(flextable)
library(officer)

vif_border <- officer::fp_border(color = "black", width = 1)

vif_footer_line <- paste(
  "Variance inflation factors (VIFs) were calculated for all fixed-effect terms in the fully adjusted Time × NfL × Biomarker interaction models shown in Table 2 using performance::check_collinearity(). All VIFs were <2.0, indicating no meaningful multicollinearity. The VIF range for remaining fixed-effect terms represents the minimum and maximum VIFs for time, the core clinical covariates (age, %VC, ΔFS, onset site, and REEC classification), and their interactions with time."
)

vif_table2_display <- Table2_vif_table %>%
  dplyr::mutate(
    dplyr::across(
      c(`NfL (main effect)`, `Biomarker (main effect)`,
        `Time × NfL`, `Time × Biomarker`, `NfL × Biomarker`,
        `Time × NfL × Biomarker`),
      ~ sprintf("%.2f", .x)
    ),
    biomarker = dplyr::if_else(
      biomarker %in% names(T1label),
      unname(unlist(T1label[biomarker])),
      biomarker
    )
  ) %>%
  dplyr::transmute(
    Biomarker = biomarker,
    NfL = `NfL (main effect)`,
    `Biomarker effect` = `Biomarker (main effect)`,
    `Time × NfL` = `Time × NfL`,
    `Time × Biomarker` = `Time × Biomarker`,
    `NfL × Biomarker` = `NfL × Biomarker`,
    `Time × NfL × Biomarker` = `Time × NfL × Biomarker`,
    `VIF range` = covariate_VIF_range
  )

vif_ft_table2 <- vif_table2_display %>%
  flextable::flextable() %>%
  flextable::add_header_row(
    values = c(
      "",
      "Main-effect VIFs",
      "Two-way interaction VIFs",
      "Three-way interaction VIF",
      "Remaining fixed-effect terms"
    ),
    colwidths = c(1, 2, 3, 1, 1)
  )%>%
  flextable::align(align = "center", part = "header") %>%
  flextable::align(align = "center", j = 2:8, part = "body") %>%
  flextable::align(align = "left", j = 1, part = "body") %>%
  flextable::bold(part = "header") %>%
  flextable::bold(j = 1, part = "body") %>%
  flextable::width(j = 1, width = 1.45) %>%
  flextable::width(j = 2:3, width = 0.75) %>%
  flextable::width(j = 4:6, width = 1.05) %>%
  flextable::width(j = 7, width = 1.35) %>%
  flextable::width(j = 8, width = 0.95) %>%
  flextable::border_remove() %>%
  flextable::hline_top(border = vif_border, part = "header") %>%
  flextable::hline_bottom(border = vif_border, part = "header") %>%
  flextable::vline(j = 6, border = officer::fp_border(width = 0.6), part = "all") %>%
  flextable::vline(j = 7, border = officer::fp_border(width = 0.6), part = "all") %>%
  flextable::add_footer_lines(values = vif_footer_line) %>%
  flextable::hline_top(border = vif_border, part = "footer") %>%
  flextable::font(fontname = "Times New Roman", part = "all") %>%
  flextable::fontsize(size = 9, part = "all") %>%
  flextable::padding(padding.top = 2, padding.bottom = 2, padding.left = 3, padding.right = 3, part = "all") %>%
  flextable::set_caption(
    caption = flextable::as_paragraph(
      flextable::as_chunk(
        "Table for Reviewer. Multicollinearity diagnostics for NfL × biomarker interaction models",
        props = officer::fp_text(font.size = 12, bold = TRUE, font.family = "Times New Roman")
      )
    ),
    word_stylename = "Table Caption",
    fp_p = officer::fp_par(text.align = "left", padding = 5)
  )

vif_ft_table2


performance::check_collinearity(
  fit_core_nfl_interaction_model_only("Ab_status_bl")
)



# Export ------------------------------------------------------------------

# flextable::save_as_docx(
#   vif_ft_table2,
#   path = "output/Tables/ReviewerTable_VIF_Table2.docx",
#   pr_section = officer::prop_section(
#     page_size = officer::page_size(orient = "landscape", width = 11, height = 8.5),
#     page_margins = officer::page_mar(bottom = 0.5, top = 0.5, left = 0.5, right = 0.5)
#   )
# )
