# ============================================================================
# NfL Interaction Models with Time*NfL*Biomarker 3-way Interactions
# Functions for analyzing NfL as effect modifier
# ============================================================================

library(lmerTest)
library(broom.mixed)
library(MuMIn)
library(dplyr)
library(tidyr)
library(purrr)
library(flextable)
library(officer)

# Formatting Utilities -----------------------------------------------------

fmt_p2 <- function(x) gtsummary::style_pvalue(x, digits = 2)
fmt_2d <- function(x) ifelse(is.na(x), NA, scales::number(x, accuracy = 0.01))
fmt_3d <- function(x) ifelse(is.na(x), NA, scales::number(x, accuracy = 0.001))
fmt_ci2 <- function(est, lcl, ucl) {
  ifelse(
    is.na(est) | is.na(lcl) | is.na(ucl),
    NA,
    sprintf("%.2f (%.2f–%.2f)", est, lcl, ucl)
  )
}

# Core Model Fitting Function ----------------------------------------------

# Fit model testing NfL as modifier of biomarker-decline association
#' Base model: Time*NfL + covariates
#' Full model: Time*NfL*Biomarker + covariates (3-way interaction)
fit_nfl_interaction <- function(biomarker_col, 
                                base_covariates = character(0),
                                interaction_covariates = character(0),
                                base_label, 
                                model_label,
                                data = df_vis,
                                cut_re = FALSE) {
  # Random effects structure
  # cut_re = TRUE: simplified random intercept/slope structure without covariance estimation
  # cut_re = FALSE: random intercept/slope structure with covariance estimation
  
  time_var <- "MonthsFromFirstVisit"
  nfl_var  <- "NfL_z"
  
  # Random effects structure
  if (cut_re) {
    re_term <- paste0("(", time_var, " || SSBS_ID)")
  } else {
    re_term <- paste0("(", time_var, " | SSBS_ID)")
  }
  
  # Check variables exist
  required_vars <- c("ALSFRSR_Total", "SSBS_ID", time_var, nfl_var, biomarker_col, 
                     base_covariates, interaction_covariates)
  missing <- setdiff(required_vars, names(data))
  if (length(missing) > 0) {
    stop("Missing variables in data: ", paste(missing, collapse = ", "))
  }
  
  # Build model formulas
  # Base model: Time*NfL + time*interaction_covariates + main_effect_covariates
  rhs_base_parts <- c(paste(time_var, "*", nfl_var))
  if (length(interaction_covariates) > 0) {
    interaction_terms <- paste(time_var, "*", interaction_covariates)
    rhs_base_parts <- c(rhs_base_parts, interaction_terms)
  }
  if (length(base_covariates) > 0) {
    rhs_base_parts <- c(rhs_base_parts, base_covariates)
  }
  rhs_base <- paste(rhs_base_parts, collapse = " + ")
  
  # Full model: Time*NfL*Biomarker + time*interaction_covariates + main_effect_covariates
  rhs_full_parts <- c(paste(time_var, "*", nfl_var, "*", biomarker_col))
  if (length(interaction_covariates) > 0) {
    interaction_terms <- paste(time_var, "*", interaction_covariates)
    rhs_full_parts <- c(rhs_full_parts, interaction_terms)
  }
  if (length(base_covariates) > 0) {
    rhs_full_parts <- c(rhs_full_parts, base_covariates)
  }
  rhs_full <- paste(rhs_full_parts, collapse = " + ")
  
  base_formula <- as.formula(paste("ALSFRSR_Total ~", rhs_base, "+", re_term))
  full_formula <- as.formula(paste("ALSFRSR_Total ~", rhs_full, "+", re_term))
  
  # Complete cases + at least 2 time points per subject
  analysis_frame <- data %>%
    dplyr::select(dplyr::all_of(required_vars)) %>%
    tidyr::drop_na() %>%
    dplyr::group_by(SSBS_ID) %>%
    dplyr::filter(dplyr::n_distinct(.data[[time_var]]) >= 2) %>%
    dplyr::ungroup()
  
  if (nrow(analysis_frame) == 0) return(NULL)
  
  ctrl <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
  
  # ML fits for AIC/BIC
  fit_base_ML <- lmerTest::lmer(base_formula, data = analysis_frame, REML = FALSE, control = ctrl)
  fit_full_ML <- lmerTest::lmer(full_formula, data = analysis_frame, REML = FALSE, control = ctrl)
  
  # REML fits for inference
  fit_base_REML <- lmerTest::lmer(base_formula, data = analysis_frame, REML = TRUE, control = ctrl)
  fit_full_REML <- lmerTest::lmer(full_formula, data = analysis_frame, REML = TRUE, control = ctrl)
  
  # Extract coefficients from REML fit
  coef_full <- broom.mixed::tidy(fit_full_REML, effects = "fixed", conf.int = TRUE)
  
  # Determine term names based on biomarker type
  is_factor_bio <- is.factor(analysis_frame[[biomarker_col]])
  
  if (is_factor_bio) {
    lvs <- levels(analysis_frame[[biomarker_col]])
    if (length(lvs) != 2) return(NULL)
    
    main_term <- paste0(time_var, ":", biomarker_col, lvs[2])
    interaction_term <- paste0(time_var, ":", nfl_var, ":", biomarker_col, lvs[2])
    nfl_term <- paste0(time_var, ":", nfl_var)
    
    display_label <- sprintf("%s (%s vs %s)", 
                             stringr::str_remove(biomarker_col, "_bl$|_z$"), 
                             lvs[2], lvs[1])
  } else {
    main_term <- paste0(time_var, ":", biomarker_col)
    interaction_term <- paste0(time_var, ":", nfl_var, ":", biomarker_col)
    nfl_term <- paste0(time_var, ":", nfl_var)
    
    display_label <- stringr::str_remove(biomarker_col, "_bl$|_z$")
  }
  
  # Extract each coefficient
  main_coef <- coef_full %>% dplyr::filter(term == main_term)
  interaction_coef <- coef_full %>% dplyr::filter(term == interaction_term)
  nfl_coef <- coef_full %>% dplyr::filter(term == nfl_term)
  
  # Handle missing terms
  if (nrow(main_coef) == 0) {
    main_coef <- tibble::tibble(estimate = NA_real_, conf.low = NA_real_, 
                                conf.high = NA_real_, p.value = NA_real_)
  }
  if (nrow(interaction_coef) == 0) {
    interaction_coef <- tibble::tibble(estimate = NA_real_, conf.low = NA_real_, 
                                       conf.high = NA_real_, p.value = NA_real_)
  }
  if (nrow(nfl_coef) == 0) {
    nfl_coef <- tibble::tibble(estimate = NA_real_, conf.low = NA_real_, 
                               conf.high = NA_real_, p.value = NA_real_)
  }
  
  # Model fit metrics
  aic_base_ML <- AIC(fit_base_ML)
  aic_full_ML <- AIC(fit_full_ML)
  
  r2_base <- suppressWarnings(MuMIn::r.squaredGLMM(fit_base_REML)[1, "R2m"])
  r2_full <- suppressWarnings(MuMIn::r.squaredGLMM(fit_full_REML)[1, "R2m"])
  
  tibble::tibble(
    biomarker_raw   = biomarker_col,
    biomarker_label = display_label,
    # Main effect (Time:Biomarker)
    main_beta       = main_coef$estimate,
    main_beta_lcl   = main_coef$conf.low,
    main_beta_ucl   = main_coef$conf.high,
    main_p_value    = main_coef$p.value,
    # Interaction effect (Time:NfL:Biomarker - 3-way interaction)
    int_beta        = interaction_coef$estimate,
    int_beta_lcl    = interaction_coef$conf.low,
    int_beta_ucl    = interaction_coef$conf.high,
    int_p_value     = interaction_coef$p.value,
    # NfL effect (Time:NfL)
    nfl_beta        = nfl_coef$estimate,
    nfl_beta_lcl    = nfl_coef$conf.low,
    nfl_beta_ucl    = nfl_coef$conf.high,
    nfl_p_value     = nfl_coef$p.value,
    # Model fit
    AIC_base        = aic_base_ML,
    AIC_full        = aic_full_ML,
    delta_AIC       = aic_full_ML - aic_base_ML,
    BIC_base        = BIC(fit_base_ML),
    BIC_full        = BIC(fit_full_ML),
    R2m_base        = r2_base,
    R2m_full        = r2_full,
    delta_R2m       = r2_full - r2_base,
    n_obs           = nrow(analysis_frame),
    n_subjects      = dplyr::n_distinct(analysis_frame$SSBS_ID),
    base_model      = base_label,
    model_label     = model_label,
    cut_re          = cut_re,
    singular_full   = lme4::isSingular(fit_full_REML)
  )
}

# Wrapper Functions --------------------------------------------------------

# Fit Age + NfL adjusted model
fit_age_nfl_interaction <- function(pred) {
  if (pred %in% c(age_nfl_main_effects, age_nfl_time_interactions, "NfL_z")) return(NULL)
  
  fit_nfl_interaction(
    biomarker_col = pred,
    base_covariates = age_nfl_main_effects,
    interaction_covariates = age_nfl_time_interactions,
    base_label = "Age+NfL",
    model_label = "Age+NfL adjusted"
  )
}

# Fit fully adjusted model with core covariates + NfL
fit_core_nfl_interaction <- function(pred) {
  # Avoid fitting for core covariates themselves
  if (pred %in% c(core_main_effects, core_time_interactions, "NfL_z")) return(NULL)
  
  fit_nfl_interaction(
    biomarker_col = pred,
    base_covariates = core_main_effects,
    interaction_covariates = core_time_interactions,
    base_label = "NfL+Core",
    model_label = "NfL+Core adjusted"
  )
}


# Table Generation Functions -----------------------------------------------

# Select specific columns from wide-format NfL interaction results
select_nfl_interaction_columns <- function(master_wide_tbl,
                                           blocks = c("Age+NfL", "NfL+Core"),
                                           metrics = c("MainBeta", "MainP", "InteractionBeta", 
                                                       "InteractionP", "NflBeta", "NflP",
                                                       "AIC", "DeltaAIC", "MarginalR2", "DeltaR2"),
                                           show_n_cols = c("subj", "obs", "both", "none")) {
  show_n_cols <- match.arg(show_n_cols)
  
  metric_map <- c(
    MainBeta        = "Main β (95% CI)",
    MainP           = "Main p-value",
    InteractionBeta = "Interaction β (95% CI)",
    InteractionP    = "Interaction p-value",
    NflBeta         = "NfL β (95% CI)",
    NflP            = "NfL p-value",
    AIC             = "AIC",
    DeltaAIC        = "ΔAIC",
    BIC             = "BIC",
    MarginalR2      = "Marginal R²",
    DeltaR2         = "ΔMarginal R²"
  )
  
  desired_cols <- "Biomarker"
  
  for (b in blocks) {
    cols_b <- paste(b, unname(metric_map[metrics]), sep = " | ")
    cols_b <- cols_b[cols_b %in% names(master_wide_tbl)]
    desired_cols <- c(desired_cols, cols_b)
    
    n_subj_col <- paste0(b, " | N (subj)")
    n_obs_col  <- paste0(b, " | N (obs)")
    
    if (show_n_cols %in% c("subj", "both") && n_subj_col %in% names(master_wide_tbl)) {
      desired_cols <- c(desired_cols, n_subj_col)
    }
    if (show_n_cols %in% c("obs", "both") && n_obs_col %in% names(master_wide_tbl)) {
      desired_cols <- c(desired_cols, n_obs_col)
    }
  }
  
  desired_cols <- desired_cols[desired_cols %in% names(master_wide_tbl)]
  
  master_wide_tbl %>%
    dplyr::select(dplyr::all_of(desired_cols))
}

# Build flextable from wide-format NfL interaction results
build_nfl_interaction_flextable <- function(master_wide_tbl,
                                            include_main_p = c("none", "all"),
                                            include_int_p = c("none", "all"),
                                            include_nfl = c("none", "all"),
                                            metrics_to_show = c("AIC", "DeltaAIC", "MarginalR2", "DeltaR2"),
                                            core_covariates = character(0),
                                            label_map = NULL) {
  include_main_p <- match.arg(include_main_p)
  include_int_p <- match.arg(include_int_p)
  include_nfl <- match.arg(include_nfl)
  
  if (!"Biomarker" %in% names(master_wide_tbl)) {
    if ("biomarker_label" %in% names(master_wide_tbl)) {
      master_wide_tbl <- dplyr::rename(master_wide_tbl, Biomarker = biomarker_label)
    } else {
      stop("build_nfl_interaction_flextable(): expected a 'Biomarker' column.")
    }
  }
  
  keep_cols <- names(master_wide_tbl)
  
  # Metric filtering
  if (!"AIC" %in% metrics_to_show) {
    keep_cols <- keep_cols[!grepl("\\| AIC$", keep_cols)]
  }
  if (!"DeltaAIC" %in% metrics_to_show) {
    keep_cols <- keep_cols[!grepl("\\| ΔAIC$", keep_cols)]
  }
  if (!"BIC" %in% metrics_to_show) {
    keep_cols <- keep_cols[!grepl("\\| BIC$", keep_cols)]
  }
  if (!"MarginalR2" %in% metrics_to_show) {
    keep_cols <- keep_cols[!grepl("\\| Marginal R²$", keep_cols)]
  }
  if (!"DeltaR2" %in% metrics_to_show) {
    keep_cols <- keep_cols[!grepl("\\| ΔMarginal R²$", keep_cols)]
  }
  
  # p-value filtering
  if (include_main_p == "none") {
    keep_cols <- keep_cols[!grepl("\\| Main p-value$", keep_cols)]
  }
  if (include_int_p == "none") {
    keep_cols <- keep_cols[!grepl("\\| Interaction p-value$", keep_cols)]
  }
  
  # NfL coefficient filtering
  if (include_nfl == "none") {
    keep_cols <- keep_cols[!grepl("\\| NfL β \\(95% CI\\)$|\\| NfL p-value$", keep_cols)]
  }
  
  master_use <- master_wide_tbl[, keep_cols, drop = FALSE]
  
  if (!is.null(label_map)) {
    transform_one <- function(b) {
      if (is.na(b) || !nzchar(b)) return(b)
      base <- sub(" \\(.*$", "", b)
      base <- sub("^log_", "", base)
      suffix <- if (grepl(" \\(", b)) sub("^[^\\(]+", "", b) else ""
      if (base %in% names(label_map)) {
        paste0(unname(unlist(label_map[base])), suffix)
      } else {
        b
      }
    }
    master_use$Biomarker <- vapply(master_use$Biomarker, transform_one, character(1))
  }
  
  # Build two-level header
  col_keys <- names(master_use)
  top_lab  <- character(length(col_keys))
  bot_lab  <- character(length(col_keys))
  
  top_lab[1] <- ""
  bot_lab[1] <- "Biomarker"
  
  if (length(col_keys) > 1) {
    for (i in 2:length(col_keys)) {
      parts <- strsplit(col_keys[i], " \\| ")[[1]]
      if (length(parts) == 2) {
        top_lab[i] <- parts[1]
        bot_lab[i] <- parts[2]
      } else {
        top_lab[i] <- col_keys[i]
        bot_lab[i] <- ""
      }
    }
  }
  
  header_df <- data.frame(
    col_keys = col_keys, 
    Block = top_lab, 
    Metric = bot_lab, 
    stringsAsFactors = FALSE
  )
  
  ft <- flextable::flextable(master_use, col_keys = col_keys)
  ft <- flextable::set_header_df(ft, mapping = header_df, key = "col_keys")
  ft <- flextable::merge_h(ft, part = "header")
  ft <- flextable::merge_v(ft, part = "header")
  ft <- flextable::align(ft, align = "center", part = "header")
  ft <- flextable::bold(ft, part = "header")
  ft <- flextable::autofit(ft)
  ft <- flextable::align(ft, align = "center", j = 2:ncol(master_use))
  ft <- flextable::bold(ft, j = 1)
  
  # Add borders
  brd <- officer::fp_border(color = "black", width = 1)
  ft  <- flextable::border_remove(ft)
  ft  <- flextable::hline(ft, i = 1, border = brd, part = "header")
  ft  <- flextable::hline(ft, i = 2, border = brd, part = "header")
  ft  <- flextable::hline_top(ft, border = brd, part = "header")
  ft  <- flextable::hline_bottom(ft, border = brd, part = "header")
  
  # Add footer notes
  abbrev_line <- "Aβ, amyloid-β; AIC, Akaike information criterion; β, regression coefficient; CI, confidence interval; ΔAIC, change in AIC; GFAP, glial fibrillary acidic protein; NfL, neurofilament light"
  model_line  <- "Linear mixed-effects models with random intercepts and slopes for time (patient as a random effect) were used, with ALSFRS-R total score as the dependent variable.
  NfL β represents the Time × NfL interaction, reflecting the association between baseline NfL (per 1-SD increase) and ALSFRS-R slope (points per month). 
  Interaction β represents the Time × NfL × Biomarker three-way interaction term, indicating how each biomarker modifies the association between NfL and ALSFRS-R slope."
                 #"Main β represents the Time:Biomarker interaction; Interaction β represents the Time:NfL:Biomarker 3-way interaction (NfL modifying biomarker effect); NfL β represents the Time:NfL interaction in the full model."
  base_line   <- "ΔAIC represents the AIC difference between the full model with biomarker-related terms, including the three-way interaction, and the corresponding base model without these terms (i.e., the Time × NfL model with the same covariates); values < −2 indicate meaningful improvement in model fit."
  
  footer_lines <- c(model_line, base_line)
  if (length(core_covariates) > 0) {
    core_line <- sprintf("All models were adjusted for the core clinical covariates listed below, including their interactions with Time. 
                         Core covariates: %s.", paste(core_covariates, collapse = ", "))
    footer_lines <- c(footer_lines, core_line, abbrev_line)
  }
  
  ft <- flextable::add_footer_lines(ft, values = footer_lines)
  ft <- flextable::hline_top(ft, border = brd, part = "footer")
  
  ft
}



# Sensitivity wrapper: Time x NfL x Age (APOE) ------------------------------------
# Age_init is excluded from interaction_covariates here so it is not
# caught by the guard clause, allowing it to be tested as the biomarker_col
fit_age_nfl_interaction_sens <- function(pred) {
  if (pred %in% c("NfL_z")) return(NULL)
  fit_nfl_interaction(
    biomarker_col          = pred,
    base_covariates        = character(0),
    interaction_covariates = character(0),  # Age_init not listed here
    base_label             = "Age+NfL",
    model_label            = "Age+NfL adjusted"
  )
}

fit_core_nfl_interaction_sens <- function(pred) {
  if (pred %in% c("NfL_z")) return(NULL)
  fit_nfl_interaction(
    biomarker_col          = pred,
    base_covariates        = character(0),          
    interaction_covariates = core_time_interactions, 
    base_label             = "NfL+Core",
    model_label            = "NfL+Core adjusted"
  )
}
