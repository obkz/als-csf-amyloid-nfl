# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
# Linear Mixed-Effects Models with Time*Biomarker Interactions
# Functions for analyzing longitudinal ALSFRS-R data
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #

library(lmerTest)
library(broom.mixed)
library(MuMIn)
library(dplyr)
library(tidyr)
library(purrr)
library(flextable)
library(officer)

# Data Preparation Functions -----------------------------------------------

# Build analysis dataset with baseline biomarkers
build_df_vis_baseline <- function(visitdf, csfdf, NotValidSample = character(),
                                  keep_multi_visit_only = TRUE,
                                  min_baseline_ALSFRS = 20) {
  baseline_labs <- csfdf %>%
    filter(ALS_label == "1") %>%
    distinct(SSBS_ID, .keep_all = TRUE) %>%
    transmute(
      SSBS_ID,
      # Core clinical features
      Age_init       = Age_init,
      VC_Percent     = VC_Percent,
      BMI_init       = BMI_init,
      deltaFS        = deltaFS,
      # Biomarkers
      Cr_bl          = Cr,
      CK_bl          = CK,
      Alb_bl         = Alb,
      CysC_bl        = cysC,
      CrCys_ratio_bl = Cr_CysC_ratio,
      NfL_bl         = NfL_csf_pgml,
      GFAP_bl        = GFAP_csf_pgml,
      pTau181_bl     = pTau181_csf,
      pTau217_bl     = pTau217_csf,
      Ab4240_bl      = Ab42_40_csf_bridged,
      Ab3840_bl      = Ab38_40_csf,
      Ab38_bl        = Ab38_csf,
      Ab40_bl        = Ab40_csf,
      Ab42_bl        = Ab42_csf,
      Ab_status_bl   = Ab_status
    ) |> 
    mutate(
      # Core clinical features scaled as z-scores for interpretability of interaction terms
      Age_init      = as.numeric(scale(Age_init)),
      VC_Percent    = as.numeric(scale(VC_Percent)),
      BMI_init      = as.numeric(scale(BMI_init)),
      deltaFS       = as.numeric(scale(deltaFS)),
      # Biomarkers scaled as z-scores for interpretability of interaction terms
      NfL_z         = as.numeric(scale(NfL_bl)),
      GFAP_z        = as.numeric(scale(GFAP_bl)),
      pTau181_z     = as.numeric(scale(pTau181_bl)),
      pTau217_z     = as.numeric(scale(pTau217_bl)),
      Ab4240_z      = as.numeric(scale(Ab4240_bl)),
      Ab3840_z      = as.numeric(scale(Ab3840_bl)),
      Ab38_z        = as.numeric(scale(Ab38_bl)),
      Ab40_z        = as.numeric(scale(Ab40_bl)),
      Ab42_z        = as.numeric(scale(Ab42_bl)),
      Cr_z          = as.numeric(scale(Cr_bl)),
      CK_z          = as.numeric(scale(CK_bl)),
      Alb_z         = as.numeric(scale(Alb_bl)),
      CrCys_ratio_z = as.numeric(scale(CrCys_ratio_bl))
    )
  
  df_vis <- visitdf %>%
    filter(!SampleID %in% NotValidSample) %>%
    filter(!is.na(ALSFRSR_Total) & ALS_label == "1" & ALSFRS_init >= min_baseline_ALSFRS) %>%
    select(
      -any_of(c(
                "Age_init","VC_Percent","BMI_init","deltaFS",
                 "Cr","CK","Alb","cysC","Cr_CysC_ratio",
                "NfL_csf_pgml","GFAP_csf_pgml",
                "pTau181_csf","pTau217_csf","Ab42_40_csf_bridged",
                "Ab_status","Ab38_40_csf"))) %>%
    left_join(ever_flags, by = "SSBS_ID") %>%
    mutate(
      MonthsFromFirstVisit = DaysFromFirstVisit_months,
      Sex = as.factor(Sex)) %>%
    left_join(baseline_labs, by = "SSBS_ID") 
  
  # Sanity check: baseline invariance within subject
  chk <- df_vis %>%
    group_by(SSBS_ID) %>%
    summarise(
      n_NfL = n_distinct(NfL_bl, na.rm = TRUE),
      n_Ab  = n_distinct(Ab4240_bl, na.rm = TRUE),
      .groups = "drop"
    )
  stopifnot(all(chk$n_NfL <= 1, na.rm = TRUE))
  stopifnot(all(chk$n_Ab  <= 1, na.rm = TRUE))
  
  if (keep_multi_visit_only) {
    vc <- df_vis %>% count(SSBS_ID, name = "valid_visit_n")
    ids <- vc %>% filter(valid_visit_n >= 2) %>% pull(SSBS_ID)
    df_vis <- df_vis %>% filter(SSBS_ID %in% ids)
  }
  
  df_vis
}

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

# Apply custom labels to biomarker names
apply_labels <- function(df, label_map) {
  if (!"Biomarker" %in% names(df)) {
    if ("biomarker_label" %in% names(df)) {
      df <- dplyr::rename(df, Biomarker = biomarker_label)
    } else {
      stop("apply_labels(): expected a 'Biomarker' or 'biomarker_label' column.")
    }
  }
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
  df$Biomarker <- vapply(df$Biomarker, transform_one, character(1))
  df
}

# Core Model Fitting Function ----------------------------------------------

# Fit linear mixed effects model with Time*Biomarker interaction
fit_lmer_interaction <- function(predictor_col, 
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
  
  if (cut_re) {
    re_term <- paste0("(", time_var, " || SSBS_ID)")
  } else {
    re_term <- paste0("(", time_var, " | SSBS_ID)")
  }
  
  # Build model formulas
  # Base model: time + main effect covariates + RE
  # Full model: time * predictor + time * interaction_covariates + main effect covariates + RE
  
  rhs_base_parts <- c(time_var)
  if (length(interaction_covariates) > 0) {
    interaction_terms <- paste(time_var, "*", interaction_covariates)
    rhs_base_parts <- c(rhs_base_parts, interaction_terms)
  }
  if (length(base_covariates) > 0) {
    rhs_base_parts <- c(rhs_base_parts, base_covariates)
  }
  rhs_base <- paste(rhs_base_parts, collapse = " + ")
  
  rhs_full_parts <- c(paste(time_var, "*", predictor_col))
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
  
  # Prepare analysis frame: complete cases with ≥2 timepoints per subject
  model_vars <- c("ALSFRSR_Total", "SSBS_ID", time_var, predictor_col, 
                  base_covariates, interaction_covariates)
  analysis_frame <- data %>%
    dplyr::select(dplyr::all_of(model_vars)) %>%
    tidyr::drop_na() %>%
    dplyr::group_by(SSBS_ID) %>%
    dplyr::filter(dplyr::n_distinct(.data[[time_var]]) >= 2) %>%
    dplyr::ungroup()
  if (nrow(analysis_frame) == 0) return(NULL)
  
  ctrl <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
  
  # Fit ML models for AIC/BIC comparison
  fit_base_ML <- lmerTest::lmer(base_formula, data = analysis_frame, REML = FALSE, control = ctrl)
  fit_full_ML <- lmerTest::lmer(full_formula, data = analysis_frame, REML = FALSE, control = ctrl)
  
  # Fit REML models for inference (beta, CI, p-value, R2)
  fit_base_REML <- lmerTest::lmer(base_formula, data = analysis_frame, REML = TRUE, control = ctrl)
  fit_full_REML <- lmerTest::lmer(full_formula, data = analysis_frame, REML = TRUE, control = ctrl)
  
  # Extract interaction term
  is_factor_pred <- is.factor(analysis_frame[[predictor_col]])
  if (is_factor_pred) {
    lvs <- levels(analysis_frame[[predictor_col]])
    if (length(lvs) != 2) return(NULL)
    interaction_term <- paste0(time_var, ":", predictor_col, lvs[2])
    display_label <- sprintf("%s (%s vs %s)",
                             stringr::str_remove(predictor_col, "_bl$|_z$"),
                             lvs[2], lvs[1])
  } else {
    interaction_term <- paste0(time_var, ":", predictor_col)
    display_label <- stringr::str_remove(predictor_col, "_bl$|_z$")
  }
  
  coef_tbl <- broom.mixed::tidy(fit_full_REML, effects = "fixed", conf.int = TRUE) %>%
    dplyr::filter(term == interaction_term)
  if (nrow(coef_tbl) == 0) return(NULL)
  
  aic_base_ML <- AIC(fit_base_ML)
  aic_full_ML <- AIC(fit_full_ML)
  
  r2_base <- suppressWarnings(MuMIn::r.squaredGLMM(fit_base_REML)[1, "R2m"])
  r2_full <- suppressWarnings(MuMIn::r.squaredGLMM(fit_full_REML)[1, "R2m"])
  
  tibble::tibble(
    biomarker_raw   = predictor_col,
    biomarker_label = display_label,
    beta            = coef_tbl$estimate,
    beta_lcl        = coef_tbl$conf.low,
    beta_ucl        = coef_tbl$conf.high,
    p_value         = coef_tbl$p.value,
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

# Wrapper Functions for Standard Adjustment Levels -------------------------

# Fit unadjusted model (time*biomarker only)
fit_unadjusted_lmer <- function(pred, interaction_covariates = character(0)) {
  fit_lmer_interaction(
    predictor_col = pred,
    base_covariates = character(0),
    interaction_covariates = character(0),  # Always empty for unadjusted
    base_label = "Null",
    model_label = "Unadjusted"
  )
}

# Fit age-adjusted model
fit_age_lmer <- function(pred, interaction_covariates = character(0)) {
  if (pred == "Age_init") return(NULL)
  main_effects <- setdiff("Age_init", interaction_covariates)
  fit_lmer_interaction(
    predictor_col = pred,
    base_covariates = main_effects,
    interaction_covariates = interaction_covariates,
    base_label = "Age",
    model_label = "Age adjusted"
  )
}

fit_age_nfl_lmer <- function(pred, interaction_covariates = character(0)) {
  if (pred %in% c("Age_init", "NfL_z")) return(NULL)
  
  interaction_covariates_use <- unique(c(interaction_covariates, "NfL_z"))
  main_effects <- setdiff(c("Age_init", "NfL_z"), interaction_covariates_use)
  
  fit_lmer_interaction(
    predictor_col = pred,
    base_covariates = main_effects,
    interaction_covariates = interaction_covariates_use,
    base_label = "Age+NfL",
    model_label = "Age+NfL adjusted"
  )
}

# Fit fully adjusted model with core covariates
fit_core_lmer <- function(pred, core_covariates, interaction_covariates = character(0)) {
  if (pred %in% core_covariates) return(NULL)
  main_effects <- setdiff(core_covariates, interaction_covariates)
  fit_lmer_interaction(
    predictor_col = pred,
    base_covariates = main_effects,
    interaction_covariates = interaction_covariates,
    base_label = "Core",
    model_label = "Core adjusted"
  )
}

# Table Generation Functions -----------------------------------------------

# Select specific columns from wide-format results table
select_lmer_columns <- function(master_wide_tbl,
                                blocks = c("Unadj", "Age", "Age+NfL", "Core"),
                                metrics = c("Beta", "p", "AIC", "DeltaAIC", "BIC", "MarginalR2", "DeltaR2"),
                                show_n_cols = c("subj", "obs", "both", "none")) {
  show_n_cols <- match.arg(show_n_cols)
  
  metric_map <- c(
    Beta      = "β (95% CI)",
    p         = "p-value",
    AIC       = "AIC",
    DeltaAIC  = "ΔAIC",
    BIC       = "BIC",
    MarginalR2 = "Marginal R²",
    DeltaR2   = "ΔMarginal R²"
  )
  
  desired_cols <- "Biomarker"
  
  for (b in blocks) {
    # Exclude delta metrics for unadjusted models
    metrics_this_block <- if (b == "Unadj") {
      setdiff(metrics, c("DeltaR2", "DeltaAIC"))
    } else {
      metrics
    }
    
    cols_b <- paste(b, unname(metric_map[metrics_this_block]), sep = " | ")
    cols_b <- cols_b[cols_b %in% names(master_wide_tbl)]
    desired_cols <- c(desired_cols, cols_b)
    
    # Add sample size columns
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

# Build flextable from wide-format results
build_lmer_flextable <- function(master_wide_tbl,
                                 include_p = c("none", "all"),
                                 metrics_to_show = c("AIC", "DeltaAIC", "BIC", "MarginalR2", "DeltaR2"),
                                 core_covariates = character(0),
                                 label_map = NULL) {
  include_p <- match.arg(include_p)
  
  if (!"Biomarker" %in% names(master_wide_tbl)) {
    if ("biomarker_label" %in% names(master_wide_tbl)) {
      master_wide_tbl <- dplyr::rename(master_wide_tbl, Biomarker = biomarker_label)
    } else {
      stop("build_lmer_flextable(): expected a 'Biomarker' column.")
    }
  }
  
  keep_cols <- names(master_wide_tbl)
  
  # Filter columns based on requested metrics
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
  if (include_p == "none") {
    keep_cols <- keep_cols[!grepl("\\| p-value$", keep_cols)]
  }
  
  master_use <- master_wide_tbl[, keep_cols, drop = FALSE]
  
  if (!is.null(label_map)) {
    master_use <- apply_labels(master_use, label_map)
  }
  
  # Build two-level header
  col_keys <- names(master_use)
  top_lab  <- character(length(col_keys))
  bot_lab  <- character(length(col_keys))
  
  top_lab[1] <- ""
  bot_lab[1] <- "Variables" #Variables/ biomarker
  
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
  model_line <- paste(
    "Linear mixed-effects models with random intercepts and slopes for time (patient as a random effect) were used, with ALSFRS-R total score as the dependent variable.",
    "β values represent the interaction between time and each variable X (per 1-SD increase for continuous variables).",
    "ΔAIC was defined as AIC [full model with time × X] − AIC [corresponding base model with the same covariates but without X]; values < −2 indicate a meaningful improvement in model fit.",
    "In all adjusted models, the variable of interest (X) and the covariates were specified with interactions with time; in unadjusted models, only the time × X term was included.",
    "N differs across models due to complete-case analysis.",
    sep = " "
  )
  footer_lines <- model_line
  if (length(core_covariates) > 0) {
    core_line <- sprintf("Core clinical features included as covariates in the fully adjusted (“Core”) models were: %s.", paste(core_covariates, collapse = ", "))
    footer_lines <- c(footer_lines, core_line)
  }
  
  ft <- flextable::add_footer_lines(ft, values = c(footer_lines,abbrev_line))
  ft <- flextable::hline_top(ft, border = brd, part = "footer")
  
  ft
}