library(dplyr)
library(purrr)
library(stringr)
library(broom)
library(survival)
library(flextable)
library(ggplot2)
library(tidyr)
library(officer)

# Utilities -----------------------------------------------------------

# Format p-values with 2 decimals (vectorized)
fmt_p2 <- function(x) {
  # style_pvalue handles <0.001 nicely; digits=2 for typical display
  gtsummary::style_pvalue(x, digits = 2)
}

# Format numeric to 2 decimals (vectorized)
fmt_2d <- function(x) {
  ifelse(is.na(x), NA, scales::number(x, accuracy = 0.01))
}

# Format numeric to 2 decimals (vectorized)
fmt_3d <- function(x) {
  ifelse(is.na(x), NA, scales::number(x, accuracy = 0.001))
}

# Compose "a (b–c)" with 2 decimals, NA-safe
fmt_ci2 <- function(est, lcl, ucl) {
  ifelse(
    is.na(est) | is.na(lcl) | is.na(ucl),
    NA,
    sprintf("%.2f (%.2f–%.2f)", est, lcl, ucl)
  )
}

fmt_ci3 <- function(est, lcl, ucl) {
  ifelse(
    is.na(est) | is.na(lcl) | is.na(ucl),
    NA,
    sprintf("%.3f (%.3f–%.3f)", est, lcl, ucl)
  )
}



# Compute Harrell's C-index (point estimate) for a fitted Cox model
c_index_point <- function(cox_fit) {
  out <- tryCatch(
    as.numeric(survival::concordance(cox_fit)$concordance),
    error = function(e) NA_real_
  )
  out
}

# Robust extractor for Harrell's C-index and its SE/CI from a coxph fit.
# 1) Prefer summary(fit)$concordance (stable: c, se).
# 2) Fallback to fit$concordance (names often: concordance, std).
# 3) Clamp CI to [0, 1] just in case of tiny rounding issues.

c_index_stats <- function(cox_fit) {
  try_summary <- tryCatch({
    cc <- summary(cox_fit)$concordance
    c_est <- unname(cc[1])
    c_se  <- unname(cc[2])
    list(c = c_est, se = c_se)
  }, error = function(e) NULL)
  
  if (!is.null(try_summary) && is.finite(try_summary$c) && is.finite(try_summary$se)) {
    c_est <- try_summary$c
    c_se  <- try_summary$se
  } else {
    # Fallback path using the model slot if available
    cc2 <- tryCatch(cox_fit$concordance, error = function(e) NULL)
    if (!is.null(cc2)) {
      # Many survival versions name them "concordance" and "std"
      c_est <- suppressWarnings(as.numeric(cc2[["concordance"]]))
      c_se  <- suppressWarnings(as.numeric(cc2[["std"]]))
    } else {
      c_est <- NA_real_; c_se <- NA_real_
    }
  }
  
  if (is.na(c_est) || is.na(c_se)) {
    return(tibble::tibble(
      cindex = NA_real_, cindex_se = NA_real_,
      c_lcl = NA_real_, c_ucl = NA_real_
    ))
  }
  
  c_lcl <- max(0, c_est - 1.96 * c_se)
  c_ucl <- min(1, c_est + 1.96 * c_se)
  
  tibble::tibble(
    cindex = c_est,
    cindex_se = c_se,
    c_lcl = c_lcl,
    c_ucl = c_ucl
  )
}


fit_and_summarize <- function(predictor_col, base_formula, full_formula, 
                              base_label, model_label, required_vars = NULL) {
  
  # Base variables always needed
  base_vars <- c("surv_time_days", "surv_event", predictor_col)
  
  # Add model-specific required variables
  if (!is.null(required_vars)) {
    base_vars <- c(base_vars, required_vars)
  }
  
  # Remove duplicates and ensure all variables exist in the dataset
  base_vars <- unique(base_vars)
  base_vars <- base_vars[base_vars %in% names(als_surv_for_screen)]
  
  # Build analysis frame with only necessary variables
  analysis_frame <- als_surv_for_screen %>%
    dplyr::select(all_of(base_vars)) %>%
    dplyr::mutate(
      !!predictor_col := {
        x <- .data[[predictor_col]]
        if (is.character(x) || is.factor(x)) {
          if (all(c("negative","positive") %in% unique(as.character(x[!is.na(x)])))) {
            factor(x, levels = c("negative","positive"))
          } else {
            factor(x)
          }
        } else x
      }
    ) %>%
    tidyr::drop_na()
  
  if (nrow(analysis_frame) == 0) {
    return(NULL)
  }
  
  # Fit Cox models
  fit_base <- survival::coxph(base_formula, data = analysis_frame)
  fit_full <- survival::coxph(full_formula, data = analysis_frame)
  
  # Extract coefficients
  tidy_full <- broom::tidy(fit_full, exponentiate = TRUE, conf.int = TRUE)
  
  coef_row <- NULL; display_label <- predictor_col
  if (is.numeric(analysis_frame[[predictor_col]])) {
    coef_row <- dplyr::filter(tidy_full, term == predictor_col)
    display_label <- stringr::str_remove(predictor_col, "^log_")
  } else if (is.factor(analysis_frame[[predictor_col]])) {
    lv <- levels(analysis_frame[[predictor_col]])
    if (length(lv) == 2) {
      term_to_pull <- paste0(predictor_col, lv[2])
      coef_row <- dplyr::filter(tidy_full, term == term_to_pull)
      display_label <- sprintf("%s (%s vs %s)",
                               stringr::str_remove(predictor_col, "^log_"), lv[2], lv[1])
    } else {
      coef_row <- tibble::tibble(estimate = NA_real_, conf.low = NA_real_,
                                 conf.high = NA_real_, p.value = NA_real_)
      display_label <- sprintf("%s (multi-level)", predictor_col)
    }
  }
  # C-index
  ci_base <- c_index_stats(fit_base)
  ci_full <- c_index_stats(fit_full)
  
  # AIC calculation
  aic_base <- AIC(fit_base)
  aic_full <- AIC(fit_full)
  
  tibble::tibble(
    biomarker_raw   = predictor_col,
    biomarker_label = display_label,
    hr              = coef_row$estimate,
    hr_lcl          = coef_row$conf.low,
    hr_ucl          = coef_row$conf.high,
    p_value         = coef_row$p.value,
    cindex_base     = ci_base$cindex,
    cindex_base_lcl = ci_base$c_lcl,
    cindex_base_ucl = ci_base$c_ucl,
    cindex_full     = ci_full$cindex,
    cindex_full_lcl = ci_full$c_lcl,
    cindex_full_ucl = ci_full$c_ucl,
    delta_cindex    = ci_full$cindex - ci_base$cindex,
    aic_base        = aic_base,
    aic_full        = aic_full,
    delta_aic       = aic_full - aic_base,
    n_analytic      = nrow(analysis_frame),
    n_events        = sum(analysis_frame$surv_event),
    base_model      = base_label,
    model_label     = model_label
  )
}



# ───────────────────────────────────────────── Wrappers for 4 model types ---
fit_unadjusted <- function(pred) {
  base_formula <- as.formula("Surv(surv_time_days, surv_event) ~ 1")
  full_formula <- as.formula(paste("Surv(surv_time_days, surv_event) ~", pred))
  fit_and_summarize(pred, base_formula, full_formula, "Null", "Unadjusted",
                    required_vars = NULL)
}


fit_age_adj <- function(pred) {
  if (pred == "Age_init") return(NULL)  # skip self-adjustment
  base_formula <- as.formula("Surv(surv_time_days, surv_event) ~ Age_init")
  full_formula <- as.formula(paste("Surv(surv_time_days, surv_event) ~ Age_init +", pred))
  fit_and_summarize(pred, base_formula, full_formula, "Age", "Age adjusted",
                    required_vars = "Age_init")
}

fit_age_NfL_adj <- function(pred) {
  if (pred %in% c("Age_init","log_NfL_csf")) return(NULL)  # skip self-adjustment
  base_formula <- as.formula("Surv(surv_time_days, surv_event) ~ Age_init +  log_NfL_csf")
  full_formula <- as.formula(paste("Surv(surv_time_days, surv_event) ~ Age_init + log_NfL_csf +", pred))
  fit_and_summarize(pred, base_formula, full_formula, "Age+NfL", "Age + NfL adjusted",
                    required_vars = c("Age_init", "log_NfL_csf"))
}


# Define Core covariates once
core_covariates <- c("Age_init", "VC_Percent", "deltaFS", "REEC_definite", "OnsetSite")

fit_core <- function(pred) {
  base_formula <- as.formula(paste("Surv(surv_time_days, surv_event) ~", paste(core_covariates, collapse = " + ")))
  full_formula <- as.formula(paste("Surv(surv_time_days, surv_event) ~", paste(c(core_covariates, pred), collapse = " + ")))
  fit_and_summarize(pred, base_formula, full_formula, "Core", "Core adjusted",
                    required_vars = core_covariates)
}


# Helper function to pretty-print one subset
make_cox_ft <- function(results_tbl, base_model_filter, title_text, show_delta = TRUE) {
  x <- results_tbl %>%
    filter(base_model == base_model_filter) %>%
    mutate(
      hr_ci_disp  = fmt_ci2(hr, hr_lcl, hr_ucl),
      p_disp      = fmt_p2(p_value),
      c_full_disp = fmt_ci3(cindex_full, cindex_full_lcl, cindex_full_ucl),
      dc_disp     = fmt_3d(delta_cindex),
      n_evt_disp  = sprintf("%d / %d", n_events, n_analytic)
    )
  
  # Column set depends on whether ΔC-index is meaningful for this base model
  if (isTRUE(show_delta)) {
    x <- x %>%
      select(
        Biomarker = biomarker_label,
        `HR (95% CI)` = hr_ci_disp,
        `p-value` = p_disp,
        `C-index (95% CI)` = c_full_disp,
        `ΔC-index` = dc_disp,
        `Events / N` = n_evt_disp
      )
  } else {
    x <- x %>%
      select(
        Biomarker = biomarker_label,
        `HR (95% CI)` = hr_ci_disp,
        `p-value` = p_disp,
        `C-index (95% CI)` = c_full_disp,
        `Events / N` = n_evt_disp
      )
  }
  
  flextable(x) |>
    autofit() |>
    align(align = "center", j = 2:ncol(x)) |>   # <-- use ncol(x) instead of last_col()
    bold(j = 1) |>
    add_header_lines(title_text) |>
    add_footer_lines(values = c(
      "Cox proportional hazards models.",
      "HRs for continuous biomarkers are per 1 log unit (where applicable).",
      "C-index 95% CIs from summary.coxph standard errors.",
      if (isTRUE(show_delta)) "ΔC-index = C(full) – C(base)." else NULL
    ))
}

make_core_na_row <- function(pred) {
  tibble::tibble(
    biomarker_raw   = pred,
    biomarker_label = pred,
    hr              = NA_real_,
    hr_lcl          = NA_real_,
    hr_ucl          = NA_real_,
    p_value         = NA_real_,
    cindex_base     = NA_real_,
    cindex_base_lcl = NA_real_,
    cindex_base_ucl = NA_real_,
    cindex_full     = NA_real_,
    cindex_full_lcl = NA_real_,
    cindex_full_ucl = NA_real_,
    aic_base        = NA_real_,    
    aic_full        = NA_real_,    
    delta_aic       = NA_real_,    
    delta_cindex    = NA_real_,
    n_analytic      = NA_integer_,
    n_events        = NA_integer_,
    base_model      = "Core",
    model_label     = "Core adjusted"
  )
}

# ──────────────────────────────────────────────────────────────
# Wrapper: Core block that yields NA if pred is already in Core
#    Otherwise, fall back to the existing fit_core(pred)
# ──────────────────────────────────────────────────────────────
fit_core_or_na <- function(pred) {
  if (pred %in% core_covariates) {
    make_core_na_row(pred)
  } else {
    fit_core(pred)
  }
}
