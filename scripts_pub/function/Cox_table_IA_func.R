# Interaction -------------------------------------------------------------

# Core interaction fitting function for Cox models
fit_and_summarize_cox_interaction <- function(predictor_col, base_formula, interaction_formula,
                                              base_label, model_label, scale_predictors = FALSE) {
  
  # Extract needed variables from formulas
  vars_base <- all.vars(base_formula)
  vars_base <- vars_base[!vars_base %in% c("surv_time_days", "surv_event")]
  vars_interaction <- all.vars(interaction_formula)
  vars_interaction <- vars_interaction[!vars_interaction %in% c("surv_time_days", "surv_event")]
  needed <- unique(c("surv_time_days", "surv_event", vars_base, vars_interaction, predictor_col))
  
  # als_surv_for_screen is already scaled 
  analysis_frame <- als_surv_for_screen %>%
    dplyr::select(dplyr::all_of(needed))
  
  # Scale numeric predictors if requested
  if (scale_predictors) {
    numeric_cols <- names(analysis_frame)[sapply(analysis_frame, is.numeric)]
    # Don't scale outcome variables or time
    cols_to_scale <- setdiff(numeric_cols, c("surv_time_days", "surv_event"))
    
    for (col in cols_to_scale) {
      if (col %in% names(analysis_frame)) {
        analysis_frame[[col]] <- as.numeric(scale(analysis_frame[[col]]))
      }
    }
  }
  
  # Handle factor conversions for variables that exist
  if ("Sex" %in% names(analysis_frame)) {
    analysis_frame$Sex <- as.factor(analysis_frame$Sex)
  }
  if (predictor_col %in% names(analysis_frame) && !is.numeric(analysis_frame[[predictor_col]])) {
    x <- analysis_frame[[predictor_col]]
    # This block assumes binary Aβ status labels ("negative", "positive")
    if (is.character(x) || is.factor(x)) {
      if (all(c("negative","positive") %in% unique(as.character(x[!is.na(x)])))) {
        analysis_frame[[predictor_col]] <- factor(x, levels = c("negative","positive"))
      } else {
        analysis_frame[[predictor_col]] <- factor(x)
      }
    }
  }
  
  analysis_frame <- tidyr::drop_na(analysis_frame)
  
  if (nrow(analysis_frame) == 0) return(NULL)
  
  # Fit Cox models
  fit_base <- survival::coxph(base_formula, data = analysis_frame)
  fit_interaction <- survival::coxph(interaction_formula, data = analysis_frame)
  
  # Extract coefficients - get raw coefficients first
  tidy_raw <- broom::tidy(fit_interaction, exponentiate = FALSE, conf.int = TRUE)
  
  # Find interaction and main effect terms
  # Note: interaction terms are extracted by exact term names from broom::tidy().
  # If interaction HRs appear as NA in another package version, check the term names in tidy_raw.
  interaction_term <- paste0("log_NfL_csf:", predictor_col)
  main_term <- predictor_col
  display_label <- predictor_col
  main_coef_row <- NULL
  interaction_coef_row <- NULL
  
  if (is.numeric(analysis_frame[[predictor_col]])) {
    display_label <- stringr::str_remove(predictor_col, "^log_")
    main_coef_row <- dplyr::filter(tidy_raw, term == main_term)
    interaction_coef_row <- dplyr::filter(tidy_raw, term == interaction_term)
    
    # Exponentiate to get HR
    if (nrow(main_coef_row) > 0) {
      main_coef_row <- main_coef_row %>%
        dplyr::mutate(
          estimate = exp(estimate),
          conf.low = exp(conf.low),
          conf.high = exp(conf.high)
        )
    }
    if (nrow(interaction_coef_row) > 0) {
      interaction_coef_row <- interaction_coef_row %>%
        dplyr::mutate(
          estimate = exp(estimate),
          conf.low = exp(conf.low),
          conf.high = exp(conf.high)
        )
    }
  } else if (is.factor(analysis_frame[[predictor_col]])) {
    lv <- levels(analysis_frame[[predictor_col]])
    if (length(lv) == 2) {
      main_term <- paste0(predictor_col, lv[2])
      interaction_term <- paste0("log_NfL_csf:", predictor_col, lv[2])
      display_label <- sprintf("%s (%s vs %s)",
                               stringr::str_remove(predictor_col, "^log_"), lv[2], lv[1])
      main_coef_row <- dplyr::filter(tidy_raw, term == main_term)
      interaction_coef_row <- dplyr::filter(tidy_raw, term == interaction_term)
      
      # Exponentiate to get HR
      if (nrow(main_coef_row) > 0) {
        main_coef_row <- main_coef_row %>%
          dplyr::mutate(
            estimate = exp(estimate),
            conf.low = exp(conf.low),
            conf.high = exp(conf.high)
          )
      }
      if (nrow(interaction_coef_row) > 0) {
        interaction_coef_row <- interaction_coef_row %>%
          dplyr::mutate(
            estimate = exp(estimate),
            conf.low = exp(conf.low),
            conf.high = exp(conf.high)
          )
      }
    } else {
      main_coef_row <- tibble::tibble(estimate = NA_real_, conf.low = NA_real_,
                                      conf.high = NA_real_, p.value = NA_real_)
      interaction_coef_row <- tibble::tibble(estimate = NA_real_, conf.low = NA_real_,
                                             conf.high = NA_real_, p.value = NA_real_)
      display_label <- sprintf("%s (multi-level)", predictor_col)
    }
  }
  
  if (is.null(main_coef_row) || nrow(main_coef_row) == 0) {
    main_coef_row <- tibble::tibble(estimate = NA_real_, conf.low = NA_real_,
                                    conf.high = NA_real_, p.value = NA_real_)
  }
  if (is.null(interaction_coef_row) || nrow(interaction_coef_row) == 0) {
    interaction_coef_row <- tibble::tibble(estimate = NA_real_, conf.low = NA_real_,
                                           conf.high = NA_real_, p.value = NA_real_)
  }
  
  # Extract NfL coefficient from interaction model
  nfl_coef_row <- dplyr::filter(tidy_raw, term == "log_NfL_csf")
  if (nrow(nfl_coef_row) == 0) {
    nfl_coef_row <- tibble::tibble(estimate = NA_real_, conf.low = NA_real_,
                                   conf.high = NA_real_, p.value = NA_real_)
  } else {
    # Exponentiate NfL coefficient
    nfl_coef_row <- nfl_coef_row %>%
      dplyr::mutate(
        estimate = exp(estimate),
        conf.low = exp(conf.low),
        conf.high = exp(conf.high)
      )
  }
  
  # C-index
  ci_base <- c_index_stats(fit_base)
  ci_interaction <- c_index_stats(fit_interaction)
  
  # AIC and BIC
  aic_base <- AIC(fit_base)
  aic_interaction <- AIC(fit_interaction)
  bic_base <- BIC(fit_base)  
  bic_interaction <- BIC(fit_interaction)
  
  # Likelihood ratio test for interaction
  lr_test <- anova(fit_base, fit_interaction, test = "Chisq")
  lr_p_value <- if (nrow(lr_test) >= 2) lr_test$`Pr(>|Chi|)`[2] else NA_real_
  
  tibble::tibble(
    biomarker_raw = predictor_col,
    biomarker_label = display_label,
    main_hr = main_coef_row$estimate,
    main_hr_lcl = main_coef_row$conf.low,
    main_hr_ucl = main_coef_row$conf.high,
    main_p_value = main_coef_row$p.value,
    interaction_hr = interaction_coef_row$estimate,
    interaction_hr_lcl = interaction_coef_row$conf.low,
    interaction_hr_ucl = interaction_coef_row$conf.high,
    interaction_p_value = interaction_coef_row$p.value,
    nfl_hr = nfl_coef_row$estimate,
    nfl_hr_lcl = nfl_coef_row$conf.low,
    nfl_hr_ucl = nfl_coef_row$conf.high,
    nfl_p_value = nfl_coef_row$p.value,
    lr_test_p = lr_p_value,
    delta_aic = aic_interaction - aic_base,  # Changed to delta
    delta_bic = bic_interaction - bic_base,  # Changed to delta
    cindex_base = ci_base$cindex,
    cindex_interaction = ci_interaction$cindex,
    delta_cindex = ci_interaction$cindex - ci_base$cindex,
    n_analytic = nrow(analysis_frame),
    n_events = sum(analysis_frame$surv_event),
    base_model = base_label,
    model_label = model_label
  )
}

# Wrapper functions for different base models
fit_nfl_age_cox_interaction <- function(pred) {
  if (pred %in% c("log_NfL_csf", "Age_init")) return(NULL)
  
  base_formula <- as.formula("Surv(surv_time_days, surv_event) ~ Age_init + log_NfL_csf")
  interaction_formula <- as.formula(paste("Surv(surv_time_days, surv_event) ~ Age_init + log_NfL_csf +", 
                                          pred, "+ log_NfL_csf:", pred))
  
  fit_and_summarize_cox_interaction(
    pred, base_formula, interaction_formula,
    "Age+NfL", "NfL*Var + Age+NfL",
    scale_predictors = FALSE  # scale_predictors is kept FALSE because variables are already pre-scaled upstream (see `als_surv_for_screen`).
  )
}

fit_nfl_core_cox_interaction <- function(pred) {
  if (pred %in% c("log_NfL_csf", core_covariates)) return(NULL)
  
  base_terms <- c(core_covariates, "log_NfL_csf")
  base_formula_str <- paste("Surv(surv_time_days, surv_event) ~", paste(base_terms, collapse = " + "))
  interaction_formula_str <- paste(base_formula_str, "+", pred, "+ log_NfL_csf:", pred)
  
  base_formula <- as.formula(base_formula_str)
  interaction_formula <- as.formula(interaction_formula_str)
  
  fit_and_summarize_cox_interaction(
    pred, base_formula, interaction_formula,
    "NfL+Core", "NfL*Var + NfL+Core",
    scale_predictors = FALSE  # scale_predictors is kept FALSE because variables are already pre-scaled upstream (see `als_surv_for_screen`).
  )
}

# NA row maker for core variables
make_core_na_row_cox_interaction <- function(pred) {
  tibble::tibble(
    biomarker_raw = pred,
    biomarker_label = pred,
    main_hr = NA_real_,
    main_hr_lcl = NA_real_,
    main_hr_ucl = NA_real_,
    main_p_value = NA_real_,
    interaction_hr = NA_real_,
    interaction_hr_lcl = NA_real_,
    interaction_hr_ucl = NA_real_,
    interaction_p_value = NA_real_,
    nfl_hr = NA_real_,
    nfl_hr_lcl = NA_real_,
    nfl_hr_ucl = NA_real_,
    nfl_p_value = NA_real_,
    lr_test_p = NA_real_,
    delta_aic = NA_real_,
    delta_bic = NA_real_,
    cindex_base = NA_real_,
    cindex_interaction = NA_real_,
    delta_cindex = NA_real_,
    n_analytic = NA_integer_,
    n_events = NA_integer_,
    base_model = "NfL+Core",
    model_label = "NfL*Var + NfL+Core"
  )
}

fit_nfl_core_or_na_cox_interaction <- function(pred) {
  if (pred %in% c("log_NfL_csf", core_covariates)) {
    make_core_na_row_cox_interaction(pred)
  } else {
    fit_nfl_core_cox_interaction(pred)
  }
}



# Helper function to cap extreme HR values
cap_hr <- function(x, max_val = 1e6) {
  ifelse(is.na(x) | !is.finite(x) | x > max_val, NA_real_, x)
}


# Column selector for Cox interaction table -----------------------------

select_cox_interaction_columns <- function(master_wide_tbl,
                                           blocks = c("Age+NfL","NfL+Core"),
                                           metrics = c("MainHR","MainP","InteractionHR","InteractionP",
                                                       "NfLHR","NfLP","LRTestP","DeltaAIC","DeltaBIC","Cindex","DeltaC"),
                                           include_N = c("none","all")) {
  include_N <- match.arg(include_N)
  
  id_col <- dplyr::case_when(
    "Biomarker" %in% names(master_wide_tbl) ~ "Biomarker",
    "biomarker_label" %in% names(master_wide_tbl) ~ "biomarker_label",
    TRUE ~ NA_character_
  )
  if (is.na(id_col)) stop("select_cox_interaction_columns(): table must have 'Biomarker' or 'biomarker_label'.")
  
  metric_map <- c(
    MainHR = "Main HR (95% CI)",
    MainP = "Main p-value",
    InteractionHR = "Interaction HR (95% CI)",
    InteractionP = "Interaction p-value",
    NfLHR = "NfL HR (95% CI)",
    NfLP = "NfL p-value",
    LRTestP = "LR test p-value",
    DeltaAIC = "ΔAIC",
    DeltaBIC = "ΔBIC",
    Cindex = "C-index",
    DeltaC = "ΔC-index"
  )
  
  desired_cols <- c(id_col)
  for (b in blocks) {
    cols_b <- paste(b, unname(metric_map[metrics]), sep = " | ")
    cols_b <- cols_b[cols_b %in% names(master_wide_tbl)]
    desired_cols <- c(desired_cols, cols_b)
    
    if (include_N == "all") {
      ncol_name <- paste0(b, " | Events/N")   # <-- changed to Events/N
      if (ncol_name %in% names(master_wide_tbl)) {
        desired_cols <- c(desired_cols, ncol_name)
      }
    }
  }
  
  out <- dplyr::select(master_wide_tbl, dplyr::all_of(desired_cols))
  if (id_col == "biomarker_label") out <- dplyr::rename(out, Biomarker = biomarker_label)
  out
}

# Add N columns
add_N_columns_cox_interaction <- function(master_wide_tbl, results_tbl, where = c("all","none")) {
  where <- match.arg(where)
  need_blocks <- switch(where, all = c("Age+NfL", "NfL+Core"), none = character(0))
  out <- master_wide_tbl
  
  for (blk in need_blocks) {
    ncol_name <- paste0(blk, " | Events/N")
    n_tbl <- results_tbl %>%
      dplyr::filter(base_model == blk) %>%
      dplyr::distinct(biomarker_raw, n_events, n_analytic) %>%
      dplyr::mutate(!!ncol_name := sprintf("%d/%d", n_events, n_analytic)) %>%
      dplyr::select(biomarker_raw, !!ncol_name)
    out <- dplyr::left_join(out, n_tbl, by = "biomarker_raw")
  }
  out
}


# Build flextable for Cox interactions ----------------------------------------

# --- PATCH 1: add include_events arg + handling ---
build_cox_interaction_flextable <- function(master_wide_tbl,
                                            include_p = c("none","all"),
                                            include_nfl = c("none","all"),
                                            include_events = c("none","all"),   # <-- added
                                            metrics_to_show = c("AIC","BIC","Cindex","DeltaC"),
                                            core_covariates = c("Age_init","VC_Percent",
                                                                "deltaFS","REEC_definite","OnsetSite"),
                                            label_map = NULL,
                                            blocks = c("Age+NfL","NfL+Core")) {
  include_p <- match.arg(include_p)
  include_nfl <- match.arg(include_nfl)
  include_events <- match.arg(include_events)   # <-- added
  
  if (!"Biomarker" %in% names(master_wide_tbl)) {
    if ("biomarker_label" %in% names(master_wide_tbl)) {
      master_wide_tbl <- dplyr::rename(master_wide_tbl, Biomarker = biomarker_label)
    } else {
      stop("build_cox_interaction_flextable(): expected a 'Biomarker' column.")
    }
  }
  
  # keep only selected block columns
  escaped_blocks <- gsub("\\+", "\\\\+", blocks)
  pattern <- paste(paste0("^", escaped_blocks, " \\|"), collapse = "|")
  keep_cols <- c("Biomarker", names(master_wide_tbl)[grepl(pattern, names(master_wide_tbl))])
  
  # filter metrics
  if (!"DeltaAIC" %in% metrics_to_show) keep_cols <- keep_cols[!grepl("\\| ΔAIC$", keep_cols)]
  if (!"DeltaBIC" %in% metrics_to_show) keep_cols <- keep_cols[!grepl("\\| ΔBIC$", keep_cols)]
  if (!"Cindex"   %in% metrics_to_show) keep_cols <- keep_cols[!grepl("\\| C-index$", keep_cols)]
  if (!"DeltaC"   %in% metrics_to_show) keep_cols <- keep_cols[!grepl("\\| ΔC-index$", keep_cols)]
  
  if (include_p == "none") {
    keep_cols <- keep_cols[!grepl("\\| Main p-value$|\\| Interaction p-value$|\\| NfL p-value$|\\| LR test p-value$", keep_cols)]
  }
  
  # Events/N on or off
  if (include_events == "none") {
    keep_cols <- keep_cols[!grepl("\\| Events/N$", keep_cols)]
  }
  
  if (include_nfl == "none") {
    keep_cols <- keep_cols[!grepl("\\| NfL HR \\(95% CI\\)$|\\| NfL p-value$", keep_cols)]
  }
  
  master_use <- master_wide_tbl[, keep_cols[keep_cols %in% names(master_wide_tbl)], drop = FALSE]
  
  if (!is.null(label_map)) master_use <- apply_labels(master_use, label_map)
  
  # two-level header
  col_keys <- names(master_use)
  top_lab <- character(length(col_keys)); bot_lab <- character(length(col_keys))
  top_lab[1] <- ""; bot_lab[1] <- "Variable"
  if (length(col_keys) > 1) {
    top_lab[-1] <- sub("^(.*?) \\| .*", "\\1", col_keys[-1])
    bot_lab[-1] <- sub("^.*? \\| (.*)$", "\\1", col_keys[-1])
  }
  header_df <- data.frame(col_keys = col_keys, Block = top_lab, Metric = bot_lab, stringsAsFactors = FALSE)
  
  ft <- flextable::flextable(master_use, col_keys = col_keys)
  ft <- flextable::set_header_df(ft, mapping = header_df, key = "col_keys")
  ft <- flextable::merge_h(ft, part = "header")
  ft <- flextable::merge_v(ft, part = "header")
  ft <- flextable::align(ft, align = "center", part = "header")
  ft <- flextable::bold(ft, part = "header")
  ft <- flextable::autofit(ft)
  ft <- flextable::align(ft, align = "center", j = 2:ncol(master_use))
  ft <- flextable::bold(ft, j = 1)
  
  brd <- officer::fp_border(color = "black", width = 1)
  ft <- flextable::border_remove(ft)
  ft <- flextable::hline(ft, i = 1, border = brd, part = "header")
  ft <- flextable::hline(ft, i = 2, border = brd, part = "header")
  ft <- flextable::hline_top(ft, border = brd, part = "header")
  ft <- flextable::hline_bottom(ft, border = brd, part = "header")
  
  abbrev_line <- "Abbreviations: HR, hazard ratio; CI, confidence interval; ΔAIC, delta Akaike information criterion; ΔBIC, delta Bayesian information criterion; C-index, Harrell's concordance index; ΔC-index, delta C-index; NfL, neurofilament light; GFAP, glial fibrillary acidic protein."
  interaction_line <- "NfL HR represents the NfL coefficient in the interaction model, and the Interaction HR represents the coefficient for the NfL × biomarker interaction term. All continuous biomarkers were log-transformed and standardized prior to modeling."
  delta_line <- "ΔAIC is defined as AIC(interaction model) − AIC(base model without the biomarker or its interaction term); negative values indicate improved model fit of the interaction model."
  #lr_line <- "LR test p-value compares interaction model vs base model without interaction."
  core_line <- sprintf("Core model includes: %s, NfL.", paste(core_covariates, collapse = ", "))
  
  ft <- flextable::add_footer_lines(ft, values = c(interaction_line, delta_line, core_line, abbrev_line))
  ft <- flextable::hline_top(ft, border = brd, part = "footer")
  
  ft
}