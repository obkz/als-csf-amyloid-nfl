# LME Interaction - Visualization functions for Time*NfL*Biomarker models
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(MuMIn)
library(ggeffects)
library(ggplot2)
library(stringr)

# Model Fitting Function ---------------------------------------------------

# Fit 3-way interaction model with optional categorization and flexible covariate interactions
fit_nfl_modifier_interaction <- function(df_vis,
                                         modifier_var,
                                         modifier_type = c("continuous", "factor"),
                                         categorize_nfl = FALSE,
                                         nfl_categorize_method = c("quartiles", "tertiles", "median_split"),
                                         categorize_modifier = FALSE,
                                         modifier_categorize_method = c("quartiles", "tertiles", "median_split"),
                                         base_covariates = c("Age_init", "VC_Percent", "deltaFS", 
                                                             "REEC_definite", "OnsetSite"),
                                         interaction_covariates = character(0),
                                         cut_re = FALSE) {
  modifier_type <- match.arg(modifier_type)
  nfl_categorize_method <- match.arg(nfl_categorize_method)
  modifier_categorize_method <- match.arg(modifier_categorize_method)
  
  # Check required variables
  stopifnot("NfL_z" %in% names(df_vis))
  if (!(modifier_var %in% names(df_vis))) {
    stop("modifier_var not found in data: ", modifier_var)
  }
  
  df_mod <- df_vis
  nfl_var <- "NfL_z"
  
  # Categorize NfL if requested
  if (categorize_nfl) {
    nfl_cat_var_name <- "NfL_z_cat"
    
    if (nfl_categorize_method == "quartiles") {
      q <- quantile(df_mod$NfL_z, probs = c(0.25, 0.75), na.rm = TRUE)
      df_mod <- df_mod %>%
        mutate(NfL_z_cat = cut(NfL_z,
                               breaks = c(-Inf, q[1], q[2], Inf),
                               labels = c("Q1", "Q2-Q3", "Q4"),
                               include.lowest = TRUE))
    } else if (nfl_categorize_method == "tertiles") {
      t <- quantile(df_mod$NfL_z, probs = c(1/3, 2/3), na.rm = TRUE)
      df_mod <- df_mod %>%
        mutate(NfL_z_cat = cut(NfL_z,
                               breaks = c(-Inf, t[1], t[2], Inf),
                               labels = c("Low", "Middle", "High"),
                               include.lowest = TRUE))
    } else if (nfl_categorize_method == "median_split") {
      med <- median(df_mod$NfL_z, na.rm = TRUE)
      df_mod <- df_mod %>%
        mutate(NfL_z_cat = cut(NfL_z,
                               breaks = c(-Inf, med, Inf),
                               labels = c("Below median", "Above median"),
                               include.lowest = TRUE))
    }
    
    nfl_var <- nfl_cat_var_name
    df_mod <- df_mod %>% mutate(!!nfl_var := as.factor(.data[[nfl_var]]))
  }
  
  # Categorize modifier if requested and it's continuous
  if (modifier_type == "continuous" && categorize_modifier) {
    cat_var_name <- paste0(modifier_var, "_cat")
    
    if (modifier_categorize_method == "quartiles") {
      q <- quantile(df_mod[[modifier_var]], probs = c(0.25, 0.75), na.rm = TRUE)
      df_mod <- df_mod %>%
        mutate(!!cat_var_name := cut(.data[[modifier_var]],
                                     breaks = c(-Inf, q[1], q[2], Inf),
                                     labels = c("Q1", "Q2-Q3", "Q4"),
                                     include.lowest = TRUE))
    } else if (modifier_categorize_method == "tertiles") {
      t <- quantile(df_mod[[modifier_var]], probs = c(1/3, 2/3), na.rm = TRUE)
      df_mod <- df_mod %>%
        mutate(!!cat_var_name := cut(.data[[modifier_var]],
                                     breaks = c(-Inf, t[1], t[2], Inf),
                                     labels = c("Low", "Middle", "High"),
                                     include.lowest = TRUE))
    } else if (modifier_categorize_method == "median_split") {
      med <- median(df_mod[[modifier_var]], na.rm = TRUE)
      df_mod <- df_mod %>%
        mutate(!!cat_var_name := cut(.data[[modifier_var]],
                                     breaks = c(-Inf, med, Inf),
                                     labels = c("Below median", "Above median"),
                                     include.lowest = TRUE))
    }
    
    modifier_var <- cat_var_name
    modifier_type <- "factor"
  }
  
  # Convert to factor if specified
  if (modifier_type == "factor") {
    df_mod <- df_mod %>% mutate(!!modifier_var := as.factor(.data[[modifier_var]]))
  }
  
  # Choose random-effects structure
  if (cut_re) {
    re_term <- "(MonthsFromFirstVisit || SSBS_ID)"
  } else {
    re_term <- "(MonthsFromFirstVisit | SSBS_ID)"
  }
  
  # Build formula with 3-way interaction and flexible covariate interactions
  # Main 3-way interaction
  fixed_parts <- paste("MonthsFromFirstVisit", "*", nfl_var, "*", modifier_var)
  
  # Add time interactions for specified covariates
  if (length(interaction_covariates) > 0) {
    interaction_terms <- paste("MonthsFromFirstVisit", "*", interaction_covariates)
    fixed_parts <- c(fixed_parts, interaction_terms)
  }
  
  # Add main effect covariates
  main_effect_covariates <- setdiff(base_covariates, interaction_covariates)
  if (length(main_effect_covariates) > 0) {
    fixed_parts <- c(fixed_parts, main_effect_covariates)
  }
  
  fixed_effects <- paste(fixed_parts, collapse = " + ")
  frm <- as.formula(paste("ALSFRSR_Total ~", fixed_effects, "+", re_term))
  
  # Fit with control options
  ctrl <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
  
  # Filter to complete cases with at least 2 time points
  model_vars <- c("ALSFRSR_Total", "SSBS_ID", "MonthsFromFirstVisit", 
                  nfl_var, modifier_var, base_covariates, interaction_covariates)
  df_analysis <- df_mod %>%
    select(all_of(model_vars)) %>%
    drop_na() %>%
    group_by(SSBS_ID) %>%
    filter(n_distinct(MonthsFromFirstVisit) >= 2) %>%
    ungroup()
  
  if (nrow(df_analysis) == 0) {
    stop("No valid data after filtering")
  }
  
  # Fit REML model for visualization
  fitted_model <- lmerTest::lmer(frm, data = df_analysis, REML = TRUE, control = ctrl)
  
  # Store metadata as attributes
  attr(fitted_model, "nfl_var") <- nfl_var
  attr(fitted_model, "nfl_was_categorized") <- categorize_nfl
  attr(fitted_model, "nfl_categorize_method") <- if (categorize_nfl) nfl_categorize_method else NULL
  
  attr(fitted_model, "original_modifier_var") <- 
    if (categorize_modifier) stringr::str_remove(modifier_var, "_cat$") else modifier_var
  attr(fitted_model, "modifier_was_categorized") <- categorize_modifier
  attr(fitted_model, "modifier_categorize_method") <- if (categorize_modifier) modifier_categorize_method else NULL
  
  attr(fitted_model, "interaction_covariates") <- interaction_covariates
  
  return(fitted_model)
}

# Plotting Function --------------------------------------------------------

# Plot predicted trajectories from 3-way interaction model
plot_nfl_modifier_interaction <- function(fm,
                                          modifier_var = NULL,
                                          modifier_type = c("continuous","factor"),
                                          level_mode = c("sd","quantile"),
                                          facet_extremes_only = TRUE,
                                          facet_keep_levels = NULL,
                                          months_range = c(0,36),
                                          months_by = 3,
                                          nfL_levels_sd = c(-1,0,+1),
                                          nfL_probs = c(.25,.5,.75),
                                          mod_levels_sd = c(-1,0,+1),
                                          mod_probs = c(.25,.5,.75),
                                          base_font_size = 20,
                                          y_limits = c(-5, 45),
                                          T1label = NULL) {
  
  modifier_type <- match.arg(modifier_type)
  level_mode    <- match.arg(level_mode)
  
  mf <- model.frame(fm)
  
  # Get NfL variable name from model
  nfl_var <- attr(fm, "nfl_var")
  if (is.null(nfl_var)) nfl_var <- "NfL_z"
  if (!(nfl_var %in% names(mf))) stop("NfL variable not found in model frame: ", nfl_var)
  
  nfl_was_categorized <- !is.null(attr(fm, "nfl_was_categorized")) && attr(fm, "nfl_was_categorized")
  
  # Auto-detect modifier variable from model if not provided
  if (is.null(modifier_var)) {
    if (!is.null(attr(fm, "original_modifier_var"))) {
      original_var <- attr(fm, "original_modifier_var")
      cat_var <- paste0(original_var, "_cat")
      if (cat_var %in% names(mf)) {
        modifier_var <- cat_var
        modifier_type <- "factor"
      } else if (original_var %in% names(mf)) {
        modifier_var <- original_var
      } else {
        stop("Could not find modifier variable in model frame")
      }
    } else {
      stop("modifier_var must be specified")
    }
  }
  
  if (!(modifier_var %in% names(mf))) stop("modifier_var not found in model frame: ", modifier_var)
  
  modifier_was_categorized <- !is.null(attr(fm, "modifier_was_categorized")) && attr(fm, "modifier_was_categorized")
  original_modifier_name <- if (modifier_was_categorized) attr(fm, "original_modifier_var") else modifier_var
  
  # Time specification
  time_spec <- paste0("MonthsFromFirstVisit [",
                      months_range[1], ":", months_range[2], " by=", months_by, "]")
  
  # NfL specification
  if (nfl_was_categorized) {
    nfL_term <- nfl_var
    nfL_labels <- levels(mf[[nfl_var]])
  } else {
    if (level_mode == "sd") {
      nfL_term   <- paste0(nfl_var, " [", paste0(nfL_levels_sd, collapse=","), "]")
      nfL_labels <- ifelse(nfL_levels_sd == 0, "0 SD",
                           paste0(ifelse(nfL_levels_sd > 0, "+", ""), nfL_levels_sd, " SD"))
    } else {
      nfL_q      <- quantile(mf[[nfl_var]], probs = nfL_probs, na.rm = TRUE)
      nfL_term   <- paste0(nfl_var, " [", paste0(sprintf("%.6f", nfL_q), collapse=","), "]")
      if (length(nfL_probs) == 3 && all(nfL_probs == c(.25,.5,.75))) {
        nfL_labels <- c("Q1 (25%)", "Q2 (50%)", "Q3 (75%)")
      } else {
        nfL_labels <- paste0("Q", seq_along(nfL_q))
      }
    }
  }
  
  # Modifier specification
  terms_vec <- NULL
  
  if (modifier_type == "factor") {
    terms_vec <- c(time_spec, nfL_term, modifier_var)
  } else {
    target_vals <- NULL
    target_labels <- NULL
    
    if (level_mode == "sd") {
      target_vals <- if (facet_extremes_only) range(mod_levels_sd) else mod_levels_sd
      target_labels <- ifelse(target_vals == 0, "0 SD",
                              paste0(ifelse(target_vals > 0, "+", ""), target_vals, " SD"))
      mod_term <- paste0(modifier_var, " [", paste0(target_vals, collapse=","), "]")
    } else {
      if (facet_extremes_only) {
        raw_q <- quantile(mf[[modifier_var]], probs = c(.25,.75), na.rm = TRUE)
        target_vals   <- raw_q
        target_labels <- c("Q1 (25%)", "Q3 (75%)")
      } else {
        q1   <- as.numeric(quantile(mf[[modifier_var]], probs = .25, na.rm = TRUE))
        med  <- as.numeric(quantile(mf[[modifier_var]], probs = .50, na.rm = TRUE))
        q3   <- as.numeric(quantile(mf[[modifier_var]], probs = .75, na.rm = TRUE))
        target_vals   <- c(q1, med, q3)
        target_labels <- c("Q1 (25%)", "Q2 (50%)", "Q3 (75%)")
      }
      mod_term <- paste0(modifier_var, " [", paste0(sprintf("%.6f", target_vals), collapse=","), "]")
    }
    terms_vec <- c(time_spec, nfL_term, mod_term)
  }
  
  # Generate predictions
  gp  <- ggpredict(fm, terms = terms_vec, ci.level = 0.95)
  dfp <- as.data.frame(gp)
  
  # Map NfL groups to labels
  nfL_groups_sorted <- unique(dfp$group)
  if (nfl_was_categorized) {
    nfL_labels <- nfL_groups_sorted
  } else {
    if (length(nfL_labels) != length(nfL_groups_sorted)) {
      nfL_labels <- paste0("NfL-L", seq_along(nfL_groups_sorted))
    }
  }
  
  dfp$nfL_group <- factor(dfp$group,
                          levels = nfL_groups_sorted,
                          labels = nfL_labels[seq_along(nfL_groups_sorted)])
  
  # Color and linetype mapping
  nfL_colors    <- setNames(c("#009E73","#56B4E9","#E69F00")  #pastel 
                            [seq_along(levels(dfp$nfL_group))],
                            levels(dfp$nfL_group))
  nfL_linetypes <- setNames(c("dashed","solid","dotdash")[seq_along(levels(dfp$nfL_group))],
                            levels(dfp$nfL_group))
  
  # Get display label for the modifier variable
  if (!is.null(T1label) && original_modifier_name %in% names(T1label)) {
    display_var_name <- T1label[[original_modifier_name]]
  } else {
    clean_var_name <- stringr::str_remove(original_modifier_name, "_bl$|_z$")
    if (!is.null(T1label) && clean_var_name %in% names(T1label)) {
      display_var_name <- T1label[[clean_var_name]]
    } else {
      display_var_name <- clean_var_name
    }
  }
  
  # Create facet labels
  if (modifier_type == "factor") {
    if (!is.null(facet_keep_levels)) {
      dfp <- dfp[dfp$facet %in% facet_keep_levels, , drop = FALSE]
    }
    dfp <- dfp[!is.na(dfp$facet), , drop = FALSE]
    lab_keys <- sort(unique(dfp$facet))
    
    if (modifier_was_categorized) {
      cat_method <- attr(fm, "modifier_categorize_method")
      if (!is.null(cat_method) && cat_method == "quartiles") {
        lab_vals <- paste0(display_var_name, ": ", lab_keys)
      } else {
        lab_vals <- paste0(display_var_name, " = ", lab_keys)
      }
    } else {
      lab_vals <- paste0(display_var_name, " = ", lab_keys)
    }
    lab_map  <- setNames(lab_vals, lab_keys)
  } else {
    suppressWarnings({
      facet_num <- as.numeric(as.character(dfp$facet))
    })
    keep <- !is.na(facet_num)
    dfp <- dfp[keep, , drop = FALSE]
    
    raw_keys <- unique(dfp$facet)
    suppressWarnings({
      raw_keys_num <- as.numeric(as.character(raw_keys))
    })
    nearest_idx_keys <- vapply(raw_keys_num, function(v) which.min(abs(v - target_vals)), integer(1))
    lab_vals <- paste0(display_var_name, " = ", target_labels[nearest_idx_keys])
    names(lab_vals) <- raw_keys
    lab_map <- lab_vals
  }
  
  # Get NfL display label
  nfl_display_name <- if (!is.null(T1label) && "NfL" %in% names(T1label)) {
    T1label[["NfL"]]
  } else {
    "NfL"
  }
  
  # Create the plot
  p <- ggplot(dfp, aes(x = x, y = predicted,
                       ymin = conf.low, ymax = conf.high,
                       color = nfL_group, fill = nfL_group, linetype = nfL_group)) +
    geom_ribbon(alpha = 0.18, colour = NA) +
    geom_line(linewidth = 1.2) +
    facet_wrap(~ facet, nrow = 1, scales = "free_y",
               labeller = as_labeller(lab_map)) +
    scale_color_manual(values = nfL_colors) +
    scale_fill_manual(values = nfL_colors) +
    scale_linetype_manual(values = nfL_linetypes) +
    scale_x_continuous(
      breaks = seq(0, months_range[2], by = 12),
      limits = months_range
    )+
    labs(
      title = paste("Predicted ALSFRS-R:", nfl_display_name, "×", display_var_name,"interaction"),
      x = "Months from baseline",
      y = "Predicted ALSFRS-R",
      color = paste(nfl_display_name, "level"),
      fill  = paste(nfl_display_name, "level"),
      linetype = paste(nfl_display_name, "level")
    ) +
    theme_classic(base_size = base_font_size) +
    coord_cartesian(ylim = y_limits) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 16, face = "bold"),
      legend.text  = element_text(size = 14),
      axis.title   = element_text(size = 16),
      axis.text    = element_text(size = 14),
      axis.title.x = element_text(margin=ggplot2::margin(t=10, b=-5)),
      axis.title.y = element_text(margin=ggplot2::margin(r=10)),
      strip.text   = element_text(size = 18, face = "bold"),
      plot.title   = element_text(size = 20, face = "bold", margin=ggplot2::margin(b=15)),
      panel.spacing = unit(1, "lines")
    )
  
  print(p)
  invisible(p)
}


# label -------------------------------------------------------------------


# Define display labels for variables
T1labelex <- list(
  Age_init = "Age at baseline",
  genetic_ALS = "Genetic status",
  genetic_label = "Genetic mutation",
  OnsetSite = "Onset site",
  VC_Percent = "Vital capacity (%)",
  deltaFS = "ΔFS",
  REEC_definite = "El Escorial criteria",
  NfL = "NfL",
  NfL_bl = "Baseline NfL",
  pTau181 = "pTau181",
  pTau181_bl = "Baseline pTau181",
  pTau217 = "pTau217", 
  pTau217_bl = "Baseline pTau217",
  Ab38 = "Aβ38",
  Ab40 = "Aβ40",
  Ab42 = "Aβ42",
  Ab4240 = "Aβ42/40",
  Ab4240_bl = "Baseline Aβ42/40",
  Ab3840 = "Aβ38/40",
  Ab3840_bl = "Baseline Aβ38/40",
  Ab_status = "Aβ status",
  Ab_status_bl = "Baseline Aβ status",
  GFAP = "GFAP",
  GFAP_bl = "Baseline GFAP",
  Cr = "Creatinine",
  Cr_bl = "Baseline creatinine",
  CK = "Creatine kinase",
  CK_bl = "Baseline CK",
  Alb = "Albumin",
  Alb_bl = "Baseline albumin",
  CrCys_ratio = "Cr/CysC ratio",
  CrCys_ratio_bl = "Baseline Cr/CysC"
)

