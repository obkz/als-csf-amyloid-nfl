# Subgroup Forest Plot with Flexible Covariate Interactions
# Time × Variable interaction effects by Aβ status subgroup and adjustment
library(dplyr)
library(forcats)
library(lmerTest)
library(broom.mixed)
library(tidyr)
library(purrr)
library(ggplot2)

# User Configuration -------------------------------------------------------

# Subgroup variable: change this to switch stratification
subgroup_var_lme_FP      <- "Ab_status_bl"     # e.g. "Ab_status_bl", "Sex"
subgroup_pos_lme_FP      <- "positive"
subgroup_neg_lme_FP      <- "negative"
subgroup_pos_disp_lme_FP <- "Aβ+"
subgroup_neg_disp_lme_FP <- "Aβ−"


# Facet pattern: which adjustment sets to display
facet_pattern_lme_FP <- c("Unadjusted", "Core-adjusted")

# Variable display pattern
variable_pattern_lme_FP <- "biomarkers" # Options: "all", "clinical", "custom"

# Predictor order (top to bottom in plot)
predictor_order_lme_FP <- c(
  #"Age", 
  "NfL", "GFAP", "pTau181", "pTau217","Aβ38/40", "Aβ42/40"
)

# Define predictors with variable names in df_vis
predictors_named_lme_FP <- c(
  NfL       = "NfL_z",
  GFAP      = "GFAP_z",
  pTau181   = "pTau181_z",
  pTau217   = "pTau217_z",
  `Aβ38/40` = "Ab3840_z",
  `Aβ42/40` = "Ab4240_z"
)

# Identify variable types
biomarkers_lme_FP <- grep("_z$", unname(predictors_named_lme_FP), value = TRUE)
clinical_vars_lme_FP <- c("Age_init", "deltaFS", "VC_Percent")
categorical_vars_lme_FP <- c("REEC_definite", "Sex", "OnsetSite")

# Define core covariates for adjustment
core_covars_all_lme_FP <- c("Age_init", "VC_Percent", "deltaFS", "REEC_binary", "OnsetSite")

core_covars_present_lme_FP <- intersect(core_covars_all_lme_FP, names(df_vis))

# Covariate Interaction Configuration --------------------------------------

interaction_covariates_config_FP <- list(
  `Unadjusted`        = character(0),  # Always no interactions for unadjusted
  `Age-adjusted`      = c("Age_init"),  
  `Age+NfL-adjusted`  = c("Age_init","NfL_z"), 
  `Core-adjusted`     = c("Age_init", "VC_Percent", "deltaFS", "REEC_binary", "OnsetSite") 
)

# Define Adjustment Sets ---------------------------------------------------

# Base covariates (main effects) for each adjustment level
adjustment_sets_base_lme_FP <- list(
  `Unadjusted`        = character(0),
  `Age-adjusted`      = intersect("Age_init", names(df_vis)),
  `Age, NfL-adjusted`  = intersect(c("Age_init", "NfL_z"), names(df_vis)),
  `Core-adjusted`     = core_covars_present_lme_FP
)


# Filter to selected facets
adjustment_sets_base_lme_FP <- adjustment_sets_base_lme_FP[facet_pattern_lme_FP]
interaction_covariates_config_FP <- interaction_covariates_config_FP[facet_pattern_lme_FP]

# Create full adjustment specifications combining base and interactions
adjustment_specs_lme_FP <- purrr::map2(
  adjustment_sets_base_lme_FP,
  interaction_covariates_config_FP,
  function(base_covars, int_covars) {
    list(
      base_covariates = base_covars,
      interaction_covariates = int_covars
    )
  }
)

# Define Subgroups ---------------------------------------------------------

# Subgroup vector (dynamically generated)
subgroups_vec_lme_FP <- c("All", subgroup_neg_disp_lme_FP, subgroup_pos_disp_lme_FP) #subgroups_vec_lme_FP <- c("All", "Aβ−", "Aβ+")


# Helper function to filter by subgroup
filter_subgroup_df_lme_FP <- function(dat, subgroup_flag) {
  if (subgroup_flag == "All")                      return(dat)
  if (subgroup_flag == subgroup_pos_disp_lme_FP)   return(dplyr::filter(dat, .data[[subgroup_var_lme_FP]] == subgroup_pos_lme_FP))
  if (subgroup_flag == subgroup_neg_disp_lme_FP)   return(dplyr::filter(dat, .data[[subgroup_var_lme_FP]] == subgroup_neg_lme_FP))
  stop("Unknown subgroup_flag")
}


# Model Fitting Function ---------------------------------------------------

# Helper function to fit lmer model with flexible covariate interactions
fit_extract_beta_single_lme_FP <- function(dat, time_col, outcome_col, focal_var, 
                                           base_covars, interaction_covars) {
  # Build formula with Time*Biomarker interaction and flexible covariate interactions
  # Random effects: (time || ID) for numerical stability
  re_term <- paste0("(", time_col, " || SSBS_ID)")
  
  # Main term: Time*focal_var
  rhs_parts <- paste(time_col, "*", focal_var)
  
  # Add time interactions for specified covariates
  if (length(interaction_covars) > 0) {
    interaction_terms <- paste(time_col, "*", interaction_covars)
    rhs_parts <- c(rhs_parts, interaction_terms)
  }
  
  # Add main effect covariates (excluding those with time interactions)
  main_effect_covars <- setdiff(base_covars, interaction_covars)
  if (length(main_effect_covars) > 0) {
    rhs_parts <- c(rhs_parts, main_effect_covars)
  }
  
  fml_txt <- paste(outcome_col, "~", paste(unique(rhs_parts), collapse = " + "), "+", re_term)
  
  tryCatch({
    # Fit model with REML
    ctrl <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
    fit_obj <- lmerTest::lmer(as.formula(fml_txt), data = dat, REML = TRUE, control = ctrl)
    
    # Extract coefficients
    td <- broom.mixed::tidy(fit_obj, effects = "fixed", conf.int = TRUE)
    
    # Identify interaction term for focal variable
    if (is.factor(dat[[focal_var]])) {
      lvs <- levels(dat[[focal_var]])
      if (length(lvs) != 2) {
        return(tibble(beta = NA_real_, LCL = NA_real_, UCL = NA_real_, p = NA_real_, comparison = NA_character_))
      }
      interaction_term <- paste0(time_col, ":", focal_var, lvs[2])
      comparison_level <- lvs[2]
    } else {
      interaction_term <- paste0(time_col, ":", focal_var)
      comparison_level <- NA_character_
    }
    
    out <- td %>% dplyr::filter(term == interaction_term) %>%
      dplyr::transmute(
        beta = estimate,
        LCL = conf.low,
        UCL = conf.high,
        p = p.value,
        comparison = comparison_level
      )
    
    if (nrow(out) == 0) {
      tibble(beta = NA_real_, LCL = NA_real_, UCL = NA_real_, p = NA_real_, comparison = NA_character_)
    } else {
      out
    }
  }, error = function(e) {
    tibble(beta = NA_real_, LCL = NA_real_, UCL = NA_real_, p = NA_real_, comparison = NA_character_)
  })
}

# Filter Predictors --------------------------------------------------------

predictors_to_use_lme_FP <- if (variable_pattern_lme_FP == "all") {
  predictors_named_lme_FP
} else if (variable_pattern_lme_FP == "biomarkers") {
  predictors_named_lme_FP[predictors_named_lme_FP %in% biomarkers_lme_FP]
} else if (variable_pattern_lme_FP == "clinical") {
  predictors_named_lme_FP[predictors_named_lme_FP %in% c(clinical_vars_lme_FP, categorical_vars_lme_FP)]
} else if (variable_pattern_lme_FP == "custom") {
  predictors_named_lme_FP[names(predictors_named_lme_FP) %in% custom_variables_lme_FP]
} else {
  predictors_named_lme_FP
}

# Create Specification Grid ------------------------------------------------

spec_grid_lme_FP <- tibble(
  predictor_label = names(predictors_to_use_lme_FP),
  predictor_var   = unname(predictors_to_use_lme_FP)
) %>%
  tidyr::crossing(
    subgroup     = subgroups_vec_lme_FP,
    adjust_label = names(adjustment_specs_lme_FP)
  ) %>%
  mutate(
    base_covars = purrr::map(adjustment_specs_lme_FP[adjust_label], "base_covariates"),
    interaction_covars = purrr::map(adjustment_specs_lme_FP[adjust_label], "interaction_covariates")
  )

# Run Models ---------------------------------------------------------------

model_results_tbl_lme_FP <- purrr::pmap_dfr(
  spec_grid_lme_FP,
  function(predictor_label, predictor_var, subgroup, adjust_label, 
           base_covars, interaction_covars) {
    # Skip if predictor is already in adjustment set
    all_adjust_vars <- c(unlist(base_covars), unlist(interaction_covars))
    if (predictor_var %in% all_adjust_vars) {
      return(tibble(
        adjust_label = adjust_label, predictor_label = predictor_label,
        predictor_var = predictor_var, subgroup = subgroup,
        n_sample = NA_integer_, n_subjects = NA_integer_,
        beta = NA_real_, LCL = NA_real_, UCL = NA_real_, 
        p = NA_real_, comparison = NA_character_
      ))
    }
    
    # Filter to subgroup
    dat_sub <- filter_subgroup_df_lme_FP(df_vis, subgroup)
    
    # Keep necessary variables
    keep_vars <- unique(c("ALSFRSR_Total", "SSBS_ID", "MonthsFromFirstVisit", 
                          predictor_var, all_adjust_vars))
    dat_sub2 <- dat_sub %>% 
      dplyr::select(dplyr::any_of(keep_vars)) %>% 
      tidyr::drop_na()
    
    # Require at least 2 time points per subject
    dat_sub2 <- dat_sub2 %>%
      dplyr::group_by(SSBS_ID) %>%
      dplyr::filter(dplyr::n_distinct(MonthsFromFirstVisit) >= 2) %>%
      dplyr::ungroup()
    
    # Check minimum sample size
    n_subj <- dplyr::n_distinct(dat_sub2$SSBS_ID)
    n_obs <- nrow(dat_sub2)
    
    if (n_subj < 20 || n_obs < 40) {
      return(tibble(
        adjust_label = adjust_label, predictor_label = predictor_label,
        predictor_var = predictor_var, subgroup = subgroup,
        n_sample = n_obs, n_subjects = n_subj,
        beta = NA_real_, LCL = NA_real_, UCL = NA_real_,
        p = NA_real_, comparison = NA_character_
      ))
    }
    
    # Fit model and extract beta
    fx <- fit_extract_beta_single_lme_FP(
      dat_sub2, 
      "MonthsFromFirstVisit", 
      "ALSFRSR_Total", 
      predictor_var, 
      unlist(base_covars),
      unlist(interaction_covars)
    )
    
    tibble(
      adjust_label = adjust_label, predictor_label = predictor_label,
      predictor_var = predictor_var, subgroup = subgroup,
      n_sample = n_obs, n_subjects = n_subj,
      beta = fx$beta, LCL = fx$LCL, UCL = fx$UCL, 
      p = fx$p, comparison = fx$comparison
    )
  }
)

# Prepare Plotting Data ----------------------------------------------------

# Apply custom order
valid_predictor_order_lme_FP <- intersect(predictor_order_lme_FP, names(predictors_to_use_lme_FP))


forest_plot_tbl_lme_FP <- model_results_tbl_lme_FP %>%
  dplyr::filter(predictor_label %in% valid_predictor_order_lme_FP) %>%
  dplyr::mutate(
    adjust_label_f = factor(adjust_label, levels = facet_pattern_lme_FP),
    subgroup_f = factor(subgroup, levels = c("All", subgroup_neg_disp_lme_FP, subgroup_pos_disp_lme_FP)), 
    #subgroup_f = factor(subgroup, levels = c("All", "Aβ−", "Aβ+")),
    predictor_f    = factor(predictor_label, levels = valid_predictor_order_lme_FP)
  ) %>%
  dplyr::distinct(adjust_label_f, predictor_f, subgroup_f, .keep_all = TRUE) %>%
  dplyr::arrange(adjust_label_f, predictor_f, subgroup_f)

# Check number of subjects per subgroup
forest_plot_tbl_lme_FP %>% print(n = Inf)


# Create y positions with spacing between predictor groups
forest_plot_tbl_lme_FP <- forest_plot_tbl_lme_FP %>%
  group_by(adjust_label_f, predictor_f) %>%
  mutate(subgroup_num = as.numeric(subgroup_f)) %>%
  ungroup() %>%
  arrange(adjust_label_f, predictor_f, subgroup_f) %>%
  group_by(adjust_label_f) %>%
  mutate(
    predictor_num = as.numeric(predictor_f),
    y_pos = max(predictor_num) * (3 + 1.8) - (subgroup_num + (predictor_num - 1) * (3 + 1.8)) + 1
  ) %>%
  ungroup()

# Create labels
forest_plot_tbl_lme_FP <- forest_plot_tbl_lme_FP %>%
  mutate(
    predictor_display_full = if_else(
      !is.na(comparison) & subgroup == "All",
      paste0(predictor_label, " (", comparison, ")"),
      as.character(predictor_label)
    ),
    predictor_display = if_else(subgroup == "All", predictor_display_full, ""),
    subgroup_display = as.character(subgroup)
  )


# Y-axis setup
y_axis_df_lme_FP <- forest_plot_tbl_lme_FP %>%
  dplyr::filter(adjust_label_f == facet_pattern_lme_FP[1], subgroup == "All") %>%
  dplyr::arrange(y_pos) %>%
  dplyr::distinct(y_pos, .keep_all = TRUE)

y_breaks_lme_FP <- y_axis_df_lme_FP$y_pos
y_labels_lme_FP <- y_axis_df_lme_FP$predictor_display_full

# X-axis labels
x_axis_labels_lme_FP <- function(x) {
  sapply(x, function(val) {
    if (abs(val) >= 1 && abs(val) == floor(abs(val))) {
      format(val, scientific = FALSE, drop0trailing = TRUE)
    } else {
      as.character(val)
    }
  })
}

# Set x-axis limits
x_breaks_lme_FP <- seq(-1.5, 1.5, by = 1)
x_limits_lme_FP <- c(-1.75, 1.75)

# 🔍 SuppleFig S10  -------------------------------------------------------

p_lme_FP <- ggplot(forest_plot_tbl_lme_FP, aes(x = beta, y = y_pos)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70", linewidth = 0.35) +
  geom_pointrange(
    data = forest_plot_tbl_lme_FP %>% dplyr::filter(!is.na(beta)),
    aes(xmin = LCL, xmax = UCL, color = subgroup_f),
    size = 0.45, linewidth = 0.6
  ) +
  scale_color_manual(
    name = "Subgroup",
    values = c("All" = "black",
               setNames("#0072B2", subgroup_neg_disp_lme_FP),
               setNames("#D55E00", subgroup_pos_disp_lme_FP)),
    labels = c("All", subgroup_neg_disp_lme_FP, subgroup_pos_disp_lme_FP)
  ) +
  facet_grid(cols = vars(adjust_label_f)) +
  scale_x_continuous(
    breaks = x_breaks_lme_FP,
    labels = x_axis_labels_lme_FP
  ) +
  scale_y_continuous(
    breaks = y_breaks_lme_FP,
    labels = y_labels_lme_FP
  ) +
  coord_cartesian(xlim = x_limits_lme_FP, clip = "off") +
  labs(
    x = "β coefficient (time × biomarker interaction)",
    y = NULL,
    #title = paste0("Time × biomarker interaction estimates by ", subgroup_var_lme_FP),
    title = "Time × biomarker interaction estimates for ALSFRS-R decline",
    #subtitle = beta_interpretation_lme_FP
  ) +
  theme_classic(base_size = 16) +
  theme(
    strip.text.x = element_text(face = "bold", size = 16),
    panel.spacing.x = unit(12, "pt"),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 15, margin = ggplot2::margin(t = 15), face = "bold"),
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 14),
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 15, face = "bold"),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    plot.margin = unit(c(12, 12, 2, 12), "pt")
  )

SuppleFig10 <- p_lme_FP
SuppleFig10


# QC ----------------------------------------------------------------------

# 0 subjects with missing beta estimates and >=20 subjects in subgroup
model_results_tbl_lme_FP %>%
  filter(is.na(beta) & !is.na(n_subjects) & n_subjects >= 20)

