# Load required packages
library(dplyr)
library(forcats)
library(survival)
library(broom)
library(tidyr)
library(purrr)
library(ggplot2)
library(flextable)
library(officer)

# Survival status convention in the original analysis dataset:
# 1 = death/tracheostomy event, 0 = censored.

# User configuration --------------------------------------------------

# Data preparation
df_surv_CoxFP <- survALS %>%
  mutate(
    event01 = case_when(
      surv_time_status %in% c("1", 1, TRUE)  ~ 1L,
      surv_time_status %in% c("0", 0, FALSE) ~ 0L,
      TRUE ~ NA_integer_
    ),
    Sex = as.factor(Sex),
  )

# Subgroup variable: change this to switch stratification
subgroup_var_CoxFP   <- "Ab_status"          # e.g. "Ab_status", "Sex"
subgroup_pos_CoxFP   <- "positive"            # label for "positive" group
subgroup_neg_CoxFP   <- "negative"            # label for "negative" group
subgroup_pos_disp_CoxFP <- "Aβ (+)"          # display label
subgroup_neg_disp_CoxFP <- "Aβ (−)"          # display label

# Scaling option
scale_clinical_CoxFP <- TRUE   # TRUE = scale all; FALSE = scale biomarkers only

# Facet pattern: which adjustment sets to display
facet_pattern_CoxFP <- c("Unadjusted", "Core-adjusted")

# Variable display pattern: "all", "biomarkers", "clinical", "custom"
variable_pattern_CoxFP <- "biomarkers"


# Predictor order (top to bottom in plot)
predictor_order_CoxFP <- c(
  #"Age", 
  "NfL", "GFAP", "pTau181", "pTau217",
  "Aβ38/40" , "Aβ42/40"
)

## Define predictors ---------------------------------

predictors_named_CoxFP <- c(
  pTau217       = "log_pTau217_csf",
  pTau181       = "log_pTau181_csf",
  GFAP          = "log_GFAP_csf",
  NfL           = "log_NfL_csf",
  Aβ42          = "log_Ab42_csf",
  Aβ40          = "log_Ab40_csf",
  Aβ38          = "log_Ab38_csf",
  `Aβ42/40`     = "log_Ab4240_csf",
  `Aβ38/40`     = "log_Ab38_40_csf",
  Cr            = "log_Cr",
  CK            = "log_CK",
  Alb           = "log_Alb",
  `Cr/CysC`     = "log_Cr_CysC_ratio",
  Age           = "Age_init",
  `ΔFS`         = "deltaFS",
  `%VC`         = "VC_Percent",
  REEC          = "REEC_definite",
  Sex           = "Sex",
  `Aβ status`    = "Ab_status",
  `Onset site`  = "OnsetSite"
)


# Identify variable types ---------------------------------------------
log_biomarkers_CoxFP <- grep("^log_", unname(predictors_named_CoxFP), value = TRUE)
clinical_vars_CoxFP <- c("Age_init", "deltaFS", "VC_Percent")
categorical_vars_CoxFP <- c("REEC_definite", "Sex", "OnsetSite")

# Determine which variables to scale
vars_to_scale_CoxFP <- if (scale_clinical_CoxFP) {
  c(log_biomarkers_CoxFP, clinical_vars_CoxFP)
} else {
  log_biomarkers_CoxFP
}

# Apply scaling
for (v_CoxFP in vars_to_scale_CoxFP) {
  if (v_CoxFP %in% names(df_surv_CoxFP) && is.numeric(df_surv_CoxFP[[v_CoxFP]])) {
    df_surv_CoxFP[[v_CoxFP]] <- as.numeric(scale(df_surv_CoxFP[[v_CoxFP]]))
  }
}

stopifnot(all(c("surv_time_days", "event01", "Ab_status") %in% names(df_surv_CoxFP)))
stopifnot(setequal(na.omit(unique(df_surv_CoxFP$event01)), c(0L, 1L)))

# Define core covariates
core_covars_all_CoxFP <- c("Age_init", "VC_Percent", "deltaFS", "REEC_definite", "OnsetSite")
core_covars_present_CoxFP <- intersect(core_covars_all_CoxFP, names(df_surv_CoxFP))

# Define adjustment sets (all possible)
adjustment_sets_list_CoxFP <- list(
  `Unadjusted`        = character(0),
  `Age-adjusted`      = intersect("Age_init", names(df_surv_CoxFP)),
  `Age, NfL-adjusted`    = intersect(c("Age_init", "log_NfL_csf"), names(df_surv_CoxFP)),
  `Core-adjusted`     = core_covars_present_CoxFP
)

# Filter to selected facets
adjustment_sets_list_CoxFP <- adjustment_sets_list_CoxFP[facet_pattern_CoxFP]

# Define subgroups
subgroups_vec_CoxFP <- c("All", subgroup_neg_disp_CoxFP, subgroup_pos_disp_CoxFP)

# Helper functions
filter_subgroup_df_CoxFP <- function(dat, subgroup_flag) {
  if (subgroup_flag == "All")                     return(dat)
  if (subgroup_flag == subgroup_pos_disp_CoxFP)   return(dplyr::filter(dat, .data[[subgroup_var_CoxFP]] == subgroup_pos_CoxFP))
  if (subgroup_flag == subgroup_neg_disp_CoxFP)   return(dplyr::filter(dat, .data[[subgroup_var_CoxFP]] == subgroup_neg_CoxFP))
  stop("Unknown subgroup_flag")
}

fit_extract_hr_single_CoxFP <- function(dat, time_col, status_col, focal_var, adjust_vars) {
  rhs_terms <- unique(c(focal_var, setdiff(adjust_vars, focal_var)))
  fml_txt   <- paste("Surv(", time_col, ",", status_col, ") ~ ", paste(rhs_terms, collapse = " + "))
  
  tryCatch({
    fit_obj   <- coxph(as.formula(fml_txt), data = dat, ties = "efron")
    td        <- broom::tidy(fit_obj, exponentiate = TRUE, conf.int = TRUE)
    
    if (is.factor(dat[[focal_var]])) {
      focal_terms <- grep(paste0("^", focal_var), td$term, value = TRUE)
      if (length(focal_terms) > 0) {
        comparison_level <- sub(paste0("^", focal_var), "", focal_terms[1])
        out <- td %>% dplyr::filter(term == focal_terms[1]) %>%
          dplyr::transmute(
            HR = estimate, LCL = conf.low, UCL = conf.high, p = p.value,
            comparison = comparison_level
          )
      } else {
        out <- tibble(HR = NA_real_, LCL = NA_real_, UCL = NA_real_, p = NA_real_, comparison = NA_character_)
      }
    } else {
      out <- td %>% dplyr::filter(term == focal_var) %>%
        dplyr::transmute(HR = estimate, LCL = conf.low, UCL = conf.high, p = p.value, comparison = NA_character_)
    }
    
    if (nrow(out) == 0) {
      tibble(HR = NA_real_, LCL = NA_real_, UCL = NA_real_, p = NA_real_, comparison = NA_character_)
    } else {
      out
    }
  }, error = function(e) {
    tibble(HR = NA_real_, LCL = NA_real_, UCL = NA_real_, p = NA_real_, comparison = NA_character_)
  })
}

# Filter predictors based on variable_pattern_CoxFP
predictors_to_use_CoxFP <- if (variable_pattern_CoxFP == "all") {
  predictors_named_CoxFP
} else if (variable_pattern_CoxFP == "biomarkers") {
  predictors_named_CoxFP[predictors_named_CoxFP %in% log_biomarkers_CoxFP]
} else if (variable_pattern_CoxFP == "clinical") {
  predictors_named_CoxFP[predictors_named_CoxFP %in% c(clinical_vars_CoxFP, categorical_vars_CoxFP)]
} else if (variable_pattern_CoxFP == "custom") {
  predictors_named_CoxFP[names(predictors_named_CoxFP) %in% custom_variables_CoxFP]
} else {
  predictors_named_CoxFP
}

# Create specification grid
spec_grid_CoxFP <- tibble(
  predictor_label = names(predictors_to_use_CoxFP),
  predictor_var   = unname(predictors_to_use_CoxFP)
) %>%
  tidyr::crossing(
    subgroup     = subgroups_vec_CoxFP,
    adjust_label = names(adjustment_sets_list_CoxFP)
  ) %>%
  mutate(adjust_vars = adjustment_sets_list_CoxFP[adjust_label])

# Run models
model_results_tbl_CoxFP <- purrr::pmap_dfr(
  spec_grid_CoxFP,
  function(predictor_label, predictor_var, subgroup, adjust_label, adjust_vars) {
    if (predictor_var %in% unlist(adjust_vars)) {
      return(tibble(
        adjust_label = adjust_label, predictor_label = predictor_label,
        predictor_var = predictor_var, subgroup = subgroup,
        n_sample = NA_integer_, HR = NA_real_, LCL = NA_real_, UCL = NA_real_, 
        p = NA_real_, comparison = NA_character_
      ))
    }
    
    dat_sub <- filter_subgroup_df_CoxFP(df_surv_CoxFP, subgroup)
    keep_vars <- unique(c("surv_time_days", "event01", predictor_var, unlist(adjust_vars)))
    dat_sub2 <- dat_sub %>% dplyr::select(dplyr::any_of(keep_vars)) %>% tidyr::drop_na()
    
    if (nrow(dat_sub2) < 20) {
      return(tibble(
        adjust_label = adjust_label, predictor_label = predictor_label,
        predictor_var = predictor_var, subgroup = subgroup,
        n_sample = nrow(dat_sub2), HR = NA_real_, LCL = NA_real_, UCL = NA_real_,
        p = NA_real_, comparison = NA_character_
      ))
    }
    
    fx <- fit_extract_hr_single_CoxFP(dat_sub2, "surv_time_days", "event01", predictor_var, unlist(adjust_vars))
    tibble(
      adjust_label = adjust_label, predictor_label = predictor_label,
      predictor_var = predictor_var, subgroup = subgroup, n_sample = nrow(dat_sub2),
      HR = fx$HR, LCL = fx$LCL, UCL = fx$UCL, p = fx$p, comparison = fx$comparison
    )
  }
)

# Apply custom order
valid_predictor_order_CoxFP <- intersect(predictor_order_CoxFP, names(predictors_to_use_CoxFP))

# Prepare plotting data
forest_plot_tbl_CoxFP <- model_results_tbl_CoxFP %>%
  dplyr::filter(predictor_label %in% valid_predictor_order_CoxFP) %>%
  dplyr::mutate(
    adjust_label_f = factor(adjust_label, levels = facet_pattern_CoxFP),
    subgroup_f     = factor(subgroup, levels = c("All", subgroup_neg_disp_CoxFP, subgroup_pos_disp_CoxFP)),  
    predictor_f    = factor(predictor_label, levels = valid_predictor_order_CoxFP)
  ) %>%
  dplyr::distinct(adjust_label_f, predictor_f, subgroup_f, .keep_all = TRUE) %>%
  dplyr::arrange(adjust_label_f, predictor_f, subgroup_f)

# Create y positions
forest_plot_tbl_CoxFP <- forest_plot_tbl_CoxFP %>%
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
forest_plot_tbl_CoxFP <- forest_plot_tbl_CoxFP %>%
  mutate(
    predictor_display_full = if_else(
      !is.na(comparison) & subgroup == "All",
      paste0(predictor_label, " (", comparison, ")"),
      as.character(predictor_label)
    ),
    predictor_display = if_else(subgroup == "All", predictor_display_full, ""),
    subgroup_display = as.character(subgroup)
  )

# Subtitle
hr_interpretation_CoxFP <- if (scale_clinical_CoxFP) {
  "Each point represents HR from separate Cox model (per 1 SD increase)"
} else {
  "Each point represents HR from separate Cox model (per 1 SD for biomarkers, per unit for clinical vars)"
}

# Y-axis setup
y_axis_df <- forest_plot_tbl_CoxFP %>%
  dplyr::filter(adjust_label_f == facet_pattern_CoxFP[1], subgroup == "All") %>%
  dplyr::arrange(y_pos) %>%
  dplyr::distinct(y_pos, .keep_all = TRUE)

y_breaks_CoxFP <- y_axis_df$y_pos
y_labels_CoxFP <- y_axis_df$predictor_display_full

# X-axis labels
x_axis_labels_CoxFP <- function(x) {
  sapply(x, function(val) {
    if (val >= 1 && val == floor(val)) {
      format(val, scientific = FALSE, drop0trailing = TRUE)
    } else {
      as.character(val)
    }
  })
}


# 📊 Fig5A ---------------------------------------------------------------

# Create plot
Fig5A_pub <- ggplot(forest_plot_tbl_CoxFP, aes(x = HR, y = y_pos)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray70", linewidth = 0.35) +
  geom_pointrange(
    data = forest_plot_tbl_CoxFP %>% dplyr::filter(!is.na(HR)),
    aes(xmin = LCL, xmax = UCL, color = subgroup_f),
    size = 0.45, linewidth = 0.6
  ) +
  scale_color_manual(
    name = "Subgroup",
    values = c("All" = "black", 
               setNames("#0072B2", subgroup_neg_disp_CoxFP), 
               setNames("#D55E00", subgroup_pos_disp_CoxFP)),
    labels = c("All", subgroup_neg_disp_CoxFP, subgroup_pos_disp_CoxFP)
  ) +
  facet_grid(cols = vars(adjust_label_f)) +
  scale_x_continuous(
    trans = "log10",
    breaks = c(0.25, 0.5, 1, 2, 4, 8, 16),
    labels = x_axis_labels_CoxFP
  ) +
  scale_y_continuous(
    breaks = y_breaks_CoxFP,
    labels = y_labels_CoxFP
  ) +
  coord_cartesian(xlim = c(0.25, 14), clip = "off") +
  labs(
    x = "Hazard ratio (HR)",
    y = NULL,
    #title = paste("Cox hazard ratios by", subgroup_var_CoxFP),
    title = "Cox hazard ratios by Aβ status",
    #subtitle = hr_interpretation_CoxFP
  ) +
  theme_classic(base_size = 16) +
  theme(
    strip.text.x = element_text(face = "bold", size = 16),
    panel.spacing.x = unit(16, "pt"),  # facet spacing
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 14, margin = ggplot2::margin(r = 10)),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 15, margin = ggplot2::margin(t = 15, b=-10), face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = ggplot2::margin(b = 10)),
    plot.subtitle = element_text(size = 14),
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 15, face = "bold"),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    plot.margin = unit(c(12, 12, 2, 12), "pt")
  )


Fig5A_pub




# Values ------------------------------------------------------------

# Prepare table data from forest plot results
table_data_CoxFP <- forest_plot_tbl_CoxFP %>%
  select(
    adjust_label_f,
    predictor_label,
    subgroup,
    n_sample,
    HR,
    LCL,
    UCL,
    p,
    comparison
  ) %>%
  arrange(adjust_label_f, predictor_label, subgroup) %>%
  mutate(
    # Format HR with confidence interval
    HR_CI = case_when(
      is.na(HR) ~ "NA",
      TRUE ~ sprintf("%.2f (%.2f–%.2f)", HR, LCL, UCL)
    ),
    # Format p-value
    p_formatted = case_when(
      is.na(p) ~ "NA",
      p < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p)
    ),
    # Create comparison label if applicable
    predictor_display = if_else(
      !is.na(comparison),
      paste0(predictor_label, " (", comparison, ")"),
      predictor_label
    )
  ) %>%
  select(
    Adjustment = adjust_label_f,
    Predictor = predictor_display,
    Subgroup = subgroup,
    N = n_sample,
    `HR (95% CI)` = HR_CI,
    `P-value` = p_formatted
  )

# Create flextable
ft_CoxFP <- flextable(table_data_CoxFP) %>%
  # Set column widths
  width(j = "Adjustment", width = 1.5) %>%
  width(j = "Predictor", width = 1.8) %>%
  width(j = "Subgroup", width = 0.8) %>%
  width(j = "N", width = 0.6) %>%
  width(j = "HR (95% CI)", width = 1.5) %>%
  width(j = "P-value", width = 0.8) %>%
  # Format header
  bold(part = "header") %>%
  align(align = "center", part = "header") %>%
  # Format body
  align(j = c("Adjustment", "Predictor", "Subgroup"), align = "left", part = "body") %>%
  align(j = c("N", "HR (95% CI)", "P-value"), align = "center", part = "body") %>%
  # Add borders
  border_outer(part = "all", border = officer::fp_border(width = 2)) %>%
  border_inner_h(part = "all", border = officer::fp_border(width = 1)) %>%
  border_inner_v(part = "all", border = officer::fp_border(width = 1)) %>%
  # Merge cells for adjustment groups
  flextable::merge_v(j = "Adjustment") %>%
  flextable::valign(j = "Adjustment", valign = "top") %>%
  # Set font
  flextable::font(fontname = "Arial", part = "all") %>%
  flextable::fontsize(size = 10, part = "all") %>%
  # Add background color to header
  flextable::bg(bg = "#D3D3D3", part = "header") %>%
  # Add alternating row colors for better readability
  flextable::theme_vanilla()


table_data_CoxFP %>% filter(Adjustment=="Core-adjusted", Predictor=="NfL")
