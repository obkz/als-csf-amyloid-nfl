library(dplyr)
library(forcats)
library(survival)
library(broom)
library(tidyr)
library(purrr)
library(ggplot2)
library(flextable)
library(officer)
library(ggtext)


# Data preparation ---------------------------------------------------------------

df_int_fp2 <- survALS |>
  mutate(
    event01       = as.integer(surv_time_status == 1),
    age_z         = as.numeric(scale(Age_init)),
    deltafs_z     = as.numeric(scale(deltaFS)),
    vc_z          = as.numeric(scale(VC_Percent)),
    log_nfl_z     = as.numeric(scale(log(NfL_csf_pgml))),
    ab_status     = Ab_status,
    surv_days     = surv_time_days,
    reec_definite = as.integer(REEC_definite),
    onset_bulbar  = as.integer(OnsetSite == "bulbar")
  ) |>
  filter(surv_days > 0) |>
  dplyr::select(surv_days, event01, age_z, deltafs_z, vc_z,
                reec_definite, onset_bulbar, log_nfl_z, ab_status)

df_ext_fp2 <- ext |>
  filter(survival_days > 0) |>
  mutate(
    event01       = as.integer(survival_status == 1),
    age_z         = as.numeric(scale(age_at_baseline)),
    deltafs_z     = as.numeric(scale(delta_fs_baseline)),
    vc_z          = as.numeric(scale(fvc_percent_predicted_baseline)),
    log_nfl_z     = as.numeric(scale(log(nfl_csf_pgml))),
    ab_status     = if_else(ab42_40_csf_ratio < 0.068, "positive", "negative"),
    surv_days     = survival_days,
    reec_definite = if_else(
      el_escorial_category_baseline == "definite", 1L, 0L, missing = NA_integer_
    ),
    onset_bulbar = case_when(
      is.na(onset_site) ~ NA_integer_,
      onset_site %in% c("bulbar", "respiratory") ~ 1L,
      TRUE ~ 0L
    )
  ) |>
  dplyr::select(surv_days, event01, age_z, deltafs_z, vc_z,
                reec_definite, onset_bulbar, log_nfl_z, ab_status)

df_pooled_fp2 <- bind_rows(
  df_int_fp2 |> mutate(cohort_fp2 = "Internal"),
  df_ext_fp2 |> mutate(cohort_fp2 = "Bologna")
) |>
  mutate(cohort_fp2 = factor(cohort_fp2, levels = c("Internal", "Bologna")))


# N labels ---------------------------------------------------------------

n_label_fp2 <- list(
  Internal = nrow(drop_na(dplyr::select(df_int_fp2, surv_days, event01, log_nfl_z))),
  Bologna = nrow(drop_na(dplyr::select(df_ext_fp2, surv_days, event01, log_nfl_z))),
  Pooled   = nrow(drop_na(dplyr::select(df_pooled_fp2, surv_days, event01, log_nfl_z)))
)

cat("N per cohort:\n")
print(n_label_fp2)

# Helper functions ------------------------------------------------

filter_subgroup_fp2 <- function(dat_fp2, subgroup_fp2) {
  if (subgroup_fp2 == "All")    return(dat_fp2)
  if (subgroup_fp2 == "Aβ (+)") return(dplyr::filter(dat_fp2, ab_status == "positive"))
  if (subgroup_fp2 == "Aβ (−)") return(dplyr::filter(dat_fp2, ab_status == "negative"))
  stop("Unknown subgroup: ", subgroup_fp2)
}

fit_nfl_hr_fp2 <- function(dat_fp2, adjust_vars_fp2, cohort_stratify_fp2 = FALSE) {
  keep_fp2 <- unique(c("surv_days", "event01", "log_nfl_z",
                       adjust_vars_fp2,
                       if (cohort_stratify_fp2) "cohort_fp2" else NULL))
  df_fp2 <- dat_fp2 |> dplyr::select(any_of(keep_fp2)) |> tidyr::drop_na()
  
  if (nrow(df_fp2) < 10) return(NULL)
  
  rhs_fp2 <- paste(c(adjust_vars_fp2, "log_nfl_z",
                     if (cohort_stratify_fp2) "strata(cohort_fp2)" else NULL),
                   collapse = " + ")
  fml_fp2 <- as.formula(paste0("Surv(surv_days, event01) ~ ", rhs_fp2))
  
  fit_fp2 <- tryCatch(coxph(fml_fp2, data = df_fp2), error = function(e) NULL)
  if (is.null(fit_fp2)) return(NULL)
  
  tidy_fp2    <- broom::tidy(fit_fp2, exponentiate = TRUE, conf.int = TRUE)
  nfl_row_fp2 <- dplyr::filter(tidy_fp2, term == "log_nfl_z")
  if (nrow(nfl_row_fp2) == 0) return(NULL)
  
  tibble::tibble(
    hr       = nfl_row_fp2$estimate,
    lcl      = nfl_row_fp2$conf.low,
    ucl      = nfl_row_fp2$conf.high,
    p        = nfl_row_fp2$p.value,
    n        = nrow(df_fp2),
    n_events = sum(df_fp2$event01)
  )
}


# Specification and model runs ------------------------------------------------

subgroups_fp2     <- c("All", "Aβ (−)", "Aβ (+)")
cohort_labels_fp2 <- c("NfL (Nagoya)", "NfL (Bologna)", "NfL (Pooled)")
facet_levels_fp2  <- c("Unadjusted", "Core-adjusted", "Core5-adjusted")
subgroup_levels_fp2 <- c("Aβ (+)", "Aβ (−)", "All")
predictor_levels_fp2 <- c("NfL (Pooled)", "NfL (Bologna)", "NfL (Nagoya)")

# Core5 onset coding was used as a supplementary harmonized approximation.
# Core3-adjusted models were treated as the primary adjusted external-validation analysis.
adjust_sets_fp2 <- list(
  "Unadjusted"     = character(0),
  "Core-adjusted"  = c("age_z", "deltafs_z", "vc_z"),
  "Core5-adjusted" = c("age_z", "deltafs_z", "vc_z", "reec_definite", "onset_bulbar")
)

spec_fp2 <- tidyr::crossing(
  subgroup_fp2     = subgroups_fp2,
  predictor_fp2    = cohort_labels_fp2,
  adjust_label_fp2 = names(adjust_sets_fp2)
)

results_fp2 <- purrr::pmap_dfr(spec_fp2, function(subgroup_fp2, predictor_fp2, adjust_label_fp2) {
  adj_vars_fp2 <- adjust_sets_fp2[[adjust_label_fp2]]
  dat_fp2 <- switch(predictor_fp2,
                    "NfL (Nagoya)"  = df_int_fp2,
                    "NfL (Bologna)" = df_ext_fp2,
                    "NfL (Pooled)"   = df_pooled_fp2
  )
  stratify_fp2 <- predictor_fp2 == "NfL (Pooled)"
  dat_sub_fp2  <- filter_subgroup_fp2(dat_fp2, subgroup_fp2)
  res_fp2      <- fit_nfl_hr_fp2(dat_sub_fp2, adj_vars_fp2, stratify_fp2)
  if (is.null(res_fp2)) return(NULL)
  res_fp2 |> mutate(
    subgroup_fp2     = subgroup_fp2,
    predictor_fp2    = predictor_fp2,
    adjust_label_fp2 = adjust_label_fp2
  )
})

cat("Results:\n")
print(results_fp2, n = Inf)

# Forest plot data preparation ------------------------------------------------

predictor_labels_fp2 <- c(
  "NfL (Pooled)"   = sprintf("**Pooled**\n\n*N* = %d", n_label_fp2$Pooled),
  "NfL (Bologna)" = sprintf("**Bologna**\n\n*N* = %d", n_label_fp2$Bologna),
  "NfL (Nagoya)"  = sprintf("**Nagoya**\n\n*N* = %d", n_label_fp2$Internal)
)

forest_fp2 <- results_fp2 |>
  mutate(
    adjust_f_fp2    = factor(adjust_label_fp2, levels = facet_levels_fp2),
    predictor_f_fp2 = factor(predictor_fp2,   levels = predictor_levels_fp2),
    subgroup_f_fp2  = factor(subgroup_fp2,     levels = subgroup_levels_fp2)
  )

subgroup_colors_fp2 <- c(
  "All"    = "black",
  "Aβ (−)" = "#0072B2",
  "Aβ (+)" = "#D55E00"
)

# Forest plot (Unadjusted + Core-adjusted only) ------------------------------------------------

Fig5_added <- ggplot(
  forest_fp2 |> dplyr::filter(adjust_label_fp2 != "Core5-adjusted"),
  aes(x = hr, y = predictor_f_fp2, color = subgroup_f_fp2)
) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray70", linewidth = 0.4) +
  geom_pointrange(
    aes(xmin = lcl, xmax = ucl),
    position = position_dodge(width = 0.55),
    size = 0.45, linewidth = 0.6
  ) +
  scale_color_manual(
    name   = "Subgroup",
    values = subgroup_colors_fp2,
    breaks = c("All", "Aβ (−)", "Aβ (+)"),
    labels = c("All", "Aβ (−)", "Aβ (+)")
  ) + 
  scale_x_continuous(
    trans  = "log10",
    breaks = c(0.5, 1, 2, 4, 8),
    labels = c("0.5", "1", "2", "4", "8")
  ) +
  scale_y_discrete(labels = predictor_labels_fp2) +
  facet_grid(cols = vars(adjust_f_fp2)) +
  coord_cartesian(xlim = c(0.3, 12), clip = "off") +
  labs(
    x     = "CSF NfL hazard ratio (HR)",
    y     = NULL,
    title = "External validation of NfL HR by Aβ status"
  ) +
  theme_classic(base_size = 15) +
  theme(
    strip.text.x       = element_text(face = "bold", size = 16),
    panel.spacing.x    = unit(16, "pt"),
    axis.ticks.y       = element_blank(),
    axis.text.y        = element_markdown(size = 14, margin = ggplot2::margin(r = 10), hjust = 1),
    axis.text.x        = element_text(size = 14),
    axis.title.x       = element_text(size = 15, face = "bold",
                                      margin = ggplot2::margin(t = 15, b = -10)),
    plot.title         = element_text(size = 20, face = "bold", hjust = 0.5,
                                      margin = ggplot2::margin(b = 10)),
    legend.position    = "bottom",
    legend.text        = element_text(size = 14),
    legend.title       = element_text(size = 15, face = "bold"),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    plot.margin        = unit(c(12, 20, 8, 12), "pt")
  )


## 📊 Figure 5C ---------------------------------------------------------------

Fig5C <- Fig5_added 

# Flextable ------------------------------------------------

table_fp2_wide <- forest_fp2 |>
  dplyr::filter(subgroup_fp2 != "All") |>
  mutate(
    cohort_display = case_when(
      predictor_fp2 == "NfL (Nagoya)"  ~ "(Nagoya)",
      predictor_fp2 == "NfL (Bologna)" ~ "(Bologna)",
      predictor_fp2 == "NfL (Pooled)"   ~ "(Pooled)"
    ),
    row_label = factor(
      paste0(cohort_display, " ", subgroup_fp2),
      levels = c(
        "(Nagoya) Aβ (−)",  "(Nagoya) Aβ (+)",
        "(Bologna) Aβ (−)", "(Bologna) Aβ (+)",
        "(Pooled) Aβ (−)",   "(Pooled) Aβ (+)"
      )
    ),
    hr_ci = sprintf("%.2f (%.2f–%.2f)", hr, lcl, ucl),
    p_fmt = case_when(p < 0.001 ~ "<0.001", TRUE ~ sprintf("%.3f", p)),
    n_evt = sprintf("%d / %d", n_events, n)
  ) |>
  dplyr::select(row_label, adjust_label_fp2, hr_ci, p_fmt, n_evt) |>
  tidyr::pivot_wider(
    names_from  = adjust_label_fp2,
    values_from = c(hr_ci, p_fmt, n_evt),
    names_sep   = "___"
  ) |>
  dplyr::select(
    row_label,
    `hr_ci___Unadjusted`,      `p_fmt___Unadjusted`,      `n_evt___Unadjusted`,
    `hr_ci___Core-adjusted`,   `p_fmt___Core-adjusted`,   `n_evt___Core-adjusted`,
    `hr_ci___Core5-adjusted`,  `p_fmt___Core5-adjusted`,  `n_evt___Core5-adjusted`
  ) |>
  arrange(row_label)

ft_fp2_CoxSub <- flextable(table_fp2_wide) |>
  set_header_labels(
    row_label                = "Cohort / Subgroup",
    `hr_ci___Unadjusted`     = "NfL HR (95% CI)",
    `p_fmt___Unadjusted`     = "p-value",
    `n_evt___Unadjusted`     = "Events / N",
    `hr_ci___Core-adjusted`  = "NfL HR (95% CI)",
    `p_fmt___Core-adjusted`  = "p-value",
    `n_evt___Core-adjusted`  = "Events / N",
    `hr_ci___Core5-adjusted` = "NfL HR (95% CI)",
    `p_fmt___Core5-adjusted` = "p-value",
    `n_evt___Core5-adjusted` = "Events / N"
  ) |>
  add_header_row(
    values    = c("", "Unadjusted",
                  "Core 3-adjusted\n(age, ΔFS, %VC)",
                  "Core 5-adjusted\n(age, ΔFS, %VC, REEC, onset)"),
    colwidths = c(1, 3, 3, 3)
  ) |>
  merge_v(part = "header") |>
  flextable::align(align = "center", part = "header") |>
  flextable::align(j = 1,    align = "left",   part = "body") |>
  flextable::align(j = 2:10, align = "center", part = "body") |>
  flextable::bold(part = "header") |>
  flextable::hline(i = c(2, 4), border = fp_border(color = "gray60", width = 1)) |>
  flextable::hline_top(border    = fp_border(color = "black", width = 1), part = "header") |>
  flextable::hline_bottom(border = fp_border(color = "black", width = 1), part = "header") |>
  flextable::hline_bottom(border = fp_border(color = "black", width = 1), part = "body") |>
  flextable::font(fontname = "Times New Roman", part = "all") |>
  flextable::fontsize(size = 10, part = "all") |>
  # title
  flextable::set_caption(caption = flextable::as_paragraph(
    flextable::as_chunk("Supplementary Table S12: External validation of Aβ status-stratified NfL prognostic values on survival.",
                        props = officer::fp_text(font.size   = 12, bold = TRUE, 
                                                 font.family = "Times New Roman"))),
    word_stylename = "Table Caption", 
    fp_p = officer::fp_par(text.align = "left", padding= 5)) |>
  flextable::add_footer_lines(
    values = c(
      "Aβ-positive status was defined using cohort-specific a priori thresholds.",
      "NfL was log-transformed and standardised within each cohort prior to modelling; hazard ratios are therefore not directly comparable across cohorts on the original scale.",
      paste0(
        "In the Bologna cohort, Core 3-adjusted models were specified as the primary adjusted analysis, ",
        "as the limited number of Aβ-positive events precluded full Core 5 adjustment without substantial risk of overfitting. ",
        "Core 5-adjusted results are shown for completeness."
      ),
      "Pooled analyses used cohort-stratified Cox models to account for between-cohort differences in baseline hazard.",
      "Abbreviations: Aβ, amyloid-beta; CI, confidence interval; ΔFS, delta functional score; HR, hazard ratio; NfL, neurofilament light chain; REEC, revised El Escorial criteria; %VC, percent predicted vital capacity."
    )
  ) |> 
  flextable::autofit()

ft_fp2_CoxSub



## 🔍 SuppleT12 ------------------------------------------------------------------

SuppleT12_pub <- ft_fp2_CoxSub

# docx
# SuppleT12_pub |> 
#   flextable::save_as_docx(
#     path = "output/Tables/SuppleT12.docx",
#     pr_section = officer::prop_section(
#       page_size = officer::page_size(orient = "landscape", width = 8.5, height = 11),
#       page_margins = officer::page_mar(bottom = 0.5, top = 0.5, left = 0.5, right = 0.5)
#     )
#   )
