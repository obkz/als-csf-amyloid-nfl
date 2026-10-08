library(flextable)
library(officer)
library(dplyr)
library(tidyr)

# Sensitivity Analysis: Covariate Handling in LME Models ----
# Compares the Time x Biomarker interaction estimate across three covariate
# specifications: (1) no additional covariates (existing fit_unadjusted_lmer),
# (2) core covariates as main effects only, no time interaction (new),
# (3) core covariates with time interactions, i.e. current "Core" model
# (existing fit_core_lmer). No changes made to fit_lmer_interaction(),
# fit_unadjusted_lmer(), or fit_core_lmer().

fit_core_mainonly_lmer <- function(pred, core_covariates = core_covariates_lmer) {
  if (pred %in% core_covariates) return(NULL)
  fit_lmer_interaction(
    predictor_col = pred,
    base_covariates = core_covariates,        # main effects only
    interaction_covariates = character(0),    # no covariate x time terms
    base_label = "Core (main effects only)",
    model_label = "Core, main effects only"
  )
}

# Run all three specifications for each Supple T8 biomarker
sensitivity_results_supp8 <- purrr::map_dfr(lmer_predictors, function(pred) {
  dplyr::bind_rows(
    fit_unadjusted_lmer(pred),
    fit_core_mainonly_lmer(pred, core_covariates = core_covariates_lmer),
    fit_core_lmer(pred, core_covariates = core_covariates_lmer,
                  interaction_covariates = interaction_covariates_config)
  )
})

sensitivity_results_supp8 %>%
  dplyr::select(biomarker_label, model_label, beta, beta_lcl, beta_ucl, p_value, n_obs, n_subjects) |> print(n=Inf)




# Sensitivity Analysis: Covariate Handling in Time x NfL x Biomarker Models ----------
# Compares the three-way interaction estimate (Time x NfL x Biomarker)
# across three covariate specifications:
#   (1) no additional covariates
#   (2) core covariates as main effects only, no time interaction
#   (3) core covariates with time interactions, i.e. the current
#       "NfL+Core adjusted" model (existing fit_core_nfl_interaction())
# No changes made to fit_nfl_interaction() or fit_core_nfl_interaction().

fit_nfl_interaction_unadjusted <- function(pred) {
  if (pred %in% "NfL_z") return(NULL)
  fit_nfl_interaction(
    biomarker_col = pred,
    base_covariates = character(0),
    interaction_covariates = character(0),
    base_label = "Null",
    model_label = "Unadjusted"
  )
}

fit_core_nfl_interaction_mainonly <- function(pred,
                                              main_effects = c(core_main_effects,
                                                               core_time_interactions)) {
  if (pred %in% c(core_main_effects, core_time_interactions, "NfL_z")) return(NULL)
  fit_nfl_interaction(
    biomarker_col = pred,
    base_covariates = main_effects,          # all core covariates as main effects
    interaction_covariates = character(0),   # no covariate x time terms
    base_label = "NfL+Core (main effects only)",
    model_label = "NfL+Core, main effects only"
  )
}

# Run all three specifications for each Table 2 biomarker
sensitivity_results_table2 <- purrr::map_dfr(nfl_interaction_predictors, function(pred) {
  dplyr::bind_rows(
    fit_nfl_interaction_unadjusted(pred),
    fit_core_nfl_interaction_mainonly(pred),
    fit_core_nfl_interaction(pred)
  )
})

sensitivity_results_table2 %>%
  dplyr::select(biomarker_label, model_label,
                int_beta, int_beta_lcl, int_beta_ucl, int_p_value,
                n_obs, n_subjects) |> print(n=Inf)


# ============================================================#
# Reviewer Table: Sensitivity of Table 2 interaction estimates
# to covariate specification
# ============================================================#

fmt_beta_ci <- function(beta, lcl, ucl) {
  sprintf("%.2f (%.2f–%.2f)", beta, lcl, ucl)
}

fmt_p_reviewer <- function(p) {
  dplyr::case_when(
    is.na(p)   ~ NA_character_,
    p < 0.001  ~ "<0.001",
    TRUE       ~ sprintf("%.3f", p)
  )
}

table2_sensitivity_display <- sensitivity_results_table2 %>%
  mutate(
    model_key = case_when(
      model_label == "Unadjusted" ~ "Unadjusted",
      model_label == "NfL+Core, main effects only" ~ "Core main effects only",
      model_label == "NfL+Core adjusted" ~ "Core + time interactions",
      TRUE ~ model_label
    ),
    beta_ci = fmt_beta_ci(
      int_beta,
      int_beta_lcl,
      int_beta_ucl
    ),
    p_disp = fmt_p_reviewer(int_p_value),
    
    # Apply manuscript labels
    biomarker_label = purrr::map2_chr(
      biomarker_raw,
      biomarker_label,
      function(raw, current_label) {
        if (raw %in% names(T1label)) {
          as.character(T1label[[raw]])
        } else {
          current_label
        }
      })
    ) %>%
  select(
    biomarker_raw,
    biomarker_label,
    model_key,
    beta_ci,
    p_disp
  ) %>%
  
  # Exclude reviewer-only analyses not shown in Table 2, if necessary
  filter(
    biomarker_raw %in% c(
      "GFAP_z",
      "pTau181_z",
      "pTau217_z",
      "Ab3840_z",
      "Ab4240_z",
      "Ab_status_bl"
    )
  ) %>%
  
  pivot_wider(
    names_from = model_key,
    values_from = c(beta_ci, p_disp),
    names_glue = "{model_key} | {.value}"
  ) %>%
  
  select(
    Biomarker = biomarker_label,
    
    `Unadjusted | beta_ci`,
    `Unadjusted | p_disp`,
    
    `Core main effects only | beta_ci`,
    `Core main effects only | p_disp`,
    
    `Core + time interactions | beta_ci`,
    `Core + time interactions | p_disp`
  )

table2_sensitivity_display


## Flextable ---------------------------------------------------------------


table2_sensitivity_display <- table2_sensitivity_display %>%
  dplyr::mutate(
    Biomarker = dplyr::case_when(
      Biomarker == "Ab3840" ~ "Aβ38/40",
      Biomarker == "Ab4240" ~ "Aβ42/40",
      grepl("^Ab_status", Biomarker) ~ "Aβ status",
      TRUE ~ Biomarker
    )
  )

reviewer_border <- officer::fp_border(
  color = "black",
  width = 1
)

sens_footer <- paste(
  "β values represent the Time × NfL × Biomarker three-way interaction.",
  "The primary model included the core clinical covariates",
  "(age, %VC, ΔFS, onset site, and REEC classification)",
  "together with their interactions with time.", 
  "All three model specifications were fitted to the same 1,421 observations from 218 participants for each biomarker."
)

ft_table2_sensitivity <- flextable::flextable(
  table2_sensitivity_display
) %>%
  
  # Lower-level column labels
  flextable::set_header_labels(
    Biomarker = "Biomarker of interest",
    
    `Unadjusted | beta_ci` = "β (95% CI)",
    `Unadjusted | p_disp` = "P-value",
    
    `Core main effects only | beta_ci` = "β (95% CI)",
    `Core main effects only | p_disp` = "P-value",
    
    `Core + time interactions | beta_ci` = "β (95% CI)",
    `Core + time interactions | p_disp` = "P-value"
  ) %>%
  
  # Upper grouped header
  flextable::add_header_row(
    values = c(
      "",
      "Unadjusted",
      "Core covariates: main effects only",
      "Core covariates + time interactions (primary)"
    ),
    colwidths = c(1, 2, 2, 2)
  ) %>%
  
  # Alignment
  flextable::align(
    align = "center",
    part = "header"
  ) %>%
  flextable::align(
    align = "left",
    j = 1,
    part = "body"
  ) %>%
  flextable::align(
    align = "center",
    j = 2:7,
    part = "body"
  ) %>%
  
  # Emphasis
  flextable::bold(part = "header") %>%
  flextable::bold(j = 1, part = "body") %>%
  
  # Borders
  flextable::border_remove() %>%
  flextable::hline_top(
    border = reviewer_border,
    part = "header"
  ) %>%
  flextable::hline_bottom(
    border = reviewer_border,
    part = "header"
  ) %>%
  
  # Slight separators between model blocks
  # flextable::vline(
  #   j = 3,
  #   border = officer::fp_border(width = 0.5),
  #   part = "all"
  # ) %>%
  # flextable::vline(
  #   j = 5,
  #   border = officer::fp_border(width = 0.5),
  #   part = "all"
  # ) %>%
  
  # Widths
  flextable::width(j = 1, width = 1.55) %>%
  flextable::width(j = c(2, 4, 6), width = 1.35) %>%
  flextable::width(j = c(3, 5, 7), width = 0.70) %>%
  
  # Footer
  flextable::add_footer_lines(
    values = sens_footer
  ) %>%
  flextable::hline_top(
    border = reviewer_border,
    part = "footer"
  ) %>%
  
  # Typography
  flextable::font(
    fontname = "Times New Roman",
    part = "all"
  ) %>%
  flextable::fontsize(
    size = 9,
    part = "all"
  ) %>%
  flextable::padding(
    padding.top = 2,
    padding.bottom = 2,
    padding.left = 3,
    padding.right = 3,
    part = "all"
  ) %>%
  
  # Caption
  flextable::set_caption(
    caption = flextable::as_paragraph(
      flextable::as_chunk(
        "Sensitivity of NfL × biomarker interaction estimates to covariate specification",
        props = officer::fp_text(
          font.size = 12,
          bold = TRUE,
          font.family = "Times New Roman"
        )
      )
    ),
    word_stylename = "Table Caption",
    fp_p = officer::fp_par(
      text.align = "left",
      padding = 5
    )
  )

ft_table2_sensitivity


sensitivity_results_table2 %>%
  dplyr::select(
    biomarker_label,
    model_label,
    n_obs,
    n_subjects
  ) %>%
  print(n = Inf)



# Sensitivity analysis for single-biomarker LME models -------------

# ---- Formatting helpers ----
fmt_beta_ci <- function(beta, lcl, ucl) {
  sprintf("%.2f (%.2f–%.2f)", beta, lcl, ucl)
}

fmt_p_reviewer <- function(p) {
  dplyr::case_when(
    is.na(p)  ~ NA_character_,
    p < 0.001 ~ "<0.001",
    TRUE      ~ sprintf("%.3f", p)
  )
}

# ---- CSF biomarkers included in Supplementary Table S8 ----
csf_biomarkers <- c(
  "NfL",
  "GFAP",
  "pTau181",
  "pTau217",
  "Ab3840",
  "Ab4240",
  "Ab_status (positive vs negative)"
)

# ---- Prepare display table ----
supp8_sensitivity_display <- sensitivity_results_supp8 %>%
  filter(biomarker_label %in% csf_biomarkers) %>%
  mutate(
    biomarker_order = match(biomarker_label, csf_biomarkers),
    
    Biomarker = case_when(
      biomarker_label == "Ab3840" ~ "Aβ38/40",
      biomarker_label == "Ab4240" ~ "Aβ42/40",
      biomarker_label == "Ab_status (positive vs negative)" ~ "Aβ status",
      TRUE ~ biomarker_label
    ),
    
    model_key = case_when(
      model_label == "Unadjusted" ~ "Unadjusted",
      model_label == "Core, main effects only" ~ "Core main effects only",
      model_label == "Core adjusted" ~ "Core + time interactions",
      TRUE ~ model_label
    ),
    
    beta_ci = fmt_beta_ci(beta, beta_lcl, beta_ucl),
    p_disp = fmt_p_reviewer(p_value)
  ) %>%
  arrange(biomarker_order) %>%
  select(
    biomarker_order,
    Biomarker,
    model_key,
    beta_ci,
    p_disp
  ) %>%
  pivot_wider(
    names_from = model_key,
    values_from = c(beta_ci, p_disp),
    names_glue = "{model_key} | {.value}"
  ) %>%
  arrange(biomarker_order) %>%
  
  # IMPORTANT: explicitly specify final order
  select(
    Biomarker,
    `Unadjusted | beta_ci`,
    `Unadjusted | p_disp`,
    `Core main effects only | beta_ci`,
    `Core main effects only | p_disp`,
    `Core + time interactions | beta_ci`,
    `Core + time interactions | p_disp`
  )

supp8_sensitivity_display


reviewer_border <- officer::fp_border(
  color = "black",
  width = 1
)

thin_border <- officer::fp_border(
  color = "black",
  width = 0.5
)

supp8_footer <- paste(
  "β values represent the Time × Biomarker interaction.",
  "All three model specifications were fitted to the same 1,421 observations",
  "from 218 participants for each biomarker.",
  "The primary model included the core clinical covariates",
  "(age, %VC, ΔFS, onset site, and REEC classification)",
  "together with their interactions with time."
)

ft_supp8_sensitivity <- flextable::flextable(
  supp8_sensitivity_display
) %>%
  
  # Lower-level headers
  flextable::set_header_labels(
    Biomarker = "Biomarker",
    
    `Unadjusted | beta_ci` = "β (95% CI)",
    `Unadjusted | p_disp` = "P-value",
    
    `Core main effects only | beta_ci` = "β (95% CI)",
    `Core main effects only | p_disp` = "P-value",
    
    `Core + time interactions | beta_ci` = "β (95% CI)",
    `Core + time interactions | p_disp` = "P-value"
  ) %>%
  
  # Upper grouped header
  flextable::add_header_row(
    values = c(
      "",
      "Unadjusted",
      "Core covariates: main effects only",
      "Core covariates + time interactions (primary)"
    ),
    colwidths = c(1, 2, 2, 2)
  ) %>%
  
  # Alignment
  flextable::align(
    align = "center",
    part = "header"
  ) %>%
  flextable::align(
    align = "left",
    j = 1,
    part = "body"
  ) %>%
  flextable::align(
    align = "center",
    j = 2:7,
    part = "body"
  ) %>%
  
  # Header emphasis
  flextable::bold(part = "header") %>%
  
  # Borders
  flextable::border_remove() %>%
  flextable::hline_top(
    border = reviewer_border,
    part = "header"
  ) %>%
  flextable::hline_bottom(
    border = reviewer_border,
    part = "header"
  ) %>%
  
  # Column widths
  flextable::width(
    j = 1,
    width = 1.30
  ) %>%
  flextable::width(
    j = c(2, 4, 6),
    width = 1.35
  ) %>%
  flextable::width(
    j = c(3, 5, 7),
    width = 0.70
  ) %>%
  
  # Footer
  flextable::add_footer_lines(
    values = supp8_footer
  ) %>%
  flextable::hline_top(
    border = reviewer_border,
    part = "footer"
  ) %>%
  flextable::align(
    align = "left",
    part = "footer"
  ) %>%
  
  # Typography
  flextable::font(
    fontname = "Times New Roman",
    part = "all"
  ) %>%
  flextable::fontsize(
    size = 9,
    part = "all"
  ) %>%
  flextable::padding(
    padding.top = 2,
    padding.bottom = 2,
    padding.left = 3,
    padding.right = 3,
    part = "all"
  ) %>%
  
  # Caption
  flextable::set_caption(
    caption = flextable::as_paragraph(
      flextable::as_chunk(
        paste0(
          "Sensitivity of single-biomarker LME estimates to covariate specification"
        ),
        props = officer::fp_text(
          font.size = 12,
          bold = TRUE,
          font.family = "Times New Roman"
        )
      )
    ),
    word_stylename = "Table Caption",
    fp_p = officer::fp_par(
      text.align = "left",
      padding = 5
    )
  )

ft_supp8_sensitivity


# Export ------------------------------------------------------------------

ft_supp8_sensitivity -> SuppleT9A_pub_added
ft_table2_sensitivity -> SuppleT9B_pub_added

SuppleT9A_pub_added
SuppleT9B_pub_added

# flextable::save_as_docx(
#   ft_table2_sensitivity,
#   path = "output/Tables/ReviewerTable_Sensitivity_Table2.docx",
#   pr_section = officer::prop_section(
#     page_size = officer::page_size(orient = "landscape", width = 11, height = 8.5),
#     page_margins = officer::page_mar(bottom = 0.5, top = 0.5, left = 0.5, right = 0.5)
#   )
# )
# 
# flextable::save_as_docx(
#   ft_supp8_sensitivity,
#   path = "output/Tables/ReviewerTable_Sensitivity_Supp8.docx",
#   pr_section = officer::prop_section(
#     page_size = officer::page_size(orient = "landscape", width = 11, height = 8.5),
#     page_margins = officer::page_mar(bottom = 0.5, top = 0.5, left = 0.5, right = 0.5)
#   )
# )

