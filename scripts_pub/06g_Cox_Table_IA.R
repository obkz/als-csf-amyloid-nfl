library(dplyr)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ #
# Cox NfL Interaction Analysis   #
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ #

source("scripts_pub/function/Cox_table_IA_func.R")

# Run interaction analysis for Cox models --------------------------

# Define interaction test variables following the slope analysis order
cox_interaction_predictors <- c(
  # Core CSF biomarkers
  "log_NfL_csf",     # Will be excluded automatically
  "log_GFAP_csf", 
  "log_pTau181_csf",
  "log_pTau217_csf",
  # Amyloid markers
  "log_Ab38_40_csf",
  "log_Ab4240_csf",
  "Ab_status")

# Remove NfL from the list
cox_interaction_predictors <- setdiff(cox_interaction_predictors, "log_NfL_csf")

# Run for biomarkers only
all_cox_interaction_results <- purrr::map_dfr(cox_interaction_predictors, function(pred) {
  dplyr::bind_rows(
    fit_nfl_age_cox_interaction(pred),
    fit_nfl_core_cox_interaction(pred)
  )
})

# Format for display ------------------------

display_cox_interaction_long <- all_cox_interaction_results %>%
  dplyr::mutate(
    # Cap extreme values before formatting
    main_hr = cap_hr(main_hr),
    main_hr_lcl = cap_hr(main_hr_lcl),
    main_hr_ucl = cap_hr(main_hr_ucl),
    interaction_hr = cap_hr(interaction_hr),
    interaction_hr_lcl = cap_hr(interaction_hr_lcl),
    interaction_hr_ucl = cap_hr(interaction_hr_ucl),
    nfl_hr = cap_hr(nfl_hr),
    nfl_hr_lcl = cap_hr(nfl_hr_lcl),
    nfl_hr_ucl = cap_hr(nfl_hr_ucl),
    # Format display strings
    main_hr_disp = fmt_ci2(main_hr, main_hr_lcl, main_hr_ucl),
    main_p_disp = fmt_p2(main_p_value),
    interaction_hr_disp = fmt_ci2(interaction_hr, interaction_hr_lcl, interaction_hr_ucl),
    interaction_p_disp = fmt_p2(interaction_p_value),
    nfl_hr_disp = fmt_ci2(nfl_hr, nfl_hr_lcl, nfl_hr_ucl),
    nfl_p_disp = fmt_p2(nfl_p_value),
    lr_p_disp = fmt_p2(lr_test_p),
    daic_disp = fmt_2d(delta_aic),
    dbic_disp = fmt_2d(delta_bic),
    c_disp = fmt_3d(cindex_interaction),
    dc_disp = dplyr::case_when(
      is.na(delta_cindex) ~ NA_character_,
      delta_cindex < 0 ~ "<0",
      abs(delta_cindex) < 0.001 ~ "~0",
      TRUE ~ fmt_3d(delta_cindex)
    ),
    model_key = dplyr::case_when(
      base_model == "Age+NfL" ~ "Age+NfL",
      base_model == "NfL+Core" ~ "NfL+Core",
      TRUE ~ base_model
    )
  )

# Pivot to wide format -----------------------
cox_interaction_long_to_wide <- display_cox_interaction_long %>%
  dplyr::select(biomarker_raw, biomarker_label, model_key,
                main_hr_disp, main_p_disp, interaction_hr_disp, interaction_p_disp,
                nfl_hr_disp, nfl_p_disp, lr_p_disp,
                daic_disp, dbic_disp, c_disp, dc_disp) %>%
  tidyr::pivot_longer(
    cols = c(main_hr_disp, main_p_disp, interaction_hr_disp, interaction_p_disp,
             nfl_hr_disp, nfl_p_disp, lr_p_disp,
             daic_disp, dbic_disp, c_disp, dc_disp),
    names_to = "metric", values_to = "value"
  ) %>%
  dplyr::mutate(metric = dplyr::recode(metric,
                                       main_hr_disp = "Main HR (95% CI)",
                                       main_p_disp = "Main p-value",
                                       interaction_hr_disp = "Interaction HR (95% CI)",
                                       interaction_p_disp = "Interaction p-value",
                                       nfl_hr_disp = "NfL HR (95% CI)",
                                       nfl_p_disp = "NfL p-value",
                                       lr_p_disp = "LR test p-value",
                                       daic_disp = "ΔAIC",
                                       dbic_disp = "ΔBIC",
                                       c_disp = "C-index",
                                       dc_disp = "ΔC-index"
  )) %>%
  dplyr::mutate(final_col = paste(model_key, metric, sep = " | "))

master_cox_interaction_wide <- cox_interaction_long_to_wide %>%
  dplyr::select(biomarker_raw, biomarker_label, final_col, value) %>%
  dplyr::distinct() %>%
  tidyr::pivot_wider(names_from = final_col, values_from = value) %>%
  dplyr::rename(Biomarker = biomarker_label)


# add Events/N columns
master_cox_interaction_wide <- add_N_columns_cox_interaction(master_cox_interaction_wide, 
                                                             all_cox_interaction_results, 
                                                             where = "all")

# Column ordering
block_order_cox_interaction <- c("Age+NfL", "NfL+Core")
metric_order_cox_interaction <- c("Main HR (95% CI)", "Main p-value", 
                                  "Interaction HR (95% CI)", "Interaction p-value",
                                  "NfL HR (95% CI)", "NfL p-value", 
                                  "LR test p-value",
                                  "ΔAIC", "ΔBIC", "C-index", "ΔC-index")

desired_cols_cox_interaction <- "Biomarker"
for (b in block_order_cox_interaction) {
  cols_b <- paste(b, metric_order_cox_interaction, sep = " | ")
  cols_b <- cols_b[cols_b %in% names(master_cox_interaction_wide)]
  desired_cols_cox_interaction <- c(desired_cols_cox_interaction, cols_b)
  
  ncol_name <- paste0(b, " | Events/N")
  if (ncol_name %in% names(master_cox_interaction_wide)) {
    desired_cols_cox_interaction <- c(desired_cols_cox_interaction, ncol_name)
  }
}

master_cox_interaction_wide <- master_cox_interaction_wide %>%
  dplyr::select(dplyr::all_of(desired_cols_cox_interaction[desired_cols_cox_interaction %in% names(master_cox_interaction_wide)]))

# Row ordering (biomarkers only - no core variables in rows)
# Define preferred order matching slope analysis
preferred_biomarker_order_cox <- c(
  "log_GFAP_csf",
  "log_pTau181_csf",
  "log_pTau217_csf",
  "log_Ab38_csf",
  "log_Ab40_csf",
  "log_Ab42_csf",
  "log_Ab38_40_csf",
  "log_Ab4240_csf",
  "Ab_status",
  "log_Cr",
  "log_CK", 
  "log_Alb",
  "log_Cr_CysC_ratio"
)

preferred_biomarker_cox_interaction <- intersect(preferred_biomarker_order_cox, cox_interaction_predictors)
unadj_p_by_label_cox_interaction <- all_cox_interaction_results %>%
  dplyr::filter(base_model == "Age+NfL") %>%
  dplyr::distinct(biomarker_raw, biomarker_label, .keep_all = TRUE) %>%
  dplyr::select(biomarker_raw, biomarker_label, interaction_p_value)

preferred_biomarker_labels_cox_interaction <- unadj_p_by_label_cox_interaction %>%
  dplyr::mutate(order = match(biomarker_raw, preferred_biomarker_cox_interaction)) %>%
  dplyr::filter(!is.na(order)) %>%
  dplyr::arrange(order) %>%
  dplyr::pull(biomarker_label)

remaining_biomarker_labels_cox_interaction <- unadj_p_by_label_cox_interaction %>%
  dplyr::filter(!(biomarker_label %in% preferred_biomarker_labels_cox_interaction)) %>%
  dplyr::arrange(is.na(interaction_p_value), interaction_p_value, biomarker_label) %>%
  dplyr::pull(biomarker_label)

row_order_cox_interaction <- c(preferred_biomarker_labels_cox_interaction, 
                               remaining_biomarker_labels_cox_interaction)

master_cox_interaction_wide <- master_cox_interaction_wide %>%
  dplyr::mutate(Biomarker = factor(Biomarker, levels = row_order_cox_interaction)) %>%
  dplyr::arrange(Biomarker) %>%
  dplyr::mutate(Biomarker = as.character(Biomarker)) %>%
  dplyr::filter(!is.na(Biomarker) & nzchar(Biomarker))



# Create final tables with flexibility ----------------------------------

Cox_test_tbl <- select_cox_interaction_columns(
  master_wide_tbl = master_cox_interaction_wide,
  blocks = c("NfL+Core"),
  metrics = c("NfLHR",#"MainHR",
              "InteractionHR","InteractionP","Cindex","DeltaC","DeltaAIC","DeltaBIC"),
  include_N = "all"
) 

# Main table with key metrics only
ft_cox_interaction_main <- build_cox_interaction_flextable(
  Cox_test_tbl,
  include_p = "all",
  include_nfl = "all",
  include_events = "all",
  metrics_to_show = c("Cindex","InteractionP","DeltaAIC"),
  label_map = T1label
)


# 🔍Supple T10 ---------------------------------------------------------------

SuppleT10_pub <- ft_cox_interaction_main %>% 
  set_caption(caption = flextable::as_paragraph(
    flextable::as_chunk("Supplementary Table S10. CSF NfL × biomarker interactions for survival in ALS",
                        props = officer::fp_text(font.size   = 12, bold = TRUE, 
                                                 font.family = "Times New Roman"))),
    word_stylename = "Table Caption", 
    fp_p = officer::fp_par(text.align = "left", padding= 5)) %>%
  flextable::fontsize(size=9, part="all") %>% 
  flextable::padding(padding.left = 0, padding.right = 0, part = "all") 


SuppleT10_pub

Cox_test_tbl %>% tibble()
