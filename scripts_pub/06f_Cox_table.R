# ================== Cox model sensitivity analyses ================== # 

# This script was used to generate Cox model sensitivity analyses for transparency purposes. The raw participant-level data are not
# publicly available due to ethical and privacy restrictions.
# Therefore, this script is not intended to be fully executable without access to the original analysis dataset.

source("scripts_pub/function/Cox_table_func.R")

# Data Prep ===============================
# Dataset expected: survALS
# Outcome columns expected: surv_time_days, surv_time_status (1 = event)

biomarker_names <- c(
  "pTau181_csf","pTau217_csf","GFAP_csf_pgml","NfL_csf_pgml",
  "Ab38_40_csf","Ab42_40_csf_bridged","Cr","CK","Alb","Cr_CysC_ratio")

# log names
ln_biomarker_names <- c(
  "log_pTau217_csf", "log_pTau181_csf",
  "log_Ab38_40_csf",  "log_Ab4240_csf", "Ab_status",
  "log_GFAP_csf", "log_NfL_csf",
  "log_Alb", "log_Cr","log_CK", "log_Cr_CysC_ratio")

  
# Preprocessing ----------------------------
# Build an analysis-ready dataset with pre-transformed biomarkers.
# :: all biomarkers = log-&-Z-scaled
# :: all core clinical predictors = Z-scaled

als_surv_for_screen <- survALS %>%
  mutate(
    surv_event = surv_time_status == "1",
    Sex        = as.factor(Sex),
    Age_init   = as.numeric(scale(Age_init)), 
    deltaFS=as.numeric(scale(deltaFS)),
    VC_Percent = as.numeric(scale(VC_Percent)),
    #all scaled
    across(all_of(ln_biomarker_names[ sapply(survALS[ln_biomarker_names], is.numeric) ]),
           ~ as.numeric(scale(.x)))
  )


# ======================== Run for all predictors ===============================

all_Cox_results <- purrr::map_dfr(ln_biomarker_names, function(pred) {
  dplyr::bind_rows(
    fit_unadjusted(pred),
    fit_age_adj(pred), 
    fit_age_NfL_adj(pred),
    fit_core(pred)
  )
})


## Cox Predictors (clinical) -----------------------------------------------

core_predictors_cox <- 
  c("Age_init", "Sex", 
    "VC_Percent", "deltaFS", "REEC_definite", "OnsetSite",
    "genetic_label",
    "BMI_init",
    "APOE_e4_carrier","APOE_e2_carrier")

all_Cox_results_corevars <- purrr::map_dfr(core_predictors_cox, function(pred) {
  dplyr::bind_rows(
    fit_unadjusted(pred),   # Null vs pred
    fit_age_adj(pred),      # Age vs Age + pred
    fit_age_NfL_adj(pred),  # Age+NfL vs Age+NfL + pred
    fit_core_or_na(pred)    # Core row NA if pred ∈ core_covariates
  )
})

all_Cox_results_combined <- dplyr::bind_rows(
  all_Cox_results,          # biomarker results (unchanged)
  all_Cox_results_corevars  # core predictors with Core=NA row
)


# Cox All table -----------------------------------------------------------

# 1) Build long display with formatted strings
display_long <- all_Cox_results_combined %>%
  mutate(
    hr_disp = ifelse(is.na(hr) | is.na(hr_lcl) | is.na(hr_ucl), NA,
                     sprintf("%.2f (%.2f–%.2f)", hr, hr_lcl, hr_ucl)),
    p_disp  = gtsummary::style_pvalue(p_value, digits = 2),
    c_disp = fmt_ci3(cindex_full, cindex_full_lcl, cindex_full_ucl),
    dc_disp = case_when(
      is.na(delta_cindex) ~ NA_character_,
      delta_cindex < 0 ~ "<0",
      abs(delta_cindex) < 0.001 ~ "~0",
      TRUE ~ fmt_3d(delta_cindex)
    ),
    daic_disp = fmt_2d(delta_aic),  
    model_key = dplyr::case_when(
      base_model == "Null"         ~ "Unadj",
      base_model == "Age"      ~ "Age",
      base_model == "Age+NfL"  ~ "Age+NfL",
      base_model == "Core"         ~ "Core",
      TRUE ~ base_model
    )
  )

# 2) Drop ΔC for Unadj entirely (no empty column in the final table)
display_long <- display_long %>%
  mutate(dc_disp = dplyr::if_else(model_key == "Unadj", NA_character_, dc_disp))

# 3) Long-to-wide with explicit column keys
long_to_wide <- display_long %>%
  dplyr::select(biomarker_label, model_key, hr_disp, p_disp, c_disp, dc_disp, daic_disp) %>%
  tidyr::pivot_longer(
    cols = c(hr_disp, p_disp, c_disp, dc_disp, daic_disp),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(metric = dplyr::recode(
    metric,
    hr_disp = "HR (95% CI)",
    p_disp  = "p-value",
    c_disp  = "C-index (95% CI)",
    dc_disp = "ΔC-index",
    daic_disp = "ΔAIC"  
  )) %>%
  dplyr::filter(!(model_key == "Unadj" & metric %in% c("ΔC-index", "ΔAIC"))) %>%   # ensure the columns do not exist
  dplyr::mutate(final_col = paste(model_key, metric, sep = " | "))

master_wide <- long_to_wide %>%
  dplyr::select(biomarker_label, final_col, value) %>%
  dplyr::distinct() %>%
  tidyr::pivot_wider(names_from = final_col, values_from = value)


block_map <- c("Null" = "Unadj", "Age" = "Age", "Age+NfL" = "Age+NfL", "Core" = "Core")

events_by_block <- all_Cox_results_combined %>%
  dplyr::filter(base_model %in% names(block_map)) %>%
  dplyr::mutate(
    block = unname(block_map[base_model]),
    eventsN = dplyr::if_else(
      is.na(n_events) | is.na(n_analytic),
      NA_character_,
      sprintf("%d / %d", n_events, n_analytic)
    )
  ) %>%
  dplyr::distinct(biomarker_label, block, eventsN) %>%
  tidyr::pivot_wider(
    id_cols = biomarker_label,
    names_from = block,
    values_from = eventsN,
    names_glue = "{block} | Events / N"
  )

# Remove any previous N-join code and use this single join
master_wide <- master_wide %>%
  dplyr::left_join(events_by_block, by = "biomarker_label")



# 4) Build desired column order: block-first (Unadj → AgeNfL → Core), each with HR → p → C → ΔC
block_order  <- c("Unadj", "Age", "Age+NfL", "Core")
metric_order <- c("HR (95% CI)", "p-value", "C-index (95% CI)", "ΔC-index", "ΔAIC")

desired_cols <- "biomarker_label"
for (b in block_order) {
  metrics_this_block <- if (b == "Unadj") metric_order[1:3] else metric_order  # Unadj excludes ΔC and ΔAIC
  # append the block metrics that actually exist
  cols_b <- paste(b, metrics_this_block, sep = " | ")
  cols_b <- cols_b[cols_b %in% names(master_wide)]
  desired_cols <- c(desired_cols, cols_b)
  
  # add Events/N if present for ANY block
  ev_col <- paste0(b, " | Events / N")
  if (ev_col %in% names(master_wide)) {
    desired_cols <- c(desired_cols, ev_col)
  }
}

desired_cols <- desired_cols[desired_cols %in% names(master_wide)]

master_wide <- master_wide %>%
  dplyr::select(dplyr::all_of(desired_cols)) %>%
  dplyr::arrange(biomarker_label) %>%
  dplyr::rename(Biomarker = biomarker_label) %>%
  dplyr::filter(!is.na(Biomarker))



# ───────────────────────────────────────────────────────────────
# Row order for Cox table (core first, then others by Unadj p)
# ───────────────────────────────────────────────────────────────

# 1) Core labels present in the Cox results (from the Unadjusted block)
core_labels_in_results <- all_Cox_results_combined %>%
  dplyr::filter(base_model == "Null", biomarker_raw %in% core_predictors_cox) %>%
  dplyr::distinct(biomarker_raw, biomarker_label)

# Preferred core order
preferred_core_order <- c("Age_init","Sex","VC_Percent","deltaFS","REEC_definite","OnsetSite")

core_order_labels <- core_labels_in_results %>%
  dplyr::mutate(order = match(biomarker_raw, preferred_core_order)) %>%
  dplyr::arrange(order) %>%
  dplyr::pull(biomarker_label)

preferred_biomarker_order <- c(
  "log_NfL_csf", "log_GFAP_csf", "log_pTau181_csf","log_pTau217_csf", 
  "log_Ab38_40_csf","log_Ab4240_csf","Ab_status", "log_Cr", "log_CK", "log_Alb", "log_Cr_CysC_ratio")

# 2) Unadjusted p-values for everything else
unadj_p_by_label <- all_Cox_results_combined %>%
  dplyr::filter(base_model == "Null") %>%
  dplyr::distinct(biomarker_raw, biomarker_label, .keep_all = TRUE) %>%
  dplyr::select(biomarker_raw, biomarker_label, p_value)

# 3) Non-core order: by p-value ascending, NA last, then alphabetical
noncore_order_labels <- unadj_p_by_label %>%
  dplyr::filter(!biomarker_label %in% core_order_labels) %>%
  dplyr::arrange(is.na(p_value), p_value, biomarker_label) %>%
  dplyr::pull(biomarker_label)

# 4) Biomarker labels that correspond to preferred_biomarker_order and exist in the data
preferred_biomarker_labels <- unadj_p_by_label %>%
  dplyr::filter(!(biomarker_raw %in% core_predictors_cox)) %>%
  dplyr::mutate(order = match(biomarker_raw, preferred_biomarker_order)) %>%
  dplyr::filter(!is.na(order)) %>%
  dplyr::arrange(order) %>%
  dplyr::pull(biomarker_label)

# 5) Remaining (not in core, not in preferred list): order by Unadj p asc, NA last, then alphabetical
remaining_biomarker_labels <- unadj_p_by_label %>%
  dplyr::filter(!(biomarker_label %in% c(core_order_labels, preferred_biomarker_labels))) %>%
  dplyr::arrange(is.na(p_value), p_value, biomarker_label) %>%
  dplyr::pull(biomarker_label)

# 6) Final row order = core (specified) → biomarkers (specified) → remaining (fallback)
row_order <- c(core_order_labels, preferred_biomarker_labels, remaining_biomarker_labels)

master_wide <- master_wide %>%
  dplyr::mutate(Biomarker = factor(Biomarker, levels = row_order)) %>%
  dplyr::arrange(Biomarker) %>%
  dplyr::mutate(Biomarker = as.character(Biomarker))



# 7) Flextable with centered numeric blocks
ft_master <- flextable::flextable(master_wide) |>
  flextable::autofit() |>
  flextable::align(align = "center", j = 2:ncol(master_wide)) |>
  flextable::bold(j = 1) |>
  flextable::add_header_lines("Single-biomarker Cox results across four model specifications") |>
  flextable::add_footer_lines(values = c(
    "Unadjusted block omits ΔC-index.",
    "p-values formatted to two decimals.",
    "C-index (95% CI) from summary.coxph.",
    "ΔC-index = C(full) – C(base) for adjusted blocks."
  ))

ft_master

# ──────────────────────────────────────────────────────────────
# Utility: Select columns for Cox model summary tables
# ──────────────────────────────────────────────────────────────
# - Works with base pipe (no RHS `{}` trick).
# - Accepts either "biomarker_label" or "Biomarker" as the ID column.
# - Skips non-existent columns safely; drops ΔC for Unadj automatically.
# ──────────────────────────────────────────────────────────────
# Cox version: focus on Events / N (not just N)
select_cox_columns <- function(master_wide_tbl,
                               blocks    = c("Unadj","Age","Age+NfL","Core"),
                               metrics   = c("HR","p","C","dC", "dAIC"),
                               include_eventsN = c("none","core_only","unadj_only","all")) {
  include_eventsN <- match.arg(include_eventsN)
  
  # Detect ID column
  id_col <- dplyr::case_when(
    "Biomarker"       %in% names(master_wide_tbl) ~ "Biomarker",
    "biomarker_label" %in% names(master_wide_tbl) ~ "biomarker_label",
    TRUE ~ NA_character_
  )
  if (is.na(id_col)) stop("select_cox_columns(): table must have 'Biomarker' or 'biomarker_label'.")
  
  # Metric key -> column suffix
  metric_map <- c(
    HR = "HR (95% CI)",
    p  = "p-value",
    C  = "C-index (95% CI)",
    dC = "ΔC-index",
    dAIC = "ΔAIC"  
  )
  
  # Build desired columns
  desired_cols <- c(id_col)
  for (b in blocks) {
    metrics_this_block <- if (identical(b, "Unadj")) setdiff(metrics, "dC") else metrics
    cols_b <- paste(b, unname(metric_map[metrics_this_block]), sep = " | ")
    cols_b <- cols_b[cols_b %in% names(master_wide_tbl)]
    desired_cols <- c(desired_cols, cols_b)
    
    # Optional Events/N columns by policy
    want_eventsN <- (include_eventsN == "all") ||
      (include_eventsN == "core_only"   && b == "Core")  ||
      (include_eventsN == "unadj_only"  && b == "Unadj")
    
    if (want_eventsN) {
      ncol_events <- paste0(b, " | Events / N")
      if (ncol_events %in% names(master_wide_tbl)) {
        desired_cols <- c(desired_cols, ncol_events)
      }
    }
  }
  
  out <- dplyr::select(master_wide_tbl, dplyr::all_of(desired_cols))
  if (id_col == "biomarker_label") out <- dplyr::rename(out, Biomarker = biomarker_label)
  out
}



# ──────────────────────────────────────────────────────────────
# Helper: strip " (lcl–ucl)" from a "x.xx (l.ll–u.uu)" string
# ──────────────────────────────────────────────────────────────
strip_ci_str <- function(x) {
  ifelse(is.na(x), NA_character_, sub("\\s*\\(.*\\)$", "", x))
}

# ──────────────────────────────────────────────────────────────
# Helper: apply a label map to the Biomarker column
# - Accepts exact matches (e.g., "NfL_csf_pgml")
# - Accepts "log_"-prefixed names (e.g., "log_NfL_csf")
# - Preserves any suffix like " (B vs A)" for factors
# ──────────────────────────────────────────────────────────────
apply_labels <- function(df, label_map) {
  if (!"Biomarker" %in% names(df)) {
    stop("apply_labels(): expected a 'Biomarker' column.")
  }
  transform_one <- function(b) {
    if (is.na(b) || !nzchar(b)) return(b)
    base <- sub(" \\(.*$", "", b)                # before any " ( ... )"
    base <- sub("^log_", "", base)               # drop log_ if present
    has_suffix <- grepl(" \\(", b)
    suffix <- if (has_suffix) sub("^[^\\(]+", "", b) else ""
    if (base %in% names(label_map)) {
      mapped <- unname(unlist(label_map[base]))
      paste0(mapped, suffix)
    } else {
      b
    }
  }
  df$Biomarker <- vapply(df$Biomarker, transform_one, character(1))
  df
}

# ──────────────────────────────────────────────────────────────
# Main: build a two-tier flextable with optional label mapping
# ──────────────────────────────────────────────────────────────

build_master_flextable <- function(master_wide_tbl,
                                   include_p = c("none","agesex_only","all"),
                                   cindex_ci = TRUE,
                                   core_covariates = c("Age_init", "VC_percent","deltaFS","REEC_definite","OnsetSite"),
                                   label_map = NULL) {
  include_p <- match.arg(include_p)
  
  # Ensure the ID column is 'Biomarker'
  if ("biomarker_label" %in% names(master_wide_tbl) && !"Biomarker" %in% names(master_wide_tbl)) {
    master_wide_tbl <- dplyr::rename(master_wide_tbl, Biomarker = biomarker_label)
  }
  if (!"Biomarker" %in% names(master_wide_tbl)) {
    stop("build_master_flextable(): expected a 'Biomarker' column.")
  }
  
  # Optionally drop p-value columns
  keep_cols <- names(master_wide_tbl)
  if (include_p == "none") {
    keep_cols <- keep_cols[!grepl("\\| p-value$", keep_cols)]
  } else if (include_p == "agesex_only") {
    keep_cols <- keep_cols[!grepl("\\| p-value$", keep_cols) | grepl("^Age \\| p-value$", keep_cols)]
  }
  master_use <- master_wide_tbl[, keep_cols, drop = FALSE]
  
  # Optionally switch C-index columns to point estimates only and rename headers
  if (!isTRUE(cindex_ci)) {
    cidx_cols <- grep("\\| C-index \\(95% CI\\)$", names(master_use), value = TRUE)
    if (length(cidx_cols)) {
      for (cc in cidx_cols) master_use[[cc]] <- strip_ci_str(master_use[[cc]])
      names(master_use) <- sub("\\| C-index \\(95% CI\\)$", "| C-index", names(master_use))
    }
  }
  
  # Apply label mapping after column filtering but before header building
  if (!is.null(label_map)) {
    master_use <- apply_labels(master_use, label_map)
  }
  
  # Build two-level header mapping (top = Block, bottom = Metric)
  col_keys <- names(master_use)
  top_lab  <- character(length(col_keys))
  bot_lab  <- character(length(col_keys))
  top_lab[1] <- "";         bot_lab[1] <- "Biomarker"
  if (length(col_keys) > 1) {
    top_lab[-1] <- sub("^(.*?) \\| .*", "\\1", col_keys[-1])
    bot_lab[-1] <- sub("^.*? \\| (.*)$", "\\1", col_keys[-1])
  }
  header_df <- data.frame(col_keys = col_keys, Block = top_lab, Metric = bot_lab, stringsAsFactors = FALSE)
  
  # Build flextable
  ft <- flextable::flextable(master_use, col_keys = col_keys)
  ft <- flextable::set_header_df(ft, mapping = header_df, key = "col_keys")
  ft <- flextable::merge_h(ft, part = "header")
  ft <- flextable::merge_v(ft, part = "header")
  ft <- flextable::align(ft, align = "center", part = "header")
  ft <- flextable::bold(ft, part = "header")
  ft <- flextable::autofit(ft)
  ft <- flextable::align(ft, align = "center", j = 2:ncol(master_use))
  ft <- flextable::bold(ft, j = 1)
  
  # Header border styling
  brd <- officer::fp_border(color = "black", width = 1)
  ft  <- flextable::border_remove(ft)
  ft  <- flextable::hline(ft, i = 1, border = brd, part = "header")
  ft  <- flextable::hline(ft, i = 2, border = brd, part = "header")
  ft  <- flextable::hline_top(ft, border = brd, part = "header")
  ft  <- flextable::hline_bottom(ft, border = brd, part = "header")
  
  # Footer lines
  
  model_line <- paste(
    "Cox proportional hazards models were used with time from baseline to the endpoint (death or tracheostomy) as the outcome.",
    "Hazard ratios (HRs) for continuous predictors correspond to a 1-SD increase in log-transformed biomarker values (where applicable);",
    "Sample sizes differ across models because complete-case analyses were performed for each combination of variables.",
    sep = " "
  )
  
  delta_aic_line <- "ΔAIC was defined as AIC [full model with X] − AIC [corresponding base model with the same covariates but without X]; values < −2 indicate a meaningful improvement in model fit."
  
  abbrev_line <- "Aβ, amyloid-β; CI, confidence interval; C-index, Harrell’s concordance index; HR, hazard ratio; GFAP, glial fibrillary acidic protein; NfL, neurofilament light"
  cindex_line <- if (isTRUE(cindex_ci)) {
    "C-index values are shown with 95% CIs; higher values indicate better discrimination."
    #ΔC-index was defined as C[full model including the biomarker] − C[corresponding base model with the same covariates but without the biomarker]; positive values indicate improved discrimination."
  } else {
    "C-index values are shown as point estimates; higher values indicate better discrimination."
    #ΔC-index was defined as C[full model including the biomarker] − C[corresponding base model with the same covariates but without the biomarker]; positive values indicate improved discrimination."
  }
  core_line <- sprintf("Core clinical features included as covariates in the fully adjusted (“Core”) models were: %s.", paste(core_covariates, collapse = ", "))
  
  ft <- flextable::add_footer_lines(ft, values = c(model_line, cindex_line, delta_aic_line,
                                                   core_line, abbrev_line))
  ft <- flextable::hline_top(ft, border = brd, part = "footer")
  
  ft
}


# RESULTS TABLES ==================================================

tbl_cox_custom <- select_cox_columns(
  master_wide_tbl = master_wide,
  blocks  = c("Unadj","Age+NfL","Core"),
  metrics = c("HR","p","C","dAIC"),  
  include_eventsN = "all"
)

ft_master_all <- build_master_flextable(
  tbl_cox_custom,
  include_p = "all",
  cindex_ci = FALSE,
  core_covariates = core_covariates,
  label_map = T1label
)

ft_master_all


# 🔍Supple T9 ---------------------------------------------------------------

SuppleT9_pub <- ft_master_all %>% 
  flextable::width(j = 1:ncol_keys(.), width = 0.55) %>% 
  flextable::width(j = 1, width = 1.5) %>%
  flextable::width(j = c(2,6,11), width = 1.2) %>%
  flextable::padding(padding.left = 0, padding.right = 0, part = "all") %>% 
  set_caption(caption = flextable::as_paragraph(
    flextable::as_chunk("Supplementary Table S9. Sensitivity analyses of baseline predictors for survival in ALS",
                        props = officer::fp_text(font.size   = 12, bold = TRUE, 
                                                 font.family = "Times New Roman"))),
    word_stylename = "Table Caption", 
    fp_p = officer::fp_par(text.align = "left", padding= 5)) %>%
  flextable::font(fontname = "Times New Roman", part = "all") %>%
  flextable::fontsize(size=9, part="all") 

SuppleT9_pub



# ===== Minimal add-on: median survival and median follow-up =====

# library(survival)
# 
# d2m <- function(d) d / 30.44
# 
# km_median <- function(fit) {
#   tb <- summary(fit)$table
#   if (is.null(dim(tb))) {
#     med <- suppressWarnings(as.numeric(tb["median"]))
#     lcl <- suppressWarnings(as.numeric(tb["0.95LCL"]))
#     ucl <- suppressWarnings(as.numeric(tb["0.95UCL"]))
#   } else {
#     med <- suppressWarnings(as.numeric(tb[,"median"]))
#     lcl <- suppressWarnings(as.numeric(tb[,"0.95LCL"]))
#     ucl <- suppressWarnings(as.numeric(tb[,"0.95UCL"]))
#   }
#   list(
#     median = ifelse(length(med) == 0, NA_real_, med),
#     lcl    = ifelse(length(lcl) == 0, NA_real_, lcl),
#     ucl    = ifelse(length(ucl) == 0, NA_real_, ucl)
#   )
# }
# 
# # Use a unique variable name to avoid namespace conflicts
# dat_surv <- als_surv_for_screen
# 
# # Median survival
# fit_surv <- survfit(Surv(surv_time_days, surv_event) ~ 1, data = dat_surv)
# ms <- km_median(fit_surv)
# 
# # Median follow-up (reverse KM)
# fit_fu <- survfit(Surv(surv_time_days, !surv_event) ~ 1, data = dat_surv)
# mf <- km_median(fit_fu)
# 
# report_line <- function(label, med, lcl, ucl) {
#   if (is.na(med)) {
#     sprintf("%s: not reached (95%% CI %.1f to %.1f months)",
#             label, d2m(lcl), d2m(ucl))
#   } else {
#     sprintf("%s: %.1f months (95%% CI %.1f to %.1f)",
#             label, d2m(med), d2m(lcl), d2m(ucl))
#   }
# }
# 
# cat(
#   report_line("Median follow up", mf$median, mf$lcl, mf$ucl), "\n",
#   report_line("Median survival",  ms$median, ms$lcl, ms$ucl), "\n",
#   sep = ""
# )
# 
