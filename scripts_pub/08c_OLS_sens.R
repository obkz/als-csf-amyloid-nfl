# ============================================================ #
# Sensitivity: slope definition
# ============================================================ #

library(dplyr)
library(tidyr)
library(purrr)

# Prepare -----------------------------------------------------------------

predictors_reg <- c(
  #CLINICAL FEATURES
  "Age_init", "Sex", "OnsetSite", "REEC_definite", "deltaFS", "VC_Percent",
  #CSF BIOMARKERS
  "NfL_csf_pgml", "GFAP_csf_pgml", "pTau181_csf", "pTau217_csf",
  "Ab38_40_csf","Ab42_40_csf_bridged", "Ab_status",
  #BLOOD BIOMARKERS
  "Cr", "CK" 
)



# ---- Make OLS slopes under arbitrary criteria ----
make_slope_data <- function(min_visits = 3, cutoff = 20) {
  
  v <- visitdf %>%
    filter(
      !visit_post_event,
      !SampleID %in% NotValidSample,
      !is.na(ALSFRSR_Total),
      ALS_label == "1"
    ) %>%
    mutate(
      DaysFromFirstVisit_months = DaysFromFirstVisit / 30.44
    ) %>%
    filter(DaysFromFirstVisit_months <= 24)
  
  # NULL = no baseline ALSFRS-R cutoff
  if (!is.null(cutoff)) {
    v <- v %>% filter(ALSFRS_init > cutoff)
  }
  
  slope_df_sens <- v %>%
    group_by(SSBS_ID) %>%
    filter(n() >= min_visits) %>%
    summarise(
      n_timepoints = n(),
      slope_total = tryCatch(
        coef(lm(ALSFRSR_Total ~ DaysFromFirstVisit_months))[2],
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    filter(is.finite(slope_total))
  
  # Join to baseline CSF/predictor data
  df_sens <- csfdf %>%
    filter(
      !SampleID %in% second,
      !SampleID %in% NotValidSample
    ) %>%
    inner_join(slope_df_sens, by = "SSBS_ID") %>%
    select(
      slope_total,
      all_of(predictors_reg)
    ) %>%
    
    # Same general scaling strategy as primary analysis
    mutate(
      across(
        where(is.numeric) & !any_of(c("slope_total", "Ab_status")),
        ~ as.numeric(scale(.x))
      )
    ) %>%
    
    # Same complete-case requirement as ML analysis
    drop_na()
  
  df_sens
}

count_slope_sample <- function(min_visits = 3, cutoff = 19) {
  
  x <- visitdf %>%
    filter(
      !visit_post_event,
      !SampleID %in% NotValidSample,
      !is.na(ALSFRSR_Total),
      ALS_label == "1"
    ) %>%
    mutate(DaysFromFirstVisit_months = DaysFromFirstVisit / 30.44) %>%
    filter(DaysFromFirstVisit_months <= 24)
  
  if (!is.null(cutoff)) {
    x <- x %>% filter(ALSFRS_init > cutoff)
  }
  
  x %>%
    count(SSBS_ID) %>%
    filter(n >= min_visits) %>%
    summarise(N = n()) %>%
    pull(N)
}


# sensitivity dataset -----------------------------------------------------


sens_datasets <- list(
  primary      = make_slope_data(min_visits = 3, cutoff = 19),
  visits2      = make_slope_data(min_visits = 2, cutoff = 19),
  no_cutoff    = make_slope_data(min_visits = 3, cutoff = NULL),
  both_relaxed = make_slope_data(min_visits = 2, cutoff = NULL)
)

purrr::imap_dfr(
  sens_datasets,
  ~ tibble(
    scenario = .y,
    N = nrow(.x),
    slope_mean = mean(.x$slope_total),
    slope_sd = sd(.x$slope_total)
  )
)

tibble(
  scenario = c("primary", "visits2", "no_cutoff", "both_relaxed"),
  N_slope = c(
    count_slope_sample(3, 19),
    count_slope_sample(2, 19),
    count_slope_sample(3, NULL),
    count_slope_sample(2, NULL)
  ),
  N_ML = c(180, 216, 180, 217)
)

isTRUE(all.equal(
  sens_datasets$primary,
  sens_datasets$no_cutoff
))

# scenario     N_slope  N_ML
# 1 primary          182   180
# 2 visits2          218   216
# 3 no_cutoff        182   180
# 4 both_relaxed     219   217




# ============================================================ #
# Reviewer #3: prediction sensitivity analysis
# Run only the models needed for comparison
# ============================================================ #

# ---- GLM specifications to retain ----

all_specs <- expand_specs4(
  core_sets,
  bio_blood,
  bio_csf,
  bio_one_any,
  bio_one_csf
)

sens_spec_names <- c(
  "Core5",
  "Core5_NfL",
  "Core5_NfL_GFAP_csf_pgml_xGFAP_csf_pgml",
  "Core5_NfL_pTau181_csf_xpTau181_csf",
  "Core5_NfL_pTau217_csf_xpTau217_csf",
  "Core5_NfL_Ab38_40_csf_xAb38_40_csf",
  "Core5_NfL_Ab42_40_csf_bridged_xAb42_40_csf_bridged"
)

# sanity check
stopifnot(all(sens_spec_names %in% names(all_specs)))

sens_specs <- all_specs[sens_spec_names]



fig4d_glm_names <- c(
  "Core1",
  "Core2",
  "Core3",
  "Core4",
  "Core5",
  "Core5_NfL",
  "Core5_NfL_GFAP_csf_pgml_xGFAP_csf_pgml",
  "Core5_NfL_pTau181_csf_xpTau181_csf",
  "Core5_NfL_pTau217_csf_xpTau217_csf",
  "Core5_NfL_Ab38_40_csf_xAb38_40_csf",
  "Core5_NfL_Ab42_40_csf_bridged_xAb42_40_csf_bridged"
)

fig4d_specs <- all_specs[fig4d_glm_names]

stopifnot(length(fig4d_specs) == 11)
names(fig4d_specs)


# Run one sensitivity scenario ---------------------

# run_slope_prediction <- function(dat,
#                                  seed = 12345,
#                                  n_repeats = 10,
#                                  K_outer = 5,
#                                  inner_k = 5) {
#   
#   set.seed(seed)
#   
#   y <- dat$slope_total
#   
#   # Same outer folds shared by all models within each scenario
#   outer_folds <- lapply(seq_len(n_repeats), function(r) {
#     caret::createFolds(
#       y,
#       k = K_outer,
#       list = TRUE,
#       returnTrain = FALSE
#     )
#   })
#   
#   # ML: same five models as primary analysis ---------
# 
#   ml_methods <- c("glm", "rf", "svmRadial", "knn", "xgb")
#   
#   ml_res <- purrr::map_dfr(ml_methods, function(m) {
#     
#     cat("\nML:", m, "\n")
#     
#     res <- nested_fit_regression_repeated(
#       data = dat,
#       outcome = "slope_total",
#       predictors = predictors_reg %>% setdiff("Ab_status"),
#       method = m,
#       K_outer = K_outer,
#       inner_k = inner_k,
#       repeats = n_repeats,
#       seed = seed,
#       scale = "auto",
#       outer_folds_repeats = outer_folds
#     )
#     
#     dplyr::bind_cols(
#       tibble::tibble(model = m),
#       res$mean_sd
#     )
#   })
#   
#   
#   # Selected GLMs only ---------------
# 
#   cat("\nSelected GLM subset models\n")
#   
#   glm_res <- run_glm_subsets(
#     data = dat,
#     outcome = "slope_total",
#     specs = sens_specs,
#     repeats = n_repeats,
#     K_outer = K_outer,
#     inner_k = inner_k,
#     seed = seed,
#     scale = "auto",
#     outer_folds_rep = outer_folds
#   )$flat_mean_sd
#   
#   
#   list(
#     N = nrow(dat),
#     ml = ml_res,
#     glm = glm_res
#   )
# }

run_fig4d_sensitivity <- function(dat,
                                  seed = 12345,
                                  n_repeats = 10,
                                  K_outer = 5,
                                  inner_k = 5) {
  
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  options(contrasts = c("contr.treatment", "contr.poly"))
  
  y <- dat$slope_total
  
  # Same procedure as primary analysis
  outer_folds <- lapply(seq_len(n_repeats), function(r) {
    caret::createFolds(
      y,
      k = K_outer,
      list = TRUE,
      returnTrain = FALSE
    )
  })
  
  # ML models shown in Fig. 4D --------------------------------------------------
  
  ml_methods <- c("glm", "rf", "svmRadial", "knn", "xgb")
  
  ml_results <- purrr::map_dfr(ml_methods, function(m) {
    
    cat("\nML:", m, "\n")
    
    res <- nested_fit_regression_repeated(
      data = dat,
      outcome = "slope_total",
      predictors = predictors_reg %>% setdiff("Ab_status"),
      method = m,
      K_outer = K_outer,
      inner_k = inner_k,
      repeats = n_repeats,
      seed = seed,
      scale = "auto",
      outer_folds_repeats = outer_folds
    )
    
    res$rep_stats %>%
      dplyr::mutate(model = m, .before = 1)
  })
  
  # GLMs shown in Fig. 4D only --------------------------------------------------
  
  glm_results <- run_glm_subsets(
    data = dat,
    outcome = "slope_total",
    specs = fig4d_specs,
    repeats = n_repeats,
    K_outer = K_outer,
    inner_k = inner_k,
    seed = seed,
    scale = "auto",
    outer_folds_rep = outer_folds
  )
  
  glm_rep <- glm_results$nested %>%
    dplyr::select(model, rep_stats) %>%
    tidyr::unnest(rep_stats)
  
  list(
    N = nrow(dat),
    ml_rep = ml_results,
    glm_rep = glm_rep
  )
}


# ============================================================ #
# RUN ---------------------------------------------------------
# ============================================================ #

# res_visits2 <- run_fig4d_sensitivity(
#   sens_datasets$visits2
# )
# 
# saveRDS(
#   res_visits2,
#   "output/reviewer5_visits2_reg_ML.rds"
# )
# 
# 
# res_both_relaxed <- run_fig4d_sensitivity(
#   sens_datasets$both_relaxed
# )
# 
# saveRDS(
#   res_both_relaxed,
#   "output/reviewer5_both_relaxed_reg_ML.rds"
# )


# ============================================================ #
# LOAD ---------------------------------------------------------
# ============================================================ #

# LOAD .rds

readRDS("output/reviewer5_visits2_reg_ML.rds") -> res_visits2
readRDS("output/reviewer5_both_relaxed_reg_ML.rds") -> res_both_relaxed


summarise_cv_r2 <- function(res) {
  
  bind_rows(
    res$glm_rep %>%
      group_by(model) %>%
      summarise(
        R2_mean = mean(R2),
        R2_sd = sd(R2),
        .groups = "drop"
      ) %>%
      mutate(type = "GLM"),
    
    res$ml_rep %>%
      group_by(model) %>%
      summarise(
        R2_mean = mean(R2),
        R2_sd = sd(R2),
        .groups = "drop"
      ) %>%
      mutate(type = "ML")
  )
}


summarise_cv_r2(res_visits2)
summarise_cv_r2(res_both_relaxed)


# original ----------------------------------------------------------------

original <- readRDS(
  "output/Slope_Models_Trained_20260619_1417.rds"
)

# sanity check
original$metadata$n_samples
# should be 180


# Models actually shown in Fig. 4D
fig4d_model_order <- c(
  "Core1",
  "Core2",
  "Core3",
  "Core4",
  "Core5",
  "Core5_NfL",
  "Core5_NfL_Ab38_40_csf_xAb38_40_csf",
  "Core5_NfL_Ab42_40_csf_bridged_xAb42_40_csf_bridged",
  "Core5_NfL_GFAP_csf_pgml_xGFAP_csf_pgml",
  "Core5_NfL_pTau181_csf_xpTau181_csf",
  "Core5_NfL_pTau217_csf_xpTau217_csf",
  "glm",
  "knn",
  "rf",
  "svmRadial",
  "xgb"
)


# GLM results
original_glm <- original$glm_results$nested %>%
  dplyr::filter(model %in% fig4d_model_order) %>%
  dplyr::select(model, rep_stats) %>%
  tidyr::unnest(rep_stats) %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(
    R2_mean = mean(R2),
    R2_sd   = sd(R2),
    .groups = "drop"
  ) %>%
  dplyr::mutate(type = "GLM")


# ML results
original_ml <- original$ml_results %>%
  dplyr::select(model, rep_stats) %>%
  tidyr::unnest(rep_stats) %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(
    R2_mean = mean(R2),
    R2_sd   = sd(R2),
    .groups = "drop"
  ) %>%
  dplyr::mutate(type = "ML")


# Same format/order as sensitivity results
original_cv_r2 <- dplyr::bind_rows(
  original_glm,
  original_ml
) %>%
  dplyr::mutate(
    model = factor(model, levels = fig4d_model_order)
  ) %>%
  dplyr::arrange(model) %>%
  dplyr::mutate(model = as.character(model))


original_cv_r2




comparison_cv_r2 <- dplyr::bind_rows(
  original_cv_r2 %>%
    dplyr::mutate(
      scenario = "Primary: ≥3 visits, ALSFRS-R ≥20",
      .before = 1
    ),
  
  summarise_cv_r2(res_visits2) %>%
    dplyr::mutate(
      scenario = "≥2 visits, ALSFRS-R ≥20",
      .before = 1
    ),
  
  summarise_cv_r2(res_both_relaxed) %>%
    dplyr::mutate(
      scenario = "≥2 visits, no ALSFRS-R cutoff",
      .before = 1
    )
)

comparison_cv_r2 |> print(n=Inf)



# Table for Reviewer ------------------------------------------------------


comparison_cv_r2_wide <- comparison_cv_r2 %>%
  dplyr::select(scenario, model, R2_mean) %>%
  tidyr::pivot_wider(
    names_from = scenario,
    values_from = R2_mean
  )

comparison_cv_r2_wide



# Sample-size information ====================================

scenario_meta <- tibble::tibble(
  scenario = c(
    "primary",
    "visits2",
    "no_cutoff",
    "both_relaxed"
  ),
  Definition = c(
    "≥3 visits, baseline ALSFRS-R ≥20 (primary)",
    "≥2 visits, baseline ALSFRS-R ≥20",
    "≥3 visits, no baseline ALSFRS-R cutoff",
    "≥2 visits, no baseline ALSFRS-R cutoff"
  ),
  N_slope = c(182, 218, 182, 219),
  N_ML    = c(180, 216, 180, 217)
)


# Combine original + sensitivity results =======================

comparison_cv_r2 <- dplyr::bind_rows(
  
  # Original primary analysis
  original_cv_r2 %>%
    dplyr::mutate(scenario = "primary"),
  
  # ≥2 visits
  summarise_cv_r2(res_visits2) %>%
    dplyr::mutate(scenario = "visits2"),
  
  # No cutoff with ≥3 visits:
  # identical sample to primary, therefore identical results
  original_cv_r2 %>%
    dplyr::mutate(scenario = "no_cutoff"),
  
  # ≥2 visits + no cutoff
  summarise_cv_r2(res_both_relaxed) %>%
    dplyr::mutate(scenario = "both_relaxed")
) %>%
  dplyr::left_join(
    scenario_meta,
    by = "scenario"
  ) %>%
  dplyr::select(
    scenario,
    Definition,
    N_slope,
    N_ML,
    model,
    R2_mean,
    R2_sd,
    type
  )

comparison_cv_r2 |> print(n=Inf)





summarise_alsfrs_sample <- function(min_visits = 3, cutoff = 19) {
  
  x <- visitdf %>%
    filter(
      !visit_post_event,
      !SampleID %in% NotValidSample,
      !is.na(ALSFRSR_Total),
      ALS_label == "1"
    ) %>%
    mutate(
      DaysFromFirstVisit_months = DaysFromFirstVisit / 30.44
    ) %>%
    filter(DaysFromFirstVisit_months <= 24)
  
  if (!is.null(cutoff)) {
    x <- x %>%
      filter(ALSFRS_init > cutoff)
  }
  
  eligible_ids <- x %>%
    count(SSBS_ID) %>%
    filter(n >= min_visits) %>%
    pull(SSBS_ID)
  
  x_eligible <- x %>%
    filter(SSBS_ID %in% eligible_ids)
  
  tibble(
    N_slope = n_distinct(x_eligible$SSBS_ID),
    
    min_baseline_ALSFRSR =
      min(x_eligible$ALSFRS_init, na.rm = TRUE),
    
    max_baseline_ALSFRSR =
      max(x_eligible$ALSFRS_init, na.rm = TRUE),
    
    min_observed_ALSFRSR =
      min(x_eligible$ALSFRSR_Total, na.rm = TRUE),
    
    max_observed_ALSFRSR =
      max(x_eligible$ALSFRSR_Total, na.rm = TRUE)
  )
}


alsfrs_sample_summary <- bind_rows(
  primary = summarise_alsfrs_sample(3, 19),
  visits2 = summarise_alsfrs_sample(2, 19),
  no_cutoff = summarise_alsfrs_sample(3, NULL),
  both_relaxed = summarise_alsfrs_sample(2, NULL),
  .id = "scenario"
)

alsfrs_sample_summary


# Effect of alternative slope definitions on analysed sample ===================
scenario_order <- c(
  "primary",
  "no_cutoff",
  "visits2",
  "both_relaxed"
)

scenario_labels <- c(
  primary =
    "≥3 visits, baseline ALSFRS-R ≥20 (primary)",
  visits2 =
    "≥2 visits, baseline ALSFRS-R ≥20",
  no_cutoff =
    "≥3 visits, no baseline ALSFRS-R cutoff",
  both_relaxed =
    "≥2 visits, no baseline ALSFRS-R cutoff"
)

table_rev_v_a <- scenario_meta %>%
  dplyr::select(
    scenario,
    N_slope,
    N_ML
  ) %>%
  dplyr::left_join(
    alsfrs_sample_summary %>%
      dplyr::select(
        scenario,
        min_baseline_ALSFRSR,
        min_observed_ALSFRSR,
        max_observed_ALSFRSR
      ),
    by = "scenario"
  ) %>%
  dplyr::mutate(
    scenario = factor(
      scenario,
      levels = scenario_order
    ),
    `Slope definition` =
      scenario_labels[as.character(scenario)],
    `ALSFRS-R range used for slope estimation` =
      paste0(
        min_observed_ALSFRSR,
        "–",
        max_observed_ALSFRSR
      )
  ) %>%
  dplyr::arrange(scenario) %>%
  dplyr::transmute(
    `Slope definition`,
    `N in prediction analysis` = N_ML,
    `Minimum baseline ALSFRS-R` = min_baseline_ALSFRSR,
    `ALSFRS-R range used for slope estimation`
  )

table_rev_v_a


# Sensitivity of cross-validated prediction performance ========================

model_order <- c(
  "Core1",
  "Core2",
  "Core3",
  "Core4",
  "Core5",
  "Core5_NfL",
  "Core5_NfL_GFAP_csf_pgml_xGFAP_csf_pgml",
  "Core5_NfL_pTau181_csf_xpTau181_csf",
  "Core5_NfL_pTau217_csf_xpTau217_csf",
  "Core5_NfL_Ab38_40_csf_xAb38_40_csf",
  "Core5_NfL_Ab42_40_csf_bridged_xAb42_40_csf_bridged",
  "glm",
  "knn",
  "rf",
  "svmRadial",
  "xgb"
)

model_labels <- c(
  Core1 =
    "Age",
  Core2 =
    "Age + %VC",
  Core3 =
    "Age + %VC + ΔFS",
  Core4 =
    "Age + %VC + ΔFS + onset site",
  Core5 =
    "Core 5: Age + %VC + ΔFS + onset site + REEC",
  Core5_NfL =
    "Core 5 + NfL",
  
  Core5_NfL_GFAP_csf_pgml_xGFAP_csf_pgml =
    "Core 5 + NfL × GFAP",
  
  Core5_NfL_pTau181_csf_xpTau181_csf =
    "Core 5 + NfL × pTau181",
  
  Core5_NfL_pTau217_csf_xpTau217_csf =
    "Core 5 + NfL × pTau217",
  
  Core5_NfL_Ab38_40_csf_xAb38_40_csf =
    "Core 5 + NfL × Aβ38/40",
  
  Core5_NfL_Ab42_40_csf_bridged_xAb42_40_csf_bridged =
    "Core 5 + NfL × Aβ42/40",
  
  glm =
    "GLM (all predictors)",
  knn =
    "KNN (all predictors)",
  rf =
    "Random forest (all predictors)",
  svmRadial =
    "SVM (all predictors)",
  xgb =
    "XGB (all predictors)"
)

scenario_labels_short <- c(
  primary =
    "Primary: ≥3 visits, ALSFRS-R ≥20",
  no_cutoff =
    "≥3 visits, no cutoff",
  visits2 =
    "≥2 visits, ALSFRS-R ≥20",
  both_relaxed =
    "≥2 visits, no cutoff"
)


table_rev_v_b <- comparison_cv_r2 %>%
  dplyr::filter(
    model %in% model_order
  ) %>%
  dplyr::mutate(
    scenario = factor(
      scenario,
      levels = scenario_order
    ),
    model = factor(
      model,
      levels = model_order
    ),
    
    `Model` =
      model_labels[as.character(model)],
    
    `Slope definition` =
      scenario_labels_short[as.character(scenario)],
    
    `Cross-validated R²` =
      sprintf(
        "%.3f ± %.3f",
        R2_mean,
        R2_sd
      )
  ) %>%
  dplyr::arrange(
    model,
    scenario
  ) %>%
  dplyr::select(
    Model,
    `Slope definition`,
    `Cross-validated R²`
  ) %>%
  tidyr::pivot_wider(
    names_from = `Slope definition`,
    values_from = `Cross-validated R²`
  )

table_rev_v_b


# flextable ---------------------------------------------------------------



ft_rev_v_a <- flextable(table_rev_v_a) %>%
  theme_booktabs() %>%
  autofit() %>%
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "all") %>%
  bold(part = "header") %>%
  valign(valign = "center", part = "all") %>%
  set_caption(
    "Effect of alternative slope definitions on the analysed sample"
  ) %>%
  fontsize(size = 9, part = "all") %>%
  padding(padding = 3, part = "all")

ft_rev_v_a



ft_rev_v_b <- flextable(table_rev_v_b) %>%
  theme_booktabs() %>%
  autofit() %>%
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "all") %>%
  bold(part = "header") %>%
  valign(valign = "center", part = "all") %>%
  set_caption(
    "Sensitivity of cross-validated ALSFRS-R slope prediction to alternative slope definitions"
  ) %>%
  add_footer_lines(
    values = paste0(
      "Values are mean cross-validated R² ± SD across 10 repeated nested five-fold ",
      "cross-validation runs. Models correspond to those presented in Fig. 4D. ",
      "The ≥3-visit analysis without the baseline ALSFRS-R criterion included the ",
      "same participants as the primary analysis and therefore yielded identical ",
      "prediction results."
    )
  ) |> 
  fontsize(size = 9, part = "all") %>%
  padding(padding = 2, part = "all")

ft_rev_v_b

ft_rev_v_b -> SuppleT10_pub_added
SuppleT10_pub_added

# 
# flextable::save_as_docx(
#   ft_rev_v_a,
#   path = "output/Tables/ReviewerTable_SlopeDef_Sample.docx",
#   pr_section = officer::prop_section(
#     page_size = officer::page_size(orient = "portrait", width = 8.5, height = 11),
#     page_margins = officer::page_mar(bottom = 0.5, top = 0.5, left = 0.5, right = 0.5)
#   )
# )
# 
# 
# flextable::save_as_docx(
#   ft_rev_v_b,
#   path = "output/Tables/ReviewerTable_SlopeDef_CV_R2.docx",
#   pr_section = officer::prop_section(
#     page_size = officer::page_size(orient = "landscape", width = 11, height = 8.5),
#     page_margins = officer::page_mar(bottom = 0.5, top = 0.5, left = 0.5, right = 0.5)
#   )
# )
