library(flextable)
library(caret)
library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)
library(ggforce)
library(future)

# OLS-based ML Regression Pipeline for ALSFRS-R Slope Prediction
# requires - df_slope_model_scaled: preprocessed dataframe with outcome and predictors

source("scripts/ml/ML_reg_func.R")

# Outcome & Predictors Definition ================================

outcome_reg <- "slope_total"  # OLS slope (negative = decline)

predictors_reg <- c(
  #CLINICAL FEATURES
  "Age_init", "Sex", "OnsetSite", "REEC_definite", "deltaFS", "VC_Percent",
  #CSF BIOMARKERS
  "NfL_csf_pgml", "GFAP_csf_pgml", "pTau181_csf", "pTau217_csf",
  "Ab38_40_csf","Ab42_40_csf_bridged", "Ab_status",
  #BLOOD BIOMARKERS
  "Cr", "CK" 
)

# Build analysis dataframe (complete cases only)
df <- df_slope_model_scaled %>%
  dplyr::select(all_of(c(outcome_reg, predictors_reg))) %>%
  tidyr::drop_na()

cat(sprintf("Final sample size for ML: N = %d\n", nrow(df)))
cat(sprintf("Outcome range: %.3f to %.3f\n", 
            min(df$slope_total), max(df$slope_total)))

outcome <- outcome_reg
predictors <- predictors_reg %>% setdiff(c("Ab_status"))  # Exclude Ab_status for ML (since Ab42/40 is included)


# Setup: Random seed and outer folds ---------------------------

RNGkind("L'Ecuyer-CMRG")
.seed_value <- 12345
set.seed(.seed_value)
options(contrasts = c("contr.treatment", "contr.poly"))

K_outer <- 5
inner_k <- 5
n_repeats <- 10

y <- df$slope_total

# Create repeated outer folds (shared across all models)
outer_folds_repeats <- lapply(seq_len(n_repeats), function(r) {
  caret::createFolds(y, k = K_outer, list = TRUE, returnTrain = FALSE)
})

# OUTER FOLDS REPEATS structure:
outer_folds_repeats

# ══════════════════════════════════════════════════════════════
# ML Methods
# ══════════════════════════════════════════════════════════════

methods <- c("glm", "rf", "svmRadial", "knn", "xgb")

# ══════════════════════════════════════════════════════════════
# !!Run!! Repeated Nested CV for ML Models -------------------------
# ══════════════════════════════════════════════════════════════

cat("\n=== Running repeated nested CV for ML models ===\n")
cat(sprintf("Configuration: %d repeats × %d outer folds × %d inner folds\n",
            n_repeats, K_outer, inner_k))
cat("This may take several minutes...\n\n")

nested_rep_tbl <- purrr::map_dfr(methods, function(m) {
  cat(sprintf("Training %s...\n", toupper(m)))
  resr <- nested_fit_regression_repeated(
    data                = df,
    outcome             = outcome,
    predictors          = predictors,
    method              = m,
    K_outer             = K_outer,
    inner_k             = inner_k,
    repeats             = n_repeats,
    seed                = .seed_value,
    scale               = "auto",
    outer_folds_repeats = outer_folds_repeats
  )
  tibble::tibble(
    model          = m,
    stats_on_mean  = list(resr$stats_on_mean),
    rep_stats      = list(resr$rep_stats),
    mean_sd        = list(resr$mean_sd),
    oof_pred_mean  = list(resr$oof_pred_mean),
    per_repeat     = list(resr$per_repeat)
  )
})


cat("\n=== ML model training complete ===\n\n")

# Results: Mean ± SD across repeats ----------------------

ml_flat_mean_sd <- nested_rep_tbl %>%
  dplyr::select(model, mean_sd) %>%
  tidyr::unnest(mean_sd) |> arrange(RMSE_mean,R2_mean)

ml_flat_mean_sd 

# sensitivity check: single eval on "mean OOF" predictions
ml_flat_on_mean <- nested_rep_tbl %>%
  dplyr::select(model, stats_on_mean) %>%
  tidyr::unnest(stats_on_mean) |> arrange(RMSE,R2)

ml_flat_on_mean


# ══════════════════════════════════════════════════════════════
# GLM Subset Analysis -----------------------------------------
# ══════════════════════════════════════════════════════════════

cat("\n=== GLM Subset Models ===\n")

# Core covariate sets (staged)
core_sets <- list(
  Core0 = c(),
  Core1 = c("Age_init"),
  Core2 = c("Age_init", "VC_Percent"),
  Core3 = c("Age_init", "VC_Percent", "deltaFS"),
  Core4 = c("Age_init", "VC_Percent", "deltaFS", "OnsetSite"),
  Core5 = c("Age_init", "VC_Percent", "deltaFS", "OnsetSite","REEC_definite")
)

# Biomarker candidates
bio_blood <- c("Cr", "CK")
bio_csf   <- c("NfL_csf_pgml", "Ab_status", "Ab38_40_csf", "Ab42_40_csf_bridged",
               "GFAP_csf_pgml", "pTau181_csf", "pTau217_csf")
bio_one_any <- c("Ab42_40_csf_bridged", "Ab38_40_csf", "GFAP_csf_pgml", 
                 "pTau181_csf", "pTau217_csf")
bio_one_csf <- c("GFAP_csf_pgml", "pTau181_csf", "pTau217_csf", "Ab38_40_csf","Ab42_40_csf_bridged")

# Expand specifications (Core + NfL + bio + interactions)
specs <- expand_specs4(core_sets, bio_blood, bio_csf, bio_one_any, bio_one_csf)


# ══════════════════════════════════════════════════════════════
# !!Run!! GLM Subset Models  -----------------------------------
# ══════════════════════════════════════════════════════════════

glm_subset_results <- run_glm_subsets(
  data                = df,
  outcome             = "slope_total",
  specs               = specs,
  repeats             = n_repeats,
  K_outer             = K_outer,
  inner_k             = inner_k,
  seed                = .seed_value,
  scale               = "auto",
  outer_folds_rep     = outer_folds_repeats
)


# ════════════════════════════════════════════════
# Save trained Results -------------------------- 
# ════════════════════════════════════════════════
# 
.timestamp_tmp <- format(Sys.time(), "%Y%m%d_%H%M")
.output_file_tmp <- paste0("output/Slope_Models_Trained_", .timestamp_tmp, ".rds")

if (!dir.exists("output")) dir.create("output")

saveRDS(
  list(
    metadata = list(
      timestamp = Sys.time(),
      seed = .seed_value,
      K_outer = K_outer,
      inner_k = inner_k,
      n_repeats = n_repeats,
      n_samples = nrow(df),
      outcome = outcome,
      predictors = predictors
    ),
    outer_folds_repeats = outer_folds_repeats,
    ml_results = nested_rep_tbl,
    glm_results = glm_subset_results,
    training_data = df
  ),
  .output_file_tmp,
  compress = "xz"
)

cat("\n=== Model Training Complete ===\n")
cat("Saved to:", .output_file_tmp, "\n")

rm(.timestamp_tmp, .output_file_tmp)

