suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(caret)
  library(pROC)
})

#LOAD
source("scripts_pub/function/ML_class_func.R")

# RUN ==============================================================

pos_class <- "ALS"
methods_to_try <- c("rf","xgb","glm","svmRadial","kknn")

# Capture arguments before running
.twostep_args <- list(
  df               = ROC_df,
  outcome          = "ALS_label",
  predictors_all   = c("Cr","CK","NfL_csf_pgml","GFAP_csf_pgml",
                       "pTau181_csf","pTau217_csf","Ab42_40_csf_bridged", "Ab38_40_csf"),
  predictors_gray  = c("Cr","CK","GFAP_csf_pgml",
                       "pTau181_csf","pTau217_csf","Ab42_40_csf_bridged", "Ab38_40_csf"),
  methods_gray     = methods_to_try,
  k                = 5,
  n_repeats        = 10,
  seed_base        = 20250912,
  seed_inner       = 12345,
  target_spec_high = 0.95,
  target_sens_low  = 0.95,
  gray_spec_targets= c(0.95),
  gray_sens_targets = c(0.95),
  nfl_baseline_sens_target = 0.95
)


rep_out <- do.call(twostep_kfold_multi_repeat, .twostep_args)



## 💾 Save trained models and results -------------------------------------

# .timestamp_save <- format(Sys.time(), "%Y%m%d_%H%M")
# .output_file_save <- paste0("output/TwoStep_Models_Trained_", .timestamp_save, ".rds")
# if (!dir.exists("output")) dir.create("output")
# 
# saveRDS(
#   list(
#     metadata = list(
#       timestamp = Sys.time(),
#       seed_base = 20250912,
#       seed_inner = 12345,
#       k_folds = 5,
#       n_repeats = 5,
#       n_samples = nrow(ROC_df),
#       n_ALS = sum(ROC_df$ALS_label == "ALS"),
#       n_nonALS = sum(ROC_df$ALS_label == "nonALS"),
#       outcome = "ALS_label",
#       predictors_all = c("Cr","CK","NfL_csf_pgml","GFAP_csf_pgml",
#                          "pTau181_csf","pTau217_csf","Ab42_40_csf_bridged", "Ab38_40_csf"),
#       predictors_gray = c("Cr","CK","GFAP_csf_pgml",
#                           "pTau181_csf","pTau217_csf","Ab42_40_csf_bridged", "Ab38_40_csf"),
#       methods_gray = methods_to_try,
#       target_spec_high = 0.95,
#       target_sens_low = 0.95,
#       gray_spec_targets = c(0.95),
#       gray_sens_targets = c(0.95),
#       nfl_baseline_sens_target = 0.95
#     ),
#     rep_out = rep_out,
#     training_data = ROC_df
#   ),
#   .output_file_save,
#   compress = "xz"
# )
# 
# cat("\n=== Two-Step Model Training Complete ===\n")
# cat("Saved to:", .output_file_save, "\n")
# cat("Key Points:\n")
# cat("  1. Two-step classification across", 5, "repeats with", 5, "-fold CV\n")
# cat("  2. NfL baseline thresholds computed per repeat\n")
# cat("  3. Gray-zone models trained on", length(methods_to_try), "methods\n")
# cat("  4. Stacked results and per-repeat details included\n\n")
# 
# rm(.timestamp_save, .output_file_save)
# 
# 
# rep_out$used_seeds
# rep_out$attempt_log
# rep_out$stacked
# 
# rep_out$runs[[2]]$compare
# #rep_out$runs[[2]]$per_method
# rep_out$runs[[2]]$meta
# 
# rep_out$runs[[2]]$per_method$rf$pooled


# ➡️ Load the trained models ---------------------------------------------------------

# Load the most recent trained model file
# latest_file <- list.files("output", 
#                           pattern = "^TwoStep_Models_Trained_.*\\.rds$", 
#                           full.names = TRUE) %>%
#   sort(decreasing = TRUE) %>%
#   head(1)
# 
# loaded <- readRDS(latest_file)
# rep_out <- loaded$rep_out
# ROC_df <- loaded$training_data
# 
# cat("Loaded:", basename(latest_file), "\n")




# 0) Choose a working variant ===============================
#    (can switch to "youden" or "sp90")

variant_pick <- "sp95"

# Helper to label method groups for plotting (optional)
method_to_group <- function(method) {
  # Simple rule: "glm" or "one:..." => GLM; others => ML
  if (grepl("^glm($|:)|^one:", method)) return("GLM")
  "ML"
}

# 1) OVERALL metrics across repeats ===============================
#    (AUC, PPV, NPV) for a chosen variant

overall_long <- rep_out$stacked$overall %>%
  # keep only the chosen variant
  filter(variant == variant_pick) %>%
  # add a group (ML/GLM)
  mutate(group = vapply(method, method_to_group, character(1))) %>%
  select(repeat_id, group, method, variant, AUC, PPV, NPV) %>%
  pivot_longer(cols = c(AUC, PPV, NPV), names_to = "metric", values_to = "value")

# Order methods by median OVERALL AUC across repeats (nice for plotting)
method_order_overall <- overall_long %>%
  filter(metric == "AUC") %>%
  group_by(method) %>%
  summarize(med = median(value, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(med)) %>%
  pull(method)

overall_long <- overall_long %>%
  mutate(method = factor(method, levels = method_order_overall),
         metric = factor(metric, levels = c("AUC","PPV","NPV")))

# Plot: box + jitter per method, facet by metric
p_overall <- ggplot(overall_long, aes(x = value, y = method, color = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  geom_jitter(height = 0.15, width = 0, alpha = 0.6, size = 1.8) +
  facet_wrap(~ metric, ncol = 3, scales = "free_x") +
  labs(title = sprintf("Overall metrics across repeats (variant = %s)", variant_pick),
       x = NULL, y = NULL, color = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())
p_overall


# 2) GRAY metrics across repeats ===============================
#    (AUC_gray, PPV_gray, NPV_gray) for the same variant
#    AUC_gray comes from stacked overall (column AUC_gray)
#    PPV_gray/NPV_gray come from stacked gray_only


# AUC_gray per repeat/method/variant
gray_auc <- rep_out$stacked$overall %>%
  filter(variant == variant_pick) %>%
  select(repeat_id, method, variant, AUC_gray)

# PPV_gray/NPV_gray per repeat/method/variant
gray_ppvnpv <- rep_out$stacked$gray_only %>%
  filter(variant == variant_pick) %>%
  select(repeat_id, method, variant, PPV_gray, NPV_gray)

# join
gray_join <- gray_auc %>%
  left_join(gray_ppvnpv, by = c("repeat_id","method","variant")) %>%
  mutate(group = vapply(method, method_to_group, character(1)))

# long for plotting
gray_long <- gray_join %>%
  pivot_longer(cols = c(AUC_gray, PPV_gray, NPV_gray),
               names_to = "metric", values_to = "value") %>%
  mutate(method = factor(method, levels = method_order_overall),
         metric = factor(metric, levels = c("AUC_gray","PPV_gray","NPV_gray"),
                         labels = c("AUC (gray)", "PPV (gray)", "NPV (gray)")))

p_gray <- ggplot(gray_long, aes(x = value, y = method, color = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  geom_jitter(height = 0.15, width = 0, alpha = 0.6, size = 1.8) +
  facet_wrap(~ metric, ncol = 3, scales = "free_x") +
  labs(title = sprintf("Gray-zone metrics across repeats (variant = %s)", variant_pick),
       x = NULL, y = NULL, color = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())
p_gray



# 3) Gray coverage across repeats ===============================
#    (proportion of gray in Step 1)
#    Buckets are shared across methods within a repeat; 
#    we can just take one method's buckets per repeat (e.g. the first one) to get the gray coverage per repeat.

gray_cov <- rep_out$stacked$buckets %>%
  filter(bucket == "gray") %>%
  transmute(repeat_id, gray_prop = prop)

p_cov <- ggplot(gray_cov, aes(x = gray_prop, y = 1)) +
  geom_boxplot(width = 0.3, outlier.shape = NA) +
  geom_jitter(height = 0.05, width = 0, alpha = 0.7, size = 1.8) +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(title = "Gray-zone coverage across repeats",
       x = "Proportion in gray") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())
p_cov


# 4) NfL baseline thresholds distribution (TRAIN-side per fold) ===============================
#    rep_out$stacked$nfl_baselines holds one row per baseline per repeat,
#    with a list-column 'th_train_all' of K thresholds (one per fold).
#    We unnest it and draw distributions for:
#      - nfl_youden
#      - nfl_sp95
#      - nfl_se95

library(ggdist) 

nfl_th <- rep_out$stacked$nfl_baselines %>%
  select(repeat_id, baseline, th_train_all) %>%
  unnest_longer(th_train_all, values_to = "threshold") %>%
  mutate(baseline = factor(baseline,
                           levels = c("nfl_sp95","nfl_youden",
                                      unique(baseline[grepl("^nfl_se", baseline)]) )))

## 📊 for reporting: NfL threshold --------------------------------------------

nfl_th %>% group_by(baseline) %>%
  summarise(
    mean = mean(threshold, na.rm=TRUE),
    sd   = sd(threshold, na.rm=TRUE),
    med = median(threshold, na.rm=TRUE),
    q25 = quantile(threshold, 0.25, na.rm=TRUE),
    q75 = quantile(threshold, 0.75, na.rm=TRUE)
  )

# Basic box + jitter
ggplot(nfl_th, aes(x = threshold, y = baseline)) +
  geom_boxplot(outlier.shape = NA, width = 0.2) +
  geom_jitter(height = 0.05, width = 0, alpha = 0.4, size = 1.2) +
  labs(title = "Distribution of NfL TRAIN thresholds across repeats (and folds)",
       x = "Threshold (NfL units on TRAIN ROC per fold)",
       y = NULL) +
  theme_bw(base_size = 16) +
  theme(panel.grid.minor = element_blank())



## 📊 for reporting: property of gray-zone cases (%) ------------------------------

gray_cov %>%
  summarise(
    mean = mean(gray_prop, na.rm=TRUE)*100,
    sd   = sd(gray_prop, na.rm=TRUE)*100,
    med = median(gray_prop, na.rm=TRUE)*100,
    q25 = quantile(gray_prop, 0.25, na.rm=TRUE)*100,
    q75 = quantile(gray_prop, 0.75, na.rm=TRUE)*100
  )



# 5) A compact summary table ===============================
#    Medians and IQR across repeats for the chosen variant

summarize_metric <- function(x) {
  c(med = median(x, na.rm = TRUE),
    q25 = quantile(x, 0.25, na.rm = TRUE),
    q75 = quantile(x, 0.75, na.rm = TRUE))
}

# overall: median & IQR per method × metric
tbl_summary_overall <- overall_long %>%
  dplyr::group_by(group, method, metric) %>%
  dplyr::summarise(
    med = stats::median(value, na.rm = TRUE),
    q25 = stats::quantile(value, 0.25, na.rm = TRUE),
    q75 = stats::quantile(value, 0.75, na.rm = TRUE),
    n   = sum(is.finite(value)),
    .groups = "drop"
  ) %>%
  dplyr::arrange(metric, dplyr::desc(med), method)

# gray: median & IQR per method × metric (AUC_gray / PPV_gray / NPV_gray)
tbl_summary_gray <- gray_long %>%
  dplyr::group_by(group, method, metric) %>%
  dplyr::summarise(
    med = stats::median(value, na.rm = TRUE),
    q25 = stats::quantile(value, 0.25, na.rm = TRUE),
    q75 = stats::quantile(value, 0.75, na.rm = TRUE),
    mean = base::mean(value, na.rm = TRUE),
    sd   = stats::sd(value, na.rm = TRUE),
    n   = sum(is.finite(value)),
    .groups = "drop"
  ) %>%
  dplyr::arrange(metric, dplyr::desc(med), method)

tbl_summary_overall
tbl_summary_gray


tbl_summary_overall <- overall_long %>%
  dplyr::group_by(group, method, metric) %>%
  dplyr::summarise(
    med = stats::median(value, na.rm = TRUE),
    q25 = stats::quantile(value, 0.25, na.rm = TRUE),
    q75 = stats::quantile(value, 0.75, na.rm = TRUE),
    mean = base::mean(value, na.rm = TRUE),
    sd   = stats::sd(value, na.rm = TRUE),
    n    = sum(is.finite(value)),
    .groups = "drop"
  ) %>%
  dplyr::arrange(metric, dplyr::desc(med), method)


# overall wide: one row per method, columns per metric with med (and IQR)
overall_wide <- tbl_summary_overall %>%
  dplyr::mutate(iqr = paste0(sprintf("%.3f", q25), "–", sprintf("%.3f", q75))) %>%
  dplyr::mutate(mediqr = paste0(sprintf("%.3f", med), " [", iqr, "]")) %>%
  dplyr::select(group, method, metric, mediqr) %>%
  tidyr::pivot_wider(names_from = metric, values_from = mediqr) %>%
  dplyr::arrange(group, method)

flextable::flextable(overall_wide)

# gray wide
gray_wide <- tbl_summary_gray %>%
  dplyr::mutate(iqr = paste0(sprintf("%.3f", q25), "–", sprintf("%.3f", q75))) %>%
  dplyr::mutate(mediqr = paste0(sprintf("%.3f", med), " [", iqr, "]")) %>%
  dplyr::select(group, method, metric, mediqr) %>%
  tidyr::pivot_wider(names_from = metric, values_from = mediqr) %>%
  dplyr::arrange(group, method)

flextable::flextable(gray_wide)


# overall wide (mean SD)

overall_wide_meanSD <- tbl_summary_overall %>%
  dplyr::mutate(mean_sd = paste0(sprintf("%.3f", mean), " (", sprintf("%.3f", sd), ")")) %>%
  dplyr::select(group, method, metric, mean_sd) %>%
  tidyr::pivot_wider(names_from = metric, values_from = mean_sd) %>%
  dplyr::arrange(group, method)

flextable::flextable(overall_wide_meanSD)

# gray wide (mean SD)

gray_wide_meanSD <- tbl_summary_gray %>%
  dplyr::mutate(mean_sd = paste0(sprintf("%.3f", mean), " (", sprintf("%.3f", sd), ")")) %>%
  dplyr::select(group, method, metric, mean_sd) %>%
  tidyr::pivot_wider(names_from = metric, values_from = mean_sd) %>%
  dplyr::arrange(group, method)

