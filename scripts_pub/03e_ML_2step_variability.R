# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# Parameters
variant_focus <- "sp95"

ml_colors2 <- c(
  "rf"        = "#C96A1B",
  "xgb"       = "#E28E2C",
  "svmRadial" = "#3F6F63",
  "kknn"      = "#7FA69A",
  "glm"       = "#4C6FA1",
  "NfL cutoff" = "gray40"
)

display_map <- c(
  rf = "RF", 
  xgb = "XGB", 
  svmRadial = "SVM", 
  kknn = "KNN", 
  glm = "GLM"
)

colors_step1 <- c(
  "Se95 (low cut-off)" = "gray90",   # light blue (rule-out)
  "Sp95 (high cut-off)" = "gray90",  # orange (rule-in)
  "Rule-out" = "#56B4E9",              # light blue
  "Gray-zone" = "gray60",              # gray
  "Rule-in" = "#E69F00",               # orange
  "PPV" = "#CC79A7",                   # purple-pink
  "NPV" = "#009E73"                    # green
)

# 1) NfL threshold variability (separate plots with individual control) ---------
nfl_thresholds <- rep_out$stacked$nfl_baselines %>%
  select(repeat_id, baseline, th_train_all) %>%
  unnest_longer(th_train_all, values_to = "threshold") %>%
  group_by(repeat_id, baseline) %>%
  summarise(threshold = mean(threshold, na.rm = TRUE), .groups = "drop")

# Se95 plot
nfl_se95 <- nfl_thresholds %>%
  filter(baseline == "nfl_se95")

p_nfl_se95 <- ggplot(nfl_se95, aes(x = threshold, y = "Se95\n(low cut-off)")) +
  geom_boxplot(outlier.shape = NA, width = 0.2, fill = colors_step1["Se95 (low cut-off)"]) +
  geom_jitter(height = 0.1, width = 0, alpha = 0.4, size = 1.0) +
  scale_x_continuous(limits = c(700, 900), breaks = seq(600, 900, 50)) +
  labs(
    title = "Step1: NfL threshold variability",
    x = NULL, y = NULL) +
  theme_classic(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 12))

# Sp95 plot
nfl_sp95 <- nfl_thresholds %>%
  filter(baseline == "nfl_sp95")

p_nfl_sp95 <- ggplot(nfl_sp95, aes(x = threshold, y = "Sp95\n(high cut-off)")) +
  geom_boxplot(outlier.shape = NA, width = 0.2, fill = colors_step1["Sp95 (high cut-off)"]) +
  geom_jitter(height = 0.1, width = 0, alpha = 0.4, size = 1.0) +
  scale_x_continuous(limits = c(5500, 7500), breaks = seq(4000, 8000, 500)) +
  labs(x = "NfL threshold (pg/mL)", y = NULL) +
  theme_classic(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 12))

# Combine threshold plots vertically
p_nfl_thresholds <- p_nfl_se95 / p_nfl_sp95 
p_nfl_thresholds

# 2) Bucket coverage - patient assignment variability ---------------------
bucket_coverage <- rep_out$stacked$buckets %>%
  mutate(prop_pct = prop * 100) %>%
  mutate(bucket = factor(bucket, 
                         levels = c("low", "gray", "high"),
                         labels = c("Rule-out", "Gray-zone", "Rule-in")))

p_bucket_coverage <- ggplot(bucket_coverage, aes(x = prop_pct, y = bucket, fill = bucket)) +
  geom_boxplot(outlier.shape = NA, width = 0.3) +
  geom_jitter(height = 0.15, width = 0, alpha = 0.4, size = 1.0) +
  scale_fill_manual(values = colors_step1, guide = "none") +
  labs(title = "Patient assignment variability",
       x = "Proportion (%)", y = NULL) +
  scale_x_continuous(limits = c(15, 70), breaks = seq(0, 100, 20)) +
  theme_classic(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 12))

p_bucket_coverage

# 3) Step1-only performance (gray-zone excluded) --------------------------

step1_performance <- lapply(seq_along(rep_out$runs), function(i) {
  any_method <- names(rep_out$runs[[i]]$per_method)[1]
  oof <- rep_out$runs[[i]]$per_method[[any_method]]$pooled$oof
  
  y_true <- oof$y_true
  bucket <- oof$bucket
  
  # Rule-in (high): PPV = ALS / (ALS + DC)
  idx_high <- bucket == "high"
  n_als_high <- sum(y_true[idx_high] == "ALS")
  n_dc_high <- sum(y_true[idx_high] == "nonALS")
  ppv_high <- if ((n_als_high + n_dc_high) > 0) {
    n_als_high / (n_als_high + n_dc_high)
  } else NA_real_
  
  # Rule-out (low): NPV = DC / (DC + ALS)
  idx_low <- bucket == "low"
  n_als_low <- sum(y_true[idx_low] == "ALS")
  n_dc_low <- sum(y_true[idx_low] == "nonALS")
  npv_low <- if ((n_dc_low + n_als_low) > 0) {
    n_dc_low / (n_dc_low + n_als_low)
  } else NA_real_
  
  tibble(
    repeat_id = i,
    PPV = ppv_high,
    NPV = npv_low
  )
})

step1_perf <- bind_rows(step1_performance) %>%
  pivot_longer(cols = c(PPV, NPV), 
               names_to = "metric", 
               values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("PPV", "NPV")))

p_step1_performance <- ggplot(step1_perf, aes(x = value, y = metric, fill = metric)) +
  geom_boxplot(outlier.shape = NA, width = 0.2) +
  geom_jitter(height = 0.1, width = 0, alpha = 0.4, size = 1.0) +
  scale_fill_manual(values = colors_step1, guide = "none") +
  labs(title = "Predictive value variability",
       x = "Predictive value outside the gray zone", y = NULL) +
  scale_x_continuous(limits = c(0.6, 1), breaks = seq(0, 1, 0.1)) +
  theme_classic(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 12))

p_step1_performance

# Combine all new plots
p_step1_summary <- (p_nfl_thresholds | p_bucket_coverage | p_step1_performance) +
  plot_annotation(
    title = "Step 1 (NfL-only) summary across 10 repeats",
    theme = theme(plot.title = element_text(size = 18, face = "bold"))
  )


p_step1_summary


# 4) Method comparison in Step 2 (overall vs. gray-zone only) -------------

# Extract data from rep_out$stacked
# Overall metrics: AUC, PPV, NPV per method × repeat
overall_data <- rep_out$stacked$overall %>%
  filter(variant == variant_focus) %>%
  select(repeat_id, method, AUC, PPV, NPV) %>%
  pivot_longer(cols = c(AUC, PPV, NPV), 
               names_to = "metric", 
               values_to = "value") %>%
  mutate(
    method_chr = as.character(method),
    method_display = dplyr::recode(method_chr, !!!display_map))


# Gray-zone metrics: AUC_gray from overall, PPV_gray/NPV_gray from gray_only
gray_auc <- rep_out$stacked$overall %>%
  filter(variant == variant_focus) %>%
  select(repeat_id, method, AUC_gray)

gray_ppv_npv <- rep_out$stacked$gray_only %>%
  filter(variant == variant_focus) %>%
  select(repeat_id, method, PPV_gray, NPV_gray)

gray_data <- gray_auc %>%
  left_join(gray_ppv_npv, by = c("repeat_id", "method")) %>%
  pivot_longer(cols = c(AUC_gray, PPV_gray, NPV_gray),
               names_to = "metric",
               values_to = "value") %>%
  mutate(metric = case_when(
    metric == "AUC_gray" ~ "AUC",
    metric == "PPV_gray" ~ "PPV",
    metric == "NPV_gray" ~ "NPV"
  )) %>%
  mutate(
    method_chr = as.character(method),
    method_display = dplyr::recode(method_chr, !!!display_map),
    method_display = factor(method_display, levels = display_map[method_order])
  )


# Method order by median overall AUC
method_order <- overall_data %>%
  filter(metric == "AUC") %>%
  group_by(method) %>%
  summarise(med_auc = median(value, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(med_auc)) %>%
  pull(method)

method_order

# Apply factor levels
overall_data <- overall_data %>%
  mutate(
    method = factor(method, levels = method_order),
    metric = factor(metric, levels = c("AUC", "PPV", "NPV")),
    method_display = factor(method_display, levels = display_map[method_order])
  )

gray_data <- gray_data %>%
  mutate(
    method = factor(method, levels = method_order),
    metric = factor(metric, levels = c("AUC", "PPV", "NPV"))
  ) 

# Separate plots for each metric to control alignment --------------


## AUC panels -------------------------

p_gray_auc <- gray_data %>%
  filter(metric == "AUC") %>%
  ggplot(aes(x = value, y = method_display, fill = method)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(height = 0.15, width = 0, alpha = 0.5, size = 1.0) +
  scale_fill_manual(values = ml_colors2, guide = "none") +
  scale_x_continuous(limits = c(0.6, 1), breaks = seq(0.5, 1, 0.1)) +
  labs(title = "Step 2: Gray-zone diagnostic performance", x = "AUC (gray)", y = NULL) +
  theme_classic(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 12))

p_overall_auc <- overall_data %>%
  filter(metric == "AUC") %>%
  ggplot(aes(x = value, y = method_display, fill = method)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(height = 0.15, width = 0, alpha = 0.5, size = 1.0) +
  scale_fill_manual(values = ml_colors2, guide = "none") +
  scale_x_continuous(limits = c(0.6, 1), breaks = seq(0.5, 1, 0.1)) +
  labs(title = "Overall diagnostic performance (Step 1 & 2)", x = "overall AUC", y = NULL) +
  theme_classic(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 12))

## PPV panels -------------------------
p_gray_ppv <- gray_data %>%
  filter(metric == "PPV") %>%
  ggplot(aes(x = value, y = method_display, fill = method)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(height = 0.15, width = 0, alpha = 0.5, size = 1.0) +
  scale_fill_manual(values = ml_colors2, guide = "none") +
  scale_x_continuous(limits = c(0.6, 1), breaks = seq(0.5, 1, 0.1)) +
  labs(x = "PPV (gray)", y = NULL) +
  theme_classic(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

p_overall_ppv <- overall_data %>%
  filter(metric == "PPV") %>%
  ggplot(aes(x = value, y = method_display, fill = method)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(height = 0.15, width = 0, alpha = 0.5, size = 1.0) +
  scale_fill_manual(values = ml_colors2, guide = "none") +
  scale_x_continuous(limits = c(0.6, 1), breaks = seq(0.5, 1, 0.1)) +
  labs(x = "overall PPV", y = NULL) +
  theme_classic(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

## NPV panels -------------------------
p_gray_npv <- gray_data %>%
  filter(metric == "NPV") %>%
  ggplot(aes(x = value, y = method_display, fill = method)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(height = 0.15, width = 0, alpha = 0.5, size = 1.0) +
  scale_fill_manual(values = ml_colors2, guide = "none") +
  scale_x_continuous(limits = c(0.3, 0.7), breaks = seq(0.3, 1, 0.1)) +
  labs(x = "NPV (gray)", y = NULL) +
  theme_classic(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

p_overall_npv <- overall_data %>%
  filter(metric == "NPV") %>%
  ggplot(aes(x = value, y = method_display, fill = method)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(height = 0.15, width = 0, alpha = 0.5, size = 1.0) +
  scale_fill_manual(values = ml_colors2, guide = "none") +
  scale_x_continuous(limits = c(0.3, 0.7), breaks = seq(0.3, 1, 0.1)) +
  labs(x = "overall NPV", y = NULL) +
  theme_classic(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# Combine with patchwork ----------------------------

p_combined <- (p_gray_auc | p_gray_ppv | p_gray_npv) /
  (p_overall_auc | p_overall_ppv | p_overall_npv) +
  plot_annotation(
    title = sprintf("Method comparison across 10 repeats (variant = %s)", variant_focus),
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

p_combined


# 🔍 Supple Fig S6 ------------------------------------------------------------

SuppleFig_Diag_2step_summary <- (p_step1_summary/p_combined) + 
  plot_layout(heights = c(1, 2)) + 
  plot_annotation(
    title = "Diagnostic performance of the two-step approach across repeated nested cross-validation",
    theme = theme(plot.title = element_text(hjust=0, size = 20, face = "bold", margin=ggplot2::margin(l=60, b=10, t=5)))
  )

SuppleFig6 <- SuppleFig_Diag_2step_summary

SuppleFig6                                  
