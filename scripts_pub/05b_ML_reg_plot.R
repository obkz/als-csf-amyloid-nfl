library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)


# Load trained results ----------------------------------------------------

#trained_results <- readRDS("output/Slope_Models_Trained_20260619_1417.rds")
nested_rep_tbl <- trained_results$ml_results
glm_subset_results <- trained_results$glm_results
metadata <- trained_results$metadata

ml_flat_mean_sd <- nested_rep_tbl %>%
  dplyr::select(model, mean_sd) %>%
  tidyr::unnest(mean_sd) |> arrange(RMSE_mean,R2_mean)

ml_flat_mean_sd 

# Compute 95% CI for GLM models (n=10 repeats)
glm_flat_mean_ci <- glm_subset_results$nested %>%
  select(model, rep_stats) %>%
  unnest(rep_stats) %>%
  group_by(model) %>%
  summarise(
    n_repeats = n(),
    R2_mean = mean(R2),
    R2_sd = sd(R2),
    R2_se = sd(R2) / sqrt(n()),
    # t-distribution with df = n-1 = 9
    R2_ci_lower = R2_mean - qt(0.975, df = n() - 1) * R2_se,
    R2_ci_upper = R2_mean + qt(0.975, df = n() - 1) * R2_se,
    RMSE_mean = mean(RMSE),
    RMSE_sd = sd(RMSE),
    .groups = "drop"
  ) %>%
  mutate(
    base_core = sub("_.*$", "", model),
    is_core_only = !grepl("_", model),
    is_additive = grepl("_NfL_", model) & !grepl("_x", model),
    is_interact = grepl("_NfL_.*_x", model)
  )

# Compute 95% CI for ML models (n=10 repeats)
ml_flat_mean_ci <- ml_flat_mean_sd %>%
  mutate(
    n_repeats = metadata$n_repeats,  # Should be 10
    R2_se = R2_sd / sqrt(n_repeats),
    R2_ci_lower = R2_mean - qt(0.975, df = n_repeats - 1) * R2_se,
    R2_ci_upper = R2_mean + qt(0.975, df = n_repeats - 1) * R2_se,
    RMSE_se = RMSE_sd / sqrt(n_repeats),
    RMSE_ci_lower = RMSE_mean - qt(0.975, df = n_repeats - 1) * RMSE_se,
    RMSE_ci_upper = RMSE_mean + qt(0.975, df = n_repeats - 1) * RMSE_se
  )


# Configuration -------------------------------------------------------

TARGET_CORE <- 5
SHOW_ML_BLOCK <- TRUE

CORE_LABELS <- c(
  "Core1" = "Age",
  "Core2" = "Age + %VC",
  "Core3" = "Core3\n(Age + %VC + ΔFS)",
  "Core4" = "Core4\n(Core3 + Onset site)",
  "Core5" = "Core5\n(Core4 + REEC)"
)

LABEL_PAD_INTERACTION <- c(
  "Aβ status" = 0, 
  "GFAP"      = 0, 
  "Aβ38/40"   = 0.0075,   "Aβ42/40"   = 0.005, 
  "pTau217"   = -0.0275,   "pTau181"   = -0.0125
)

LABEL_PAD_ML <- c("GLM (all)" = 0.0175, "RF"= -0.1, "XGB"= -0.12, 
                  "SVM" = -0.0175,"KNN" = -0.00375)

# Helper functions
shorten_alpha <- function(x) {
  x <- gsub("^GFAP_csf_pgml$", "GFAP", x)
  x <- gsub("^pTau181_csf$",  "pTau181", x)
  x <- gsub("^pTau217_csf$",  "pTau217", x)
  x <- gsub("^Ab_status$",    "Aβ status", x)
  x <- gsub("^Ab38_40_csf$",  "Aβ38/40", x)
  x <- gsub("^Ab42_40_csf_bridged$", "Aβ42/40", x)
  x
}

offset_map_fixed <- function(x, width = 0.7) {
  u <- unique(as.character(x))
  if (length(u) <= 1) return(setNames(0, u))
  setNames(seq(-width/2, width/2, length.out = length(u)), u)
}

get_extra_pad <- function(alpha_short, y_lo, manual_vec, step = 0.012) {
  man <- as.numeric(manual_vec[as.character(alpha_short)])
  auto <- {
    ord <- order(y_lo, decreasing = FALSE)
    k <- seq_along(y_lo) - 1
    ex <- numeric(length(y_lo))
    ex[ord] <- (k %% 4) * step
    ex
  }
  ifelse(is.na(man), auto, man)
}

pretty_ml_label <- function(x) {
  dplyr::recode(x,
                "glm" = "GLM (all)",
                "knn" = "KNN",
                "svmRadial" = "SVM",
                "rf" = "RF",
                "xgb" = "XGB",
                .default = x)
}

# Prepare data for plotting ------------------------------------------------

# Parse GLM data with CI
df0 <- glm_flat_mean_ci

# Build core chain
all_cores <- paste0("Core", 1:5)
core_chain <- all_cores[1:TARGET_CORE]
core_chain <- intersect(core_chain, unique(df0$base_core))

# Core-only points
df_core_points <- df0 %>%
  filter(is_core_only, base_core %in% core_chain) %>%
  arrange(base_core) %>%
  select(model, base_core, R2_mean, R2_sd, R2_se, R2_ci_lower, R2_ci_upper)

# Target Core + NfL point
target_nfl_model <- paste0("Core", TARGET_CORE, "_NfL")
df_target_nfl <- df0 %>%
  filter(model == target_nfl_model) %>%
  mutate(base_core = paste0("Core", TARGET_CORE)) %>%
  select(model, base_core, R2_mean, R2_sd, R2_se, R2_ci_lower, R2_ci_upper)

# Combine left chain
nfl_stage_name <- paste0("Core", TARGET_CORE, " + NfL")
df_line <- bind_rows(df_core_points, df_target_nfl) %>%
  mutate(stage = factor(
    c(core_chain, if (nrow(df_target_nfl)) nfl_stage_name else NULL),
    levels = c(core_chain, nfl_stage_name)
  ))

# Interaction models
df_int <- df0 %>%
  filter(grepl(paste0("^Core", TARGET_CORE, "_NfL_"), model) & 
           grepl("_x", model)) %>%
  mutate(
    stage = factor(paste0("Core", TARGET_CORE, " + NfL × α (interaction)"),
                   levels = c(core_chain, nfl_stage_name,
                              paste0("Core", TARGET_CORE, " + NfL × α (interaction)"))),
    alpha_raw = sub(paste0("^Core", TARGET_CORE, "_NfL_"), "", model),
    alpha_raw = sub("^(.+?)_x.*$", "\\1", alpha_raw),
    alpha_short = shorten_alpha(alpha_raw)
  )

# Stage levels
stage_lvls <- c(core_chain, nfl_stage_name,
                paste0("Core", TARGET_CORE, " + NfL × α (interaction)"))

x_int_base <- length(core_chain) + 2

if (SHOW_ML_BLOCK) {
  x_ml_anchor <- length(core_chain) + 3
  stage_lvls <- c(stage_lvls, "ML")
} else {
  x_ml_anchor <- NULL
}

# Desired order for alpha labels
desired_order_short <- c("GFAP","pTau181","pTau217","Aβ38/40","Aβ42/40")

# Compute offsets
if (nrow(df_int) > 0) {
  off_int <- offset_map_fixed(factor(df_int$alpha_short, 
                                     levels = desired_order_short), 
                              width = 0.7)
}

# Position interaction cloud
label_pad_base <- 0.06

if (nrow(df_int) > 0) {
  df_int_sorted <- df_int %>%
    arrange(R2_mean)
  
  off_int <- offset_map_fixed(
    factor(df_int_sorted$alpha_short, levels = unique(df_int_sorted$alpha_short)), 
    width = 0.7
  )
  
  df_int_pos <- df_int_sorted %>%
    mutate(
      x_num = x_int_base + unname(off_int[as.character(alpha_short)]),
      y_lo = R2_ci_lower,
      extra_pad = get_extra_pad(alpha_short, y_lo, LABEL_PAD_INTERACTION, 
                                step = 0.012),
      y_label = y_lo - (label_pad_base + extra_pad)
    )
} else {
  df_int_pos <- data.frame()
}

# Position left chain
df_left_pos <- df_line %>%
  mutate(
    x_num = as.numeric(factor(stage, levels = stage_lvls)),
    type = ifelse(as.character(stage) %in% core_chain, "Core-only", 
                  paste0("Core + NfL"))
  )

# Core + NfL band for shading
band_core_nfl <- if (nrow(df_target_nfl)) {
  list(ymin = df_target_nfl$R2_ci_lower,
       ymax = df_target_nfl$R2_ci_upper)
} else list(ymin = NA, ymax = NA)

# ML block preparation
if (SHOW_ML_BLOCK) {
  
  #ml_res_sorted <- ml_flat_mean_ci %>% arrange(desc(R2_mean))
  ml_res_sorted <- ml_flat_mean_ci %>% arrange(R2_mean)
  
  ml_labels_pretty <- pretty_ml_label(ml_res_sorted$model)
  ml_offset_map <- offset_map_fixed(ml_labels_pretty, width = 0.7)
  
  ml_plot_data <- ml_res_sorted %>%
    transmute(
      ml_label = pretty_ml_label(model),
      R2_mean, R2_sd, R2_se, R2_ci_lower, R2_ci_upper,
      x_num = x_ml_anchor + unname(ml_offset_map[pretty_ml_label(model)]),
      y_lo = R2_ci_lower,
      y_hi = R2_ci_upper
    ) %>%
    mutate(
      extra_pad = get_extra_pad(ml_label, y_lo, LABEL_PAD_ML, step = 0.012),
      y_label = y_lo - (0.04 + extra_pad)
    )
} else {
  ml_plot_data <- data.frame()
}

# X-axis labels
x_breaks <- df_left_pos$x_num
x_labels <- sapply(as.character(df_line$stage), function(stage) {
  if (stage %in% names(CORE_LABELS)) {
    CORE_LABELS[stage]
  } else {
    stage
  }
})

x_breaks <- c(x_breaks, x_int_base)
x_labels <- c(x_labels, 
              paste0("Core", TARGET_CORE, " + NfL × α (interaction)"))

if (SHOW_ML_BLOCK) {
  x_breaks <- c(x_breaks, x_ml_anchor)
  x_labels <- c(x_labels, "ML")
}

x_angle <- 20

# Plot with 95% CI ------------------------------------------------------
# The error bars summarise variability across repeated CV runs.
# They should not be interpreted as participant-level sampling confidence intervals.

Fig4D <- ggplot() +
  geom_line(data = df_left_pos, aes(x = x_num, y = R2_mean, group = 1),
            linewidth = 1.2, color = "black") +
  geom_point(data = df_left_pos, aes(x = x_num, y = R2_mean, shape = type),
             size = 3.8, color = "black", fill = "white", stroke = 1) +
  geom_errorbar(data = df_left_pos,
                aes(x = x_num, 
                    ymin = pmax(0, R2_ci_lower), 
                    ymax = pmin(1, R2_ci_upper)),
                width = 0.12, linewidth = 0.7, color = "black") +
  {if (nrow(df_int_pos) > 0) list(
    geom_point(data = df_int_pos, aes(x = x_num, y = R2_mean),
               size = 3.2, color = "#F58518"),
    geom_errorbar(data = df_int_pos,
                  aes(x = x_num, y = R2_mean,
                      ymin = pmax(0, R2_ci_lower), 
                      ymax = pmin(1, R2_ci_upper)),
                  width = 0.05, linewidth = 0.6, color = "#F58518", 
                  alpha = 0.8),
    geom_segment(data = df_int_pos, 
                 aes(x = x_num, xend = x_num, y = y_lo, yend = y_label),
                 color = "#C8660E", linetype = "dotted", linewidth = 0.5, 
                 alpha = 0.75),
    geom_text(data = df_int_pos, 
              aes(x = x_num, y = y_label, label = alpha_short),
              color = "#C8660E", angle = x_angle, hjust = 0.5, vjust = 1, 
              size = 5)
  )} +
  {if (SHOW_ML_BLOCK && nrow(ml_plot_data) > 0) list(
    geom_point(data = ml_plot_data, aes(x = x_num, y = R2_mean), 
               size = 3.2, color = "#6C5CE7"),
    geom_errorbar(data = ml_plot_data,
                  aes(x = x_num, y = R2_mean, ymin = y_lo, ymax = y_hi),
                  width = 0.05, linewidth = 0.6, color = "#6C5CE7", 
                  alpha = 0.85),
    geom_segment(data = ml_plot_data, 
                 aes(x = x_num, xend = x_num, y = y_lo, yend = y_label),
                 color = "#5A4FD6", linetype = "dotted", linewidth = 0.5, 
                 alpha = 0.75),
    geom_text(data = ml_plot_data, 
              aes(x = x_num, y = y_label, label = ml_label),
              color = "#5A4FD6", angle = x_angle, hjust = 0.5, vjust = 1, 
              size = 5)
  )} +
  scale_x_continuous(
    breaks = x_breaks,
    labels = x_labels,
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  scale_y_continuous(
    breaks = seq(0, 0.3, by = 0.05),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  labs(
    x = NULL, 
    y = "Cross-validated R²",
    title = "Predictive performance for ALSFRS-R slope in nested cross-validation"
    #subtitle = "Error bars represent 95% confidence intervals (n = 10 repeats, df = 9)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 14, hjust = 0.5, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16, margin = ggplot2::margin(r = 15)),
    plot.title = element_text(face = "bold", size = 20),
    plot.subtitle = element_text(size = 13),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.margin = ggplot2::margin(10, 10, 10, 10)
  )


# 📊 Figure 4D ------------------------------------------------------------------
# Note: These error bars/CIs reflect variability across repeated CV runs,
# not participant-level sampling confidence intervals.

Fig4D


# Sanity check ------------------------------------------------------------

nested_rep_tbl %>%
  dplyr::mutate(
    r2_cor = purrr::map_dbl(oof_pred_mean, ~ cor(.x, df$slope_total)^2),
    r2_oos = purrr::map_dbl(oof_pred_mean, ~ {
      rss <- sum((df$slope_total - .x)^2)
      tss <- sum((df$slope_total - mean(df$slope_total))^2)
      1 - rss/tss
    }),
    bias   = purrr::map_dbl(oof_pred_mean, ~ mean(.x - df$slope_total)),
    slope  = purrr::map_dbl(oof_pred_mean, ~ coef(lm(df$slope_total ~ .x))[2])
  ) %>%
  dplyr::select(model, r2_cor, r2_oos, bias, slope)

