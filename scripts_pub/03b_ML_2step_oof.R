suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(ggpattern)
})


# Utilities: OOF-based metrics and PPV-based representative selection --------

# OOF overall prediction for a two-step workflow (high→pos, low→neg, gray→final_pred)
oof_overall_pred <- function(kf_result, variant = "sp95") {
  y_true    <- kf_result$pooled$oof$y_true
  bucket    <- kf_result$pooled$oof$bucket
  final_g   <- kf_result$pooled$oof$final_pred[[variant]]
  
  pos <- levels(y_true)[2]; neg <- levels(y_true)[1]
  pred <- rep(NA_character_, length(y_true))
  pred[bucket == "high"] <- pos
  pred[bucket == "low"]  <- neg
  pred[bucket == "gray"] <- as.character(final_g[bucket == "gray"])
  factor(pred, levels = c(neg, pos))
}

# Calculate PPV/NPV from OOF predictions
oof_ppv_npv <- function(kf_result, variant = "sp95") {
  y_true <- kf_result$pooled$oof$y_true
  pred   <- oof_overall_pred(kf_result, variant)
  pos <- levels(y_true)[2]; neg <- levels(y_true)[1]
  TP <- sum(y_true == pos & pred == pos)
  FP <- sum(y_true == neg & pred == pos)
  TN <- sum(y_true == neg & pred == neg)
  FN <- sum(y_true == pos & pred == neg)
  tibble(TP = TP, FP = FP, TN = TN, FN = FN,
         PPV = if ((TP+FP)>0) TP/(TP+FP) else NA_real_,
         NPV = if ((TN+FN)>0) TN/(TN+FN) else NA_real_)
}

# Select representative run based on PPV/NPV median proximity
select_representative_run_by <- function(rep_out,
                                         method  = "rf",
                                         variant = "sp95",
                                         metric  = c("PPV","NPV")) {
  metric <- match.arg(metric)
  rows <- lapply(seq_along(rep_out$runs), function(i) {
    kf <- rep_out$runs[[i]]$per_method[[method]]
    mn <- oof_ppv_npv(kf, variant = variant)
    tibble(run_id = i, seed = rep_out$used_seeds[i],
           PPV = mn$PPV, NPV = mn$NPV)
  })
  df <- bind_rows(rows)
  target <- if (metric == "PPV") df$PPV else df$NPV
  med <- median(target, na.rm = TRUE)
  best_idx <- which.min(abs(target - med))
  best <- df[best_idx, , drop = FALSE]
  message(sprintf("Selected run %d (seed=%d) by %s at %s: value=%.3f (median=%.3f)",
                  best$run_id, best$seed, metric, variant, target[best_idx], med))
  list(run_id = best$run_id,
       seed   = best$seed,
       kf     = rep_out$runs[[best$run_id]]$per_method[[method]],
       table  = df,
       metric = metric,
       variant = variant)
}


# Waterfall data builders (OOF-based, separated "All" row) ------------

# Extract Step1 waterfall data, separated into "All" and classification buckets
extract_twostep_step1_oof <- function(kf_result, variant = "sp95") {
  y_true <- kf_result$pooled$oof$y_true
  bucket <- kf_result$pooled$oof$bucket
  
  pos_class <- levels(y_true)[2]
  neg_class <- levels(y_true)[1]
  
  # "All" row
  all_row <- tibble(
    stage = "Step 1",
    group = sprintf("All\n(n=%d)", length(y_true)),
    n_total = length(y_true),
    n_als = sum(y_true == pos_class),
    n_dc = sum(y_true == neg_class)
  )
  
  # Classification buckets only
  buckets_rows <- bind_rows(
    tibble(
      stage = "Step 1",
      group = sprintf("95%% Sp.\nRule-in\n(n=%d)", sum(bucket == "high")),
      n_total = sum(bucket == "high"),
      n_als = sum(y_true == pos_class & bucket == "high"),
      n_dc = sum(y_true == neg_class & bucket == "high"),
      PPV = n_als / n_total
    ),
    tibble(
      stage = "Step 1",
      group = sprintf("Gray-zone\n(n=%d)", sum(bucket == "gray")),
      n_total = sum(bucket == "gray"),
      n_als = sum(y_true == pos_class & bucket == "gray"),
      n_dc = sum(y_true == neg_class & bucket == "gray")
    ),
    tibble(
      stage = "Step 1",
      group = sprintf("95%% Se.\nRule-out\n(n=%d)", sum(bucket == "low")),
      n_total = sum(bucket == "low"),
      n_als = sum(y_true == pos_class & bucket == "low"),
      n_dc = sum(y_true == neg_class & bucket == "low"),
      NPV = n_dc / n_total
    )
  ) |>
    mutate(
      prop_als = n_als / n_total,
      prop_dc  = n_dc  / n_total,
      stage = factor(stage),
      group = factor(group, levels = group)
    )
  
  all_row <- all_row |>
    mutate(
      prop_als = n_als / n_total,
      prop_dc  = n_dc  / n_total,
      stage = factor(stage),
      group = factor(group, levels = group)
    )
  
  list(all = all_row, buckets = buckets_rows)
}

# Extract Step2 waterfall data (gray zone resolution only, no "Gray zone patients" row)
extract_twostep_step2_oof <- function(kf_result, variant = "sp95") {
  y_true <- kf_result$pooled$oof$y_true
  bucket <- kf_result$pooled$oof$bucket
  final_pred <- kf_result$pooled$oof$final_pred[[variant]]
  
  pos_class <- levels(y_true)[2]
  neg_class <- levels(y_true)[1]
  
  gray_idx <- bucket == "gray"
  
  gray_to_pos <- gray_idx & final_pred == pos_class
  gray_to_neg <- gray_idx & final_pred == neg_class
  
  bind_rows(
    tibble(
      stage = "Step 2",
      group = sprintf("Additional tests (+)\n(n=%d)", sum(gray_to_pos)),
      n_total = sum(gray_to_pos),
      n_als = sum(y_true[gray_to_pos] == pos_class),
      n_dc = sum(y_true[gray_to_pos] == neg_class)
    ),
    tibble(
      stage = "Step 2",
      group = sprintf("Additional tests (-)\n(n=%d)", sum(gray_to_neg)),
      n_total = sum(gray_to_neg),
      n_als = sum(y_true[gray_to_neg] == pos_class),
      n_dc = sum(y_true[gray_to_neg] == neg_class)
    )
  ) |>
    mutate(
      prop_als = n_als / n_total,
      prop_dc  = n_dc  / n_total,
      stage = factor(stage),
      group = factor(group, levels = group)
    )
}

# Extract overall waterfall (combined rule-in, careful follow, rule-out), separated "All"
extract_overall_oof <- function(kf_result, variant = "sp95") {
  y_true     <- kf_result$pooled$oof$y_true
  bucket     <- as.character(kf_result$pooled$oof$bucket)
  final_pred <- as.character(kf_result$pooled$oof$final_pred[[variant]])
  
  lv  <- levels(y_true)
  pos <- if ("ALS" %in% lv) "ALS" else lv[2]
  neg <- setdiff(lv, pos)[1]
  
  idx_high     <- bucket == "high"
  idx_low      <- bucket == "low"
  idx_gray     <- bucket == "gray"
  idx_gray_pos <- idx_gray & final_pred == pos
  idx_gray_neg <- idx_gray & final_pred == neg
  
  idx_combined_ri <- idx_high | idx_gray_pos
  idx_careful     <- idx_gray_neg
  idx_ruleout     <- idx_low
  
  # "All" row
  all_row <- tibble(
    stage   = "Overall",
    group   = sprintf("All\n(n=%d)", length(y_true)),
    n_total = length(y_true),
    n_als   = sum(y_true == pos),
    n_dc    = sum(y_true == neg)
  )
  
  # Classification rows only
  build_row <- function(lbl, idx) {
    tibble(
      stage   = "Overall",
      group   = sprintf("%s\n(n=%d)", lbl, sum(idx, na.rm = TRUE)),
      n_total = sum(idx, na.rm = TRUE),
      n_als   = sum(y_true[idx] == pos, na.rm = TRUE),
      n_dc    = sum(y_true[idx] == neg, na.rm = TRUE)
    )
  }
  
  class_rows <- bind_rows(
    build_row("Rule-in", idx_combined_ri),
    build_row("Still in gray-zone",   idx_careful),
    build_row("Rule-out",         idx_ruleout)
  ) |>
    mutate(
      prop_als = n_als / n_total,
      prop_dc  = n_dc  / n_total,
      stage = factor(stage),
      group = factor(group, levels = group)
    )
  
  all_row <- all_row |>
    mutate(
      prop_als = n_als / n_total,
      prop_dc  = n_dc  / n_total,
      stage = factor(stage),
      group = factor(group, levels = group)
    )
  
  list(all = all_row, classification = class_rows)
}

# Extract NfL one-step waterfall (fold-wise thresholds), separated "All"
extract_nfl_onestep_oof <- function(kf_result,
                                    variant = c("sp95","se95","youden"),
                                    outcome = "ALS_label",
                                    nfl_var = "NfL_csf_pgml") {
  variant <- match.arg(variant)
  base_tbl <- kf_result$pooled$nfl_baselines
  row_idx <- switch(variant,
                    sp95   = which(base_tbl$baseline == "nfl_sp95"),
                    se95   = which(base_tbl$baseline == "nfl_se95"),
                    youden = which(base_tbl$baseline == "nfl_youden"))
  th_list <- base_tbl$th_train_all[[row_idx]]
  
  df_clean <- kf_result$data$df_clean
  folds    <- kf_result$data$folds
  pos <- "ALS"; neg <- "nonALS"
  truth <- factor(df_clean[[outcome]], levels = c(neg, pos))
  
  # Apply fold-wise thresholds
  pred  <- rep(NA_character_, nrow(df_clean))
  for (i in seq_along(folds)) {
    idx <- folds[[i]]; th <- th_list[i]
    pred[idx] <- ifelse(df_clean[[nfl_var]][idx] >= th, pos, neg)
  }
  pred <- factor(pred, levels = c(neg, pos))
  
  TP <- sum(truth == pos & pred == pos)
  FP <- sum(truth == neg & pred == pos)
  TN <- sum(truth == neg & pred == neg)
  FN <- sum(truth == pos & pred == neg)
  
  # "All" row
  all_row <- tibble(
    stage = "NfL only",
    group = sprintf("All\n(n=%d)", length(truth)),
    n_total = length(truth),
    n_als   = sum(truth == pos),
    n_dc    = sum(truth == neg),
    prop_als = n_als / n_total,
    prop_dc  = n_dc  / n_total
  )
  
  # Classification rows only
  class_rows <- tibble(
    stage = "NfL only",
    group = c(
      sprintf("Rule-in\n(n=%d)", TP + FP),
      sprintf("Rule-out\n(n=%d)", FN + TN)
    ),
    n_total = c(TP + FP, FN + TN),
    n_als   = c(TP, FN),
    n_dc    = c(FP, TN),
    prop_als = n_als / n_total,
    prop_dc  = n_dc  / n_total,
    PPV = c(TP/(TP+FP), NA_real_),
    NPV = c(NA_real_, TN/(TN+FN))
  ) |>
    mutate(stage = factor(stage),
           group = factor(group, levels = group))
  
  all_row <- all_row |>
    mutate(stage = factor(stage),
           group = factor(group, levels = group))
  
  list(all = all_row, classification = class_rows)
}


# PPV-based representative selection and plot generation --------------

# source("scripts_pub/03a_ML_2step_train.R")

# Select representative run based on PPV median (Sp95)
selected_ppv <- select_representative_run_by(rep_out, method = "rf", variant = "sp95", metric = "PPV")
kf_selected <- selected_ppv$kf

kf_selected$pooled$step1_thresholds
selected_ppv$run_id
selected_ppv$seed

# Display OOF PPV/NPV for selected run
print(oof_ppv_npv(kf_selected, variant = "sp95"))

## Data extraction for plots (OOF-based) -------------------------------

## Extract Step1/2 data 
step1_data <- extract_twostep_step1_oof(kf_selected, variant = "sp95")
step2_data <- extract_twostep_step2_oof(kf_selected, variant = "sp95")

## Extract Overall data 
overall_data <- extract_overall_oof(kf_selected, variant = "sp95")

## Extract NfL one-step data (Sp/Se)
nfl_sp95_data <- extract_nfl_onestep_oof(kf_selected, variant = "sp95")
nfl_se95_data <- extract_nfl_onestep_oof(kf_selected, variant = "se95")



