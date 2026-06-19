source("scripts_pub/03b_ML_2step_oof.R")

# Functions -------------------------------------------------------------------

# Pick a representative repeat by minimizing the L2 distance to per-metric medians aggregated ACROSS methods (method-agnostic).
# Metrics used: overall AUC, gray AUC, gray PPV, gray NPV, gray coverage.
pick_representative_repeat_allmethods <- function(rep_out, variant = "sp95") {
  runs <- rep_out$runs
  stopifnot(length(runs) > 0)
  
  rows <- lapply(seq_along(runs), function(i) {
    per_method <- runs[[i]]$per_method
    # collect metrics for each method in this repeat
    ms <- lapply(names(per_method), function(m) {
      pm <- per_method[[m]]
      vo <- pm$pooled$variants_overall
      go <- pm$pooled$gray_overall
      bc <- pm$pooled$bucket_counts
      
      vo <- vo[vo$variant == variant, , drop = FALSE]
      go <- go[go$variant == variant, , drop = FALSE]
      if (!"PPV_gray" %in% names(go) && "PPV" %in% names(go)) {
        names(go)[names(go)=="PPV"]  <- "PPV_gray"
        names(go)[names(go)=="NPV"]  <- "NPV_gray"
        names(go)[names(go)=="Sens"] <- "Sens_gray"
        names(go)[names(go)=="Spec"] <- "Spec_gray"
      }
      cov_g <- bc$prop[bc$bucket=="gray"]
      
      data.frame(
        repeat_id = i, method = m,
        AUC = vo$AUC[1], AUC_gray = vo$AUC_gray[1],
        PPV_gray = go$PPV_gray[1], NPV_gray = go$NPV_gray[1],
        gray_cov = cov_g,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, ms)
  })
  long <- dplyr::bind_rows(rows)
  
  # median per metric across all methods & repeats
  med <- long |>
    tidyr::pivot_longer(cols = c(AUC, AUC_gray, PPV_gray, NPV_gray, gray_cov),
                        names_to = "metric", values_to = "value") |>
    dplyr::group_by(metric) |>
    dplyr::summarise(med = median(value, na.rm = TRUE), .groups = "drop")
  
  # distance per (repeat, method), then average across methods
  long2 <- long |>
    tidyr::pivot_longer(cols = c(AUC, AUC_gray, PPV_gray, NPV_gray, gray_cov),
                        names_to = "metric", values_to = "value") |>
    dplyr::left_join(med, by = "metric") |>
    dplyr::mutate(sq = (value - med)^2) |>
    dplyr::group_by(repeat_id, method) |>
    dplyr::summarise(dist = sqrt(sum(sq, na.rm = TRUE)), .groups = "drop") |>
    dplyr::group_by(repeat_id) |>
    dplyr::summarise(mean_dist = mean(dist, na.rm = TRUE), .groups = "drop")
  
  long2$repeat_id[which.min(long2$mean_dist)]
}

# --- pick a "representative repeat" close to per-metric medians ---
# Inputs:
#   rep_out      : returned by twostep_kfold_multi_repeat()
#   variant      : "sp95" / "sp90" / "youden" etc.
#   anchor_method: choose a method to anchor the representativeness (e.g., "rf")
# Returns: integer repeat_id
pick_representative_repeat <- function(rep_out, variant = "sp95", anchor_method = "rf") {
  runs <- rep_out$runs
  stopifnot(length(runs) > 0, anchor_method %in% names(runs[[1]]$per_method))
  
  # gather per-repeat metrics for the anchor method
  rows <- lapply(seq_along(runs), function(i) {
    pm <- runs[[i]]$per_method[[anchor_method]]
    vo <- pm$pooled$variants_overall
    go <- pm$pooled$gray_overall
    bc <- pm$pooled$bucket_counts
    # overall
    vo <- vo[vo$variant == variant, , drop = FALSE]
    # gray
    go <- go[go$variant == variant, , drop = FALSE]
    # normalize gray names if needed
    if (!"PPV_gray" %in% names(go) && "PPV" %in% names(go)) {
      names(go)[names(go)=="PPV"] <- "PPV_gray"
      names(go)[names(go)=="NPV"] <- "NPV_gray"
      names(go)[names(go)=="Sens"] <- "Sens_gray"
      names(go)[names(go)=="Spec"] <- "Spec_gray"
    }
    cov_g <- bc$prop[bc$bucket=="gray"]
    data.frame(
      repeat_id = i,
      AUC       = vo$AUC[1],
      AUC_gray  = vo$AUC_gray[1],
      PPV_gray  = go$PPV_gray[1],
      NPV_gray  = go$NPV_gray[1],
      gray_cov  = cov_g,
      stringsAsFactors = FALSE
    )
  })
  df <- dplyr::bind_rows(rows)
  
  # median vector
  med <- vapply(df[, c("AUC","AUC_gray","PPV_gray","NPV_gray","gray_cov")], median, numeric(1), na.rm = TRUE)
  
  # L2 distance to medians
  dist <- apply(df[, names(med)], 1, function(x) sqrt(sum((x - med)^2, na.rm = TRUE)))
  df$dist_to_med <- dist
  
  # argmin
  df$repeat_id[which.min(df$dist_to_med)]
}


# Convert one repeat-run into a km-like structure (compare tables)
as_km_like <- function(one_run) {
  per_method <- one_run$per_method
  stopifnot(is.list(per_method))
  
  # 1) overall
  overall_rows <- lapply(names(per_method), function(m) {
    ro <- per_method[[m]]$pooled$variants_overall
    ro$method <- m
    ro[, c("method","variant","AUC","AUC_gray","Sens","Spec","PPV","NPV")]
  })
  overall_cmp <- dplyr::bind_rows(overall_rows) %>%
    dplyr::arrange(variant, dplyr::desc(AUC), dplyr::desc(Sens))
  
  # 2) gray-only (confusion-based)
  gray_rows <- lapply(names(per_method), function(m) {
    rg <- per_method[[m]]$pooled$gray_overall
    # normalize if needed
    if (!"PPV_gray" %in% names(rg) && "PPV" %in% names(rg)) {
      names(rg)[names(rg)=="PPV"] <- "PPV_gray"
      names(rg)[names(rg)=="NPV"] <- "NPV_gray"
      names(rg)[names(rg)=="Sens"] <- "Sens_gray"
      names(rg)[names(rg)=="Spec"] <- "Spec_gray"
    }
    rg$method <- m
    rg[, c("method","variant","gray_n","TP","TN","FP","FN","Sens_gray","Spec_gray","PPV_gray","NPV_gray")]
  })
  gray_cmp <- dplyr::bind_rows(gray_rows)
  
  # 3) bucket coverage (identical across methods within a repeat)
  any_m <- names(per_method)[1]
  bucket_counts <- per_method[[any_m]]$pooled$bucket_counts
  
  # 4) NfL baselines (with OOF AUC + threshold summaries)
  nfl_baselines <- per_method[[any_m]]$pooled$nfl_baselines
  nfl_roc       <- per_method[[any_m]]$pooled$nfl_roc
  
  list(
    meta = list(methods_gray = names(per_method)),
    per_method = per_method,
    compare = list(
      overall   = overall_cmp,
      gray_only = gray_cmp,
      buckets   = bucket_counts,
      nfl       = list(baselines = nfl_baselines, roc = nfl_roc)
    )
  )
}


# Build SHAP "handles" for gray-only models from one repeat-run.
# Each handle matches compute_shap_for_handle_cls() input:
#   model, method, predictors, train_idx, test_idx, rep_id, fold_id, model_name, pos_class
build_gray_shap_handles_from_rep_run <- function(one_run, rep_id, pos_class = "ALS") {
  out <- list()
  k <- 0L
  for (m in names(one_run$per_method)) {
    folds <- one_run$per_method[[m]]$per_fold
    if (!length(folds)) next
    
    for (f in seq_along(folds)) {
      ts2 <- folds[[f]]
      
      mdl   <- ts2$gray_model$fit
      preds <- ts2$data$predictors_gray
      
      # NOTE: twostep_run_on_split() uses clean_side() which created row_id for train/test
      # tr_ids <- ts2$data$train_df$row_id
      # te_ids <- ts2$data$test_df$row_id
      
      # --- gray only ---
      buck <- ts2$step1$test_step1$bucket
      gray_idx_local <- which(buck == "gray")
      if (!length(gray_idx_local)) next  # if no gray in this fold, skip
      
      # tr_ids <- ts2$data$train_df$row_id
      # te_ids <- ts2$data$test_df$row_id[gray_idx_local]
      
      # Use global row indices from the outer CV fold structure. --------------------#
      # The row_id inside ts2$data$train_df and ts2$data$test_df is local to each split.
      global_test_idx <- one_run$per_method[[m]]$data$folds[[f]]
      global_train_idx <- setdiff(
        seq_len(nrow(one_run$per_method[[m]]$data$df_clean)),
        global_test_idx
      )
      
      df_global <- one_run$per_method[[m]]$data$df_clean
      
      cutoff_low  <- ts2$step1$cutoff_low
      cutoff_high <- ts2$step1$cutoff_high
      
      gray_train_idx_global <- global_train_idx[
        df_global$NfL_csf_pgml[global_train_idx] > cutoff_low &
          df_global$NfL_csf_pgml[global_train_idx] < cutoff_high
      ]
      
      tr_ids <- gray_train_idx_global
      te_ids <- global_test_idx[gray_idx_local]
      # --- End of global index mapping ----------------------------------------------#
      
      k <- k + 1L
      out[[k]] <- list(
        model       = mdl,
        method      = m,
        predictors  = preds,
        train_idx   = tr_ids,
        test_idx    = te_ids,
        rep_id      = as.character(rep_id),
        fold_id     = as.character(f),
        model_name  = paste0("gray_", m),
        pos_class   = pos_class
      )
    }
  }
  out
}

# RUN ---------------------------------------------------------------------

rep_id <- selected_ppv$run_id
run_rep <- rep_out$runs[[rep_id]]
shap_handles <- build_gray_shap_handles_from_rep_run(run_rep, rep_id, pos_class = "ALS")


# XGB gray SHAP
TwostepXGBshap <- plot_shap_for_repeat_cls_from_handles(
  shap_handles = shap_handles,
  df = run_rep$per_method[[1]]$data$df_clean %||% rep_out$runs[[rep_id]]$per_method[[1]]$data$df_clean,
  method = "gray_xgb",   
  rep_id = rep_id,
  pos_class = "ALS",
  nsim = 80, bg_n = 200, test_subsample = 200
)

TwostepXGBshap

# RF gray SHAP
TwostepRFshap <- plot_shap_for_repeat_cls_from_handles(
  shap_handles = shap_handles,
  df = run_rep$per_method[[1]]$data$df_clean %||% rep_out$runs[[rep_id]]$per_method[[1]]$data$df_clean,
  method = "gray_rf",   
  rep_id = rep_id,
  pos_class = "ALS",
  nsim = 80, bg_n = 200, test_subsample = 200
)

TwostepRFshap | TwostepXGBshap


# 📊 Figure 3D ------------------------------------------------------------------

Fig3D <- TwostepRFshap +
  scale_x_discrete(labels = function(x) ifelse(x %in% names(variable_label_map), variable_label_map[x], x))+
  theme_bw(base_size = 16) +
  labs(
    title = "Feature contributions within the NfL gray-zone (RF model)",
    y = "SHAP value (impact on prediction)"
  ) +
  theme(
    plot.title   = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 14),
    axis.title.x  = element_text(margin=ggplot2::margin(t=15,b=-10)),
    axis.text.y  = element_text(size = 14),
    axis.text.x  = element_text(size = 14),
    legend.title = element_text(size = 12, margin=ggplot2::margin(r=25,b=15)),
    legend.text  = element_text(size = 10),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width  = unit(1, "cm"),
    legend.position = "bottom"
  )

Fig3D



# QC/ Verify test set ALS labels ------------------------------
 
# df_check5 <- run_rep$per_method$rf$per_fold[[5]]$data$test_df
# df_check4 <- run_rep$per_method$rf$per_fold[[4]]$data$test_df
# df_check3 <- run_rep$per_method$rf$per_fold[[3]]$data$test_df
# df_check2 <- run_rep$per_method$rf$per_fold[[2]]$data$test_df
# df_check1 <- run_rep$per_method$rf$per_fold[[1]]$data$test_df
# 
# bind_rows(
#   table(df_check1$ALS_label),
#   table(df_check2$ALS_label),
#   table(df_check3$ALS_label),
#   table(df_check4$ALS_label),
#   table(df_check5$ALS_label)
# )
