# ========================================================= #
# Two-step diagnostic pipeline (single 80/20 split)
# ========================================================== #

source("scripts_pub/function/ML_class_func.R")

twostep_run <- function(df,
                        outcome = "ALS_label",
                        predictors_all,
                        predictors_gray,
                        method_gray = "rf",
                        train_frac = 0.80,
                        seed_split = 12345,
                        seed_inner = 12345,
                        target_spec_high = 0.95,
                        target_sens_low  = 0.95,
                        gray_spec_targets = c(0.95, 0.90),
                        gray_sens_targets = NULL) {  
  # --------------------- small internal helpers ------------------------
  
  # Pick threshold at (or near) a target specificity on a pROC::roc
  pick_threshold_by_spec <- function(roc_obj, target_spec = 0.95,
                                     policy = c("at_least","nearest","at_most")) {
    policy <- match.arg(policy)
    tab <- pROC::coords(roc_obj, x="all",
                        ret=c("threshold","specificity","sensitivity"),
                        transpose = FALSE)
    tab <- as.data.frame(tab)
    suppressWarnings(tab <- tab[order(tab$threshold), ])
    if (policy == "at_least") {
      cand <- tab[tab$specificity >= target_spec, , drop = FALSE]
      if (nrow(cand) == 0L) return(NA_real_)
      return(cand$threshold[which.max(cand$sensitivity)])
    }
    if (policy == "at_most") {
      cand <- tab[tab$specificity <= target_spec, , drop = FALSE]
      if (nrow(cand) == 0L) return(NA_real_)
      return(cand$threshold[which.max(cand$specificity)])
    }
    idx <- which.min(abs(tab$specificity - target_spec))
    return(tab$threshold[idx])
  }
  pick_threshold_by_sens <- function(roc_obj, target_sens = 0.95,
                                     policy = c("at_least","nearest","at_most")) {
    policy <- match.arg(policy)
    tab <- pROC::coords(roc_obj, x="all",
                        ret=c("threshold","specificity","sensitivity"),
                        transpose = FALSE)
    tab <- as.data.frame(tab)
    suppressWarnings(tab <- tab[order(tab$threshold), ])
    if (policy == "at_least") {
      cand <- tab[tab$sensitivity >= target_sens, , drop = FALSE]
      if (nrow(cand) == 0L) return(NA_real_)
      return(cand$threshold[which.max(cand$specificity)])
    }
    if (policy == "at_most") {
      cand <- tab[tab$sensitivity <= target_sens, , drop = FALSE]
      if (nrow(cand) == 0L) return(NA_real_)
      return(cand$threshold[which.max(cand$sensitivity)])
    }
    idx <- which.min(abs(tab$sensitivity - target_sens))
    return(tab$threshold[idx])
  }
  
  compute_by_bucket_tbl <- function(y_true, bucket, final_pred,
                                    note_lr_plus_high = NULL,
                                    note_lr_minus_low = NULL) {
    stopifnot(length(y_true) == length(final_pred), length(bucket) == length(y_true))
    neg <- levels(y_true)[1]; pos <- levels(y_true)[2]
    idx_g <- which(bucket == "gray"); idx_l <- which(bucket == "low"); idx_h <- which(bucket == "high")
    # gray confusion
    if (length(idx_g) > 0) {
      TPg <- sum(y_true[idx_g] == pos & final_pred[idx_g] == pos)
      TNg <- sum(y_true[idx_g] == neg & final_pred[idx_g] == neg)
      FPg <- sum(y_true[idx_g] == neg & final_pred[idx_g] == pos)
      FNg <- sum(y_true[idx_g] == pos & final_pred[idx_g] == neg)
      sens_g <- ifelse((TPg+FNg)>0, TPg/(TPg+FNg), NA_real_)
      spec_g <- ifelse((TNg+FPg)>0, TNg/(TNg+FPg), NA_real_)
      PPV_g  <- ifelse((TPg+FPg)>0, TPg/(TPg+FPg), NA_real_)
      NPV_g  <- ifelse((TNg+FNg)>0, TNg/(TNg+FNg), NA_real_)
    } else {
      TPg <- TNg <- FPg <- FNg <- NA_integer_
      sens_g <- spec_g <- PPV_g <- NPV_g <- NA_real_
    }
    band_metrics <- function(idx, band_call) {
      if (length(idx) == 0) return(list(n=0L, TP=NA, TN=NA, FP=NA, FN=NA, PPV=NA_real_, NPV=NA_real_))
      yb <- y_true[idx]
      pb <- factor(rep(band_call, length(idx)), levels = levels(y_true))
      TPb <- sum(yb == pos & pb == pos); TNb <- sum(yb == neg & pb == neg)
      FPb <- sum(yb == neg & pb == pos); FNb <- sum(yb == pos & pb == neg)
      list(
        n   = length(idx),
        TP  = TPb, TN = TNb, FP = FPb, FN = FNb,
        PPV = ifelse((TPb+FPb)>0, TPb/(TPb+FPb), NA_real_),
        NPV = ifelse((TNb+FNb)>0, TNb/(TNb+FNb), NA_real_)
      )
    }
    bm_low  <- band_metrics(idx_l, levels(y_true)[1]) # neg
    bm_high <- band_metrics(idx_h, levels(y_true)[2]) # pos
    tibble::tibble(
      bucket = c("low","gray","high"),
      n      = c(bm_low$n, length(idx_g), bm_high$n),
      PPV    = c(bm_low$PPV,  PPV_g,  bm_high$PPV),
      NPV    = c(bm_low$NPV,  NPV_g,  bm_high$NPV),
      Sens   = c(NA_real_,    sens_g, NA_real_),
      Spec   = c(NA_real_,    spec_g, NA_real_),
      TP     = c(bm_low$TP,  TPg,  bm_high$TP),
      TN     = c(bm_low$TN,  TNg,  bm_high$TN),
      FP     = c(bm_low$FP,  FPg,  bm_high$FP),
      FN     = c(bm_low$FN,  FNg,  bm_high$FN),
      LR_note = c(note_lr_minus_low %||% NA_character_,
                  NA_character_,
                  note_lr_plus_high %||% NA_character_)
    )
  }
  
  `%||%` <- function(x, y) if (is.null(x)) y else x
  neg <- paste0("non", pos_class); pos <- pos_class
  
  # --------------------- 1) clean + split ---------------------------------
  df_clean <- df %>%
    dplyr::select(all_of(c(outcome, predictors_all))) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(
      row_id    = dplyr::row_number(),
      !!outcome := factor(ifelse(.data[[outcome]] == 1, pos, neg),
                          levels = c(neg, pos))
    )
  set.seed(seed_split)
  train_idx <- caret::createDataPartition(df_clean[[outcome]], p = train_frac, list = FALSE)
  train_df  <- df_clean[train_idx, ]
  test_df   <- df_clean[-train_idx, ]
  stopifnot(nrow(train_df) + nrow(test_df) == nrow(df_clean))
  
  # --------------------- 2) TRAIN ROC for NfL (direction safe) ------------
  roc_auto <- pROC::roc(
    response  = factor(train_df[[outcome]], levels = c(neg, pos)),
    predictor = train_df$NfL_csf_pgml,
    direction = "auto"
  )
  if (roc_auto$direction != "<") {
    roc_nfl <- pROC::roc(
      response  = factor(train_df[[outcome]], levels = c(neg, pos)),
      predictor = -train_df$NfL_csf_pgml,
      direction = "<"
    )
    nfl_sign <- -1
  } else {
    roc_nfl  <- roc_auto
    nfl_sign <-  1
  }
  # High-specificity side (Spec >= target_spec_high), tie-break by Sens
  roc_tbl <- tibble::tibble(
    threshold   = roc_nfl$thresholds,
    sensitivity = roc_nfl$sensitivities,
    specificity = roc_nfl$specificities
  ) %>% dplyr::filter(is.finite(threshold))
  
  pick_target <- function(tbl, target_col, target_val, tie_col) {
    hit <- tbl %>%
      dplyr::filter(.data[[target_col]] >= target_val) %>%
      dplyr::arrange(dplyr::desc(.data[[tie_col]])) %>%
      dplyr::slice_head(n = 1)
    if (nrow(hit) == 0) {
      hit <- tbl %>%
        dplyr::mutate(gap = abs(.data[[target_col]] - target_val)) %>%
        dplyr::arrange(gap, dplyr::desc(.data[[tie_col]])) %>%
        dplyr::slice_head(n = 1)
    }
    hit
  }
  pick_high <- pick_target(roc_tbl, "specificity", target_spec_high, "sensitivity")
  pick_low  <- pick_target(roc_tbl, "sensitivity", target_sens_low,  "specificity")
  cutoff_high <- nfl_sign * pick_high$threshold
  cutoff_low  <- nfl_sign * pick_low$threshold
  auc_train_nfl <- as.numeric(pROC::auc(roc_nfl))
  
  # --------------------- 3) gray-only ML on TRAIN -------------------------
  gray_zone_train <- train_df %>%
    dplyr::filter(NfL_csf_pgml > cutoff_low, NfL_csf_pgml < cutoff_high)
  stopifnot(length(unique(gray_zone_train[[outcome]])) == 2)
  
  set.seed(seed_inner)
  fit_gray <- fit_model(
    train_data = gray_zone_train,
    test_data  = gray_zone_train,   # placeholder; model selection uses TRAIN OOF
    outcome    = outcome,
    predictors = predictors_gray,
    method     = method_gray
  )
  th_star_youden <- fit_gray$th_star
  roc_oof_gray   <- fit_gray$roc_oof
  
  # Gray threshold variants on OOF ROC: Sp95 / Sp90 (nearest policy fallback)
  th_sp_vec <- vapply(gray_spec_targets, function(sp) {
    th <- suppressWarnings(pROC::coords(roc_oof_gray, x = sp,
                                        input = "specificity", ret = "threshold"))
    if (is.na(th)) th <- pick_threshold_by_spec(roc_oof_gray, target_spec = sp, policy = "nearest")
    as.numeric(th)
  }, numeric(1))
  names(th_sp_vec) <- paste0("sp", round(100*gray_spec_targets))
  
  th_se_vec <- NULL
  if (!is.null(gray_sens_targets) && length(gray_sens_targets) > 0L) {
    th_se_vec <- vapply(gray_sens_targets, function(se) {
      th <- suppressWarnings(pROC::coords(roc_oof_gray, x = se,
                                          input = "sensitivity", ret = "threshold"))
      if (is.na(th)) th <- pick_threshold_by_sens(roc_oof_gray, target_sens = se, policy = "nearest")
      as.numeric(th)
    }, numeric(1))
    names(th_se_vec) <- paste0("se", round(100*gray_sens_targets))
  }
  
  # --------------------- 4) apply to TEST (step1 & step2) -----------------
  test_step1 <- test_df %>%
    dplyr::mutate(
      step1 = dplyr::case_when(
        NfL_csf_pgml <= cutoff_low  ~ neg,
        NfL_csf_pgml >= cutoff_high ~ pos,
        TRUE                        ~ "gray"
      ),
      bucket = dplyr::case_when(
        NfL_csf_pgml <= cutoff_low  ~ "low",
        NfL_csf_pgml >= cutoff_high ~ "high",
        TRUE                        ~ "gray"
      ),
      bucket = factor(bucket, levels = c("low","gray","high"))
    )
  is_gray_test <- (test_step1$step1 == "gray")
  y_true <- factor(test_step1[[outcome]], levels = c(neg, pos))
  
  # Probabilities for gray; outside gray use hard 0/1 as continuous score
  prob_final <- rep(NA_real_, nrow(test_step1))
  prob_final[!is_gray_test & test_step1$step1 == neg] <- 0
  prob_final[!is_gray_test & test_step1$step1 == pos] <- 1
  if (any(is_gray_test)) {
    prob_gray <- fit_gray$predict_fun(test_step1[is_gray_test, , drop = FALSE])
    prob_final[is_gray_test] <- prob_gray
  } else {
    prob_gray <- numeric(0)
  }
  
  # Final labels per variant
  final_pred_youden <- test_step1$step1
  if (any(is_gray_test)) {
    final_pred_youden[is_gray_test] <- ifelse(prob_final[is_gray_test] >= th_star_youden, pos, neg)
  }
  
  final_pred_list <- list(youden = final_pred_youden)
  for (nm in names(th_sp_vec)) {
    fp <- test_step1$step1
    if (any(is_gray_test)) {
      fp[is_gray_test] <- ifelse(prob_final[is_gray_test] >= th_sp_vec[[nm]], pos, neg)
    }
    final_pred_list[[nm]] <- fp
  }
  
  if (!is.null(th_se_vec)) {
    for (nm in names(th_se_vec)) {
      fp <- test_step1$step1
      if (any(is_gray_test)) {
        fp[is_gray_test] <- ifelse(prob_final[is_gray_test] >= th_se_vec[[nm]], pos, neg)
      }
      final_pred_list[[nm]] <- fp
    }
  }
  
  # --------------------- 5) overall metrics + by-bucket -------------------
  roc_2step <- pROC::roc(response = y_true, predictor = prob_final,
                         levels = c(neg, pos), direction = "<")
  auc_2step <- as.numeric(pROC::auc(roc_2step))
  
  overall_from_pred <- function(pred) {
    TP <- sum(y_true == pos & pred == pos)
    TN <- sum(y_true == neg & pred == neg)
    FP <- sum(y_true == neg & pred == pos)
    FN <- sum(y_true == pos & pred == neg)
    list(
      TP=TP, TN=TN, FP=FP, FN=FN,
      Sens = TP / (TP + FN),
      Spec = TN / (TN + FP),
      PPV  = ifelse((TP+FP)>0, TP/(TP+FP), NA_real_),
      NPV  = ifelse((TN+FN)>0, TN/(TN+FN), NA_real_)
    )
  }
  overall_tbl <- tibble::tibble(
    variant = names(final_pred_list),
    AUC     = rep(auc_2step, length(final_pred_list)),
    Sens    = NA_real_, Spec = NA_real_, PPV = NA_real_, NPV = NA_real_
  )
  for (i in seq_along(final_pred_list)) {
    m <- overall_from_pred(factor(final_pred_list[[i]], levels = c(neg, pos)))
    overall_tbl$Sens[i] <- m$Sens; overall_tbl$Spec[i] <- m$Spec
    overall_tbl$PPV[i]  <- m$PPV;  overall_tbl$NPV[i]  <- m$NPV
  }
  
  n_pos_tot <- sum(y_true == pos); n_neg_tot <- sum(y_true == neg)
  p_high_given_pos <- ifelse(n_pos_tot>0, sum(y_true==pos & test_step1$bucket=="high")/n_pos_tot, NA_real_)
  p_high_given_neg <- ifelse(n_neg_tot>0, sum(y_true==neg & test_step1$bucket=="high")/n_neg_tot, NA_real_)
  LRp_high <- ifelse(is.finite(p_high_given_pos) & is.finite(p_high_given_neg) & p_high_given_neg>0,
                     p_high_given_pos / p_high_given_neg, NA_real_)
  p_notlow_given_pos <- ifelse(n_pos_tot>0, sum(y_true==pos & test_step1$bucket!="low")/n_pos_tot, NA_real_)
  p_low_given_neg    <- ifelse(n_neg_tot>0, sum(y_true==neg & test_step1$bucket=="low")/n_neg_tot, NA_real_)
  LRm_low <- ifelse(is.finite(p_notlow_given_pos) & is.finite(p_low_given_neg) & p_low_given_neg>0,
                    p_notlow_given_pos / p_low_given_neg, NA_real_)
  lr_plus_note  <- sprintf("LR+_high=%.3f (=%s/%s)", LRp_high, "P(high|ALS)","P(high|nonALS)")
  lr_minus_note <- sprintf("LR−_low=%.3f (=%s/%s)",  LRm_low,  "P(not low|ALS)","P(low|nonALS)")
  
  variants_long <- dplyr::bind_rows(lapply(names(final_pred_list), function(vn) {
    compute_by_bucket_tbl(
      y_true, test_step1$bucket,
      factor(final_pred_list[[vn]], levels = c(neg, pos)),
      note_lr_plus_high = lr_plus_note,
      note_lr_minus_low = lr_minus_note
    ) %>% dplyr::mutate(variant = vn, .before = bucket)
  }))
  variants_wide <- variants_long %>%
    tidyr::pivot_wider(
      id_cols = bucket,
      names_from = variant,
      values_from = c(n, PPV, NPV, Sens, Spec, TP, TN, FP, FN, LR_note),
      names_sep = "."
    )
  
  run_time <- Sys.time()
  ts2 <- list(
    meta = list(
      timestamp   = run_time,
      seed_split  = seed_split,
      seed_inner  = seed_inner,
      pos_class   = pos_class,
      neg_class   = neg,
      R_version   = R.version.string
    ),
    data = list(
      df_clean      = df_clean,
      train_idx     = train_idx,
      test_idx      = setdiff(seq_len(nrow(df_clean)), train_idx),
      train_df      = train_df,
      test_df       = test_df,
      predictors_all  = predictors_all,
      predictors_gray = predictors_gray
    ),
    step1 = list(
      cutoff_low     = cutoff_low,
      cutoff_high    = cutoff_high,
      auc_train_nfl  = auc_train_nfl,
      test_step1     = test_step1,
      gray_n_test    = sum(test_step1$bucket == "gray"),
      total_test     = nrow(test_step1),
      coverage1      = mean(test_step1$bucket != "gray"),
      sens1          = { 
        keep <- which(test_step1$bucket != "gray")
        yb   <- factor(test_step1[[outcome]][keep], levels = c(neg, pos))
        pb   <- factor(test_step1$step1[keep], levels = c(neg, pos))
        TP1 <- sum(yb==pos & pb==pos); FN1 <- sum(yb==pos & pb==neg)
        ifelse((TP1+FN1)>0, TP1/(TP1+FN1), NA_real_)
      },
      spec1          = {
        keep <- which(test_step1$bucket != "gray")
        yb   <- factor(test_step1[[outcome]][keep], levels = c(neg, pos))
        pb   <- factor(test_step1$step1[keep], levels = c(neg, pos))
        TN1 <- sum(yb==neg & pb==neg); FP1 <- sum(yb==neg & pb==pos)
        ifelse((TN1+FP1)>0, TN1/(TN1+FP1), NA_real_)
      }
    ),
    gray_model = list(
      method        = method_gray,
      fit           = fit_gray$model,
      th_star_gray  = th_star_youden,
      oof_auc       = if (!is.null(roc_oof_gray)) as.numeric(pROC::auc(roc_oof_gray)) else NA_real_,
      oof_df        = if (!is.null(roc_oof_gray)) {
        tibble::tibble(
          ALS_label = roc_oof_gray$response,
          prob      = roc_oof_gray$predictor
        )
      } else NULL,
      th_variants   = c(youden = th_star_youden, th_sp_vec, th_se_vec) 
    ),
    test_2step = list(
      prob_final   = prob_final,
      final_pred   = final_pred_list$youden,
      y_true       = y_true,
      auc_2step    = auc_2step,
      variants_overall = overall_tbl,
      roc_obj      = roc_2step
    ),
    baselines = list(
      nfl_only = NULL,
      step1_only = list(
        coverage = mean(test_step1$bucket != "gray"),
        sens     = { 
          keep <- which(test_step1$bucket != "gray")
          yb   <- factor(test_step1[[outcome]][keep], levels = c(neg, pos))
          pb   <- factor(test_step1$step1[keep], levels = c(neg, pos))
          TP1 <- sum(yb==pos & pb==pos); FN1 <- sum(yb==pos & pb==neg)
          ifelse((TP1+FN1)>0, TP1/(TP1+FN1), NA_real_)
        },
        spec     = {
          keep <- which(test_step1$bucket != "gray")
          yb   <- factor(test_step1[[outcome]][keep], levels = c(neg, pos))
          pb   <- factor(test_step1$step1[keep], levels = c(neg, pos))
          TN1 <- sum(yb==neg & pb==neg); FP1 <- sum(yb==neg & pb==pos)
          ifelse((TN1+FP1)>0, TN1/(TN1+FP1), NA_real_)
        }
      )
    ),
    plots = list(),
    by_bucket = list(
      variants      = variants_long,
      variants_wide = variants_wide
    )
  )
  class(ts2) <- c("twostep_fit", class(ts2))
  
  print.twostep_fit <- function(x, ...) {
    cat("Two-step diagnostic run\n")
    cat(" - Time:    ", as.character(x$meta$timestamp), "\n", sep = "")
    cat(" - Pos/Neg: ", x$meta$pos_class, " / ", x$meta$neg_class, "\n", sep = "")
    cat(" - TRAIN NfL AUC: ", sprintf("%.3f", x$step1$auc_train_nfl), "\n", sep = "")
    if (!is.null(x$test_2step$auc_2step)) {
      sens <- x$test_2step$variants_overall$Sens[x$test_2step$variants_overall$variant=="youden"] %||% NA_real_
      spec <- x$test_2step$variants_overall$Spec[x$test_2step$variants_overall$variant=="youden"] %||% NA_real_
      cat(" - TEST  2-step AUC/Sens/Spec: ",
          sprintf("%.3f / %.3f / %.3f", x$test_2step$auc_2step, sens, spec),
          "\n", sep = "")
    }
    cat(" - Gray% (TEST): ", if (!is.null(x$step1$total_test)) {
      sprintf("%.1f%%", 100 * x$step1$gray_n_test / x$step1$total_test)
    } else "NA", "\n", sep = "")
    invisible(x)
  }
  ts2
}

# ====================================================================== #
# Two-step diagnostic pipeline on an explicit TRAIN/TEST split
# ====================================================================== #

twostep_run_on_split <- function(train_df,
                                 test_df,
                                 outcome,
                                 predictors_all,
                                 predictors_gray,
                                 method_gray = "rf",
                                 seed_inner = 12345,
                                 target_spec_high = 0.95,
                                 target_sens_low  = 0.95,
                                 gray_spec_targets = c(0.95, 0.90),
                                 gray_sens_targets = NULL) {  # ← added
  `%||%` <- function(x, y) if (is.null(x)) y else x
  neg <- paste0("non", pos_class); pos <- pos_class
  
  ## Added: ensure character vectors (2026-01-08)
  predictors_all  <- as.character(unlist(predictors_all, use.names = FALSE))
  predictors_gray <- as.character(unlist(predictors_gray, use.names = FALSE))
  
  # Validate all columns exist and are numeric
  stopifnot(all(predictors_all %in% names(train_df)))
  stopifnot(all(predictors_all %in% names(test_df)))
  
  # ---- 0) ensure clean columns / factors (both sides) ------------------
  clean_side <- function(df) {
    df %>%
      dplyr::select(all_of(c(outcome, predictors_all))) %>%
      tidyr::drop_na() %>%
      dplyr::mutate(
        row_id    = dplyr::row_number(),
        !!outcome := factor(ifelse(.data[[outcome]] %in% c(1, pos), pos, neg),
                            levels = c(neg, pos))
      )
  }
  train_df <- clean_side(train_df)
  test_df  <- clean_side(test_df)
  
  # ---- 1) TRAIN ROC for NfL (direction safe) --------------------------
  roc_auto <- pROC::roc(
    response  = factor(train_df[[outcome]], levels = c(neg, pos)),
    predictor = train_df$NfL_csf_pgml,
    direction = "auto"
  )
  if (roc_auto$direction != "<") {
    roc_nfl <- pROC::roc(
      response  = factor(train_df[[outcome]], levels = c(neg, pos)),
      predictor = -train_df$NfL_csf_pgml,
      direction = "<"
    )
    nfl_sign <- -1
  } else {
    roc_nfl  <- roc_auto
    nfl_sign <-  1
  }
  roc_tbl <- tibble::tibble(
    threshold   = roc_nfl$thresholds,
    sensitivity = roc_nfl$sensitivities,
    specificity = roc_nfl$specificities
  ) %>% dplyr::filter(is.finite(threshold))
  
  pick_target <- function(tbl, target_col, target_val, tie_col) {
    hit <- tbl %>%
      dplyr::filter(.data[[target_col]] >= target_val) %>%
      dplyr::arrange(dplyr::desc(.data[[tie_col]])) %>%
      dplyr::slice_head(n = 1)
    if (nrow(hit) == 0) {
      hit <- tbl %>%
        dplyr::mutate(gap = abs(.data[[target_col]] - target_val)) %>%
        dplyr::arrange(gap, dplyr::desc(.data[[tie_col]])) %>%
        dplyr::slice_head(n = 1)
    }
    hit
  }
  pick_high   <- pick_target(roc_tbl, "specificity", target_spec_high, "sensitivity")
  pick_low    <- pick_target(roc_tbl, "sensitivity", target_sens_low,  "specificity")
  cutoff_high <- nfl_sign * pick_high$threshold
  cutoff_low  <- nfl_sign * pick_low$threshold
  auc_train_nfl <- as.numeric(pROC::auc(roc_nfl))
  
  # ---- 2) gray-only ML on TRAIN ---------------------------------------
  gray_zone_train <- train_df %>%
    dplyr::filter(NfL_csf_pgml > cutoff_low, NfL_csf_pgml < cutoff_high)
  stopifnot(length(unique(gray_zone_train[[outcome]])) == 2)
  
  set.seed(seed_inner)
  fit_gray <- fit_model(
    train_data = gray_zone_train,
    test_data  = gray_zone_train,   # OOF inside
    outcome    = outcome,
    predictors = predictors_gray,
    method     = method_gray
  )
  th_star_youden <- fit_gray$th_star
  roc_oof_gray   <- fit_gray$roc_oof
  
  pick_threshold_by_spec <- function(roc_obj, target_spec = 0.95,
                                     policy = c("at_least","nearest","at_most")) {
    policy <- match.arg(policy)
    tab <- pROC::coords(roc_obj, x="all",
                        ret=c("threshold","specificity","sensitivity"),
                        transpose = FALSE)
    tab <- as.data.frame(tab)
    suppressWarnings(tab <- tab[order(tab$threshold), ])
    if (policy == "at_least") {
      cand <- tab[tab$specificity >= target_spec, , drop = FALSE]
      if (nrow(cand) == 0L) return(NA_real_)
      return(cand$threshold[which.max(cand$sensitivity)])
    }
    if (policy == "at_most") {
      cand <- tab[tab$specificity <= target_spec, , drop = FALSE]
      if (nrow(cand) == 0L) return(NA_real_)
      return(cand$threshold[which.max(cand$specificity)])
    }
    idx <- which.min(abs(tab$specificity - target_spec))
    return(tab$threshold[idx])
  }
  pick_threshold_by_sens <- function(roc_obj, target_sens = 0.95,
                                     policy = c("at_least","nearest","at_most")) {
    policy <- match.arg(policy)
    tab <- pROC::coords(roc_obj, x="all",
                        ret=c("threshold","specificity","sensitivity"),
                        transpose = FALSE)
    tab <- as.data.frame(tab)
    suppressWarnings(tab <- tab[order(tab$threshold), ])
    if (policy == "at_least") {
      cand <- tab[tab$sensitivity >= target_sens, , drop = FALSE]
      if (nrow(cand) == 0L) return(NA_real_)
      return(cand$threshold[which.max(cand$specificity)])
    }
    if (policy == "at_most") {
      cand <- tab[tab$sensitivity <= target_sens, , drop = FALSE]
      if (nrow(cand) == 0L) return(NA_real_)
      return(cand$threshold[which.max(cand$sensitivity)])
    }
    idx <- which.min(abs(tab$sensitivity - target_sens))
    return(tab$threshold[idx])
  }
  th_sp_vec <- vapply(gray_spec_targets, function(sp) {
    th <- suppressWarnings(pROC::coords(roc_oof_gray, x = sp,
                                        input = "specificity", ret = "threshold"))
    if (is.na(th)) th <- pick_threshold_by_spec(roc_oof_gray, target_spec = sp, policy = "nearest")
    as.numeric(th)
  }, numeric(1))
  names(th_sp_vec) <- paste0("sp", round(100*gray_spec_targets))
  
  th_se_vec <- NULL
  if (!is.null(gray_sens_targets) && length(gray_sens_targets) > 0L) {
    th_se_vec <- vapply(gray_sens_targets, function(se) {
      th <- suppressWarnings(pROC::coords(roc_oof_gray, x = se,
                                          input = "sensitivity", ret = "threshold"))
      if (is.na(th)) th <- pick_threshold_by_sens(roc_oof_gray, target_sens = se, policy = "nearest")
      as.numeric(th)
    }, numeric(1))
    names(th_se_vec) <- paste0("se", round(100*gray_sens_targets))
  }
  
  # ---- 3) apply to TEST ------------------------------------------------
  test_step1 <- test_df %>%
    dplyr::mutate(
      step1 = dplyr::case_when(
        NfL_csf_pgml <= cutoff_low  ~ neg,
        NfL_csf_pgml >= cutoff_high ~ pos,
        TRUE                        ~ "gray"
      ),
      bucket = dplyr::case_when(
        NfL_csf_pgml <= cutoff_low  ~ "low",
        NfL_csf_pgml >= cutoff_high ~ "high",
        TRUE                        ~ "gray"
      ),
      bucket = factor(bucket, levels = c("low","gray","high"))
    )
  is_gray_test <- (test_step1$step1 == "gray")
  y_true <- factor(test_step1[[outcome]], levels = c(neg, pos))
  
  # continuous score for overall ROC
  prob_final <- rep(NA_real_, nrow(test_step1))
  prob_final[!is_gray_test & test_step1$step1 == neg] <- 0
  prob_final[!is_gray_test & test_step1$step1 == pos] <- 1
  prob_gray <- numeric(0)
  if (any(is_gray_test)) {
    need <- sum(is_gray_test)
    expected <- if (inherits(fit_gray$model, "train") && !is.null(fit_gray$model$xNames)) {
      fit_gray$model$xNames
    } else {
      predictors_gray
    }
    Xg <- as.data.frame(test_df[is_gray_test, expected, drop = FALSE])
    for (j in seq_along(Xg)) {
      if (!is.numeric(Xg[[j]])) Xg[[j]] <- suppressWarnings(as.numeric(Xg[[j]]))
    }
    prob_gray <- NULL
    if (inherits(fit_gray$model, "train")) {
      pr_df <- try(predict(fit_gray$model, newdata = Xg, type = "prob"), silent = TRUE)
      if (!inherits(pr_df, "try-error") && is.data.frame(pr_df) && pos %in% colnames(pr_df)) {
        v <- as.numeric(pr_df[[pos]])
        if (!all(is.na(v))) prob_gray <- v
      }
      if (is.null(prob_gray)) {
        raw_pred <- try(predict(fit_gray$model, newdata = Xg, type = "raw"), silent = TRUE)
        if (!inherits(raw_pred, "try-error")) {
          raw_pred <- as.character(raw_pred)
          prob_gray <- ifelse(raw_pred == pos, 0.75, 0.25)
        }
      }
      if ((is.null(prob_gray) || all(is.na(prob_gray))) &&
          inherits(fit_gray$model$finalModel, "ksvm")) {
        X_pp <- Xg
        if (!is.null(fit_gray$model$preProcess)) {
          X_pp <- try(predict(fit_gray$model$preProcess, X_pp), silent = TRUE)
          if (inherits(X_pp, "try-error")) X_pp <- Xg
        }
        pr <- try(kernlab::predict(fit_gray$model$finalModel,
                                   newdata = as.matrix(X_pp), type = "probabilities"),
                  silent = TRUE)
        if (!inherits(pr, "try-error")) {
          pr <- as.matrix(pr)
          pos_col <- which(colnames(pr) %in% pos)
          if (length(pos_col) == 1L) prob_gray <- as.numeric(pr[, pos_col])
        }
        if (is.null(prob_gray) || all(is.na(prob_gray))) {
          dv <- try(kernlab::predict(fit_gray$model$finalModel,
                                     newdata = as.matrix(X_pp), type = "decision"),
                    silent = TRUE)
          if (!inherits(dv, "try-error")) {
            v <- as.numeric(dv)
            r <- rank(v, ties.method = "average")
            prob_gray <- (r - min(r)) / max(1, (max(r) - min(r)))
          }
        }
      }
      
    } else if (inherits(fit_gray$model, "xgb.Booster")) {
      prob_gray <- suppressWarnings(fit_gray$predict_fun(Xg))
    } else {
      prob_gray <- suppressWarnings(fit_gray$predict_fun(Xg))
    }
    if (is.null(prob_gray) || length(prob_gray) == 0L || all(is.na(prob_gray))) {
      prob_gray <- rep(0.5, need)
    }
    if (length(prob_gray) != need) prob_gray <- rep(0.5, need)
    prob_gray[!is.finite(prob_gray)] <- 0.5
    prob_final[is_gray_test] <- as.numeric(prob_gray)
  }
  
  final_pred_youden <- test_step1$step1
  if (any(is_gray_test)) {
    final_pred_youden[is_gray_test] <- ifelse(prob_final[is_gray_test] >= th_star_youden, pos, neg)
  }
  final_pred_list <- list(youden = final_pred_youden)
  for (nm in names(th_sp_vec)) {
    fp <- test_step1$step1
    if (any(is_gray_test)) {
      fp[is_gray_test] <- ifelse(prob_final[is_gray_test] >= th_sp_vec[[nm]], pos, neg)
    }
    final_pred_list[[nm]] <- fp
  }
  
  if (!is.null(th_se_vec)) {
    for (nm in names(th_se_vec)) {
      fp <- test_step1$step1
      if (any(is_gray_test)) {
        fp[is_gray_test] <- ifelse(prob_final[is_gray_test] >= th_se_vec[[nm]], pos, neg)
      }
      final_pred_list[[nm]] <- fp
    }
  }
  
  roc_2step <- pROC::roc(response = y_true, predictor = prob_final,
                         levels = c(neg, pos), direction = "<")
  auc_2step <- as.numeric(pROC::auc(roc_2step))
  
  run_time <- Sys.time()
  ts2 <- list(
    meta = list(
      timestamp   = run_time,
      seed_inner  = seed_inner,
      pos_class   = pos_class,
      neg_class   = neg,
      R_version   = R.version.string
    ),
    data = list(
      train_df      = train_df,
      test_df       = test_df,
      predictors_all  = predictors_all,
      predictors_gray = predictors_gray
    ),
    step1 = list(
      cutoff_low     = cutoff_low,
      cutoff_high    = cutoff_high,
      auc_train_nfl  = auc_train_nfl,
      test_step1     = test_step1,
      gray_n_test    = sum(test_step1$bucket == "gray"),
      total_test     = nrow(test_step1),
      coverage1      = mean(test_step1$bucket != "gray")
    ),
    gray_model = list(
      method        = method_gray,
      fit           = fit_gray$model,
      th_star_gray  = th_star_youden,
      oof_auc       = if (!is.null(roc_oof_gray)) as.numeric(pROC::auc(roc_oof_gray)) else NA_real_,
      th_variants   = c(youden = th_star_youden, th_sp_vec, th_se_vec)  
    ),
    test_2step = list(
      prob_final   = prob_final,
      final_pred   = final_pred_list$youden,
      y_true       = y_true,
      auc_2step    = auc_2step,
      variants_pred= final_pred_list,   # keep per-variant final labels
      roc_obj      = roc_2step
    )
  )
  class(ts2) <- c("twostep_fit", class(ts2))
  ts2
}


# ====================================================================== #
# Two-step diagnostic pipeline with outer K-fold CV (OOF pooling)
# ====================================================================== #

twostep_kfold <- function(df,
                          outcome = "ALS_label",
                          predictors_all,
                          predictors_gray,
                          method_gray = "rf",
                          k = 5,
                          seed_folds = 12345,
                          seed_inner = 12345,
                          target_spec_high = 0.95,
                          target_sens_low  = 0.95,
                          gray_spec_targets = c(0.95, 0.90),
                          nfl_baseline_sens_target = NULL,  # e.g., 0.95 to also add Se95
                          gray_sens_targets = NULL) {      
  `%||%` <- function(x, y) if (is.null(x)) y else x
  neg <- paste0("non", pos_class); pos <- pos_class
  
  # ---- 0) Clean once on FULL data --------------------------------------
  df_clean <- df %>%
    dplyr::select(all_of(c(outcome, predictors_all))) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(
      row_id    = dplyr::row_number(),
      !!outcome := factor(ifelse(.data[[outcome]] %in% c(1, pos), pos, neg),
                          levels = c(neg, pos))
    )
  
  # ---- 1) Stratified K folds -------------------------------------------
  set.seed(seed_folds)
  folds <- caret::createFolds(df_clean[[outcome]], k = k, returnTrain = FALSE)
  
  # ---- 2) OOF containers (2-step) --------------------------------------
  n <- nrow(df_clean)
  oof_prob_final <- rep(NA_real_, n)
  oof_y_true     <- df_clean[[outcome]]
  oof_bucket_chr <- rep(NA_character_, n)
  
  # Variants for gray decision (Youden + Sp targets)
  variant_names  <- c("youden", paste0("sp", round(100*gray_spec_targets)))
  if (!is.null(gray_sens_targets) && length(gray_sens_targets) > 0L) {
    variant_names <- c(variant_names, paste0("se", round(100*gray_sens_targets)))  
  }
  oof_final_pred <- lapply(
    variant_names,
    function(.) factor(rep(NA_character_, n), levels = c(neg, pos))
  )
  names(oof_final_pred) <- variant_names
  
  # ---- 3) NFL-only baselines OOF containers ----------------------------
  nfl_pred_youden <- factor(rep(NA_character_, n), levels = c(neg, pos))
  nfl_pred_sp95   <- factor(rep(NA_character_, n), levels = c(neg, pos))
  nfl_pred_sens95 <- if (!is.null(nfl_baseline_sens_target))
    factor(rep(NA_character_, n), levels = c(neg, pos)) else NULL
  
  oof_nfl_score <- rep(NA_real_, n)
  th_youden_vec <- rep(NA_real_, length(folds))
  th_sp95_vec   <- rep(NA_real_, length(folds))
  th_se_vec     <- if (!is.null(nfl_baseline_sens_target)) rep(NA_real_, length(folds)) else NULL
  step1_cut_low_vec  <- rep(NA_real_, length(folds))
  step1_cut_high_vec <- rep(NA_real_, length(folds))
  per_fold_ts2 <- vector("list", length(folds))
  
  # ---- 4) Fold loop -----------------------------------------------------
  for (i in seq_along(folds)) {
    idx_test  <- folds[[i]]
    idx_train <- setdiff(seq_len(nrow(df_clean)), idx_test)
    tr <- df_clean[idx_train, ]
    te <- df_clean[idx_test,  ]
    
    ts2 <- twostep_run_on_split(
      train_df           = tr,
      test_df            = te,
      outcome            = outcome,
      predictors_all     = predictors_all,
      predictors_gray    = predictors_gray,
      method_gray        = method_gray,
      seed_inner         = seed_inner,
      target_spec_high   = target_spec_high,
      target_sens_low    = target_sens_low,
      gray_spec_targets  = gray_spec_targets,
      gray_sens_targets  = gray_sens_targets  
    )
    per_fold_ts2[[i]] <- ts2
    
    step1_cut_low_vec[i]  <- ts2$step1$cutoff_low
    step1_cut_high_vec[i] <- ts2$step1$cutoff_high
    
    oof_prob_final[idx_test] <- ts2$test_2step$prob_final
    oof_bucket_chr[idx_test] <- as.character(ts2$step1$test_step1$bucket)
    
    for (vn in names(ts2$test_2step$variants_pred)) {
      oof_final_pred[[vn]][idx_test] <-
        factor(ts2$test_2step$variants_pred[[vn]], levels = c(neg, pos))
    }
    
    oof_nfl_score[idx_test] <- te$NfL_csf_pgml
    
    roc_nfl_tr <- pROC::roc(tr[[outcome]], tr$NfL_csf_pgml,
                            levels = c(neg, pos), direction = "<")
    
    th_youden <- as.numeric(pROC::coords(roc_nfl_tr, x = "best",
                                         best.method = "youden",
                                         ret = "threshold", transpose = FALSE))
    th_youden_vec[i] <- th_youden
    nfl_pred_youden[idx_test] <- factor(
      ifelse(te$NfL_csf_pgml >= th_youden, pos, neg), levels = c(neg, pos)
    )
    
    th_sp95 <- suppressWarnings(pROC::coords(
      roc_nfl_tr, x = 0.95, input = "specificity", ret = "threshold"
    ))
    if (is.na(th_sp95)) {
      tbl <- pROC::coords(
        roc_nfl_tr, x = "all",
        ret = c("threshold","specificity","sensitivity"),
        transpose = FALSE
      ) |> as.data.frame()
      cand <- tbl[tbl$specificity >= 0.95, , drop = FALSE]
      if (nrow(cand) == 0L) {
        cand <- tbl[order(abs(tbl$specificity - 0.95), -tbl$sensitivity), ][1, , drop = FALSE]
      } else {
        cand <- cand[order(-cand$sensitivity), ][1, , drop = FALSE]
      }
      th_sp95 <- cand$threshold[1]
    }
    th_sp95_vec[i] <- th_sp95
    nfl_pred_sp95[idx_test] <- factor(
      ifelse(te$NfL_csf_pgml >= th_sp95, pos, neg), levels = c(neg, pos)
    )
    
    if (!is.null(nfl_baseline_sens_target)) {
      th_se <- suppressWarnings(pROC::coords(
        roc_nfl_tr, x = nfl_baseline_sens_target,
        input = "sensitivity", ret = "threshold"
      ))
      if (is.na(th_se)) {
        tbl <- pROC::coords(
          roc_nfl_tr, x = "all",
          ret = c("threshold","specificity","sensitivity"),
          transpose = FALSE
        ) |> as.data.frame()
        cand <- tbl[tbl$sensitivity >= nfl_baseline_sens_target, , drop = FALSE]
        if (nrow(cand) == 0L) {
          cand <- tbl[order(abs(tbl$sensitivity - nfl_baseline_sens_target), -tbl$specificity), ][1, , drop = FALSE]
        } else {
          cand <- cand[order(-cand$specificity), ][1, , drop = FALSE]
        }
        th_se <- cand$threshold[1]
      }
      th_se_vec[i] <- th_se
      nfl_pred_sens95[idx_test] <- factor(
        ifelse(te$NfL_csf_pgml >= th_se, pos, neg), levels = c(neg, pos)
      )
    }
  } # end folds
  
  # ---- 5) Pooled OOF metrics -------------------------------------------
  roc_oof <- pROC::roc(response = oof_y_true, predictor = oof_prob_final,
                       levels = c(neg, pos), direction = "<")
  auc_oof <- as.numeric(pROC::auc(roc_oof))
  
  bucket_factor <- factor(oof_bucket_chr, levels = c("low","gray","high"))
  gray_idx <- which(bucket_factor == "gray")
  if (length(gray_idx) >= 5L) {
    roc_gray <- pROC::roc(response = oof_y_true[gray_idx],
                          predictor = oof_prob_final[gray_idx],
                          levels = c(neg, pos), direction = "<")
    auc_gray <- as.numeric(pROC::auc(roc_gray))
  } else {
    roc_gray <- NULL
    auc_gray <- NA_real_
  }
  
  overall_from_pred <- function(pred) {
    TP <- sum(oof_y_true == pos & pred == pos, na.rm = TRUE)
    TN <- sum(oof_y_true == neg & pred == neg, na.rm = TRUE)
    FP <- sum(oof_y_true == neg & pred == pos, na.rm = TRUE)
    FN <- sum(oof_y_true == pos & pred == neg, na.rm = TRUE)
    list(TP=TP, TN=TN, FP=FP, FN=FN,
         Sens = TP/(TP+FN),
         Spec = TN/(TN+FP),
         PPV  = ifelse((TP+FP)>0, TP/(TP+FP), NA_real_),
         NPV  = ifelse((TN+FN)>0, TN/(TN+FN), NA_real_))
  }
  
  variants_overall <- tibble::tibble(
    variant = names(oof_final_pred),
    AUC     = rep(auc_oof, length(oof_final_pred)),
    AUC_gray= rep(auc_gray, length(oof_final_pred)),
    Sens    = NA_real_, Spec = NA_real_, PPV = NA_real_, NPV = NA_real_
  )
  for (i in seq_along(oof_final_pred)) {
    m <- overall_from_pred(oof_final_pred[[i]])
    variants_overall$Sens[i] <- m$Sens
    variants_overall$Spec[i] <- m$Spec
    variants_overall$PPV[i]  <- m$PPV
    variants_overall$NPV[i]  <- m$NPV
  }
  
  bucket_counts <- as.data.frame(table(bucket_factor), stringsAsFactors = FALSE)
  names(bucket_counts) <- c("bucket","n")
  bucket_counts$prop <- bucket_counts$n / sum(bucket_counts$n)
  
  gray_tables <- lapply(names(oof_final_pred), function(vn) {
    pred  <- oof_final_pred[[vn]][gray_idx]
    truth <- oof_y_true[gray_idx]
    TP <- sum(truth == pos & pred == pos, na.rm = TRUE)
    TN <- sum(truth == neg & pred == neg, na.rm = TRUE)
    FP <- sum(truth == neg & pred == pos, na.rm = TRUE)
    FN <- sum(truth == pos & pred == neg, na.rm = TRUE)
    tibble::tibble(
      variant = vn,
      gray_n  = length(gray_idx),
      TP = TP, TN = TN, FP = FP, FN = FN,
      Sens_gray = ifelse((TP+FN)>0, TP/(TP+FN), NA_real_),
      Spec_gray = ifelse((TN+FP)>0, TN/(TN+FP), NA_real_),
      PPV_gray  = ifelse((TP+FP)>0, TP/(TP+FP), NA_real_),
      NPV_gray  = ifelse((TN+FN)>0, TN/(TN+FN), NA_real_)
    )
  })
  gray_overall <- dplyr::bind_rows(gray_tables)
  
  # ---- 6) NFL-only OOF ROC + baseline table -----------------------------
  roc_oof_nfl <- pROC::roc(oof_y_true, oof_nfl_score,
                           levels = c(neg, pos), direction = "<")
  auc_oof_nfl <- as.numeric(pROC::auc(roc_oof_nfl))
  
  make_row <- function(name, pred, th_vec) {
    m <- overall_from_pred(pred)
    tibble::tibble(
      baseline = name,
      Sens = m$Sens, Spec = m$Spec, PPV = m$PPV, NPV = m$NPV,
      AUC_oof = auc_oof_nfl,
      # th_train_median = stats::median(th_vec, na.rm = TRUE),
      # th_train_mean   = mean(th_vec, na.rm = TRUE),
      # th_train_sd     = stats::sd(th_vec, na.rm = TRUE),
      # th_train_all    = list(th_vec)
      th_train_median = stats::median(as.numeric(th_vec), na.rm = TRUE),  
      th_train_mean   = mean(as.numeric(th_vec), na.rm = TRUE),           
      th_train_sd     = stats::sd(as.numeric(th_vec), na.rm = TRUE),      
      th_train_all    = list(as.numeric(th_vec))                          
    )
  }
  baselines_tbl <- dplyr::bind_rows(
    make_row("nfl_youden", nfl_pred_youden, th_youden_vec),
    make_row("nfl_sp95",   nfl_pred_sp95,   th_sp95_vec),
    if (!is.null(nfl_pred_sens95)) make_row(
      sprintf("nfl_se%d", round(100*nfl_baseline_sens_target)),
      nfl_pred_sens95, th_se_vec
    )
  )
  
  out <- list(
    meta = list(
      k           = k,
      seed_folds  = seed_folds,
      seed_inner  = seed_inner,
      pos_class   = pos_class,
      neg_class   = neg
    ),
    data = list(
      df_clean = df_clean,
      folds    = folds
    ),
    per_fold = per_fold_ts2,
    pooled = list(
      oof = list(
        y_true       = oof_y_true,
        prob_final   = oof_prob_final,
        final_pred   = oof_final_pred,
        bucket       = bucket_factor,
        roc_oof      = roc_oof,
        auc_oof      = auc_oof,
        roc_gray     = roc_gray,
        auc_gray     = auc_gray
      ),
      bucket_counts    = bucket_counts,
      gray_overall     = gray_overall,
      variants_overall = variants_overall,
      nfl_baselines    = baselines_tbl,
      nfl_roc          = list(roc_oof = roc_oof_nfl, auc_oof = auc_oof_nfl),
      step1_thresholds = list(
        low  = step1_cut_low_vec,
        high = step1_cut_high_vec
      ),
      nfl_thresholds = list(
        youden = th_youden_vec,
        sp95   = th_sp95_vec,
        se95   = th_se_vec
      )
    )
  )
  class(out) <- c("twostep_kfold", class(out))
  out
}

# ====================================================================== #
# Multi-method outer-CV runner (wraps twostep_kfold for each method)
# ====================================================================== #

twostep_kfold_multi <- function(df,
                                outcome = "ALS_label",
                                predictors_all,
                                predictors_gray,
                                methods_gray = c("rf","xgb","glm","svmRadial","kknn"),
                                k = 5,
                                seed_folds = 12345,
                                seed_inner = 12345,
                                target_spec_high = 0.95,
                                target_sens_low  = 0.95,
                                gray_spec_targets = c(0.95, 0.90),
                                nfl_baseline_sens_target = 0.95,
                                gray_sens_targets = NULL) {  
  # Run the outer CV for each gray method
  results <- vector("list", length(methods_gray))
  names(results) <- methods_gray
  
  for (m in methods_gray) {
    message(sprintf("[twostep_kfold_multi] Running method_gray='%s' ...", m))
    results[[m]] <- twostep_kfold(
      df               = df,
      outcome          = outcome,
      predictors_all   = predictors_all,
      predictors_gray  = predictors_gray,
      method_gray      = m,
      k                = k,
      seed_folds       = seed_folds,
      seed_inner       = seed_inner,
      target_spec_high = target_spec_high,
      target_sens_low  = target_sens_low,
      gray_spec_targets= gray_spec_targets,
      nfl_baseline_sens_target = nfl_baseline_sens_target,
      gray_sens_targets = gray_sens_targets   
    )
  }
  
  # ---- Build comparison tables (as before) ------------------------------
  overall_rows <- lapply(names(results), function(m) {
    ro <- results[[m]]$pooled$variants_overall %>%
      dplyr::mutate(method = m)
    
    # gray side  PPV_gray / NPV_gray 
    rg <- results[[m]]$pooled$gray_overall %>%
      dplyr::transmute(
        variant, gray_n,
        Sens_gray, Spec_gray, PPV_gray, NPV_gray
      )
    
    ro %>%
      dplyr::left_join(rg, by = "variant") %>%
      dplyr::select(
        method, variant,
        AUC, AUC_gray,
        Sens, Spec, PPV, NPV,
        PPV_gray, NPV_gray, Sens_gray, Spec_gray, gray_n
      )
  })
  
  overall_cmp <- dplyr::bind_rows(overall_rows) %>%
    dplyr::arrange(variant, dplyr::desc(AUC), dplyr::desc(Sens))
  
  gray_rows <- lapply(names(results), function(m) {
    rg <- results[[m]]$pooled$gray_overall
    rg$method <- m
    rg[, c("method","variant","gray_n","TP","TN","FP","FN","Sens_gray","Spec_gray","PPV_gray","NPV_gray")]
  })
  gray_cmp <- dplyr::bind_rows(gray_rows) %>%
    dplyr::arrange(variant, dplyr::desc(Sens_gray), dplyr::desc(PPV_gray))
  
  any_method <- methods_gray[[1]]
  bucket_counts <- results[[any_method]]$pooled$bucket_counts
  nfl_baselines <- results[[any_method]]$pooled$nfl_baselines
  nfl_roc       <- results[[any_method]]$pooled$nfl_roc
  
  # ---- Step1 thresholds (low/high) per method/fold (tidy) ----------
  step1_th_tbl <- lapply(names(results), function(m) {
    sts <- results[[m]]$pooled$step1_thresholds
    tibble::tibble(
      method = m,
      fold   = seq_along(sts$low),
      cutoff_low  = as.numeric(sts$low),
      cutoff_high = as.numeric(sts$high)
    )
  }) |> dplyr::bind_rows()
  
  # ---- NFL baseline TRAIN thresholds per method/fold (tidy) --------
  nfl_th_tbl <- lapply(names(results), function(m) {
    nth <- results[[m]]$pooled$nfl_thresholds
    tibble::tibble(
      method = m,
      fold   = seq_along(nth$youden),
      youden = as.numeric(nth$youden),
      sp95   = as.numeric(nth$sp95),
      se95   = if (!is.null(nth$se95)) as.numeric(nth$se95) else NA_real_
    )
  }) |> dplyr::bind_rows()
  
  out <- list(
    meta = list(
      methods_gray = methods_gray,
      k            = k,
      seed_folds   = seed_folds,
      seed_inner   = seed_inner,
      gray_spec_targets = gray_spec_targets,
      outcome      = outcome
    ),
    per_method = results,           # each is a full twostep_kfold object
    compare = list(
      overall   = overall_cmp,
      gray_only = gray_cmp,
      buckets   = bucket_counts,
      nfl       = list(baselines = nfl_baselines, roc = nfl_roc),
      # NEW: handy tables for variability plots
      step1_thresholds = step1_th_tbl,   # → violin/box/jitter across folds
      nfl_thresholds   = nfl_th_tbl      # → also available if needed
    )
  )
  class(out) <- c("twostep_kfold_multi", class(out))
  out
}


# ====================================================================== #
# Repeated outer-CV runner that GUARANTEES n_repeats successful runs
# ====================================================================== #
twostep_kfold_multi_repeat <- function(
    df,
    outcome = "ALS_label",
    predictors_all,
    predictors_gray,
    methods_gray = c("rf","xgb","glm","svmRadial","kknn"),
    k = 5,
    n_repeats = 10,
    seed_base = 20250913,
    seed_inner = 12345,
    target_spec_high = 0.95,
    target_sens_low  = 0.95,
    gray_spec_targets = c(0.95), #c(0.95,0.90)
    gray_sens_targets = NULL,
    nfl_baseline_sens_target = 0.95,
    max_attempts = 100
) {
  # --- small validator: check a kfold_multi object looks usable ---
  is_valid_km <- function(km) {
    if (is.null(km) || !is.list(km)) return(FALSE)
    req <- c("per_method","compare")
    if (!all(req %in% names(km))) return(FALSE)
    ok_methods <- length(km$per_method) > 0
    ok_overall <- is.data.frame(km$compare$overall) && nrow(km$compare$overall) > 0
    ok_gray    <- is.data.frame(km$compare$gray_only) && nrow(km$compare$gray_only) > 0
    ok_methods && ok_overall && ok_gray
  }
  
  runs <- list()
  used_seeds <- integer(0)
  
  attempt_log <- tibble::tibble(
    attempt = integer(),
    seed_folds = integer(),
    status = character(),
    error = character()
  )
  
  successes <- 0L
  attempts  <- 0L
  seed_cur  <- seed_base
  
  while (successes < n_repeats) {
    attempts <- attempts + 1L
    if (attempts > max_attempts) {
      warning(sprintf("Stopped after %d attempts; collected %d/%d successes.",
                      max_attempts, successes, n_repeats))
      break
    }
    
    msg <- sprintf("[attempt %d] running twostep_kfold_multi (seed_folds=%d)...",
                   attempts, seed_cur)
    message(msg)
    
    km_try <- tryCatch({
      twostep_kfold_multi(
        df               = df,
        outcome          = outcome,
        predictors_all   = predictors_all,
        predictors_gray  = predictors_gray,
        methods_gray     = methods_gray,
        k                = k,
        seed_folds       = seed_cur,
        seed_inner       = seed_inner,
        target_spec_high = target_spec_high,
        target_sens_low  = target_sens_low,
        gray_spec_targets= gray_spec_targets,
        gray_sens_targets= gray_sens_targets, #add
        nfl_baseline_sens_target = nfl_baseline_sens_target
      )
    }, error = function(e) {
      structure(list(.error = e), class = "try-error")
    })
    
    if (inherits(km_try, "try-error") || !is_valid_km(km_try)) {
      err_msg <- if (inherits(km_try, "try-error")) km_try$.error$message else "invalid object"
      attempt_log <- dplyr::bind_rows(
        attempt_log,
        tibble::tibble(
          attempt   = attempts,
          seed_folds= seed_cur,
          status    = "fail",
          error     = err_msg
        )
      )
      # advance seed for next attempt
      seed_cur <- seed_cur + 1L
      next
    }
    
    # success: store run
    successes <- successes + 1L
    runs[[successes]] <- km_try
    used_seeds <- c(used_seeds, seed_cur)
    attempt_log <- dplyr::bind_rows(
      attempt_log,
      tibble::tibble(
        attempt   = attempts,
        seed_folds= seed_cur,
        status    = "ok",
        error     = ""
      )
    )
    # advance seed for next attempt as well
    seed_cur <- seed_cur + 1L
  }
  
  if (length(runs) == 0L) {
    stop("No successful repeats were collected. Check data and hyperparameters.")
  }
  
  # --- build stacked tables for convenience (with repeat_id) ---
  stack_overall <- dplyr::bind_rows(lapply(seq_along(runs), function(i) {
    runs[[i]]$compare$overall %>%
      dplyr::mutate(repeat_id = i, .before = 1L)
  }))
  stack_gray <- dplyr::bind_rows(lapply(seq_along(runs), function(i) {
    runs[[i]]$compare$gray_only %>%
      dplyr::mutate(repeat_id = i, .before = 1L)
  }))
  stack_buckets <- dplyr::bind_rows(lapply(seq_along(runs), function(i) {
    runs[[i]]$compare$buckets %>%
      dplyr::mutate(repeat_id = i, .before = 1L)
  }))
  # NfL thresholds distribution per repeat (from the first method in that repeat)
  stack_nfl_th <- dplyr::bind_rows(lapply(seq_along(runs), function(i) {
    any_method <- runs[[i]]$meta$methods_gray[[1]]
    runs[[i]]$per_method[[any_method]]$pooled$nfl_baselines %>%
      dplyr::mutate(repeat_id = i, .before = 1L)
  }))
  
  out <- list(
    runs = runs,
    used_seeds = used_seeds,
    attempt_log = attempt_log,
    stacked = list(
      overall   = stack_overall,
      gray_only = stack_gray,
      buckets   = stack_buckets,
      nfl_baselines = stack_nfl_th
    )
  )
  class(out) <- c("twostep_kfold_multi_repeat", class(out))
  out
}