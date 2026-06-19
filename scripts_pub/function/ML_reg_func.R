library(dplyr)
library(caret)
library(pROC)
library(xgboost)
library(fastshap)
library(ggplot2)
library(ggbeeswarm)  
library(ggforce)     
library(tidyr)
library(viridisLite)


# Function ----------------------------------------------------------------

common_ctrl_reg <- trainControl(
  method = "cv",
  number = 5,             # 5 means 5-fold cross-validation (i.e., each test is ~20% of the data)
  savePredictions = "final",
  allowParallel = FALSE,  # FALSE is default for small data # (to avoid parallel overhead); set TRUE for larger data sets
)


fmt_mean_sd <- function(m, s, digits = 3, na_char = "-") {
  ifelse(is.na(m) | is.na(s),
         na_char,
         sprintf(paste0("%.", digits, "f ± %.", digits, "f"), m, s))
}


#── Utility: Automatically determine whether scaling is needed for each model ─────────────
needs_scaling <- function(method) {
  # Distance-based or kernel-based methods require scaling
  method %in% c("knn", "svmRadial")
}



# Prepare design matrices with dummy encoding (train-fitted, leakage-safe)
prepare_design_matrix <- function(train_data, test_data, predictors) {
  # 1) Select predictor frame
  Xtr_raw <- train_data[, predictors, drop = FALSE]
  Xte_raw <- test_data[,  predictors, drop = FALSE]
  
  # 2) Learn dummy mapping on TRAIN only
  dmy <- caret::dummyVars(~ ., data = Xtr_raw, fullRank = TRUE)
  
  # 3) Apply mapping to train/test
  Xtr <- as.data.frame(predict(dmy, newdata = Xtr_raw))
  Xte <- as.data.frame(predict(dmy, newdata = Xte_raw))
  
  # 4) Ensure numeric matrix-like frames
  stopifnot(all(sapply(Xtr, is.numeric)), all(sapply(Xte, is.numeric)))
  
  list(Xtr = Xtr, Xte = Xte, dummy = dmy, features = colnames(Xtr))
}


# Regression fitter with safe dummy encoding for categorical predictors

fit_model_reg <- function(train_data, test_data, outcome, predictors,
                          method, seed = NULL, scale = c("auto","yes","no"),
                          tr_control = common_ctrl_reg) {
  if (!is.null(seed)) set.seed(seed)
  
  scale <- match.arg(scale)
  needs_scaling <- function(m) m %in% c("knn","svmRadial","glm")
  
  # ---------- build design matrices ----------
  if (is.character(predictors)) {
    dm <- prepare_design_matrix(train_data, test_data, predictors)
    Xtr <- dm$Xtr; Xte <- dm$Xte
    feat_names <- dm$features
    ytr <- train_data[[outcome]]
    yte <- test_data[[outcome]]
    
    sp <- NULL
    if (scale == "auto") scale <- if (needs_scaling(method)) "yes" else "no"
    if (scale == "yes") {
      m <- sapply(Xtr, mean, na.rm = TRUE)
      s <- sapply(Xtr, sd,   na.rm = TRUE); s[s == 0 | is.na(s)] <- 1
      Xtr <- as.data.frame(sweep(sweep(Xtr, 2, m, "-"), 2, s, "/"))
      Xte <- as.data.frame(sweep(sweep(Xte, 2, m, "-"), 2, s, "/"))
      sp <- list(mean = m, sd = s)
    }
    Xtr_for_shap <- Xtr
    fmla   <- as.formula(paste(outcome, "~", paste(feat_names, collapse = " + ")))
    train_x <- cbind(setNames(data.frame(ytr), outcome), Xtr)
    test_x  <- cbind(setNames(data.frame(yte), outcome), Xte)
    
  } else if (inherits(predictors, "formula")) {
    fmla    <- predictors
    train_x <- train_data
    test_x  <- test_data
    feat_names <- all.vars(fmla)[-1]
    sp <- NULL
    
    # <<< ensure scaling even in formula branch for SVM/kNN
    preproc_vec <- NULL
    if (scale == "auto") scale <- if (needs_scaling(method)) "yes" else "no"
    if (scale == "yes") preproc_vec <- c("center","scale")
  } else {
    stop("predictors must be either a character vector or a formula.")
  }
  
  # ---------- XGBoost (custom small-data tuning) ----------
  ### XGB -------------------------------------------------------
  if (method == "xgb" && is.character(predictors)) {
    if (!is.null(seed)) set.seed(seed)
    
    dtrain <- xgboost::xgb.DMatrix(as.matrix(Xtr), label = ytr)
    
    # ---- Random search configuration ----
    n_trials     <- 24        # number of random hyperparameter sets to try
    max_time_sec <- 120       # time limit (seconds) for search
    nrounds_cap  <- 2000      # upper bound for boosting rounds
    esr          <- 50        # early stopping rounds
    nfold_inner  <- if (!is.null(tr_control$number)) tr_control$number else 5
    
    t0 <- proc.time()[[3]]
    
    # Helper functions for sampling distributions
    r_unif      <- function(a, b) runif(1, a, b)
    r_logunif   <- function(a, b) exp(runif(1, log(a), log(b)))
    r_discrete  <- function(v)   sample(v, 1)
    
    best <- NULL; best_rmse <- Inf; best_nrounds <- 100
    
    for (t in seq_len(n_trials)) {
      # ---- Sample one candidate hyperparameter set ----
      params <- list(
        objective        = "reg:squarederror",
        eval_metric      = "rmse",
        nthread          = 1,
        max_depth        = r_discrete(2:5),
        eta              = r_logunif(0.02, 0.3),       # learning rate (log-uniform)
        min_child_weight = r_discrete(c(1, 2, 3, 5, 7)),
        subsample        = r_unif(0.6, 0.95),
        colsample_bytree = r_unif(0.6, 0.95),
        lambda           = r_logunif(1e-4, 10),        # L2 regularization
        alpha            = r_logunif(1e-4, 1)          # L1 regularization
      )
      
      # Cross-validation with early stopping
      cv <- xgboost::xgb.cv(
        params = params, data = dtrain,
        nrounds = nrounds_cap, nfold = nfold_inner,
        early_stopping_rounds = esr, verbose = 0
      )
      
      rmse_mean <- min(cv$evaluation_log$test_rmse_mean)
      if (rmse_mean < best_rmse) {
        best_rmse   <- rmse_mean
        best        <- params
        best_nrounds <- cv$best_iteration
      }
      
      # ---- Time-based stopping condition ----
      if ((proc.time()[[3]] - t0) > max_time_sec) break
    }
    
    # Fallback if no configuration was successful
    if (is.null(best)) {
      best <- list(objective="reg:squarederror", eval_metric="rmse",
                   nthread=1, max_depth=3, eta=0.1, min_child_weight=3,
                   subsample=0.8, colsample_bytree=0.8, lambda=1, alpha=0)
      best_nrounds <- 200
    }
    
    # Final model training with the best hyperparameters
    mdl  <- xgboost::xgb.train(params = best, data = dtrain,
                               nrounds = best_nrounds, verbose = 0)
    pred <- as.numeric(predict(mdl, xgboost::xgb.DMatrix(as.matrix(Xte))))
    
    # Variable importance (optional)
    VIP <- tryCatch({
      imp <- xgboost::xgb.importance(feature_names = feat_names, model = mdl)
      if (!NROW(imp)) return(NULL)
      df <- data.frame(Feature = imp$Feature, Overall = imp$Gain, row.names = imp$Feature)
      df <- df[, "Overall", drop = FALSE]
      structure(list(importance = df, model = "xgb", calledFrom = "varImp"),
                class = "varImp.train")
    }, error = function(e) NULL)
    
    return(list(
      model = mdl,
      pred  = pred,
      vip   = VIP,
      scaled = (scale == "yes"),
      scale_mean = if (scale == "yes") sp$mean else NULL,
      scale_sd   = if (scale == "yes") sp$sd   else NULL,
      dummy = dm$dummy,
      features = feat_names,
      train_X_for_shap = Xtr_for_shap
    ))
  }
  
  
  # ---------- caret branch (glm, rf, knn, svmRadial, ...) ----------
  ctrl <- if (is.null(tr_control)) common_ctrl_reg else tr_control
  
  # ----- method-specific tuneGrid / preProcess (conservative) -----
  tune_grid <- NULL
  pre_proc  <- NULL
  p <- length(feat_names)
  n_tr <- nrow(train_x)
  
  if (method == "rf") {
    # ranger: robust defaults for small n
    method <- "ranger"
    pre_proc <- "zv"  # trees do not need scaling
    mtry_vals <- unique(pmax(1, round(c(sqrt(p)*0.5, sqrt(p), sqrt(p)*1.5))))
    tune_grid <- expand.grid(
      mtry = mtry_vals,
      splitrule = c("variance"),
      min.node.size = c(1, 5, 10)
    )
    
    # Save training design matrix for SHAP (RF doesn't scale, so just save Xtr)
    # if (is.character(predictors)) {
    #   Xtr_for_shap <- Xtr
    # }
    if (!exists("Xtr_for_shap")) {
      Xtr_for_shap <- Xtr  # RF doesn't scale, so use Xtr directly
    }
    
  } else if (method == "knn") {
    # Use kknn backend with conservative settings (no PCA).
    method <- "kknn"
    
    # Preprocess: handle skew/outliers, then center/scale, drop zero-variance.
    # No PCA to avoid breaking neighborhood geometry on small n.
    #pre_proc <- c("YeoJohnson", "center", "scale", "zv")
    # (manual scaling already applied to Xtr upstream)
    pre_proc <- c("YeoJohnson", "zv")
    
    # Set a sensible upper bound for k: at most ~30% of training size, odd only.
    n_tr <- nrow(train_x)
    k_max_cap <- max(7, min(45, floor(0.3 * n_tr)))
    k_grid <- seq(5, k_max_cap, by = 2)
    
    # Distance: try Manhattan (1) and Euclidean (2).
    # Kernel: rectangular (unweighted) is most stable on small n.
    tune_grid <- expand.grid(
      kmax     = k_grid,
      distance = c(1, 2),
      kernel   = c("rectangular")
    )
    
  } else if (method == "svmRadial") {
    # build a sigma grid around sigest on the TRAIN predictors
    #pre_proc <- c("center","scale","zv")
    pre_proc <- c("zv")
    # extract numeric matrix after potential preProcess (caret will re-do internally;
    # here we only need a quick estimate for sigma scale)
    X_for_sigma <- model.matrix(fmla, data = train_x)[, -1, drop = FALSE]
    # safe guard against constant columns
    csd <- apply(X_for_sigma, 2, sd); csd[csd == 0 | is.na(csd)] <- 1
    X_for_sigma <- sweep(X_for_sigma, 2, csd, "/")
    
    sig <- tryCatch(kernlab::sigest(X_for_sigma, frac = 1)[1], error = function(e) NA_real_)
    if (!is.finite(sig) || sig <= 0) sig <- 0.1
    sigma_grid <- sig * c(0.25, 0.5, 1, 2, 4)
    
    Cs <- 2 ^ seq(-4, 8, by = 1)
    tune_grid <- expand.grid(C = Cs, sigma = sigma_grid)
    
  } else if (method == "svmLinear") {
    pre_proc <- c("center","scale","zv")
    Cs <- 2 ^ seq(-6, 8, by = 1)
    tune_grid <- expand.grid(C = Cs)
    
  } else if (method == "glm") {
    pre_proc <- NULL
    tune_grid <- NULL
  }
  
  mdl <- caret::train(
    fmla, data = train_x,
    method = method,
    metric = "RMSE",
    trControl = ctrl,
    tuneGrid  = tune_grid,
    preProcess = pre_proc
  )
  
  
  
  pred <- as.numeric(predict(mdl, newdata = test_x))
  VIP  <- tryCatch(caret::varImp(mdl), error = function(e) NULL)
  
  # FIXED: Save train_X_for_shap for SHAP calculations
  # For caret branch, extract it from the fitted data
  if (exists("Xtr_for_shap")) {
    # Already defined in is.character(predictors) branch
    train_X_for_shap_final <- Xtr_for_shap
  } else {
    # Extract from train_x (remove outcome column)
    train_X_for_shap_final <- train_x[, feat_names, drop = FALSE]
  }
  
  list(
    model = mdl,
    pred  = pred,
    vip   = VIP,
    #scaled = (scale == "yes" || (!is.null(pre_proc))),
    scaled = (scale == "yes"),
    scale_mean = if (exists("sp") && !is.null(sp)) sp$mean else NULL,
    scale_sd   = if (exists("sp") && !is.null(sp)) sp$sd   else NULL,
    #dummy = if (exists("dm")) dm$dummy else NULL,
    dummy = if (exists("dm") && !is.null(dm)) dm$dummy else NULL,
    features = feat_names,
    #train_X_for_shap = if (exists("Xtr_for_shap")) Xtr_for_shap else NULL
    #train_X_for_shap = if (exists("Xtr_for_shap") && !is.null(Xtr_for_shap)) Xtr_for_shap else Xtr  # Fallback to Xtr
    train_X_for_shap = train_X_for_shap_final
  )
}



# Repeated nested CV wrapper ------------------------------------------------------


# 1. Compute adjusted R²
compute_adjR2 <- function(y_true, y_pred, p) {
  n <- length(y_true)
  rss <- sum((y_true - y_pred)^2)
  tss <- sum((y_true - mean(y_true))^2)
  
  r2 <- 1 - rss/tss
  r2_adj <- 1 - (1 - r2) * (n - 1) / (n - p - 1)
  
  tibble::tibble(R2 = r2, R2_adj = r2_adj)
}

# 2. Append adjusted R² to metrics
add_adjR2_to_metrics <- function(metrics_tbl, y_true, y_pred, predictors) {
  p <- if (inherits(predictors, "formula")) {
    length(attr(terms(predictors), "term.labels"))
  } else {
    length(predictors)
  }
  
  adj <- compute_adjR2(y_true, y_pred, p)
  dplyr::bind_cols(metrics_tbl, adj %>% dplyr::select(R2_adj))
}

# 3. Repeated Nested CV (with adjR²)

nested_fit_regression_repeated <- function(data, outcome, predictors, method,
                                           K_outer = 5, inner_k = 5,
                                           repeats = 5, seed = 12345, scale = "auto",
                                           outer_folds_repeats = NULL,
                                           show_progress = TRUE) {
  n <- nrow(data); y <- data[[outcome]]
  
  per_rep <- vector("list", repeats)
  oof_mat <- matrix(NA_real_, nrow = n, ncol = repeats)
  
  if (show_progress && interactive()) {
    cat(sprintf("[Repeated nested] %s | repeats=%d, outer=%d, inner=%d\n",
                method, repeats, K_outer, inner_k))
  }
  
  for (r in seq_len(repeats)) {
    if (show_progress && interactive()) cat(sprintf("\r  repeat %d/%d ...", r, repeats))
    
    ofolds <- if (is.null(outer_folds_repeats)) NULL else outer_folds_repeats[[r]]
    
    res <- nested_fit_regression(
      data        = data,
      outcome     = outcome,
      predictors  = predictors,
      method      = method,
      K_outer     = K_outer,
      inner_k     = inner_k,
      seed        = seed + 1000L * r,
      scale       = scale,
      outer_folds = ofolds
    )
    
    per_rep[[r]] <- res
    oof_mat[, r] <- res$oof_pred
  }
  
  if (show_progress && interactive()) cat("\r  repeats done        \n")
  
  # (1) metrics per repeat (with adjR²)
  rep_stats <- purrr::map_dfr(seq_len(repeats), function(r) {
    s <- extract_metrics_reg(y, oof_mat[, r])
    s <- add_adjR2_to_metrics(s, y, oof_mat[, r], predictors)
    dplyr::mutate(s, repeat_id = r, .before = 1)
  })
  
  # (2) metrics on averaged OOF prediction across repeats
  oof_mean <- rowMeans(oof_mat, na.rm = TRUE)
  stats_on_mean <- extract_metrics_reg(y, oof_mean)
  stats_on_mean <- add_adjR2_to_metrics(stats_on_mean, y, oof_mean, predictors)
  
  # (3) mean ± sd across repeats
  mean_sd <- rep_stats %>%
    dplyr::summarise(
      RMSE_mean = mean(RMSE), RMSE_sd = sd(RMSE),
      MAE_mean  = mean(MAE),  MAE_sd  = sd(MAE),
      R2_mean   = mean(R2),   R2_sd   = sd(R2),
      R2_adj_mean = mean(R2_adj), R2_adj_sd = sd(R2_adj),
      Bias_mean = mean(Bias), Bias_sd = sd(Bias),
      MAPE_mean = mean(MAPE), MAPE_sd = sd(MAPE)
    )
  
  list(
    per_repeat     = per_rep,      # list of results (each from nested_fit_regression)
    oof_pred_mat   = oof_mat,      # n x repeats matrix
    oof_pred_mean  = oof_mean,     # averaged OOF prediction (length n)
    rep_stats      = rep_stats,    # tibble per repeat
    stats_on_mean  = stats_on_mean,# tibble on averaged OOF prediction
    mean_sd        = mean_sd       # tibble mean±sd across repeats
  )
}



#── Regression metrics ────────────────────────────────────────────────────────
extract_metrics_reg <- function(truth, pred) {
  stopifnot(length(truth) == length(pred))
  resid <- pred - truth
  tibble::tibble(
    RMSE = sqrt(mean(resid^2)),
    MAE  = mean(abs(resid)),
    R2cor= cor(truth, pred, use = "complete.obs")^2,
    R2   = 1 - sum(resid^2) / sum((truth - mean(truth))^2),
    Bias = mean(resid),
    MAPE = mean(abs(resid) / pmax(abs(truth), .Machine$double.eps)) * 100
  )
}


# GLM subset specifications generator ---------------------------------------------

expand_specs4 <- function(core_sets, bio_blood, bio_csf, bio_one_any, bio_one_csf) {
  specs <- list()
  
  for (core_name in names(core_sets)) {
    core_vars <- core_sets[[core_name]]
    
    # (0) Core-only
    specs[[core_name]] <- core_vars
    
    # (1) Core + NfL
    specs[[paste0(core_name, "_NfL")]] <- c(core_vars, "NfL_csf_pgml")
    
    # (2) Core + NfL + one biomarker (additive)
    for (bio in bio_one_any) {
      base_nm   <- paste0(core_name, "_NfL_", bio)      # e.g., Core1_NfL_pTau181_csf
      base_covs <- c(core_vars, "NfL_csf_pgml", bio)
      specs[[base_nm]] <- base_covs
      
      # (3) Interaction (bio-only): add exactly NfL : bio
      nm_int <- paste0(base_nm, "_x", bio)              # e.g., Core1_NfL_pTau181_csf_xpTau181_csf
      rhs <- paste(
        paste(core_vars, collapse = " + "),
        "+ NfL_csf_pgml +", bio, "+ NfL_csf_pgml :", bio
      )
      specs[[nm_int]] <- rhs                             # store as a single formula string
    }
  }
  specs
}

# Main function: run GLM subsets with repeated nested CV -----------------------

run_glm_subsets <- function(data, outcome, specs,
                            K_outer = 5, inner_k = 5, repeats = 5,
                            seed = 12345, scale = "auto",
                            outer_folds_rep = outer_reps) {
  
  # guess ID col from data
  guess_id_col <- function(df) {
    cands <- c("id","ID","subject_id","SubjectID","patient_id","PatientID","pid")
    cands[cands %in% names(df)][1]
  }
  id_col <- guess_id_col(data)
  id_vec <- if (!is.na(id_col)) data[[id_col]] else seq_len(nrow(data))
  
  # Iterate over specs and run repeated nested CV per model
  results_nested <- purrr::imap_dfr(specs, function(preds, spec_name) {
    
    # (A) Additive model: predictors as char vector
    if (is.character(preds) && length(preds) > 1) {
      res <- nested_fit_regression_repeated(
        data       = data,
        outcome    = outcome,
        predictors = preds,
        method     = "glm",
        K_outer    = K_outer,
        inner_k    = inner_k,
        repeats    = repeats,
        seed       = seed,
        scale      = scale,
        outer_folds_repeats = outer_folds_rep
      )
    } else {
      # (B) Interaction model: single string → formula
      if (is.character(preds) && length(preds) == 1) {
        preds <- as.formula(paste(outcome, "~", preds))
      }
      res <- nested_fit_regression_repeated(
        data       = data,
        outcome    = outcome,
        predictors = preds,
        method     = "glm",
        K_outer    = K_outer,
        inner_k    = inner_k,
        repeats    = repeats,
        seed       = seed,
        scale      = scale,
        outer_folds_repeats = outer_folds_rep
      )
    }
    
    # coerce oof_pred_mean into (id,y,yhat)
    oof_raw <- res$oof_pred_mean
    yhat <- if (is.data.frame(oof_raw) && ncol(oof_raw) == 1) oof_raw[[1]] else as.numeric(oof_raw)
    stopifnot(length(yhat) == nrow(data))
    oof_tbl <- tibble::tibble(
      id   = as.character(id_vec),
      y    = data[[outcome]],
      yhat = yhat
    )
    
    tibble::tibble(
      model         = spec_name,
      stats_on_mean = list(res$stats_on_mean),
      rep_stats     = list(res$rep_stats),
      mean_sd       = list(res$mean_sd),
      oof_pred_mean = list(oof_tbl)  
    )
  })
  
  # Flat table 1: mean ± sd
  flat_mean_sd <- results_nested %>%
    dplyr::select(model, mean_sd) %>%
    tidyr::unnest(mean_sd) %>%
    dplyr::relocate(model, RMSE_mean, RMSE_sd, MAE_mean, MAE_sd,
                    R2_mean, R2_sd, R2_adj_mean, R2_adj_sd,
                    Bias_mean, Bias_sd, .before = tidyselect::last_col())
  
  # Flat table 2: on mean preds
  flat_on_mean <- results_nested %>%
    dplyr::select(model, stats_on_mean) %>%
    tidyr::unnest(stats_on_mean) %>%
    dplyr::relocate(model, RMSE, MAE, R2, R2_adj, Bias, MAPE)
  
  list(
    nested       = results_nested,
    flat_mean_sd = flat_mean_sd,
    flat_on_mean = flat_on_mean
  )
}


# SHAP --------------------------------------------------------------------


# Build test design matrix matching training preprocessing
build_X_for_shap <- function(fit_obj, test_data) {
  # Extract features from test data (raw)
  X <- test_data[, setdiff(colnames(test_data), "slope_total"), drop = FALSE]
  
  # Apply dummy encoding using stored dummy object
  if (!is.null(fit_obj$dummy)) {
    # Use caret's predict.dummyVars to transform test data
    X_dummy <- predict(fit_obj$dummy, newdata = X)
    X <- as.data.frame(X_dummy)
  }
  
  # Apply scaling if needed
  if (!is.null(fit_obj$scaled) && fit_obj$scaled && 
      !is.null(fit_obj$scale_mean) && !is.null(fit_obj$scale_sd)) {
    numeric_cols <- intersect(names(fit_obj$scale_mean), colnames(X))
    for (col in numeric_cols) {
      if (!is.na(fit_obj$scale_sd[[col]]) && fit_obj$scale_sd[[col]] > 0) {
        X[[col]] <- (X[[col]] - fit_obj$scale_mean[[col]]) / fit_obj$scale_sd[[col]]
      }
    }
  }
  
  # Ensure column order matches training features
  X <- X[, fit_obj$features, drop = FALSE]
  
  return(as.data.frame(X))
}



# --- Nested CV runner (keeps fits & test sets for pooled SHAP) ---
nested_fit_regression <- function(data, outcome, predictors, method,
                                  K_outer = 5, inner_k = 5,
                                  seed = 12345, scale = "auto",
                                  outer_folds = NULL) {
  set.seed(seed)
  n <- nrow(data)
  y <- data[[outcome]]
  
  # If folds are not provided, create them here
  if (is.null(outer_folds)) {
    outer_folds <- caret::createFolds(y, k = K_outer, list = TRUE, returnTrain = FALSE)
  } else {
    # Validate provided folds
    stopifnot(is.list(outer_folds), length(outer_folds) == K_outer)
  }
  
  # Inner CV control (for hyperparameter tuning)
  inner_ctrl <- caret::trainControl(method = "cv", number = inner_k, savePredictions = "final")
  
  oof_pred   <- rep(NA_real_, n)
  fold_stats <- vector("list", K_outer)
  fit_list   <- vector("list", K_outer)
  test_list  <- vector("list", K_outer)
  
  for (i in seq_len(K_outer)) {
    te_idx <- outer_folds[[i]]
    tr_idx <- setdiff(seq_len(n), te_idx)
    tr <- data[tr_idx, , drop = FALSE]
    te <- data[te_idx,  , drop = FALSE]
    
    fit <- fit_model_reg(
      train_data = tr, test_data = te,
      outcome = outcome, predictors = predictors,
      method = method, seed = seed + i, scale = scale,
      tr_control = inner_ctrl
    )
    
    oof_pred[te_idx] <- fit$pred
    fold_stats[[i]]  <- extract_metrics_reg(truth = te[[outcome]], pred = fit$pred) |>
      dplyr::mutate(outer_fold = i, .before = 1)
    fit_list[[i]]    <- fit
    test_list[[i]]   <- te
  }
  
  list(
    oof_pred   = oof_pred,
    oof_stats  = extract_metrics_reg(truth = data[[outcome]], pred = oof_pred),
    fold_stats = dplyr::bind_rows(fold_stats),
    fit_list   = fit_list,
    test_list  = test_list
  )
}
