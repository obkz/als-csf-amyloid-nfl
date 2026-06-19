#── PACKAGES ────────────────────────────────────────────────────────────────
library(dplyr)
library(caret)
library(pROC)
library(xgboost)
library(fastshap)
library(ggplot2)
library(ggbeeswarm)   # for the beeswarm geom
library(ggforce)      # for geom_sina
library(tidyr)
library(viridisLite)
library(cowplot)
library(patchwork)
library(flextable)

# ------ COMMON TRAIN CONTROL -----------------------------------------------
common_ctrl <- trainControl(
  method          = "cv",           # k‐fold CV
  number          = 5,              # 5 to 10 folds
  classProbs      = TRUE,           # needed for AUC
  summaryFunction = twoClassSummary,# measure ROC, sens, spec
  savePredictions = "final"
)

# positive‐class label (must match the factor levels)
pos_class <- "ALS"

.make_predict_fun <- function(model, method, pos, predictors, scaler = NULL) {
  force(method); force(predictors); force(scaler)
  if (method == "one") {
    return(function(newdata) {
      x_raw <- newdata[[predictors]]
      mu <- scaler$mean; sdv <- scaler$sd
      x <- (x_raw - mu) / sdv
      predict(model, newdata = data.frame(x = x), type = "response")
    })
  }
  if (method == "xgb") {
    return(function(newdata) {
      X <- data.matrix(newdata[, predictors, drop = FALSE])
      dnew <- xgboost::xgb.DMatrix(X, missing = NA_real_)
      as.numeric(predict(model, dnew))
    })
  }
  # caret 系（rf / glm / svmRadial / knn など）
  return(function(newdata) {
    as.numeric(predict(model, newdata = newdata, type = "prob")[, pos])
  })
}


# 1) FIT MODEL -----------------------------------------------

fit_model <- function(train_data,
                      test_data,
                      outcome,       #ALS_label
                      predictors,    #list of predictors
                      method
                      ) {
  
  # ---- define class labels ONCE (used by all branches) ----
  
  # normalize predictors to a plain character vector
  predictors <- as.character(unlist(predictors, use.names = FALSE))
  stopifnot(length(predictors) >= 1L)
  stopifnot(all(predictors %in% names(train_data)))
  
  pos <- pos_class                         # e.g., "ALS"
  #neg <- paste0("non", pos_class)         # e.g., "nonALS"
  lv <- levels(factor(c(train_data[[outcome]], test_data[[outcome]])))
  stopifnot(length(lv) == 2, pos %in% lv)
  neg <- setdiff(lv, pos)
  train_data[[outcome]] <- factor(train_data[[outcome]], levels = c(pos, neg))
  test_data [[outcome]] <- factor(test_data [[outcome]], levels = c(pos, neg))

  train_data <- as.data.frame(train_data)
  test_data  <- as.data.frame(test_data)
  Xtr <- train_data[, predictors, drop = FALSE]
  Xte <- test_data [, predictors, drop = FALSE]

  stopifnot(all(vapply(train_data[, predictors, drop=FALSE], is.numeric, logical(1))))
  stopifnot(all(vapply(test_data [, predictors, drop=FALSE], is.numeric, logical(1))))

  # build formula once
  fmla <- as.formula(paste(outcome, "~", paste(predictors, collapse = " + ")))

  # ---------------- ONE-VARIABLE GLM BRANCH ----------------
  if (method == "one") {
    # validate: single predictor, numeric
    if (length(predictors) != 1) stop("`method='one'` expects exactly 1 predictor.")
    xtr_raw <- train_data[[predictors]]
    xte_raw <- test_data [[predictors]]
    if (!is.numeric(xtr_raw)) stop("Predictor for `one` must be numeric.")

    # scale by training stats
    mu <- mean(xtr_raw, na.rm = TRUE)
    sdv <- stats::sd(xtr_raw, na.rm = TRUE); if (isTRUE(all.equal(sdv, 0))) sdv <- 1
    xtr <- (xtr_raw - mu) / sdv
    xte <- (xte_raw - mu) / sdv

    # inner CV (OOF) to fix threshold (Youden) without touching outer test
    y_tr_fac <- factor(train_data[[outcome]], levels = c(neg, pos))
    folds_in <- caret::createFolds(y_tr_fac, k = 5, returnTrain = FALSE)
    oof_prob <- numeric(length(xtr)); oof_prob[] <- NA_real_

    for (idx_te_in in folds_in) {
      idx_tr_in <- setdiff(seq_along(xtr), idx_te_in)
      df_in <- data.frame(y = as.integer(y_tr_fac[idx_tr_in] == pos),
                          x = xtr[idx_tr_in])
      mod_in <- glm(y ~ x, family = binomial(), data = df_in)
      oof_prob[idx_te_in] <- predict(mod_in,
                                     newdata = data.frame(x = xtr[idx_te_in]),
                                     type = "response")
    }

    # ROC on OOF (training-only) → fixed cutoff
    roc_oof <- pROC::roc(response  = y_tr_fac,
                         predictor = oof_prob,
                         levels    = c(neg, pos),
                         direction = "<")
    th_star <- as.numeric(pROC::coords(roc_oof, x = "best",
                                       best.method = "youden",
                                       ret = "threshold", transpose = FALSE))

    # fit final GLM on full training (same scaling)
    df_full <- data.frame(y = as.integer(y_tr_fac == pos), x = xtr)
    mod <- glm(y ~ x, family = binomial(), data = df_full)

    # predict on outer test
    y_te_fac <- factor(test_data[[outcome]], levels = c(neg, pos))
    prob <- predict(mod, newdata = data.frame(x = xte), type = "response")

    # outer test ROC (no threshold tuning)
    roc_obj <- pROC::roc(response  = y_te_fac,
                         predictor = prob,
                         levels    = c(neg, pos),
                         direction = "<")

    # simple varImp-like (absolute coef for the single variable)
    vip_df <- data.frame(Overall = abs(coef(mod)["x"]), row.names = predictors)
    vip_like <- list(importance = vip_df, model = "one", calledFrom = "varImp")
    class(vip_like) <- "varImp.train"

    # predict function for SHAP later
    predict_fun <- .make_predict_fun(mod, "one", pos, predictors, scaler = list(mean = mu, sd = sdv))

    return(list(
      method ="one",
      model = mod, roc = roc_obj, vip = vip_like,
      th_star = th_star, roc_oof = roc_oof,
      scaler = list(mean = mu, sd = sdv),
      test_prob = prob,
      predict_fun = predict_fun
      )
           )
  }

  # XGB pattern ---------------------------------------------------------------------
  if (method == "xgb") {
    #────────────────────────────────────────────────────────────
    # Train-time model selection uses xgb.cv on TRAIN only
    # to pick the optimal number of boosting rounds (best_iteration)
    # with early stopping. External TEST is untouched -> no leakage.
    #────────────────────────────────────────────────────────────

    # Ensure pure numeric matrices for XGBoost
    Xtr_df <- train_data[, predictors, drop = FALSE]
    Xte_df <- test_data [, predictors, drop = FALSE]

    # Defensive: abort if any list-columns slipped in
    if (any(vapply(Xtr_df, is.list, logical(1))) || any(vapply(Xte_df, is.list, logical(1)))) {
      stop("XGBoost branch: predictors produced list-columns. ",
           "Please unnest/convert those columns before modelling.")
    }

    Xtr <- data.matrix(Xtr_df)
    Xte <- data.matrix(Xte_df)

    # Build DMatrices; be explicit about missing
    dtrain <- xgboost::xgb.DMatrix(
      data  = Xtr,
      label = as.integer(train_data[[outcome]] == pos_class),
      missing = NA_real_
    )
    dtest <- xgboost::xgb.DMatrix(
      data = Xte,
      missing = NA_real_
    )

    # XGBoost parameters (safe defaults). Add nthread if you want parallelism control.
    params <- list(
      objective   = "binary:logistic",
      eval_metric = "auc"
      # , nthread = parallel::detectCores()  # uncomment if you want to pin threads
    )

    # Cross-validated training (OOF predictions are returned via prediction=TRUE)
    cv_res <- xgboost::xgb.cv(
      params                = params,
      data                  = dtrain,
      nrounds               = 2000,       # high cap; early stopping stops earlier
      nfold                 = 5,          # aligned with caret CV
      stratified            = TRUE,
      prediction            = TRUE,       # OOF probabilities
      early_stopping_rounds = 20,
      verbose               = 0
    )

    # Robust best iteration: prefer cv_res$best_iteration, else infer from eval log
    best_iter <- cv_res$best_iteration
    if (is.null(best_iter) || is.na(best_iter)) {
      el <- cv_res$evaluation_log
      # test_auc_mean column name can vary across versions; be defensive
      auc_col <- intersect(colnames(el), c("test_auc_mean", "test-auc-mean"))
      best_iter <- if (length(auc_col) == 1) which.max(el[[auc_col]]) else nrow(el)
    }

    # OOF predictions & training-side ROC (for threshold fixing)
    oof_prob  <- cv_res$pred
    oof_label <- factor(train_data[[outcome]], levels = c(neg, pos_class))

    roc_oof <- pROC::roc(
      response  = oof_label,
      predictor = oof_prob,
      levels    = c(neg, pos_class),
      direction = "<"
    )
    th_star <- as.numeric(pROC::coords(
      roc_oof, x = "best", best.method = "youden",
      ret = "threshold", transpose = FALSE
    ))

    # Final training on full TRAIN using the selected number of rounds
    xgb_mod <- xgboost::xgb.train(
      params        = params,
      data          = dtrain,
      nrounds       = best_iter,
      print_every_n = 10,
      verbose       = 0
    )

    # Predict probabilities on TEST
    prob <- predict(xgb_mod, dtest)

    # External test ROC
    roc_obj <- pROC::roc(
      response  = test_data[[outcome]],
      predictor = prob,
      levels    = c(neg,pos), #c(paste0("non", pos_class), pos_class),
      direction = "<"
    )

    # Variable importance (Gain). Use matrix colnames in case predictors were transformed.

    if (ncol(Xtr) == 1L) {
      vip_df <- data.frame(Overall = 1, row.names = colnames(Xtr))
    } else {
      vip_matrix <- tryCatch(
        xgboost::xgb.importance(feature_names = colnames(Xtr), model = xgb_mod),
        error = function(e) NULL
      )

      if (is.null(vip_matrix) || NROW(vip_matrix) == 0) {
        vip_df <- data.frame(Overall = NA_real_, row.names = colnames(Xtr))
      } else {
        vip_df0 <- as.data.frame(vip_matrix, stringsAsFactors = FALSE)

        if (!("Feature" %in% names(vip_df0))) {
          vip_df0$Feature <- colnames(Xtr)[seq_len(NROW(vip_df0))]
        }
        if (!("Gain" %in% names(vip_df0))) {
          vip_df0$Gain <- NA_real_
        }

        vip_df <- vip_df0[, c("Feature","Gain"), drop = FALSE]
        names(vip_df) <- c("Feature","Overall")
        rownames(vip_df) <- vip_df$Feature
        vip_df$Feature <- NULL
      }
    }

    vip_like_varImp <- list(
      importance = vip_df,
      model      = "xgb",
      calledFrom = "varImp"
    )
    class(vip_like_varImp) <- "varImp.train"

    predict_fun <- .make_predict_fun(xgb_mod, "xgb", pos, predictors, scaler = NULL)

    return(list(
      model      = xgb_mod,
      roc        = roc_obj,
      vip        = vip_like_varImp,
      predict_fun= predict_fun,
      method     = "xgb",
      th_star    = th_star,     # fixed on TRAIN OOF (Youden)
      roc_oof    = roc_oof,     # training-side ROC
      test_prob  = prob
    ))
  }

  else {

    # -------------------- BRANCH: caret models (rf/glm/svmRadial/kknn/...) ---

    # Use x/y interface so caret only sees our cleaned numeric predictor matrix.
    # This avoids any formula/model.matrix surprises with list-like columns elsewhere.

    # 1) Build clean numeric X for TRAIN/TEST over *only* `predictors`
    Xtr_df <- as.data.frame(train_data[, predictors, drop = FALSE])
    Xte_df <- as.data.frame(test_data [, predictors, drop = FALSE])

    # hard-coerce every column to plain numeric atomic vectors (no list/matrix)
    fix_num <- function(x) {
      if (is.data.frame(x)) { if (ncol(x)==1L) x <- x[[1]] else stop("data.frame col >1") }
      if (is.matrix(x)) x <- x[,1]
      if (is.list(x))  x <- vapply(x, function(z) as.numeric(z[[1]]), numeric(1))
      if (is.factor(x)) x <- as.character(x)
      if (is.character(x)) x <- suppressWarnings(as.numeric(x))
      x[!is.finite(x)] <- NA_real_
      as.numeric(x)
    }
    Xtr_df[] <- lapply(Xtr_df, fix_num)
    Xte_df[] <- lapply(Xte_df, fix_num)

    # 2) Target as factor with correct order
    
    # y_tr <- factor(train_data[[outcome]], levels = c(neg, pos))
    # y_te <- factor(test_data [[outcome]], levels = c(neg, pos))
    
    extract_outcome <- function(df, col) {
      x <- df[[col]]
      if (is.list(x) && !is.factor(x)) x <- unlist(x, use.names = FALSE)
      if (is.matrix(x)) x <- x[, 1]
      if (is.data.frame(x)) x <- x[[1]]
      factor(as.character(x), levels = c(neg, pos))
    }
    
    y_tr <- extract_outcome(train_data, outcome)
    y_te <- extract_outcome(test_data, outcome)
    
    stopifnot(is.factor(y_tr), is.atomic(y_tr), length(y_tr) == nrow(train_data))
    stopifnot(is.factor(y_te), is.atomic(y_te), length(y_te) == nrow(test_data))


    # 3) Train via caret using x/y (preProcess only touches X we pass)

    tune_grid <- NULL
    if (length(predictors) == 1L && method == "rf") {
      tune_grid <- data.frame(mtry = 1L)
    }
    
    #=== KKNN specific tuning grid (larger k, distance, kernel) ===
    if (method == "kknn") {
      tune_grid <- expand.grid(
        kmax    = c(15, 25, 35, 45, 60),   # big k
        distance= c(1, 2),                 # Minkowski(=1,2)
        kernel  = c("optimal","triangular","epanechnikov")  
      )
    }

    caret_mod <- caret::train(
      x          = Xtr_df,
      y          = y_tr,
      method     = method,                # "rf", "svmRadial", "kknn", "glm", ...
      metric     = "ROC",
      trControl  = common_ctrl,           # twoClassSummary / savePredictions="final"
      preProcess = c("center", "scale"),  # keep as you had it
      tuneGrid   = tune_grid              # for single-predictor RF edge case
    )

    # 4) OOF ROC from caret's saved predictions
    pred_df <- caret_mod$pred
    bt <- caret_mod$bestTune
    if (!is.null(bt) && nrow(bt) > 0 && !is.null(pred_df)) {
      for (nm in names(bt)) pred_df <- pred_df[pred_df[[nm]] == bt[[nm]], , drop = FALSE]
    }
    
    pred_oof <- pred_df[[pos_class]]
    
    # Safe extraction handling different return types
    if (is.data.frame(pred_oof)) {
      pred_oof <- pred_oof[[1]]
    } else if (is.matrix(pred_oof)) {
      pred_oof <- pred_oof[, 1]
    } else if (is.list(pred_oof)) {
      pred_oof <- unlist(pred_oof, recursive = TRUE, use.names = FALSE)
    }
    
    pred_oof <- as.vector(pred_oof)
    pred_oof <- as.numeric(pred_oof)
    storage.mode(pred_oof) <- "double"
    
    if (!is.numeric(pred_oof) || !is.atomic(pred_oof) || is.matrix(pred_oof)) {
      stop("OOF predictor is not an atomic numeric vector (pos_class column).")
    }

    roc_oof <- pROC::roc(
      response  = factor(pred_df$obs, levels = c(neg, pos_class)),
      predictor = pred_df[[pos_class]],
      levels    = c(neg, pos_class),
      direction = "<"
    )
    th_star <- as.numeric(pROC::coords(
      roc_oof, x = "best", best.method = "youden",
      ret = "threshold", transpose = FALSE
    ))

    # 5) External TEST prediction (use our clean TEST X)
    
    pred_result <- predict(caret_mod, newdata = Xte_df, type = "prob")
    
    # Safe extraction handling different return types
    if (is.data.frame(pred_result)) {
      prob <- pred_result[[pos_class]]
    } else if (is.matrix(pred_result)) {
      prob <- pred_result[, pos_class]
    } else if (is.list(pred_result)) {
      prob <- pred_result[[pos_class]]
    } else {
      prob <- pred_result
    }
    
    # Triple ensure it's numeric
    prob <- as.vector(prob)
    prob <- as.numeric(prob)
    storage.mode(prob) <- "double"
    
    # Safety check
    if (!is.numeric(prob) || !is.atomic(prob) || is.matrix(prob)) {
      stop("Failed to convert prediction to atomic numeric vector")
    }
    
    roc_obj <- pROC::roc(
      response  = y_te,
      predictor = prob,
      levels    = c(neg, pos_class),
      direction = "<"
    )

    # 6) varImp policy unchanged
    if (method %in% c("knn","kknn","svmRadial")) {
      VIP_caret <- NULL
    } else {
      VIP_caret <- caret::varImp(caret_mod)
    }

    predict_fun <- .make_predict_fun(caret_mod, method, pos, predictors, scaler = NULL)

    return(list(
      method     = method,
      model      = caret_mod,
      roc        = roc_obj,
      vip        = VIP_caret,
      predict_fun= predict_fun,
      th_star    = th_star,
      roc_oof    = roc_oof,
      test_prob  = prob
    ))

  }
}


# 2) EXTRACT METRICS -----------------------------------------------

extract_metrics <- function(roc_obj, cutoff, boot_n=2000, stratified=TRUE,
                            calcCI = TRUE) {
  if (missing(cutoff) || is.null(cutoff) || !is.finite(cutoff)) {
    stop("`cutoff` must be provided and finite. Decide it on training/CV (e.g., OOF).")
  }
  
  # Positive / negative class labels
  neg_class <- roc_obj$levels[1]
  pos_class <- roc_obj$levels[2]
  
  # AUC and its 95% CI (threshold-independent)
  auc_val <- pROC::auc(roc_obj)
  
  # Operating point metrics at the provided cutoff
  coords_res <- pROC::coords(
    roc_obj,
    x         = cutoff,
    input     = "threshold",
    ret       = c("sensitivity", "specificity", "tp", "tn", "fn", "fp", "accuracy", "threshold"),
    transpose = FALSE
  )
  
  # Extract confusion matrix counts
  TP <- as.numeric(coords_res["tp"])
  TN <- as.numeric(coords_res["tn"])
  FP <- as.numeric(coords_res["fp"])
  FN <- as.numeric(coords_res["fn"])
  
  conf_mat <- matrix(
    c(TP, FP, FN, TN),
    nrow = 2, byrow = TRUE,
    dimnames = list(
      Prediction = c(pos_class, neg_class),
      Reference  = c(pos_class, neg_class)
    )
  )
  
  # Default: no CIs (used for nested CV aggregation)
  AUC_lower <- AUC_upper <- Sens_lower <- Sens_upper <- Spec_lower <- Spec_upper <- NA_real_
  
  if (calcCI) {
    # AUC CI (bootstrap, stratified if specified)
    auc_ci <- pROC::ci.auc(
      roc_obj,
      method           = "bootstrap",
      boot.n           = boot_n,
      boot.stratified  = stratified,
      conf.level       = 0.95
    )
    AUC_lower <- as.numeric(auc_ci["2.5%"])
    AUC_upper <- as.numeric(auc_ci["97.5%"])
    
    # Sensitivity CI at fixed cutoff
    sens_ci <- pROC::ci.se(
      roc_obj,
      specificities    = as.numeric(coords_res["specificity"]),
      boot.n           = boot_n,
      boot.stratified  = stratified
    )
    Sens_lower <- as.numeric(sens_ci[1, "2.5%"])
    Sens_upper <- as.numeric(sens_ci[1, "97.5%"])
    
    # Specificity CI at fixed cutoff
    spec_ci <- pROC::ci.sp(
      roc_obj,
      sensitivities    = as.numeric(coords_res["sensitivity"]),
      boot.n           = boot_n,
      boot.stratified  = stratified
    )
    Spec_lower <- as.numeric(spec_ci[1, "2.5%"])
    Spec_upper <- as.numeric(spec_ci[1, "97.5%"])
  }
  
  # Output metrics table
  metrics <- tibble::tibble(
    AUC        = as.numeric(auc_val),
    AUC_lower  = AUC_lower,
    AUC_upper  = AUC_upper,
    optCut_off = as.numeric(cutoff),
    Sens       = as.numeric(coords_res["sensitivity"]),
    Sens_lower = Sens_lower,
    Sens_upper = Sens_upper,
    Spec       = as.numeric(coords_res["specificity"]),
    Spec_lower = Spec_lower,
    Spec_upper = Spec_upper,
    accuracy   = as.numeric(coords_res["accuracy"]),
    TP = TP, TN = TN, FP = FP, FN = FN  # raw counts for later aggregation
  )
  
  return(list(metrics = metrics, matrix = conf_mat))
}






#── 3) confusion matrix ──────────────────────────────────────────────────────

# Suppose `result$matrix` is the 2×2 matrix:
#
#              Reference
# Prediction  ALS nonALS
#     ALS     42      8
#     nonALS   5     30
#

plotConfMatrix <- function(conf_matrix, title = NULL) {
  require(ggplot2)
  
  # 1) Convert the confusion matrix (2×2) into a tidy data frame
  cm_df <- as.data.frame(as.table(conf_matrix))
  colnames(cm_df) <- c("Prediction", "Reference", "Count")
  
  # 2) Add a fill color: gray for TP/TN (where Prediction == Reference), white otherwise
  cm_df$fill_color <- ifelse(
    as.character(cm_df$Prediction) == as.character(cm_df$Reference),
    "gray80",  # true positives / true negatives
    "white"    # false positives / false negatives
  )
  
  # 3) Build the ggplot
  p <- ggplot(cm_df, aes(x = Reference, y = Prediction, fill = fill_color)) +
    geom_tile(color = "black") +              # Draw tiles with black borders
    geom_text(aes(label = Count), size = 8) + # Overlay the counts in large text
    scale_fill_identity() +                   # Use the fill_color column directly (no legend)
    scale_x_discrete(position = "top") +      # Draw the x-axis on top
    labs(                                     # Axis labels and optional title
      x     = "Reference",
      y     = "Prediction",
      title = title
    ) +
    theme_minimal(base_size = 14) +           # Minimal theme as base
    theme(                                    # Remove ticks, axis lines, and grid lines
      axis.ticks         = element_blank(),
      axis.line          = element_blank(),
      panel.grid.major   = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.title.y       = element_text(angle = 90, vjust = 0.5, hjust = 0.5), # Rotate Y-axis title and center it
      axis.text.y        = element_text(angle = 90, vjust = 0.5, hjust = 0.5), # Rotate Y-axis labels and center them
      axis.text.x.top    = element_text(vjust = -5),       # Nudge X-axis labels downward (they're on top)
      axis.text.y.left   = element_text(vjust = -5)        # Nudge Y-axis labels inward if needed
    ) +
    coord_fixed(ratio = 1)                    # Fix aspect ratio so each cell is square
  
  # Print the plot
  print(p)
  # Return the ggplot object (invisible) for further customization if desired
  invisible(p)
}

#── 4)   PLOT SHAP ─────────────────────────────────────────────────────────────
# Given a trained model and test_data (predictors only),
# computes fastshap::explain() and returns a ggplot of
# feature‐wise SHAP beeswarm (geom_sina).
SHAPlot <- function(model, 
                    method,
                    test_data,
                    predictors,
                    seed=12345,
                    title="SHAP Summary Plot") {
  # wrapper for caret models:
  predfun <- function(object, newdata) {
    predict(object, newdata=newdata, type="prob")[, pos_class]
  }
  
  X <- test_data[, predictors]
  if (method=="xgb") {
    # for xgb.Booster
    predfun <- function(object, newdata) {
      predict(object, xgb.DMatrix(as.matrix(newdata)))
    }
  }
  
  # set seed
  set.seed(seed)
  
  # 1) compute SHAP contributions
  shap_vals <- fastshap::explain(
    object       = model,     # trained model
    X            = X,         # data.frame of test set predictors
    pred_wrapper = predfun,   # wrapper to return class probabilities for SHAP calculation
    nsim         = 50         # number of Monte Carlo simulations per observation
  )
  
  # 2) pivot to long form
  # Each row will represent a single SHAP value for one feature in one observation
  shap_long <- as.data.frame(shap_vals) %>%
    mutate(obs = row_number()) %>%   # Add row index to identify observation (=Patient)
    pivot_longer(
      cols      = -obs,              # Pivot all columns except observation index
      names_to  = "feature",         # Feature name becomes a column
      values_to = "phi"              # SHAP value for that feature
    ) %>%
    left_join(
      X %>% mutate(obs = row_number()) %>%
        pivot_longer(-obs, names_to="feature", values_to="feature_value"),
      by = c("obs","feature")
    )
  
  # 3) color by feature’s rank within each feature
  shap_long <- shap_long %>%
    group_by(feature) %>%
    # Normalize feature_value by feature
    mutate(
      rankfvalue = rank(feature_value),
      stdfvalue = (rankfvalue - min(rankfvalue)) / 
        (max(rankfvalue) - min(rankfvalue))
    ) %>% ungroup()
  
  #arrange
  mean_phi <- shap_long %>%
    group_by(feature) %>%
    summarise(mean_value = mean(abs(phi)), .groups = "drop") %>% 
    arrange(mean_value)
  
  shap_long$feature <- factor(shap_long$feature, levels = mean_phi$feature)
  
  # 4) plot
  ggplot(shap_long, aes(x = feature, y = phi, color = stdfvalue)) +
    geom_sina(size = 1.5, alpha = 0.8, scale="area") +
    coord_flip() +
    geom_text(
      data = mean_phi,
      aes(x = feature, y = -Inf, label = sprintf("%.3f", mean_value)),
      hjust = -0.2, fontface = "bold", size = 3,
      inherit.aes = FALSE
    ) +
    scale_color_gradient(low = "#FFCC33", high = "#6600CC", 
                         name = "Feature value", breaks = c(0, 1), labels = c("Low", "High")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
      title = title,
      x = NULL,
      y = "SHAP value (impact on model output)"
    ) +
    theme_bw(base_size = 13) +
    theme(axis.ticks.y = element_blank(), legend.position = "bottom")
  
}


# SUMMARY PLOT -----------------------------------------------------------------

make_model_metrics_plot <- function(all_metrics,
                                    model_name_map = c(
                                      knn = "KNN",
                                      kknn = "KNN",
                                      svmRadial = "SVM",
                                      xgb = "XGB",
                                      rf = "RF",
                                      glm = "GLM",
                                      NfLonly = "NfL alone"
                                    ),
                                    model_order = c("GLM","KNN","SVM","XGB","RF","NfL alone"),
                                    metrics_to_show = c("AUC","Sens","Spec"),
                                    dashed_at = c(AUC = 0.5, Sens = 0.5, Spec = 0.5),
                                    sort_by_auc = c("none","desc","asc"),
                                    colormap = NULL) {
  # --- 0) Validate required columns in the input data frame ----------------
  # The function assumes that point estimates and their 95% CIs are present.
  req_cols <- c("model",
                "AUC","AUC_lower","AUC_upper",
                "Sens","Sens_lower","Sens_upper",
                "Spec","Spec_lower","Spec_upper")
  missing <- setdiff(req_cols, names(all_metrics))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }
  
  # --- 1) Recode model names and ensure factor ordering -------------------
  # Convert algorithm identifiers to human-readable labels where applicable.
  df <- all_metrics %>%
    dplyr::mutate(
      model = as.character(.data$model),
      model = dplyr::if_else(.data$model %in% names(model_name_map),
                             unname(model_name_map[.data$model]),
                             .data$model),
      model = factor(.data$model, levels = model_order)
    ) %>% dplyr::arrange(model)
  
  # --- 2) Reshape into long format for plotting ---------------------------
  # Stack AUC, Sens, and Spec into a single "metric" column with estimates and CIs.
  metrics_long <- dplyr::bind_rows(
    df %>% dplyr::transmute(model, metric = "AUC",  est = .data$AUC,  lower = .data$AUC_lower,  upper = .data$AUC_upper),
    df %>% dplyr::transmute(model, metric = "Sens", est = .data$Sens, lower = .data$Sens_lower, upper = .data$Sens_upper),
    df %>% dplyr::transmute(model, metric = "Spec", est = .data$Spec, lower = .data$Spec_lower, upper = .data$Spec_upper)
  ) %>%
    dplyr::filter(metric %in% metrics_to_show) %>%
    dplyr::mutate(metric = factor(metric, levels = metrics_to_show))
  
  # --- 3) Create a forestplot-style text table ----------------------------
  # Format each metric as "estimate (lower–upper)" and pivot to a wide layout.
  tbl_text <- metrics_long %>%
    dplyr::mutate(ci = sprintf("%.2f (%.2f–%.2f)", est, lower, upper)) %>%
    dplyr::select(model, metric, ci) %>%
    tidyr::pivot_wider(names_from = metric, values_from = ci) %>% 
    dplyr::mutate(
      Model = factor(as.character(model), levels = model_order)
    ) %>%
    dplyr::arrange(desc(Model)) %>%          
    dplyr::select(-model) %>%
    dplyr::mutate(Model = as.character(Model)) %>%   # for rbind
    dplyr::select(Model, everything())
    
  # Prepend a header row (character matrix) for use with forestplot-like APIs.
  table_text <- rbind(
    c("Model",
      if ("AUC"  %in% metrics_to_show) "AUC (95% CI)"  else NULL,
      if ("Sens" %in% metrics_to_show) "Sens (95% CI)" else NULL,
      if ("Spec" %in% metrics_to_show) "Spec (95% CI)" else NULL),
    as.matrix(tbl_text)
  )
  
  # --- 4) Prepare reference (dashed) lines per metric ---------------------
  # These serve as visual thresholds (e.g., 0.5 for AUC/Sens/Spec).
  vline_df <- tibble::tibble(
    metric = names(dashed_at),
    vline  = as.numeric(dashed_at)
  ) %>% dplyr::filter(metric %in% metrics_to_show)
  
  # --- 5) Build the visualization ----------------------------------------
  # Points represent estimates; horizontal error bars represent 95% CIs.
  p <- ggplot2::ggplot(metrics_long, ggplot2::aes(x = est, y = model, color = model)) +
    geom_point(size = 3) +
    geom_errorbar(ggplot2::aes(xmin = lower, xmax = upper), height = 0.2) +
    facet_wrap(~ metric, nrow = 1) +
    geom_text(ggplot2::aes(label = sprintf("%.2f", est)),
                       hjust = 0.5, vjust = 0.5, nudge_x = 0.0, nudge_y = 0.25,
                       size = 5.0, show.legend = FALSE) +
    geom_vline(data = vline_df,
                        ggplot2::aes(xintercept = vline),
                        linetype = "dashed", color = "grey70") +
    scale_x_continuous(breaks = c(0.5, 0.75, 1.0), limits = c(0.6, 1.0),labels = c("0.5", "0.75", "1.0"),
                       )+
    labs(x = NULL, y = NULL, color = "Model") +
    theme_bw(base_size = 18) +
    theme(legend.position = "none") 

  if (!is.null(colormap)) {
    p <- p + ggplot2::scale_color_manual(values = colormap)
  } else {
    p <- p + ggplot2::scale_color_brewer(palette = "Set2")
  }
  
  # --- 6) Return all useful objects ---------------------------------------
  # Returning multiple objects allows downstream customization if needed.
  list(
    plot = p,
    table_text = table_text,
    metrics_long = metrics_long,
    tbl_text = tbl_text
  )
}




# Nested CV ------------------------------------------------------------


ci_t_per_group <- function(x) {
  n <- length(x); m <- mean(x); s <- sd(x)
  tcrit <- qt(0.975, df = n - 1)
  c(mean = m, lower = m - tcrit * s / sqrt(n), upper = m + tcrit * s / sqrt(n))
}
clamp01 <- function(x) pmax(0, pmin(1, x))


run_on_outer_fold_metrics <- function(test_idx, fold_id, rep_id) {
  train_idx  <- setdiff(seq_len(nrow(df)), test_idx)
  train_data <- df[train_idx, ]
  test_data  <- df[test_idx,  ]
  
  purrr::imap_dfr(models, function(spec, model_name) {
    message(sprintf(">> METRICS %s / %s / %s", rep_id, fold_id, model_name))
    fit <- fit_model(train_data, test_data, "ALS_label", spec$preds, spec$method)
    met <- extract_metrics(fit$roc, cutoff = fit$th_star, calcCI = FALSE)$metrics
    dplyr::mutate(met, model = model_name, outer_fold = fold_id, iter = rep_id)
  }, .progress = TRUE)
}


# Compute metrics -------------------------------------------------------------------

# ---- Compute metrics once with a controllable threshold rule ----
# threshold_rule:
#   - "youden"      : optimize Youden index on the given data (default)
#   - "fixed"       : use a fixed threshold passed via 'target_value' (numeric)
#   - "target_spec" : choose threshold achieving closest specificity to 'target_value' (0-1)
#   - "target_sens" : choose threshold achieving closest sensitivity to 'target_value' (0-1)
compute_metrics_once <- function(df,
                                 neg_class = "nonALS", pos_class = "ALS",
                                 threshold_rule = c("youden","fixed","target_spec","target_sens"),
                                 target_value = NULL) {
  threshold_rule <- match.arg(threshold_rule)
  df <- tibble::as_tibble(df)
  df$y_true  <- factor(df$y_true, levels = c(neg_class, pos_class))
  df$y_score <- as.numeric(df$y_score)
  
  roc_obj <- pROC::roc(
    response  = df$y_true,
    predictor = df$y_score,
    levels    = c(neg_class, pos_class),
    direction = "<"
  )
  auc_val <- as.numeric(pROC::auc(roc_obj))
  
  # --- decide threshold according to the rule ---
  decide_threshold <- function() {
    if (threshold_rule == "youden") {
      best <- pROC::coords(
        roc_obj, x = "best", best.method = "youden",
        ret = c("threshold","sensitivity","specificity"), transpose = FALSE
      )
      best <- as.data.frame(best)
      if (nrow(best) > 1L) {
        best$balance <- abs(best$sensitivity - best$specificity)
        best <- best[order(best$balance), , drop = FALSE]
        if (nrow(best) > 1L && best$balance[1] == best$balance[2]) {
          med_score <- stats::median(df$y_score, na.rm = TRUE)
          d <- abs(best$threshold - med_score)
          best <- best[order(d), , drop = FALSE]
        }
      }
      return(as.numeric(best$threshold[1]))
    } else if (threshold_rule == "fixed") {
      if (is.null(target_value) || !is.numeric(target_value) || length(target_value) != 1L) {
        stop("For threshold_rule='fixed', provide numeric target_value as the threshold.")
      }
      return(as.numeric(target_value))
    } else if (threshold_rule == "target_spec") {
      if (is.null(target_value)) stop("Provide target_value in [0,1] for 'target_spec'.")
      cr <- pROC::coords(
        roc_obj, x = target_value, input = "specificity",
        ret = c("threshold","sensitivity","specificity"), transpose = FALSE
      )
      return(as.numeric(cr["threshold"]))
    } else if (threshold_rule == "target_sens") {
      if (is.null(target_value)) stop("Provide target_value in [0,1] for 'target_sens'.")
      cr <- pROC::coords(
        roc_obj, x = target_value, input = "sensitivity",
        ret = c("threshold","sensitivity","specificity"), transpose = FALSE
      )
      return(as.numeric(cr["threshold"]))
    }
  }
  
  th <- decide_threshold()
  
  # Sens/Spec at that chosen threshold
  ss <- pROC::coords(roc_obj, x = th, input = "threshold",
                     ret = c("sensitivity","specificity"), transpose = FALSE)
  sen <- as.numeric(ss["sensitivity"])
  spe <- as.numeric(ss["specificity"])
  
  pred_pos <- df$y_score >= th
  tp <- sum(df$y_true == pos_class & pred_pos)
  tn <- sum(df$y_true == neg_class & !pred_pos)
  fp <- sum(df$y_true == neg_class & pred_pos)
  fn <- sum(df$y_true == pos_class & !pred_pos)
  acc <- (tp + tn) / nrow(df)
  
  tibble::tibble(
    AUC        = auc_val,
    optCut_off = th,
    Sens       = sen,
    Spec       = spe,
    accuracy   = acc,
    TP = tp, TN = tn, FP = fp, FN = fn
  )
}


# ---- Patient-level bootstrap (aligned with the same threshold rule) ----
bootstrap_metrics <- function(df, B = 2000, stratified = TRUE, seed = 12345,
                              neg_class = "nonALS", pos_class = "ALS",
                              show_progress = TRUE,
                              threshold_rule = c("youden","fixed","target_spec","target_sens"),
                              target_value = NULL) {
  threshold_rule <- match.arg(threshold_rule)
  
  set.seed(seed)
  df <- tibble::as_tibble(df)
  df$y_true  <- factor(df$y_true, levels = c(neg_class, pos_class))
  df$y_score <- as.numeric(df$y_score)
  
  ids_all <- unique(df$row_id)
  ids_pos <- unique(df$row_id[df$y_true == pos_class])
  ids_neg <- unique(df$row_id[df$y_true == neg_class])
  n_pos   <- length(ids_pos); n_neg <- length(ids_neg)
  
  if (show_progress) {
    pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
  }
  
  out_list <- vector("list", B)
  for (b in seq_len(B)) {
    if (stratified) {
      samp_ids <- c(sample(ids_pos, size = n_pos, replace = TRUE),
                    sample(ids_neg, size = n_neg, replace = TRUE))
    } else {
      samp_ids <- sample(ids_all, size = length(ids_all), replace = TRUE)
    }
    
    #df_b <- dplyr::semi_join(df, tibble::tibble(row_id = samp_ids), by = "row_id")
    
    df_b <- tibble::tibble(
      row_id = samp_ids,
      .boot_id = seq_along(samp_ids)
    ) %>%
      dplyr::left_join(df, by = "row_id")
    
    # NOTE: if threshold_rule == "fixed" and target_value is NULL, we should pass
    # the full-sample fixed threshold. Handle that upstream (see summarise_with_boot).
    out_list[[b]] <- compute_metrics_once(
      df_b, neg_class = neg_class, pos_class = pos_class,
      threshold_rule = threshold_rule, target_value = target_value
    )
    
    if (show_progress) utils::setTxtProgressBar(pb, b)
  }
  
  boot_stats <- dplyr::bind_rows(out_list)
  
  list(
    AUC_ci  = stats::quantile(boot_stats$AUC,        probs = c(0.025, 0.975), na.rm = TRUE),
    Sens_ci = stats::quantile(boot_stats$Sens,       probs = c(0.025, 0.975), na.rm = TRUE),
    Spec_ci = stats::quantile(boot_stats$Spec,       probs = c(0.025, 0.975), na.rm = TRUE),
    Cut_ci  = stats::quantile(boot_stats$optCut_off, probs = c(0.025, 0.975), na.rm = TRUE)
  )
}


# ---- Wrapper that keeps point estimate & CI on the *same* rule ----
summarise_with_boot <- function(df_group, B = 2000, stratified = TRUE, seed = 12345,
                                neg_class = "nonALS", pos_class = "ALS",
                                show_progress = TRUE,
                                threshold_rule = c("youden","fixed","target_spec","target_sens"),
                                target_value = NULL) {
  threshold_rule <- match.arg(threshold_rule)
  
  # If using fixed threshold but value is not provided, compute it once on full data
  fixed_th <- NULL
  if (threshold_rule == "fixed") {
    if (is.null(target_value)) {
      # derive the "fixed" threshold from full data via Youden by default
      tmp <- compute_metrics_once(df_group, neg_class, pos_class, threshold_rule = "youden")
      fixed_th <- tmp$optCut_off
    } else {
      fixed_th <- target_value
    }
  }
  
  point <- compute_metrics_once(
    df_group, neg_class = neg_class, pos_class = pos_class,
    threshold_rule = threshold_rule,
    target_value   = if (threshold_rule == "fixed") fixed_th else target_value
  )
  
  ci <- bootstrap_metrics(
    df_group, B = B, stratified = stratified, seed = seed,
    neg_class = neg_class, pos_class = pos_class,
    show_progress = show_progress,
    threshold_rule = threshold_rule,
    target_value   = if (threshold_rule == "fixed") fixed_th else target_value
  )
  
  dplyr::mutate(
    point,
    AUC_lower    = unname(ci$AUC_ci[1]),
    AUC_upper    = unname(ci$AUC_ci[2]),
    Sens_lower   = unname(ci$Sens_ci[1]),
    Sens_upper   = unname(ci$Sens_ci[2]),
    Spec_lower   = unname(ci$Spec_ci[1]),
    Spec_upper   = unname(ci$Spec_ci[2]),
    optCut_lower = unname(ci$Cut_ci[1]),
    optCut_upper = unname(ci$Cut_ci[2])
  )
}


# Evaluate all (model × iter) groups with a top-level progress bar,
# and show an inner progress bar for the bootstrap inside each group.
evaluate_rep_level_with_progress <- function(oof_preds,
                                             B = 2000, stratified = TRUE, seed = 12345,
                                             neg_class = "nonALS", pos_class = "ALS",
                                             threshold_rule = "youden") {
  # Identify unique groups to iterate over.
  groups <- oof_preds %>%
    dplyr::distinct(model, iter) %>%
    dplyr::arrange(model, iter)
  
  # Top-level progress bar across groups.
  pb <- utils::txtProgressBar(min = 0, max = nrow(groups), style = 3)
  on.exit(try(close(pb), silent = TRUE), add = TRUE)
  
  res_list <- vector("list", nrow(groups))
  
  for (i in seq_len(nrow(groups))) {
    g  <- groups[i, ]
    df <- oof_preds %>%
      dplyr::filter(model == g$model, iter == g$iter) %>%
      dplyr::select(row_id, y_true, y_score)
    
    # Inner bootstrap will also show its own progress bar.
    out <- summarise_with_boot(
      df,
      B = B, stratified = stratified, seed = seed,
      neg_class = neg_class, pos_class = pos_class,threshold_rule = threshold_rule,
      show_progress = TRUE
    ) %>%
      dplyr::mutate(model = g$model, iter = g$iter, .before = dplyr::everything())
    
    res_list[[i]] <- out
    utils::setTxtProgressBar(pb, i)
  }
  
  dplyr::bind_rows(res_list)
}


# Evaluate the repeat-averaged ensemble per model, with a progress bar per model
# and an inner progress bar for the bootstrap.
evaluate_ensemble_with_progress <- function(oof_preds,
                                            B = 2000, stratified = TRUE, seed = 20240909,
                                            neg_class = "nonALS", pos_class = "ALS",
                                            threshold_rule = "youden") {
  # Prepare the per-patient repeat-averaged score per model.
  y_by_id <- oof_preds %>% dplyr::distinct(row_id, y_true)
  oof_mean <- oof_preds %>%
    dplyr::group_by(model, row_id) %>%
    dplyr::summarise(y_score = mean(as.numeric(y_score), na.rm = TRUE), .groups = "drop") %>%
    dplyr::left_join(y_by_id, by = "row_id")
  
  models <- oof_mean %>%
    dplyr::distinct(model) %>%
    dplyr::arrange(model)
  
  # Top-level progress bar across models.
  pb <- utils::txtProgressBar(min = 0, max = nrow(models), style = 3)
  on.exit(try(close(pb), silent = TRUE), add = TRUE)
  
  res_list <- vector("list", nrow(models))
  
  for (i in seq_len(nrow(models))) {
    m  <- models$model[i]
    df <- oof_mean %>%
      dplyr::filter(model == m) %>%
      dplyr::select(row_id, y_true, y_score)
    
    # Inner bootstrap will also show its own progress bar.
    out <- summarise_with_boot(
      df,
      B = B, stratified = stratified, seed = seed,
      neg_class = neg_class, pos_class = pos_class,threshold_rule = threshold_rule,
      show_progress = TRUE
    ) %>%
      dplyr::mutate(model = m, .before = dplyr::everything())
    
    res_list[[i]] <- out
    utils::setTxtProgressBar(pb, i)
  }
  
  dplyr::bind_rows(res_list)
}



# ROC for each rep --------------------------------------------------------

# Plot overlaid ROC curves for all models within a given repeat.
# Returns a list: $plot (ggplot), $auc_table (per-model AUC), $roc_list (pROC roc objects)
plot_repeat_roc <- function(oof_preds, rep_id,
                            neg_class = "nonALS", pos_class = "ALS",
                            smooth = FALSE, legacy_axes = TRUE, line_size = 0.9,
                            colormap = NULL) {
  model_labels <- c(
    "NfLonly" = "NfL alone"
  )
  
  df_rep <- oof_preds %>%
    dplyr::filter(iter == rep_id) %>%
    dplyr::mutate(
      y_true  = factor(y_true, levels = c(neg_class, pos_class)),
      y_score = as.numeric(y_score)
    )
  
  models <- sort(unique(df_rep$model))
  roc_list <- list()
  auc_df <- dplyr::tibble()
  
  for (m in models) {
    df_m <- df_rep %>% dplyr::filter(model == m)
    if (length(unique(df_m$y_true)) < 2) next
    
    r <- pROC::roc(response = df_m$y_true, predictor = df_m$y_score,
                   levels = c(neg_class, pos_class), direction = "<")
    if (isTRUE(smooth)) r <- pROC::smooth.roc(r)
    auc_val <- as.numeric(pROC::auc(r))
    
    pretty_m <- ifelse(m %in% names(model_labels), model_labels[[m]], m)
    nm <- sprintf("%s (AUC=%.2f)", pretty_m, auc_val)
    
    roc_list[[nm]] <- r
    auc_df <- dplyr::bind_rows(auc_df, tibble::tibble(model = pretty_m, AUC = auc_val, label = nm))
  }
  
  auc_df <- auc_df %>% arrange(desc(AUC))
  ordered_labels <- auc_df$label
  
  gp <- pROC::ggroc(roc_list, legacy.axes = legacy_axes, size = line_size) +
    ggplot2::labs(
      title = sprintf("ROC curves — %s", rep_id),
      x = if (legacy_axes) "1 - Specificity" else "False Positive Rate",
      y = "Sensitivity",
      color = "Model"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = element_text(face = "bold")) 

  
    if (!is.null(colormap)) {
      base_names <- sub(" \\(AUC=.*\\)$", "", ordered_labels)   # "GLM (AUC=...)" -> "GLM"
      values_vec <- unname(colormap[base_names])                # color vector
      # fallback
      if (any(is.na(values_vec))) {
        warning("Some models missing in colormap; using default palette for those.")
        gp <- gp + ggplot2::scale_color_manual(breaks = ordered_labels, values = values_vec, na.translate = FALSE)
      } else {
        gp <- gp + ggplot2::scale_color_manual(breaks = ordered_labels, values = values_vec)
      }
    } else {
      gp <- gp + ggplot2::scale_color_brewer(palette = "Set2")
    }
  
  list(plot = gp, auc_table = auc_df, roc_list = roc_list)
}



# Plot overlaid ROC curves for the repeat-averaged ensemble per model.

# ---- Pretty ensemble ROC (mean over repeats) ----
plot_ensemble_roc <- function(
    oof_mean,
    colormap = NULL,          # option
    neg_class = "nonALS", pos_class = "ALS",
    smooth = FALSE, legacy_axes = TRUE, line_size = 0.9,
    title = "Ensemble ROC (mean over repeats)",
    model_labels = c(                      # option
      "NfLonly" = "NfL alone"
    ),
    show_auc_in_legend = TRUE              # TRUE
) {
  models  <- sort(unique(oof_mean$model))
  roc_list <- list()
  auc_df   <- dplyr::tibble()
  
  for (m in models) {
    df_m <- oof_mean %>% dplyr::filter(model == m)
    if (length(unique(df_m$y_true)) < 2) next
    
    r <- pROC::roc(
      response  = factor(df_m$y_true, levels = c(neg_class, pos_class)),
      predictor = as.numeric(df_m$y_score),
      direction = "<"
    )
    if (isTRUE(smooth)) r <- pROC::smooth.roc(r)
    
    auc_val  <- as.numeric(pROC::auc(r))
    pretty_m <- if (m %in% names(model_labels)) model_labels[[m]] else m
    label    <- if (show_auc_in_legend) sprintf("%s (AUC=%.2f)", pretty_m, auc_val) else pretty_m
    
    roc_list[[label]] <- r
    auc_df <- dplyr::bind_rows(
      auc_df,
      tibble::tibble(model = m, pretty_model = pretty_m, AUC = auc_val, label = label)
    )
  }
  
  if (length(roc_list) == 0L) stop("No valid ROC curves to plot for ensemble.")
  
  # Legend order = AUC descending
  auc_df <- auc_df %>% dplyr::arrange(dplyr::desc(AUC))
  ordered_labels <- auc_df$label
  
  gp <- pROC::ggroc(roc_list, legacy.axes = legacy_axes, size = line_size) +
    ggplot2::labs(
      title = title,
      x = if (legacy_axes) "1 - Specificity" else "False Positive Rate",
      y = "Sensitivity",
      color = "Model"
    ) +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold")) 
  
  if (!is.null(colormap)) {
    base_names <- sub(" \\(AUC=.*\\)$", "", ordered_labels)
    values_vec <- unname(colormap[base_names])
    if (any(is.na(values_vec))) {
      warning("Some models missing in colormap; using default palette for those.")
      gp <- gp + ggplot2::scale_color_manual(breaks = ordered_labels, values = values_vec, na.translate = FALSE)
    } else {
      gp <- gp + ggplot2::scale_color_manual(breaks = ordered_labels, values = values_vec)
    }
  } else {
    gp <- gp + ggplot2::scale_color_brewer(palette = "Set2")
  }
  
  list(plot = gp, auc_table = auc_df, roc_list = roc_list)
}


# ---- Choose representative repeat across ALL models ----
# oof_summary: tibble with columns [model, iter, AUC] 
choose_rep_for_all_models <- function(
    oof_summary,
    metric = "AUC",
    tie_break = c("closest_to_median_then_min_iter", "closest_to_median_then_random")
) {
  tie_break <- match.arg(tie_break)
  
  rep_stats <- oof_summary %>%
    dplyr::group_by(iter) %>%
    dplyr::summarise(mean_metric = mean(.data[[metric]], na.rm = TRUE), .groups = "drop")
  
  med_val <- stats::median(rep_stats$mean_metric, na.rm = TRUE)
  
  rep_stats <- rep_stats %>%
    dplyr::mutate(dist = abs(mean_metric - med_val)) %>%
    dplyr::arrange(dist)
  
  # tie break
  top_dist <- rep_stats$dist[1]
  tied <- rep_stats %>% dplyr::filter(dist == top_dist)
  
  chosen_iter <- if (nrow(tied) == 1) {
    tied$iter[1]
  } else if (tie_break == "closest_to_median_then_min_iter") {
    # iter = "Rep1","Rep2" => choose minimum numeric part
    parse_num <- function(x) suppressWarnings(as.numeric(gsub("[^0-9]+", "", x)))
    tied$iter[order(ifelse(is.na(parse_num(tied$iter)), Inf, parse_num(tied$iter)),
                    tied$iter)][1]
  } else { # "closest_to_median_then_random"
    sample(tied$iter, 1)
  }
  
  chosen_iter
}


choose_rep_for_all_models_2 <- function(
    oof_summary,
    metric = "AUC",
    baseline_model = "NfLonly",
    baseline_window = c("within_eps", "within_iqr", "within_quantiles"),
    eps = 0.02,
    q = c(0.25, 0.75),
    tie_break = c("closest_to_median_then_min_iter", "closest_to_median_then_random"),
    return_details = FALSE,
    seed = 12345
) {
  baseline_window <- match.arg(baseline_window)
  tie_break       <- match.arg(tie_break)
  
  parse_num <- function(x) suppressWarnings(as.numeric(gsub("[^0-9]+", "", x)))
  
  # ────────────────────────────────────────────────
  # Stage 0: basic checks
  # ────────────────────────────────────────────────
  if (!all(c("model", "iter", metric) %in% names(oof_summary))) {
    stop("oof_summary must contain columns: model, iter, and the specified metric.")
  }
  if (!baseline_model %in% unique(oof_summary$model)) {
    stop("baseline_model not found in oof_summary$model.")
  }
  
  # ────────────────────────────────────────────────
  # Stage 1: identify eligible repeats where baseline model is 'typical'
  # ────────────────────────────────────────────────
  base_tbl <- oof_summary %>%
    dplyr::filter(.data$model == baseline_model) %>%
    dplyr::select(iter, base_metric = dplyr::all_of(metric))
  
  base_med <- stats::median(base_tbl$base_metric, na.rm = TRUE)
  
  eligible_iters <- switch(
    baseline_window,
    within_eps = {
      base_tbl %>%
        dplyr::mutate(dist_base = abs(base_metric - base_med)) %>%
        dplyr::filter(dist_base <= eps) %>%
        dplyr::pull(iter)
    },
    within_iqr = {
      qs <- stats::quantile(base_tbl$base_metric, probs = c(0.25, 0.75), na.rm = TRUE)
      base_tbl %>%
        dplyr::filter(base_metric >= qs[1], base_metric <= qs[2]) %>%
        dplyr::pull(iter)
    },
    within_quantiles = {
      qs <- stats::quantile(base_tbl$base_metric, probs = q, na.rm = TRUE)
      base_tbl %>%
        dplyr::filter(base_metric >= qs[1], base_metric <= qs[2]) %>%
        dplyr::pull(iter)
    }
  )
  
  # Fallback: if filter is too strict, use all iters (but keep baseline distances for transparency)
  if (length(eligible_iters) == 0) {
    eligible_iters <- unique(oof_summary$iter)
  }
  
  # ────────────────────────────────────────────────
  # Stage 2: among eligible repeats, choose the one closest to the median mean performance
  # ────────────────────────────────────────────────
  rep_stats <- oof_summary %>%
    dplyr::filter(.data$iter %in% eligible_iters) %>%
    dplyr::group_by(iter) %>%
    dplyr::summarise(mean_metric = mean(.data[[metric]], na.rm = TRUE), .groups = "drop")
  
  med_val <- stats::median(rep_stats$mean_metric, na.rm = TRUE)
  
  rep_stats <- rep_stats %>%
    dplyr::mutate(dist_mean = abs(mean_metric - med_val)) %>%
    dplyr::arrange(dist_mean)
  
  top_dist <- rep_stats$dist_mean[1]
  tied <- rep_stats %>% dplyr::filter(dist_mean == top_dist)
  
  chosen_iter <- if (nrow(tied) == 1) {
    tied$iter[1]
  } else if (tie_break == "closest_to_median_then_min_iter") {
    tied$iter[order(
      ifelse(is.na(parse_num(tied$iter)), Inf, parse_num(tied$iter)),
      tied$iter
    )][1]
  } else {
    set.seed(seed)
    sample(tied$iter, 1)
  }
  
  if (!isTRUE(return_details)) return(chosen_iter)
  
  # ────────────────────────────────────────────────
  # return diagnostics
  # ────────────────────────────────────────────────
  base_diag <- base_tbl %>%
    dplyr::mutate(
      base_median = base_med,
      dist_base = abs(base_metric - base_med),
      eligible = iter %in% eligible_iters
    )
  
  rep_diag <- rep_stats %>%
    dplyr::mutate(
      overall_median_mean = med_val,
      chosen = iter == chosen_iter
    )
  
  list(
    chosen_iter = chosen_iter,
    baseline_model = baseline_model,
    baseline_window = baseline_window,
    eps = eps,
    quantiles = q,
    baseline_diagnostics = base_diag,
    repeat_diagnostics = rep_diag
  )
}



# Main --------------------------------------------------------------------


# ————————————————————————————————————————————————————————————————
# RUN one outer fold once per model
#   - Keeps the original OOF structure identical
#   - ALSO returns a lightweight handle for SHAP later (no retraining)
# Dependencies: uses global df, models, pos_class, fit_model
# ————————————————————————————————————————————————————————————————
run_on_outer_fold_oof_and_handle <- function(test_idx, fold_id, rep_id, outcome = "ALS_label") {
  train_idx  <- setdiff(seq_len(nrow(df)), test_idx)
  train_data <- df[train_idx, , drop=FALSE]
  test_data  <- df[test_idx , , drop=FALSE]
  
  train_data <- as.data.frame(train_data)
  test_data  <- as.data.frame(test_data)
  
  # iterate models exactly as before
  purrr::imap(models, function(spec, model_name) {
    message(sprintf(">> OOF %s / %s / %s", rep_id, fold_id, model_name))
    
    # normalize predictor names (character vector)
    preds <- as.character(unlist(spec$preds, use.names = FALSE))
    
    # fit once (inner-CV happens inside fit_model)
    #fit <- fit_model(train_data, test_data, outcome, preds, spec$method)
    
    fit <- tryCatch(
      fit_model(train_data, test_data, outcome, preds, spec$method),
      error = function(e) {
        # Error handling: save debug info to RDS
        saveRDS(
          list(
            msg = conditionMessage(e),
            rep_id = rep_id, fold_id = fold_id, model_name = model_name, method = spec$method,
            preds = preds, preds_class = class(preds), preds_typeof = typeof(preds),
            # train/test indices
            train_col_class = sapply(train_data[, preds, drop = FALSE], function(x) paste(class(x), collapse="/")),
            train_col_typeof= sapply(train_data[, preds, drop = FALSE], typeof),
            test_col_class  = sapply(test_data [, preds, drop = FALSE], function(x) paste(class(x), collapse="/")),
            test_col_typeof = sapply(test_data [, preds, drop = FALSE], typeof)
          ),
          file = sprintf("debug_fit_fail_%s_%s_%s.rds", rep_id, fold_id, model_name)
        )
        stop(e)
      }
    )
    
    # ---- Safe extraction of test_prob ----
    test_prob_raw <- fit$test_prob
    
    # Ensure it's a plain numeric vector
    if (is.list(test_prob_raw) && !is.data.frame(test_prob_raw)) {
      test_prob_raw <- unlist(test_prob_raw, recursive = TRUE, use.names = FALSE)
    }
    if (is.matrix(test_prob_raw)) {
      test_prob_raw <- as.vector(test_prob_raw)
    }
    if (is.data.frame(test_prob_raw)) {
      test_prob_raw <- test_prob_raw[[1]]
    }
    
    # Final conversion
    test_prob <- as.numeric(test_prob_raw)
    storage.mode(test_prob) <- "double"
    
    # Verify
    if (!is.numeric(test_prob) || !is.atomic(test_prob) || length(test_prob) != nrow(test_data)) {
      stop(sprintf("test_prob has invalid structure: class=%s, length=%d, expected=%d", 
                   class(test_prob), length(test_prob), nrow(test_data)))
    }
    
    
    # ---- (A) OOF tibble: identical columns as the current pipeline ----
    oof_tbl <- tibble::tibble(
      row_id     = test_idx,                            # original df row-index
      iter       = rep_id,                              # e.g., "Rep1"
      outer_fold = fold_id,                             # e.g., "Fold1"
      model      = model_name,                          # "RF","XGB","KNN","SVM","GLM","NfLonly"
      y_true     = factor(test_data[[outcome]], levels = c(paste0("non", pos_class), pos_class)),
      y_score    = as.numeric(fit$test_prob)            # predicted prob on outer test
    )
    
    # ---- (B) Per-fold handle for SHAP (small; no SHAP is computed here) ----
    # Note: we do NOT store train/test data here to keep it light.
    handle <- list(
      model      = fit$model,                           # trained model object
      method     = spec$method,                         # "rf","xgb","svmRadial","kknn","one",...
      predictors = preds,                               # the exact feature set used
      pos_class  = pos_class,
      train_idx  = train_idx,                           # indices only (to rebuild bg later if needed)
      test_idx   = test_idx,                            # indices only (to know target rows)
      rep_id     = rep_id,
      fold_id    = fold_id,
      model_name = model_name
    )
    
    # ---- (C) Per-repeat scaffolding: store per fold per model objects ----
    # This is what "per_repeat" will keep. It mirrors: repeats -> folds -> models
    per_fold_unit <- list(
      model_name = model_name,
      handle     = handle,
      # OPTIONAL: stash extras that make later SHAP faster/reproducible.
      # Keeping them NULL keeps memory footprint low; you can add if you like.
      # e.g., "train_X_for_shap", "preprocess_params", etc.
      extras     = NULL
    )
    
    list(oof = oof_tbl, unit = per_fold_unit, handle = handle)
  })
}


# ————————————————————————————————————————————————————————————————
# Driver that runs ALL repeats/folds once and returns:
#   - oof_preds    : same tibble you already use downstream
#   - per_repeat   : nested list [RepX][FoldY][[ per-model units ]]
#   - shap_handles : flat list of handles for bulk SHAP (optional)
# ————————————————————————————————————————————————————————————————


run_all_oof_collect <- function(outer_folds_list, outcome = "ALS_label") {
  
  cat("=== Processing", length(outer_folds_list), "repeats sequentially ===\n\n")
  
  all_results <- vector("list", length(outer_folds_list))
  
  for (r in seq_along(outer_folds_list)) {
    rep_label <- paste0("Rep", r)
    cat("--- Processing", rep_label, "---\n")
    
    # CRITICAL: Clean up state before each repeat
    if (r > 1) {
      # Force garbage collection
      gc(verbose = FALSE)
      
      # Clear caret's internal cache (if exists)
      if (exists(".Random.seed", envir = .GlobalEnv)) {
        rm(.Random.seed, envir = .GlobalEnv)
      }
      
      # Reset seed
      set.seed(.seed_value + r)
      
      # Small delay to allow R to stabilize
      Sys.sleep(0.1)
    }
    
    folds_one_rep <- outer_folds_list[[r]]
    
    rep_results <- tryCatch({
      purrr::imap(folds_one_rep, function(test_idx, fold_id) {
        fold_label <- as.character(fold_id)
        run_on_outer_fold_oof_and_handle(
          test_idx = test_idx, 
          fold_id = fold_label, 
          rep_id = rep_label, 
          outcome = outcome
        )
      })
    }, error = function(e) {
      cat("ERROR in", rep_label, ":", conditionMessage(e), "\n")
      NULL
    })
    
    if (!is.null(rep_results)) {
      all_results[[r]] <- rep_results
      
      # Count predictions
      n_pred <- sum(sapply(rep_results, function(fold) {
        sum(sapply(fold, function(model) nrow(model$oof)))
      }))
      
      # Count SHAP handles
      n_shap <- sum(sapply(rep_results, function(fold) length(fold)))
      
      cat("Completed", rep_label, ":", n_pred, "OOF predictions,", n_shap, "SHAP handles\n\n")
    } else {
      cat("Completed", rep_label, ": 0 OOF predictions, 0 SHAP handles\n\n")
    }
  }
  
  # Flatten results
  oof_list <- list()
  per_repeat <- list()
  flat_handles <- list()
  
  for (r in seq_along(all_results)) {
    if (is.null(all_results[[r]])) next
    
    rep_label <- paste0("Rep", r)
    per_repeat[[rep_label]] <- list()
    
    for (k in seq_along(all_results[[r]])) {
      fold_label <- names(all_results[[r]])[k] %||% as.character(k)
      per_repeat[[rep_label]][[fold_label]] <- list()
      
      per_model <- all_results[[r]][[k]]
      for (m in seq_along(per_model)) {
        oof_list[[length(oof_list) + 1L]] <- per_model[[m]]$oof
        per_repeat[[rep_label]][[fold_label]][[m]] <- per_model[[m]]$unit
        flat_handles[[length(flat_handles) + 1L]] <- per_model[[m]]$handle
      }
    }
  }
  
  list(
    oof_preds = dplyr::bind_rows(oof_list),
    per_repeat = per_repeat,
    shap_handles = flat_handles
  )
}


# nested SHAP --------------------------------------------------------------------


# Small helper: list available repeats per model (for sanity checks / UI)
list_available_repeats <- function(shap_handles, model) {
  idx <- vctrs::vec_detect(vctrs::field(shap_handles, "model_name") %||% 
                             purrr::map_chr(shap_handles, ~ .x$model_name %||% NA_character_),
                           model)
  reps <- unique(purrr::map_int(shap_handles[idx], ~ .x$rep_id))
  sort(reps)
}

# Helpers (kept minimal and inline) ---------------

# Normalize method names to a canonical key (vectorized, NA-safe)
.normalize_method <- function(m) {
  m <- tolower(as.character(m))
  m <- gsub("\\s+", "", m)
  dplyr::case_when(
    m %in% c("svmradial","svm")          ~ "svm",
    m %in% c("randomforest","rf")        ~ "rf",
    m %in% c("xgboost","xgb","gbtree")   ~ "xgb",
    m %in% c("glmnet","glm")             ~ "glm",
    m %in% c("kknn","knn")               ~ "knn",
    TRUE                                 ~ m
  )
}

# Coerce repeat ids like "Rep1"/"1" → integer 1 (vectorized)
.coerce_rep_int <- function(x) {
  x_chr <- as.character(x)
  suppressWarnings(as.integer(gsub("\\D+", "", x_chr)))
}

# Build a normalized availability table (model/rep/fold) from handles
.inspect_handles_normalized <- function(shap_handles) {
  tibble::tibble(
    idx      = seq_along(shap_handles),
    model    = .normalize_method(purrr::map_chr(shap_handles, ~ .x$model_name %||% .x$method %||% NA_character_)),
    rep      = .coerce_rep_int(purrr::map_chr(shap_handles, ~ .x$rep_id %||% NA_character_)),
    fold     = suppressWarnings(as.integer(gsub("\\D+", "", purrr::map_chr(shap_handles, ~ .x$fold_id %||% NA_character_))))
  )
}

# Core: compute SHAP rows for ONE handle (classification) --------------------------------

# Returns tibble: model, iter, outer_fold, feature, shap_value, feature_value

compute_shap_for_handle_cls <- function(handle,
                                        data = NULL, df = NULL,
                                        # SHAP parameters (defaults are reasonable for a quick approximation)
                                        nsim = 80, bg_n = 200, test_subsample = 200,
                                        pos_class = NULL, ...) {
  # Accept both `data=` and `df=`; prefer `df` if both provided.
  if (is.null(df)) df <- data
  if (is.null(df)) stop("compute_shap_for_handle_cls: either `data` or `df` must be provided.")
  
  # Bind-rows safe empty schema for errors/skip
  empty_schema <- tibble::tibble(
    model = character(), iter = character(), outer_fold = character(),
    feature = character(), shap_value = numeric(), feature_value = numeric()
  )
  
  # Validate required fields in handle
  need <- c("model","method","predictors","train_idx","test_idx","rep_id","fold_id","model_name")
  miss <- setdiff(need, names(handle))
  if (length(miss)) {
    message(sprintf("!! SHAP skipped (missing fields in handle: %s)", paste(miss, collapse = ", ")))
    return(empty_schema)
  }
  
  # Defaults
  if (is.null(pos_class) && !is.null(handle$pos_class)) pos_class <- as.character(handle$pos_class)
  
  # Unpack
  mdl        <- handle$model
  method     <- tolower(as.character(handle$method))
  preds      <- as.character(unlist(handle$predictors, use.names = FALSE))
  train_idx  <- as.integer(handle$train_idx)
  test_idx   <- as.integer(handle$test_idx)
  rep_id     <- as.character(handle$rep_id)
  fold_id    <- as.character(handle$fold_id)
  model_name <- tolower(as.character(handle$model_name))
  
  # Clip indices to range
  n_all <- nrow(df)
  train_idx <- intersect(train_idx, seq_len(n_all))
  test_idx  <- intersect(test_idx,  seq_len(n_all))
  if (!length(train_idx) || !length(test_idx)) {
    message("!! SHAP skipped (empty train/test indices)")
    return(empty_schema)
  }
  
  # Build background/test frames over the same predictor columns
  X_bg <- as.data.frame(df[train_idx, preds, drop = FALSE])
  X_te <- as.data.frame(df[test_idx,  preds, drop = FALSE])
  
  # Optional subsampling for speed
  if (is.finite(test_subsample) && nrow(X_te) > test_subsample) {
    set.seed(12345 + suppressWarnings(as.integer(gsub("\\D","", fold_id))))
    X_te <- X_te[sample.int(nrow(X_te), test_subsample), , drop = FALSE]
  }
  if (nrow(X_bg) > bg_n) {
    set.seed(54321 + suppressWarnings(as.integer(gsub("\\D","", fold_id))))
    X_bg <- X_bg[sample.int(nrow(X_bg), bg_n), , drop = FALSE]
  }
  
  # Positive-class probability wrapper (caret/xgb/generic)
  predict_pos_prob_cls <- function(fit_obj, newdata, method, pos_class) {
    if (inherits(fit_obj, "xgb.Booster")) {
      return(as.numeric(predict(fit_obj, xgboost::xgb.DMatrix(as.matrix(newdata)))))
    }
    if (inherits(fit_obj, "train")) {
      probs <- try(predict(fit_obj, newdata = newdata, type = "prob"), silent = TRUE)
      if (!inherits(probs, "try-error") && is.data.frame(probs) && !is.null(pos_class) &&
          pos_class %in% colnames(probs)) {
        return(as.numeric(probs[[pos_class]]))
      }
      return(as.numeric(predict(fit_obj, newdata = newdata, type = "raw")))
    }
    probs <- try(predict(fit_obj, newdata = newdata, type = "prob"), silent = TRUE)
    if (!inherits(probs, "try-error") && is.data.frame(probs) && !is.null(pos_class) &&
        pos_class %in% colnames(probs)) {
      return(as.numeric(probs[[pos_class]]))
    }
    rawp <- try(predict(fit_obj, newdata = newdata, type = "response"), silent = TRUE)
    if (!inherits(rawp, "try-error")) return(as.numeric(rawp))
    as.numeric(predict(fit_obj, newdata = newdata))
  }
  
  # Compute SHAP (TreeSHAP for xgb; fastshap otherwise)
  out <- tryCatch({
    if (inherits(mdl, "xgb.Booster")) {
      # Exact TreeSHAP
      dX <- xgboost::xgb.DMatrix(as.matrix(X_te))
      colnames(dX) <- colnames(X_te)
      contrib <- predict(mdl, dX, predcontrib = TRUE)
      contrib <- as.data.frame(contrib)
      bias_cols <- grep("BIAS", colnames(contrib), ignore.case = TRUE)
      if (length(bias_cols)) {
        contrib <- contrib[, -bias_cols, drop = FALSE]
      } else if (ncol(contrib) == ncol(X_te) + 1) {
        contrib <- contrib[, seq_len(ncol(X_te)), drop = FALSE]
      }
      colnames(contrib) <- colnames(X_te)
      
      shap_df <- tibble::as_tibble(contrib)
      shap_df$.row_id <- seq_len(nrow(shap_df))
      feat_vals <- tibble::as_tibble(X_te); feat_vals$.row_id <- seq_len(nrow(X_te))
      
      long <- shap_df |>
        tidyr::pivot_longer(-.row_id, names_to = "feature", values_to = "shap_value") |>
        dplyr::left_join(
          feat_vals |>
            tidyr::pivot_longer(-.row_id, names_to = "feature", values_to = "feature_value"),
          by = c(".row_id","feature")
        ) |>
        dplyr::mutate(model = model_name, iter = rep_id, outer_fold = fold_id) |>
        dplyr::select(model, iter, outer_fold, feature, shap_value, feature_value)
      
      long$shap_value    <- as.numeric(long$shap_value)
      long$feature_value <- suppressWarnings(as.numeric(long$feature_value))
      tibble::as_tibble(long)
      
    } else {
      # fastshap (model-agnostic)
      pred_wrapper <- function(object, newdata) {
        predict_pos_prob_cls(object, newdata, method = method, pos_class = pos_class)
      }
      set.seed(12345 + suppressWarnings(as.integer(gsub("\\D","", fold_id))))
      ex <- fastshap::explain(
        object       = mdl,
        X            = X_bg,
        pred_wrapper = pred_wrapper,
        nsim         = nsim,
        newdata      = X_te
      )
      ex <- as.matrix(ex); storage.mode(ex) <- "double"
      shap_df <- tibble::as_tibble(ex)
      shap_df$.row_id <- seq_len(nrow(shap_df))
      feat_vals <- tibble::as_tibble(X_te); feat_vals$.row_id <- seq_len(nrow(X_te))
      
      long <- shap_df |>
        tidyr::pivot_longer(-.row_id, names_to = "feature", values_to = "shap_value") |>
        dplyr::left_join(
          feat_vals |>
            tidyr::pivot_longer(-.row_id, names_to = "feature", values_to = "feature_value"),
          by = c(".row_id","feature")
        ) |>
        dplyr::mutate(model = model_name, iter = rep_id, outer_fold = fold_id) |>
        dplyr::select(model, iter, outer_fold, feature, shap_value, feature_value)
      
      long$shap_value    <- as.numeric(long$shap_value)
      long$feature_value <- suppressWarnings(as.numeric(long$feature_value))
      tibble::as_tibble(long)
    }
  }, error = function(e) {
    message(sprintf("!! SHAP failed (rep=%s fold=%s model=%s): %s",
                    rep_id, fold_id, model_name, conditionMessage(e)))
    empty_schema
  })
  
  out
}

# —————————————————————————————————————————————————————————————————————
# Beeswarm from shap_handles for a chosen repeat (pooled over folds)
#   * Adds per-feature color normalization (rank → [0,1]) used previously
#   * Lightweight progress messages + robust errors
# —————————————————————————————————————————————————————————————————————
plot_shap_for_repeat_cls_from_handles <- function(
    shap_handles, df, method,
    rep_id, pos_class,
    nsim = 80, bg_n = 200, test_subsample = 200,
    title = NULL
) {
  # Normalize keys on BOTH sides
  method_norm <- .normalize_method(method)
  rep_norm    <- .coerce_rep_int(rep_id)
  
  # Inspect and filter handles for (method, rep)
  avail <- .inspect_handles_normalized(shap_handles)
  hits <- avail |>
    dplyr::filter(model == method_norm, rep == rep_norm) |>
    dplyr::arrange(fold)
  
  if (nrow(hits) == 0L) {
    stop(sprintf(
      "No handles for method=%s, rep=%s. Seen models=%s; reps for %s=%s",
      method_norm, rep_norm,
      paste0(sort(unique(avail$model)), collapse = ","),
      method_norm, paste0(sort(unique(avail$rep[avail$model == method_norm])), collapse = ",")
    ))
  }
  
  # Iterate folds and collect SHAP rows
  chunks <- vector("list", nrow(hits))
  for (i in seq_len(nrow(hits))) {
    h <- shap_handles[[ hits$idx[i] ]]
    message(sprintf("[SHAP] model=%s rep=%d fold=%d ...", method_norm, rep_norm, hits$fold[i]))
    chunks[[i]] <- tryCatch({
      compute_shap_for_handle_cls(
        handle = h, df = df, pos_class = pos_class,
        nsim = nsim, bg_n = bg_n, test_subsample = test_subsample
      )
    }, error = function(e) {
      message(sprintf("!! SHAP skipped for fold=%d (%s)", hits$fold[i], e$message))
      tibble::tibble(model = character(), iter = character(), outer_fold = character(),
                     feature = character(), shap_value = numeric(), feature_value = numeric())
    })
  }
  shap_long <- dplyr::bind_rows(chunks)
  if (!nrow(shap_long)) stop("SHAP produced no rows after pooling folds.")
  
  # ---- Per-feature color normalization (rank → [0,1]; constant → 0) ----
  shap_long <- shap_long %>%
    dplyr::group_by(feature) %>%
    dplyr::mutate(
      rankfvalue = rank(feature_value, ties.method = "average", na.last = "keep"),
      rng_min    = suppressWarnings(min(rankfvalue, na.rm = TRUE)),
      rng_max    = suppressWarnings(max(rankfvalue, na.rm = TRUE)),
      stdfvalue  = ifelse(
        is.finite(rng_min) & is.finite(rng_max) & (rng_max > rng_min),
        (rankfvalue - rng_min) / (rng_max - rng_min),
        0
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-rankfvalue, -rng_min, -rng_max)
  
  # Order features by mean(|SHAP|)
  mean_phi <- shap_long |>
    dplyr::group_by(feature) |>
    dplyr::summarise(mean_value = mean(abs(shap_value), na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(mean_value)
  shap_long$feature <- factor(shap_long$feature, levels = mean_phi$feature)
  
  # Title
  if (is.null(title)) {
    title <- sprintf("SHAP beeswarm (pooled over folds) — %s, repeat=%d", method_norm, rep_norm)
  }
  
  # Plot (color = stdfvalue in [0,1])
  ggplot2::ggplot(shap_long, ggplot2::aes(x = feature, y = shap_value, color = stdfvalue)) +
    ggforce::geom_sina(size = 1.3, alpha = 0.75, scale = "area") +
    ggplot2::coord_flip() +
    ggplot2::geom_text(
      data = mean_phi,
      ggplot2::aes(x = feature, y = -Inf, label = sprintf("%.3f", mean_value)),
      hjust = -0.2, fontface = "bold", size = 3, inherit.aes = FALSE
    ) +
    ggplot2::scale_color_gradient(name = "Feature value (within-feature rank)",
                                  limits = c(0, 1),
                                  low = "#FFCC33", high = "#6600CC",
                                  oob = scales::squish) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::labs(title = title, x = NULL, y = "SHAP value (impact on prediction)") +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(axis.ticks.y = ggplot2::element_blank(),
                   legend.position = "bottom")
}


# ---- Mean |SHAP| ranking from shap_handles across repeats & folds (classification) ----

compute_mean_shap_ranking_cls_from_handles <- function(
    shap_handles, df, method,
    repeats_to_use = NULL, pos_class = "ALS",
    seed = 12345, 
    # SHAP parameters (defaults are reasonable for a quick approximation)
    nsim = 30, bg_n = 80, test_subsample = 120, topk = 30
) {
  # pick handles for the target method (case-insensitive)
  mh <- shap_handles[ purrr::map_lgl(shap_handles, ~ tolower(.x$model_name) == tolower(method)) ]
  if (!length(mh)) stop(sprintf("No handles for method=%s", method))
  
  # which reps
  reps_all <- sort(unique(purrr::map_int(mh, ~ as.integer(gsub("\\D+","", .x$rep_id)))))
  reps     <- if (is.null(repeats_to_use)) reps_all else intersect(repeats_to_use, reps_all)
  if (!length(reps)) stop("No matching repeats in handles for the requested filter.")
  
  # count tasks for progress
  Ks <- vapply(reps, function(r) sum(purrr::map_int(mh, ~ as.integer(gsub("\\D+","", .x$rep_id)) == r)), integer(1))
  total_tasks <- sum(Ks)
  cat(sprintf(">> SHAP ranking(handles) %s | reps=%s | nsim=%d, bg_n=%d, test_subsample=%d\n",
              tolower(method), paste(reps, collapse=","), nsim, bg_n, test_subsample))
  
  chunks <- list(); task_counter <- 0L
  
  for (r in reps) {
    mh_r <- mh[purrr::map_lgl(mh, ~ as.integer(gsub("\\D+","", .x$rep_id)) == r)]
    K <- length(mh_r)
    for (i in seq_len(K)) {
      task_counter <- task_counter + 1L
      pct <- floor(100 * (task_counter - 1) / max(1L, total_tasks))
      cat(sprintf("\r   - [%d/%d] rep=%d fold=%s ... %3d%%",
                  task_counter, total_tasks, r, as.character(mh_r[[i]]$fold_id), pct))
      utils::flush.console()
      
      piece <- tryCatch({
        compute_shap_for_handle_cls(
          handle = mh_r[[i]], df = df, pos_class = pos_class,
          nsim = nsim, bg_n = bg_n, test_subsample = test_subsample
        )
      }, error = function(e) {
        cat(sprintf("\n!! SHAP failed (rep=%d, fold=%s): %s\n", r, as.character(mh_r[[i]]$fold_id), conditionMessage(e)))
        NULL
      })
      
      if (!is.null(piece) && nrow(piece)) {
        piece$repeat_id <- r
        chunks[[length(chunks) + 1L]] <- piece
      }
    }
  }
  cat(sprintf("\r   - done %d/%d ... 100%%\n", total_tasks, total_tasks)); utils::flush.console()
  
  if (!length(chunks)) stop("No SHAP values collected for ranking.")
  all_shap <- dplyr::bind_rows(chunks)
  
  # rank by mean(|SHAP|)
  rank_tbl <- all_shap %>%
    dplyr::group_by(feature) %>%
    dplyr::summarise(mean_abs_shap = mean(abs(shap_value), na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_abs_shap)) %>%
    dplyr::slice_head(n = topk)
  
  rank_tbl
}

