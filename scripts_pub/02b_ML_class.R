library(flextable)

source("scripts_pub/function/ML_class_func.R")

# Data prep ---------------------------------------------------------------

predictors <- c("NfL_csf_pgml", "GFAP_csf_pgml", "pTau181_csf", "pTau217_csf", 
                "Ab38_40_csf", "Ab42_40_csf_bridged","Cr", "CK")

predictors2 <- c("NfL_csf_pgml")

flatten_numeric <- function(col) {
  if (is.list(col)) {
    col <- vapply(col, function(x) if (length(x)==1) as.numeric(x) else NA_real_, numeric(1))
  }
  if (is.matrix(col) && ncol(col)==1) col <- as.numeric(col)
  if (is.factor(col)) col <- as.character(col)
  if (is.character(col)) col <- suppressWarnings(as.numeric(col))
  col[!is.finite(col)] <- NA_real_
  as.numeric(col)
}

df <- ROC_df %>%
  dplyr::select(ALS_label, dplyr::all_of(predictors)) %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(predictors), flatten_numeric)) %>%
  tidyr::drop_na()  
df$ALS_label <- factor(ifelse(df$ALS_label == 1, "ALS", "nonALS"))


df %>% nrow()
df %>% filter(ALS_label=="ALS") %>% nrow()
df %>% filter(ALS_label=="nonALS") %>% nrow()


# SETTINGS ---------------------------------------------------------
.seed_value <- 12345
set.seed(.seed_value)
outer_k <- 5        # number of outer folds; 5 means each outer test is ~20% (fold 1 uses fold 2-5 for training)
repeats <- 10       # number of repetitions of the entire outer 5-fold split
pos_class <- "ALS"    

# Models
models <- list(
  GLM  = list(method = "glm",       preds = predictors),
  RF   = list(method = "rf",        preds = predictors),
  SVM  = list(method = "svmRadial", preds = predictors),
  KNN  = list(method = "kknn",       preds = predictors),
  XGB  = list(method = "xgb",       preds = predictors),
  NfLonly  = list(method = "one",   preds = predictors2)  # NfL only model (for comparison)
)

## ── OUTER FOLDS ──────────────────────────
# make the k-fold (outer_k) * shared by all models (glm, rf, svm...)
outer_folds_list <- lapply(seq_len(repeats), function(r) {
  caret::createFolds(df$ALS_label, k = outer_k, returnTrain = FALSE)
})


# RUN: repeated nest-CV（repeats × outer_k folds）---------------------------

## (!) RUN the nested OOF + SHAP handle collection (across all outer folds & repeats)

res_all <- run_all_oof_collect(outer_folds_list)  # trains once per model/fold

# keep using your existing downstream code:
oof_preds <- res_all$oof_preds


# ---- OOF pred value vs True label / per each repeat ----------
# wait for 6 models (RF, XGN, NfL alone,...) * k repeats (6 models × k repeats)
rep_level_eval <- evaluate_rep_level_with_progress(
  oof_preds,
  B = 1000, stratified = TRUE, seed = .seed_value,
  neg_class = "nonALS", pos_class = "ALS",
  threshold_rule = "fixed" # "youden", "fixed"
) 

rep_level_eval 



# 💾 Save trained models and results -----------------------------------------

# .timestamp_tmp <- format(Sys.time(), "%Y%m%d_%H%M")
# .output_file_tmp <- paste0("output/Classification_Models_Trained_", .timestamp_tmp, ".rds")
# 
# if (!dir.exists("output")) dir.create("output")
# 
# saveRDS(
#   list(
#     metadata = list(
#       timestamp = Sys.time(),
#       seed = .seed_value,
#       outer_k = outer_k,
#       repeats = repeats,
#       n_samples = nrow(df),
#       n_ALS = sum(df$ALS_label == "ALS"),
#       n_nonALS = sum(df$ALS_label == "nonALS"),
#       outcome = "ALS_label",
#       predictors = predictors,
#       pos_class = pos_class,
#       models = names(models)
#     ),
#     outer_folds_list = outer_folds_list,
#     res_all = res_all,
#     rep_level_eval = rep_level_eval,
#     training_data = df
#   ),
#   .output_file_tmp,
#   compress = "xz"
# )
# 
# cat("\n=== Model Training Complete ===\n")
# cat("Saved to:", .output_file_tmp, "\n")
# cat("Key Points:\n")
# cat("  1. Out-of-fold predictions across", repeats, "repeats\n")
# cat("  2. SHAP values stored for RF/XGB/KNN/SVM models\n")
# cat("  3. All models trained on identical CV folds\n")
# cat("  4. Bootstrap evaluation (B=1000) included\n\n")
# 
# rm(.timestamp_tmp, .output_file_tmp)


# Load --------------------------------------------------------------------

# Load previously trained models
# .saved_file <- "output/Classification_Models_Trained_20260108_1452.rds"
# saved_results <- readRDS(.saved_file)
# 
# # Extract essential components
# outer_folds_list <- saved_results$outer_folds_list
# res_all <- saved_results$res_all
# rep_level_eval <- saved_results$rep_level_eval
# df <- saved_results$training_data
# 
# # Derive other variables from res_all
# oof_preds <- res_all$oof_preds
# shap_handles <- res_all$shap_handles
# per_repeat <- res_all$per_repeat
# 
# cat("Loaded results from:", saved_results$metadata$timestamp, "\n")
# cat("Models included:", paste(saved_results$metadata$models, collapse = ", "), "\n")
# 
# rm(.saved_file)


# Analysis ----------------------------------------------------------------

## choose representative repeat
represent_rep_id <- choose_rep_for_all_models(rep_level_eval, metric = "AUC")
represent_rep_id


ml_colors <- c(
  "RF"        = "#C96A1B",  # dark orange
  "XGB"       = "#E28E2C",  # mid orange
  "SVM"       = "#3F6F63",  # darker green-gray
  "KNN"       = "#7FA69A",  # mid green-gray
  "GLM"       = "#4C6FA1",  # muted blue
  "NfL alone" = "gray25"
)


ROC_ML_1step <- oof_preds %>% plot_repeat_roc(rep_id = represent_rep_id, colormap = ml_colors)

ROC_ML_1step_metrics <- rep_level_eval %>% filter(iter == represent_rep_id) %>% 
  make_model_metrics_plot(colormap = ml_colors,
                          model_order = c("KNN","GLM","NfL alone","SVM","XGB","RF") 
  )


# ---- ROC based on Mean Predicted Value ------------------------------------

# Compute mean predicted value per model/row across all repeats
oof_mean <- oof_preds %>%
  dplyr::group_by(model, row_id) %>%
  dplyr::summarise(
    y_score = mean(as.numeric(y_score), na.rm = TRUE),
    y_true  = dplyr::first(y_true),  # same for all rows
    .groups = "drop"
  )

# Evaluate AUC based on mean predicted value (across all repeats)
ensemble_eval <- evaluate_ensemble_with_progress(
  oof_preds,
  B = 2000, stratified = TRUE, seed = .seed_value,
  neg_class = "nonALS", pos_class = "ALS",
  threshold_rule = "fixed"
)

ensemble_eval 

ROC_ML_1step_mean_Pub <- oof_mean %>% plot_ensemble_roc(colormap = ml_colors,
                                                        title="Ensemble OOF ROC (based on mean predicted value of all repeats)")

ROC_ML_1step_mean_Pub <- ROC_ML_1step_mean_Pub$plot+
  labs(
    title = "Out-of-fold classification performance across models",
    x = "1 - Specificity",
    y = "Sensitivity") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  coord_equal() +  # square plotting area
  theme_classic(base_size = 15) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
  theme(
    legend.position = c(0.75,0.25),
    #legend.position = "right",
    legend.text  = element_text(size = 11),
    legend.title = element_text(size = 12),
    axis.title.x = element_text(margin=ggplot2::margin(t=10)),
    axis.title.y = element_text(margin=ggplot2::margin(r=10)),
    plot.title = element_text(hjust = 0, size = 20, face = "bold", margin=ggplot2::margin(b=15))
  )

ROC_ML_1step_mean_Pub

ROC_ML_1step_mean_metrics <- ensemble_eval %>% 
  make_model_metrics_plot(colormap = ml_colors, 
                          model_order = c("NfL alone","GLM","KNN","SVM","XGB","RF"))

ROC_ML_1step_mean_metrics$plot <-ROC_ML_1step_mean_metrics$plot+
  theme(plot.margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0))


# 📊 Fogire 2CD --------------------------------------------------------------

Fig2CD <- (ROC_ML_1step_mean_Pub|ROC_ML_1step_mean_metrics$plot) + plot_layout(widths = c(1, 0.8))
Fig2CD


# Nested CV results ----------------------------------------------------------

rep_level_eval 
model_levels_vec <- c("NfL alone", "GLM","KNN","SVM","XGB","RF")  


# 🔍 SuppleFig5 : Across 10 repeats -------------------------------------------------------

SuppleFig5 <- ggplot(rep_level_eval %>%
                       mutate(model = dplyr::recode(model, "NfLonly" = "NfL alone"),
                              model = factor(model, levels = model_levels_vec),
                              iter = factor(iter, levels = paste0("Rep", 1:10)),
                              is_representative = (iter == represent_rep_id)),
                     aes(x = model, y = AUC)) +
  geom_boxplot(aes(fill = model), outlier.shape = NA, alpha = 0.8, width = 0.35,
               show.legend = FALSE) +
  geom_jitter(size = 0.5, alpha = 0.5, color = "gray20", width=0.2) +
  scale_fill_manual(values = ml_colors, breaks = model_levels_vec) +
  labs(title = "Variability in diagnostic performance across repeated nested cross-validation",
       x = "Model", y = "Cross-validated AUC") +
  scale_y_continuous(limits = c(0.75, 1), breaks = seq(0.75, 1, 0.05)) +
  coord_flip()+
  theme_bw(base_size = 18) +
  theme(legend.position = "right",
        plot.title = element_text(size=20, face="bold", margin=ggplot2::margin(b=15)),
        axis.title.x = element_text(size=16, margin=ggplot2::margin(t=15)),
        axis.title.y = element_text(size=16, margin=ggplot2::margin(r=5))
  ) 

SuppleFig5


# Nested SHAP results -----------------------------------------------------

per_repeat    <- res_all$per_repeat
shap_handles  <- res_all$shap_handles

p_xgb_class <- plot_shap_for_repeat_cls_from_handles(
  shap_handles, df, method = "xgb",
  rep_id = represent_rep_id,                 # representative repeat
  pos_class = "ALS",
  nsim = 80, bg_n = 200, test_subsample = 200
)

p_xgb_class

p_rf_class <- plot_shap_for_repeat_cls_from_handles(
  shap_handles, df, method = "rf",
  rep_id = represent_rep_id,                 # representative repeat
  pos_class = "ALS",
  nsim = 80, bg_n = 200, test_subsample = 200
)

Fig2E_pub <- p_rf_class +
  scale_x_discrete(labels = function(x) ifelse(x %in% names(variable_label_map), variable_label_map[x], x))+
  theme_bw(base_size = 16) +
  labs(
    title = "SHAP summary plot for RF model",
    y = "SHAP value (impact on model output)"
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

Fig2E_pub

# mean |SHAP| -------------------------------------------------------------

# per-model ranking from shap_handles
rank_xgb_cls <- compute_mean_shap_ranking_cls_from_handles(shap_handles, df, "xgb", repeats_to_use = 1:10,
                                                           nsim = 80, bg_n=200, test_subsample = 200) %>% mutate(model="xgb")
rank_rf_cls  <- compute_mean_shap_ranking_cls_from_handles(shap_handles, df, "rf",  repeats_to_use = 1:10,
                                                           nsim = 80, bg_n=200, test_subsample = 200) %>% mutate(model="rf")
rank_knn_cls  <- compute_mean_shap_ranking_cls_from_handles(shap_handles, df, "knn",  repeats_to_use = 1:10,
                                                            nsim = 80, bg_n=200, test_subsample = 200) %>% mutate(model="knn")
rank_svm_cls  <- compute_mean_shap_ranking_cls_from_handles(shap_handles, df, "svm",  repeats_to_use = 1:10,
                                                            nsim = 80, bg_n=200, test_subsample = 200) %>% mutate(model="svm")


# 💾 Save SHAP calc --------------------------------------------------------------------
# .timestamp_shap <- format(Sys.time(), "%Y%m%d_%H%M")
# .output_file_shap <- paste0("output/SHAP_Rankings_Publication_", .timestamp_shap, ".rds")
# if (!dir.exists("output")) dir.create("output")
# 
# saveRDS(
#   list(
#     metadata = list(
#       timestamp = Sys.time(),
#       nsim = 80,
#       bg_n = 200,
#       test_subsample = "all",
#       repeats_used = 1:10,
#       topk = 30,
#       models = c("xgb", "rf", "knn", "svm"),
#       n_features = length(unique(rank_all_cls$feature)),
#       total_tasks = 10 * 5 * 4  # 10 repeats × 5 folds × 4 models
#     ),
#     rank_xgb = rank_xgb_cls,
#     rank_rf = rank_rf_cls,
#     rank_knn = rank_knn_cls,
#     rank_svm = rank_svm_cls,
#     rank_all = rank_all_cls
#   ),
#   .output_file_shap,
#   compress = "xz"
# )
# 
# cat("\n=== SHAP Rankings Saved ===\n")
# cat("Saved to:", .output_file_shap, "\n")
# cat("Key Points:\n")
# cat("  1. SHAP values computed across 10 repeats (50 folds per model)\n")
# cat("  2. All test samples used (no subsampling)\n")
# cat("  3. Models: XGB, RF, KNN, SVM\n")
# cat("  4. Mean |SHAP| rankings for top", 30, "features\n\n")
# 
# rm(.timestamp_shap, .output_file_shap)

# ─────────────────────────────────────────────────────────────────────────
# To load later:
# ─────────────────────────────────────────────────────────────────────────
shap_results <- readRDS("output/SHAP_Rankings_Publication_20260108_1617.rds")
rank_xgb_cls <- shap_results$rank_xgb
rank_rf_cls <- shap_results$rank_rf
rank_knn_cls <- shap_results$rank_knn
rank_svm_cls <- shap_results$rank_svm
cat("Loaded SHAP results from:", format(shap_results$metadata$timestamp), "\n")


# ---- Feature importance barplot -------------------------------------------

# Bind and normalize to shares
rank_all_cls <- dplyr::bind_rows(rank_xgb_cls, rank_rf_cls, rank_knn_cls, rank_svm_cls) %>%
  mutate(feature = as.character(feature),
         model   = as.character(model)) %>%
  complete(feature, model,
           fill = list(mean_abs_shap = 0))  # if a feature wasn't in topk for a model

rank_all_cls

# By default use raw mean(|SHAP|) averaged across models (like your regression code)
order_tbl_raw <- rank_all_cls %>%
  group_by(feature) %>%
  summarise(order_key = mean(mean_abs_shap, na.rm = TRUE), .groups = "drop") %>%
  arrange(order_key)


# Within-model % normalization for cross-model comparability ----
# Convert each model's bars to percentage contribution within that model.
rank_plot_norm <- rank_all_cls %>%
  group_by(model) %>%
  mutate(rel_norm = mean_abs_shap / sum(mean_abs_shap, na.rm = TRUE)) %>%
  ungroup()

# Recompute feature order using across-model mean of normalized values
order_tbl_norm <- rank_plot_norm %>%
  group_by(feature) %>%
  summarise(order_key = mean(rel_norm, na.rm = TRUE), .groups = "drop") %>%
  arrange(order_key)

rank_plot_norm <- rank_plot_norm %>%
  mutate(feature = factor(feature, levels = order_tbl_norm$feature),
         model   = factor(model, levels = c("knn","svm","xgb","rf")),
         model = dplyr::recode(model,
                               "rf"  = "RF",
                               "xgb" = "XGB",
                               "svm" = "SVM",
                               "knn" = "KNN"))



# 📊 Figure 2EF --------------------------------------------------------------

Fig2F_pub <- ggplot(rank_plot_norm,
                    aes(x = feature, y = rel_norm, fill = model)) +
  geom_col(position = "dodge", color="gray95", linewidth = 0.2, alpha=0.8) +
  scale_x_discrete(labels = function(x) ifelse(x %in% names(variable_label_map), variable_label_map[x], x))+
  coord_flip() +
  labs(title = "Normalized |SHAP| across repeats",
       x = NULL, y = "Normalized |SHAP|") +
  theme_bw(base_size = 16) +
  theme(
    plot.title   = element_text(size = 18, face = "bold", margin=ggplot2::margin(b=15)),
    axis.title = element_text(size = 14),
    axis.title.x  = element_text(size=16, margin=ggplot2::margin(t=15,b=-10)),
    axis.text.y  = element_text(size = 14),
    axis.text.x  = element_text(size = 14),
    legend.text  = element_text(size = 12),
  )+
  scale_fill_manual(values = ml_colors) +
  theme(legend.position = "bottom")


Fig2EF <- Fig2E_pub|Fig2F_pub

Fig2EF


# QC/ Raw bars ----------------------------

# rank_plot_raw <- rank_all_cls %>%
#   mutate(feature = factor(feature, levels = order_tbl_raw$feature),
#          model   = factor(model, levels = c("knn","svm","rf","xgb")))
# 
# ggplot(rank_plot_raw,
#        aes(x = feature, y = mean_abs_shap, fill = model)) +
#   geom_col(position = "dodge") +
#   coord_flip() +
#   labs(title = "Mean |SHAP| ranking across all repeats",
#        x = NULL, y = "Mean |SHAP| (model scale)") +
#   theme_bw(base_size = 13) +
#   scale_fill_brewer(palette = "Set1") +
#   theme(legend.position = "bottom")
