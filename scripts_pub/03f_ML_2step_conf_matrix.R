# Confusion matrices for the two-step framework (Supplementary Figure S7)

library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)
library(purrr)

pos     <- "ALS"
neg     <- "nonALS"
neg_lab <- "DC"

# Extract OOF components from the representative run
.cm2_y_true     <- kf_selected$pooled$oof$y_true
.cm2_bucket     <- as.character(kf_selected$pooled$oof$bucket)
.cm2_final_pred <- as.character(kf_selected$pooled$oof$final_pred[["sp95"]])
.cm2_folds      <- kf_selected$data$folds
.cm2_df         <- kf_selected$data$df_clean
.cm2_n_total    <- length(.cm2_y_true)

# Helper: NfL one-step predictions using fold-wise thresholds
.cm2_get_nfl_pred <- function(variant = c("se95", "sp95")) {
  variant  <- match.arg(variant)
  base_tbl <- kf_selected$pooled$nfl_baselines
  row_idx  <- which(base_tbl$baseline == paste0("nfl_", variant))
  th_list  <- base_tbl$th_train_all[[row_idx]]
  pred     <- rep(NA_character_, nrow(.cm2_df))
  for (i in seq_along(.cm2_folds)) {
    idx       <- .cm2_folds[[i]]
    pred[idx] <- ifelse(.cm2_df$NfL_csf_pgml[idx] >= th_list[i], pos, neg)
  }
  factor(pred, levels = c(pos, neg))
}

# Helper: build tidy 2x2 confusion matrix tibble (ALS top, DC bottom)
.cm2_build_tbl_2x2 <- function(y_true, y_pred) {
  pred_levels <- c(pos, neg_lab)
  ref_levels  <- c(pos, neg_lab)
  y_true_lab  <- as.character(y_true); y_true_lab[y_true_lab == neg] <- neg_lab
  y_pred_lab  <- as.character(y_pred); y_pred_lab[y_pred_lab == neg] <- neg_lab
  expand.grid(
    Prediction = factor(pred_levels, levels = pred_levels),
    Reference  = factor(ref_levels,  levels = ref_levels),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      Count    = mapply(function(p, r) sum(y_true_lab == r & y_pred_lab == p),
                        as.character(Prediction), as.character(Reference)),
      diagonal = as.character(Prediction) == as.character(Reference)
    ) %>%
    dplyr::group_by(Prediction) %>%
    dplyr::mutate(n_row = sum(Count), Prop = Count / n_row) %>%
    dplyr::ungroup()
}

# Helper: compute PPV/NPV from 2x2 tibble
.cm2_metrics_2x2 <- function(tbl) {
  tp <- tbl$Count[tbl$Prediction == pos     & tbl$Reference == pos]
  tn <- tbl$Count[tbl$Prediction == neg_lab & tbl$Reference == neg_lab]
  fp <- tbl$Count[tbl$Prediction == pos     & tbl$Reference == neg_lab]
  fn <- tbl$Count[tbl$Prediction == neg_lab & tbl$Reference == pos]
  list(sens = tp/(tp+fn), spec = tn/(tn+fp),
       ppv  = tp/(tp+fp), npv  = tn/(tn+fn))
}

.cm2_ml_2x2 <- function(tbl) {
  m  <- .cm2_metrics_2x2(tbl)
  tp <- tbl$Count[tbl$Prediction == pos     & tbl$Reference == pos]
  tn <- tbl$Count[tbl$Prediction == neg_lab & tbl$Reference == neg_lab]
  fp <- tbl$Count[tbl$Prediction == pos     & tbl$Reference == neg_lab]
  fn <- tbl$Count[tbl$Prediction == neg_lab & tbl$Reference == pos]
  n_als <- tp + fn
  n_dc  <- tn + fp
  tibble::tibble(
    Prediction = factor(c(pos, neg_lab), levels = c(pos, neg_lab)),
    label      = c(sprintf("PPV = %.2f", m$ppv),
                   sprintf("NPV = %.2f", m$npv))
  )
}


# Helper: plot 2x2 confusion matrix
.cm2_ggplot_2x2 <- function(tbl, metric_labels, title = NULL, base_size = 14) {
  pred_levels <- c(pos, neg_lab)
  
  # Column totals for FP/FN rate labels
  .cm2_col_totals <- tbl %>%
    dplyr::group_by(Reference) %>%
    dplyr::summarise(n_col = sum(Count), .groups = "drop")
  fp_count <- tbl$Count[tbl$Prediction == pos     & tbl$Reference == neg_lab]
  fn_count <- tbl$Count[tbl$Prediction == neg_lab & tbl$Reference == pos]
  n_als    <- sum(tbl$Count[tbl$Reference == pos])
  n_dc     <- sum(tbl$Count[tbl$Reference == neg_lab])
  
  .cm2_col_labels <- tibble::tibble(
    Reference = factor(c(pos, neg_lab), levels = c(pos, neg_lab)),
    label     = c(sprintf("FN = %.1f%%", fn_count/n_als*100),
                  sprintf("FP = %.1f%%", fp_count/n_dc*100))
  )
  
  ggplot(tbl, aes(x = Reference, y = Prediction)) +
    geom_tile(aes(fill = diagonal), color = "gray30", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%d\n(%.0f%%)", Count, Prop * 100)),
              size = 5, fontface = "bold", lineheight = 1.1) +
    geom_text(data = metric_labels,
              aes(x = 2.65, y = Prediction, label = label),
              inherit.aes = FALSE,
              hjust = 0, size = 4, color = "gray20") +
    geom_text(data = .cm2_col_labels,
              aes(x = Reference, y = 0.35, label = label),
              inherit.aes = FALSE,
              size = 4, color = "gray20") +
    scale_fill_manual(values = c("TRUE" = "gray80", "FALSE" = "white"),
                      guide  = "none") +
    scale_x_discrete(position = "top", expand = expansion(add = c(0.5, 1.4))) +
    scale_y_discrete(limits = rev(pred_levels), expand = expansion(add = c(0.7, 0.5))) +
    labs(title = title, x = "Final diagnosis", y = "Prediction") +
    theme_minimal(base_size = base_size) +
    theme(
      axis.ticks      = element_blank(),
      panel.grid      = element_blank(),
      axis.text.y     = element_text(size = 12),
      axis.text.x.top = element_text(size = 12, vjust = -0.5),
      plot.title      = element_text(size = 16, hjust = 0.5),
      plot.margin     = margin(t = 15, r = 0, b = 15, l = 0),
      axis.title.x    = element_text(hjust = 0.3)
    ) +
    coord_fixed(clip="off")
}
# Helper: build tidy 2xN tibble with DC display label
.cm2_build_matrix_tbl <- function(pred_group, pred_levels) {
  y_true_lab <- as.character(.cm2_y_true)
  y_true_lab[y_true_lab == neg] <- neg_lab
  ref_levels <- c(pos, neg_lab)
  expand.grid(
    Prediction = factor(pred_levels, levels = pred_levels),
    Reference  = factor(ref_levels,  levels = ref_levels),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      Count = mapply(function(p, r) {
        sum(pred_group == p & y_true_lab == r, na.rm = TRUE)
      }, as.character(Prediction), as.character(Reference))
    ) %>%
    dplyr::group_by(Prediction) %>%
    dplyr::mutate(n_total = sum(Count), Prop = Count / n_total) %>%
    dplyr::ungroup()
}

# Helper: metric labels for 2xN matrix
.cm2_build_metric_labels <- function(tbl, pred_levels) {
  purrr::map_dfr(pred_levels, function(p) {
    sub       <- dplyr::filter(tbl, Prediction == p)
    n         <- sum(sub$Count)
    pct_total <- n / .cm2_n_total * 100
    als_pct   <- sub$Count[sub$Reference == pos] / n * 100
    label <- if (p == "Rule-in") {
      sprintf("PPV = %.2f", sub$Count[sub$Reference == pos] / n)
    } else if (p == "Rule-out") {
      sprintf("NPV = %.2f", sub$Count[sub$Reference == neg_lab] / n)
    } else {
      # sprintf("n = %d (%.0f%%)\nALS: %.0f%%", n, pct_total, als_pct)
      sprintf("n = %d \n (%.0f%%)", n, pct_total)
    }
    tibble::tibble(Prediction = factor(p, levels = pred_levels), label = label)
  })
}

# Helper: plot 2xN confusion matrix
.cm2_ggplot_matrix <- function(tbl, correct_cells, metric_labels,
                               title = NULL, base_size = 14) {
  pred_levels <- levels(tbl$Prediction)
  tbl <- tbl %>%
    dplyr::mutate(
      correct = mapply(function(p, r) {
        any(correct_cells$Prediction == p & correct_cells$Reference == r)
      }, as.character(Prediction), as.character(Reference))
    )
  
  #tp, fn
  fn_count <- tbl$Count[tbl$Prediction == "Rule-out" & tbl$Reference == pos]
  fp_count <- tbl$Count[tbl$Prediction == "Rule-in"  & tbl$Reference == neg_lab]
  n_als    <- sum(tbl$Count[tbl$Reference == pos])
  n_dc     <- sum(tbl$Count[tbl$Reference == neg_lab])
  
  .cm2_col_labels <- tibble::tibble(
    Reference = factor(c(pos, neg_lab), levels = c(pos, neg_lab)),
    label     = c(sprintf("FN = %.1f%%", fn_count/n_als*100),
                  sprintf("FP = %.1f%%", fp_count/n_dc*100))
  )
  
  ggplot(tbl, aes(x = Reference, y = Prediction)) +
    geom_tile(aes(fill = correct), color = "gray30", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%d\n(%.0f%%)", Count, Prop * 100)),
              size = 5, fontface = "bold", lineheight = 1.1) +
    geom_text(data = metric_labels,
              aes(x = 2.65, y = Prediction, label = label),
              inherit.aes = FALSE,
              hjust = 0, size = 4, color = "gray20", lineheight = 1.2) +
    geom_text(data = .cm2_col_labels,
              aes(x = Reference, y = 0.35, label = label),
              inherit.aes = FALSE,
              size = 4, color = "gray20") +
    scale_fill_manual(values = c("TRUE" = "gray80", "FALSE" = "white"),
                      guide  = "none") +
    scale_x_discrete(position = "top", expand = expansion(add = c(0.5, 1.4))) +
    scale_y_discrete(limits = rev(pred_levels)) +
    labs(title = title, x = "Final diagnosis", y = "Prediction") +
    theme_minimal(base_size = base_size) +
    theme(
      axis.ticks      = element_blank(),
      panel.grid      = element_blank(),
      axis.text.y     = element_text(size = 12),
      axis.text.x.top = element_text(size = 12, vjust = -0.5),
      plot.title      = element_text(size = 16, hjust = 0.5),
      plot.margin     = margin(t = 15, r = 40, b = 25, l = 0),
      axis.title.x    = element_text(hjust = 0.3)
    ) +
    coord_fixed(clip="off")
}

# Panel A-left: NfL 95% Sensitivity cut-off (2x2) -------------------------
.cm2_pred_se95 <- .cm2_get_nfl_pred("se95")
.cm2_tbl_se95  <- .cm2_build_tbl_2x2(as.character(.cm2_y_true),
                                     as.character(.cm2_pred_se95))
.cm2_panel_A_left <- .cm2_ggplot_2x2(
  .cm2_tbl_se95,
  metric_labels = .cm2_ml_2x2(.cm2_tbl_se95),
  title = "95% Sensitivity cut-off"
)

# Panel A-right: NfL 95% Specificity cut-off (2x2) ------------------------
.cm2_pred_sp95 <- .cm2_get_nfl_pred("sp95")
.cm2_tbl_sp95  <- .cm2_build_tbl_2x2(as.character(.cm2_y_true),
                                     as.character(.cm2_pred_sp95))
.cm2_panel_A_right <- .cm2_ggplot_2x2(
  .cm2_tbl_sp95,
  metric_labels = .cm2_ml_2x2(.cm2_tbl_sp95),
  title = "95% Specificity cut-off"
)

# Panel B-left: Step 1 (2x3) ----------------------------------------------
.cm2_pred_step1   <- dplyr::case_when(
  .cm2_bucket == "high" ~ "Rule-in",
  .cm2_bucket == "gray" ~ "Gray-zone",
  .cm2_bucket == "low"  ~ "Rule-out"
)
.cm2_levels_step1  <- c("Rule-in", "Gray-zone", "Rule-out")
.cm2_tbl_step1     <- .cm2_build_matrix_tbl(.cm2_pred_step1, .cm2_levels_step1)
.cm2_metrics_step1 <- .cm2_build_metric_labels(.cm2_tbl_step1, .cm2_levels_step1)
.cm2_correct_step1 <- tibble::tibble(Prediction = c("Rule-in", "Rule-out"),
                                     Reference  = c(pos, neg_lab))
.cm2_panel_B_left <- .cm2_ggplot_matrix(
  .cm2_tbl_step1,
  correct_cells = .cm2_correct_step1,
  metric_labels = .cm2_metrics_step1,
  title         = "Step 1: applying two NfL cut-offs"
)

# Panel B-right: Two-step overall (2x3) -----------------------------------
.cm2_pred_overall   <- dplyr::case_when(
  .cm2_bucket == "high" | (.cm2_bucket == "gray" & .cm2_final_pred == pos) ~ "Rule-in",
  .cm2_bucket == "gray" & .cm2_final_pred == neg                           ~ "Still in gray-zone",
  .cm2_bucket == "low"                                                     ~ "Rule-out"
)
.cm2_levels_overall  <- c("Rule-in", "Still in gray-zone", "Rule-out")
.cm2_tbl_overall     <- .cm2_build_matrix_tbl(.cm2_pred_overall, .cm2_levels_overall)
.cm2_metrics_overall <- .cm2_build_metric_labels(.cm2_tbl_overall, .cm2_levels_overall)
.cm2_correct_overall <- tibble::tibble(Prediction = c("Rule-in", "Rule-out"),
                                       Reference  = c(pos, neg_lab))
.cm2_panel_B_right <- .cm2_ggplot_matrix(
  .cm2_tbl_overall,
  correct_cells = .cm2_correct_overall,
  metric_labels = .cm2_metrics_overall,
  title         = "Step 2: gray-zone refinement using RF"
)


# Assemble Panel A block (title + two 2x2) --------------------------------
.cm2_A_title <- ggdraw() +
  draw_label("One-step dichotomisation using a single NfL cut-off",
             fontface = "bold", size = 20, hjust = 0, vjust = 0.5,
             x = 0.05)

.cm2_A_panels <- plot_grid(
  .cm2_panel_A_left, .cm2_panel_A_right,
  nrow = 1
)

.cm2_A_block <- plot_grid(
  .cm2_A_title, .cm2_A_panels,
  ncol        = 1,
  rel_heights = c(0.1, 1)
)

# Assemble Panel B block (title + two 2x3) --------------------------------
.cm2_B_title <- ggdraw() +
  draw_label("Two-step approach",
             fontface = "bold", size = 20, hjust = 0, vjust = 0.5,
             x = 0.05)

.cm2_B_panels <- plot_grid(
  .cm2_panel_B_left, .cm2_panel_B_right,
  nrow = 1
)

.cm2_B_block <- plot_grid(
  .cm2_B_title, .cm2_B_panels,
  ncol        = 1,
  rel_heights = c(0.1, 1)
)

# 🔍 SuppleFig7  -------------------------------------
SuppleFig7_added <- plot_grid(
  .cm2_A_block,
  .cm2_B_block,
  ncol        = 1,
  rel_heights = c(1, 1.5),
  labels      = c("A", "B"),
  label_size  = 22,
  label_x     = 0,
  label_y     = 1.02
)+theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))


SuppleFig7_added

#ggsave("output/supple_figs/SuppleFig7_added.pdf", SuppleFig7_added, device = cairo_pdf, width = 320, height = 280, units = "mm")

