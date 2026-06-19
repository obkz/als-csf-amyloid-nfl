# Visualization function for regression slope prediction
# With prediction error band (tolerance zone)
library(ggplot2)
library(dplyr)
library(patchwork)

# Helper: calculate metrics for annotation
calc_metrics_for_plot <- function(truth, pred) {
  residual <- truth - pred
  rss <- sum(residual^2)
  tss <- sum((truth - mean(truth))^2)
  r2 <- 1 - rss/tss
  rmse <- sqrt(mean(residual^2))
  mae <- mean(abs(residual))
  bias <- mean(residual)
  
  list(R2 = r2, RMSE = rmse, MAE = mae, Bias = bias)
}

# Predicted vs observed scatter plot with flexible color options
plot_predobs_reg_slope <- function(truth, pred, model_name = "Model",
                                   point_size = 2, point_alpha = 0.9,
                                   show_loess = FALSE,
                                   xlim = NULL, ylim = NULL,
                                   show_tolerance_band = FALSE,
                                   tolerance_width = 0.5,
                                   title = NULL,
                                   color_by = c("abs_residual", "none", "residual"),
                                   color_palette = c("viridis", "gray", "blue_red", "custom"),
                                   color_low = NULL,
                                   color_high = NULL,
                                   color_mid = NULL,
                                   color_single = "steelblue",
                                   legend_name = "Absolute error in\nslope prediction") {
  
  stopifnot(length(truth) == length(pred))
  
  color_by <- match.arg(color_by)
  color_palette <- match.arg(color_palette)
  
  df_plot <- data.frame(
    observed = truth,
    predicted = pred,
    residual = truth - pred,
    abs_residual = abs(truth - pred)
  )
  
  metrics <- calc_metrics_for_plot(truth, pred)
  
  anno_text <- sprintf(
    "R² = %.3f\nRMSE = %.3f\n",#MAE = %.3f\nBias = %.3f",
    metrics$R2, metrics$RMSE#, metrics$MAE, metrics$Bias
  )
  
  # Axis limits
  if (is.null(xlim) || is.null(ylim)) {
    lims <- range(c(df_plot$observed, df_plot$predicted))
    if (is.null(xlim)) xlim <- lims
    if (is.null(ylim)) ylim <- lims
  }
  
  # Default title
  if (is.null(title)) {
    title <- paste0("Predicted vs Observed Slope: ", model_name)
  }
  
  # Base plot
  p <- ggplot(df_plot, aes(x = predicted, y = observed))
  
  # Add tolerance band (y = x ± tolerance_width)
  if (show_tolerance_band) {
    p <- p + geom_ribbon(
      data = data.frame(x = seq(xlim[1], xlim[2], length.out = 100)),
      aes(x = x, ymin = x - tolerance_width, ymax = x + tolerance_width),
      fill = "gray80", alpha = 0.3, inherit.aes = FALSE
    )
  }
  
  p <- p + geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
                       color = "gray50", linewidth = 0.5)
  
  # Add points with color mapping
  if (color_by == "none") {
    p <- p + geom_point(color = color_single, size = point_size, alpha = point_alpha)
    
  } else if (color_by == "abs_residual") {
    p <- p + geom_point(aes(color = abs_residual), size = point_size, alpha = point_alpha)
    
    if (color_palette == "viridis") {
      p <- p + scale_color_viridis_c(option = "viridis", direction = -1, name = legend_name)
    } else if (color_palette == "gray") {
      low_col <- if (!is.null(color_low)) color_low else "gray90"
      high_col <- if (!is.null(color_high)) color_high else "gray10"
      p <- p + scale_color_gradient(low = low_col, high = high_col, name = legend_name)
    } else if (color_palette == "blue_red") {
      low_col <- if (!is.null(color_low)) color_low else "#3498DB"
      high_col <- if (!is.null(color_high)) color_high else "#E74C3C"
      p <- p + scale_color_gradient(low = low_col, high = high_col, name = legend_name)
    } else if (color_palette == "custom") {
      if (is.null(color_low) || is.null(color_high)) {
        stop("For 'custom' palette, please specify both color_low and color_high")
      }
      p <- p + scale_color_gradient(low = color_low, high = color_high, name = legend_name)
    }
    
  } else if (color_by == "residual") {
    p <- p + geom_point(aes(color = residual), size = point_size, alpha = point_alpha)
    
    if (color_palette == "blue_red") {
      low_col <- if (!is.null(color_low)) color_low else "blue"
      mid_col <- if (!is.null(color_mid)) color_mid else "white"
      high_col <- if (!is.null(color_high)) color_high else "red"
      p <- p + scale_color_gradient2(low = low_col, mid = mid_col, high = high_col, 
                                     midpoint = 0, name = legend_name)
    } else if (color_palette == "custom") {
      if (is.null(color_low) || is.null(color_mid) || is.null(color_high)) {
        stop("For diverging 'custom' palette, please specify color_low, color_mid, and color_high")
      }
      p <- p + scale_color_gradient2(low = color_low, mid = color_mid, high = high_col, 
                                     midpoint = 0, name = legend_name)
    } else {
      p <- p + scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                     midpoint = 0, name = legend_name)
    }
  }
  
  # Add LOESS if requested
  if (show_loess) {
    p <- p + geom_smooth(method = "loess", se = TRUE, color = "#E74C3C", 
                         linewidth = 1, alpha = 0.2)
  }
  
  p <- p +
    annotate("text", x = xlim[1], y = ylim[2], 
             label = anno_text, hjust = 0, vjust = 1, 
             size = 4, fontface = "bold") +
    coord_equal(xlim = xlim, ylim = ylim) +
    labs(
      title = title,
      x = "Predicted slope (points/month)",
      y = "Observed slope (points/month)",
      caption = if (show_tolerance_band) {
        if (show_loess) {
          sprintf("dashed line: perfect prediction (y=x)\nGray band: ±%.1f points/month tolerance\nRed curve: LOESS smoothing", tolerance_width)
        } else {
          sprintf("dashed line: perfect prediction (y=x)\nGray band: ±%.1f points/month tolerance", tolerance_width)
        }
      } else {
        if (show_loess) {
          "dashed line: perfect prediction (y=x)\nRed curve: LOESS smoothing"
        } else {
          "dashed line: perfect prediction (y=x)"
        }
      }
    ) +
    theme_classic(base_size = 16) +
    theme(
      legend.position = if (color_by == "none") "none" else "bottom",
      plot.caption = element_text(hjust = 0, size = 12, color = "gray40"),
      legend.title = element_text(size = 14),
      plot.title = element_text(size=18, hjust = 0, face = "bold")
    )
  
  p
}


# Load --------------------------------------------------------------------

# trained models and results from repeated nested CV (from 05a_ML_reg_train.R)

#trained_results <- readRDS("output/XXX.rds")
#nested_rep_tbl <- trained_results$ml_results
#glm_subset_results <- trained_results$glm_results
#metadata <- trained_results$metadata

# df <- as.data.frame(trained_results$training_data)
# 
# cat("=== Loaded Trained Models ===\n")
# cat("Training timestamp:", format(.trained_tmp$metadata$timestamp, "%Y-%m-%d %H:%M:%S"), "\n")
# cat("Sample size:", .trained_tmp$metadata$n_samples, "\n")
# cat("Random seed:", .trained_tmp$metadata$seed, "\n")
# cat("CV configuration:", .trained_tmp$metadata$n_repeats, "repeats ×", 
#     .trained_tmp$metadata$K_outer, "outer folds ×", 
#     .trained_tmp$metadata$inner_k, "inner folds\n")
# cat("Outcome:", .trained_tmp$metadata$outcome, "\n")
# cat("Number of predictors:", length(.trained_tmp$metadata$predictors), "\n\n")
# cat("Predictors used:", paste(.trained_tmp$metadata$predictors, collapse = ", "), "\n\n")


# RF ----------------------------------------------------------------------

SuppleFig_slope_RF <- plot_predobs_reg_slope(
  model_name = "RF",
  pred = nested_rep_tbl$oof_pred_mean[[which(nested_rep_tbl$model == "rf")]],
  truth = df$slope_total,
  title = "(A) Cross-validated prediction of ALSFRS-R slopes by RF model",
  legend_name = "Absolute Error\n(points/month)",
  color_palette = "custom",
  xlim = c(-3.5, 0.5), ylim = c(-3.5, 0.5),
  color_low  = "gray80",color_high = "black"
)

SuppleFig_slope_RF


# 🔍 Supple Fig.11 ---------------------------------------------------------------

SuppleFigS11A <- SuppleFig_slope_RF

SuppleFigS11A
