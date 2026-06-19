# Load required library
library(gghalves)
library(patchwork)


plot_biomarker_comparison <- function(
    data,
    group_var = "ALS_label",
    biomarkers = NULL,
    y_log = FALSE,
    biomarker_labels = NULL,
    title = "CSF Biomarker Comparison",
    ncol = NULL,
    add_stats = TRUE,
    stat_method = "fdr",
    show_xlab = TRUE,
    show_ylab = TRUE,  
    show_n = TRUE,  # NEW: control whether to show sample sizes
    show_points = FALSE,
    show_violin = TRUE,
    trim_outliers = TRUE,
    trim_percentile = 0.9999,
    color_palette = c("#E64B35", "#D4AC0D", "#4DBBD5", "#00A087", "#3C5488", 
                      "#F39B7F", "#8491B4", "#91D1C2", "#DC0000"),
    pairwise_comparisons = NULL,
    return_stats = FALSE,
    custom_scales = NULL 
) {
  
  # Default labels if not provided
  if (is.null(biomarker_labels)) {
    biomarker_labels <- c(
      Cr = "Blood Cr (mg/dL)",
      CK = "Blood CK (U/L)",
      Ab38_csf = "Aβ38 (pg/mL)",
      Ab40_csf = "Aβ40 (pg/mL)",
      Ab42_csf = "Aβ42 (pg/mL)",
      Ab42_40_csf_bridged = "Aβ42/40 ratio",
      Ab38_40_csf = "Aβ38/40 ratio",
      Ab42_38_csf = "Aβ42/38 ratio",
      pTau181_csf = "pTau181 (pg/mL)",
      pTau217_csf = "pTau217 (pg/mL)",
      pTau181_ab40 = "pTau181/Aβ40 ratio",
      pTau217_ab40 = "pTau217/Aβ40 ratio",
      GFAP_csf_pgml = "GFAP (pg/mL)",
      NfL_csf_pgml = "NfL (pg/mL)",
      GFAP_csf_NTKbr = "GFAP (NTK-br)",
      NfL_csf_NTKbr = "NfL (NTK-br)",
      NfL_GFAP_csf = "NfL/GFAP ratio",
      GFAP_ab40_csf = "GFAP/Aβ40 ratio",
      log_Ab4240_csf = "Aβ42/40",
      log_Ab3840_csf = "Aβ38/40",
      log_pTau181_csf = "pTau181",
      log_pTau217_csf = "pTau217",
      log_GFAP_csf   = "GFAP",
      log_NfL_csf    = "NfL"
    )
  }
  
  # Reorder factor levels for specific groupings
  if (group_var == "Apos_label") {
    data <- data %>%
      mutate(
        !!sym(group_var) := factor(
          .data[[group_var]],
          levels = c(
            "ALS (Aβ-)", "DC (Aβ-)",
            "ALS (Aβ+)", "DC (Aβ+)"
          )
        )
      )
  }
  if (group_var == "Apos_label2") {
    data <- data %>%
      mutate(
        !!sym(group_var) := factor(
          .data[[group_var]],
          levels = c(
            "ALS (Aβ-)", "mimics (Aβ-)", "ND (Aβ-)",
            "ALS (Aβ+)", "mimics (Aβ+)", "ND (Aβ+)"
          )
        )
      )
  }
  
  # Prepare data for plotting
  plot_data <- data %>%
    filter(!is.na(.data[[group_var]])) %>%
    select(all_of(c(group_var, biomarkers))) %>%
    pivot_longer(
      cols = all_of(biomarkers),
      names_to = "Biomarker",
      values_to = "Value"
    ) %>%
    filter(!is.na(Value)) %>%
    mutate(
      Biomarker_label = biomarker_labels[Biomarker],
      Biomarker_label = factor(Biomarker_label, 
                               levels = biomarker_labels[biomarkers])
    )
  
  # Calculate sample sizes for each group AND biomarker combination
  group_counts_by_biomarker <- plot_data %>%
    group_by(Biomarker_label, .data[[group_var]]) %>%
    summarise(n = n(), .groups = "drop")
  
  # Filter outliers if requested
  if (trim_outliers) {
    plot_data <- plot_data %>%
      group_by(Biomarker_label) %>%
      mutate(
        lower_limit = quantile(Value, 1 - trim_percentile, na.rm = TRUE),
        upper_limit = quantile(Value, trim_percentile, na.rm = TRUE)
      ) %>%
      filter(Value >= lower_limit & Value <= upper_limit) %>%
      select(-lower_limit, -upper_limit) %>%
      ungroup()
  }
  
  # Statistical analysis if requested
  sig_data <- NULL
  stats_results <- NULL
  
  if (add_stats) {
    if (!is.null(pairwise_comparisons)) {
      
      pairwise_results <- list()
      
      for (biomarker in biomarkers) {
        biomarker_label <- biomarker_labels[biomarker]
        
        for (i in seq_along(pairwise_comparisons)) {
          comparison <- pairwise_comparisons[[i]]
          group1 <- comparison[1]
          group2 <- comparison[2]
          
          comparison_data <- data %>%
            filter(.data[[group_var]] %in% c(group1, group2)) %>%
            select(all_of(c(group_var, biomarker))) %>%
            filter(!is.na(.data[[biomarker]]))%>%
            # ADDED: Drop unused levels to prevent wilcox.test errors
            droplevels()
          
          if (nrow(comparison_data) > 0) {
            test_result <- wilcox.test(
              as.formula(paste(biomarker, "~", group_var)),
              data = comparison_data
            )
            
            pairwise_results[[length(pairwise_results) + 1]] <- data.frame(
              Biomarker = biomarker,
              Biomarker_label = biomarker_label,
              group1 = group1,
              group2 = group2,
              p.value = test_result$p.value,
              comparison_id = i
            )
          }
        }
      }
      
      pairwise_df <- bind_rows(pairwise_results)
      
      if (nrow(pairwise_df) > 0) {
        pairwise_df <- pairwise_df %>%
          mutate(
            q.value = p.adjust(p.value, method = stat_method),
            significance = case_when(
              q.value < 0.001 ~ "***",
              q.value < 0.01 ~ "**",
              q.value < 0.05 ~ "*",
              TRUE ~ ""
            )
          )
        
        stats_results <- pairwise_df
        
        pairwise_df_sig <- pairwise_df %>%
          filter(significance != "")
        
        if (nrow(pairwise_df_sig) > 0) {
          group_levels <- levels(data[[group_var]])
          if (is.null(group_levels)) {
            group_levels <- unique(data[[group_var]])
          }
          
          sig_data <- pairwise_df_sig %>%
            left_join(
              plot_data %>%
                group_by(Biomarker_label) %>%
                summarise(max_value = max(Value, na.rm = TRUE), .groups = "drop"),
              by = "Biomarker_label",
              relationship = "many-to-one"
            ) %>%
            group_by(Biomarker_label) %>%
            mutate(
              x1 = match(group1, group_levels),
              x2 = match(group2, group_levels),
              y_pos = max_value * (1.1 + 0.2 * (comparison_id - 1))
            ) %>%
            select(Biomarker_label, group1, group2, significance, x1, x2, y_pos) %>%
            distinct()
        }
      }
      
    } else {
      stat_data <- data %>%
        select(all_of(c(group_var, biomarkers)))
      
      label_list <- as.list(biomarker_labels[biomarkers])
      names(label_list) <- biomarkers
      
      stat_table <- stat_data %>%
        tbl_summary(
          by = all_of(group_var),
          missing = "no",
          label = label_list
        ) %>%
        add_p() %>%
        add_q(method = stat_method)
      
      q_values <- stat_table$table_body %>%
        select(label, q.value) %>%
        mutate(
          significance = case_when(
            q.value < 0.001 ~ "***",
            q.value < 0.01 ~ "**",
            q.value < 0.05 ~ "*",
            TRUE ~ ""
          )
        )
      
      stats_results <- q_values
      
      plot_data <- plot_data %>%
        left_join(
          q_values %>% select(label, significance),
          by = c("Biomarker_label" = "label")
        )
      
      sig_data <- plot_data %>%
        filter(significance != "") %>%
        group_by(Biomarker, Biomarker_label, significance) %>%
        summarise(
          y_pos = max(Value, na.rm = TRUE) * 1.0,
          .groups = "drop"
        ) %>%
        mutate(
          x1 = 1,
          x2 = n_distinct(data[[group_var]])
        )
    }
  }
  
  if (is.null(ncol)) {
    ncol <- length(biomarkers)
  }
  
  p <- ggplot(
    plot_data,
    aes(x = .data[[group_var]], y = Value, fill = .data[[group_var]])
  ) +
    geom_boxplot(
      width = 0.25,
      outlier.shape = if (show_points) NA else 16,  
      outlier.color = "gray10",
      outlier.size = 1,
      outlier.alpha = 0.5,
      alpha = 0.7
    )
  
  if (show_points) {
    p <- p +
      geom_jitter(
        width = 0.05,
        size = 0.125,
        color = "gray10",
        alpha = 0.3
      )
  }
  
  if (show_violin) {
    p <- p +
      geom_half_violin(
        side = "r",
        width = 0.75,
        trim = TRUE,
        alpha = 0.5,
        color = NA,
        position = position_nudge(x = 0.25)
      )
  }
  
  p <- p +
    facet_wrap(
      ~ Biomarker_label,
      scales = "free_y",
      ncol = ncol
    )
  
  if (!is.null(custom_scales)) {
    if (requireNamespace("ggh4x", quietly = TRUE)) {
      
      scale_list <- list()
      
      for (biomarker in names(custom_scales)) {
        biomarker_label <- biomarker_labels[biomarker]
        scale_info <- custom_scales[[biomarker]]
        
        if (y_log) {
          scale_list[[biomarker_label]] <- ggplot2::scale_y_log10(
            limits = scale_info$limits,
            breaks = scale_info$breaks,
            labels = scale_info$labels,
            expand = if(!is.null(scale_info$expand)) scale_info$expand else ggplot2::expansion(mult = c(0.05, 0.20))
          )
        } else {
          scale_list[[biomarker_label]] <- ggplot2::scale_y_continuous(
            limits = scale_info$limits,
            breaks = scale_info$breaks,
            labels = scale_info$labels,
            expand = if(!is.null(scale_info$expand)) scale_info$expand else ggplot2::expansion(mult = c(0.05, 0.20))
          )
        }
      }
      
      p <- p + ggh4x::facetted_pos_scales(y = scale_list)
      
    } else {
      warning("Package 'ggh4x' is required for custom_scales. Install with: install.packages('ggh4x')")
    }
  } else {
    if (y_log) {
      p <- p + ggplot2::scale_y_log10(
        expand = ggplot2::expansion(mult = c(0.05, 0.20)),
        labels = function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
      )
    } else {
      p <- p + ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0.05, 0.20))
      )
    }
  }
  
  # Create x-axis labels with biomarker-specific sample sizes
  # Create x-axis labels with biomarker-specific sample sizes
  if (show_xlab && show_n) {
    # Create a named vector of labels for each facet and group combination
    label_data <- group_counts_by_biomarker %>%
      mutate(
        label_with_n = paste0(.data[[group_var]], "\nn=", n)
      )
    
    # Create a custom labeller function
    x_labels <- label_data %>%
      select(Biomarker_label, !!sym(group_var), label_with_n) %>%
      split(.$Biomarker_label) %>%
      lapply(function(df) {
        setNames(df$label_with_n, df[[group_var]])
      })
    
    # Apply labels using ggh4x if available, otherwise use a simpler approach
    if (requireNamespace("ggh4x", quietly = TRUE)) {
      p <- p + ggh4x::facet_wrap2(
        ~ Biomarker_label,
        scales = "free_y",
        ncol = ncol,
        strip = ggh4x::strip_vanilla()
      ) +
        ggh4x::facetted_pos_scales(
          x = lapply(x_labels, function(labs) {
            scale_x_discrete(labels = labs)
          })
        )
    } else {
      # Fallback: add sample size as text below x-axis
      p <- p +
        geom_text(
          data = group_counts_by_biomarker,
          aes(x = .data[[group_var]], label = paste0("n=", n)),
          y = 0,
          vjust = 2,
          size = 4.5,
          color = "gray30",
          inherit.aes = FALSE
        )
    }
  } else if (show_xlab) {
    p <- p +
      scale_x_discrete(labels = function(x) x)
  } else {
    p <- p +
      scale_x_discrete(labels = NULL)
  }
  
  p <- p +
    scale_fill_manual(values = color_palette) +
    labs(
      title = title,
      x = NULL,
      y = if(show_ylab) "Biomarker level" else NULL
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      text = element_text(size = 18),
      axis.text = element_text(size = 16),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16,
                                 margin = ggplot2::margin(t = 5, r = 0, b = 0, l = 0)),
      axis.title = element_text(size = 20),
      axis.title.x = element_blank(),
      axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 15, b = 0, l = 0)),
      strip.text = element_text(size = 18, margin = ggplot2::margin(b = 10, t = 5)),
      plot.title = element_text(size = 22, face = "bold"),
      panel.spacing = ggplot2::unit(1.5, "lines")
    )
  
  if (!is.null(sig_data) && nrow(sig_data) > 0) {
    
    if (y_log) {
      whisker_ratio <- 0.85  
    } else {
      whisker_ratio <- 0.98
    }
    
    p <- p +
      geom_segment(
        data = sig_data,
        aes(x = x1, xend = x2, y = y_pos, yend = y_pos),
        inherit.aes = FALSE,
        linewidth = 0.5
      ) +
      geom_segment(
        data = sig_data,
        aes(x = x1, xend = x1, y = y_pos * whisker_ratio, yend = y_pos),
        inherit.aes = FALSE,
        linewidth = 0.5
      ) +
      geom_segment(
        data = sig_data,
        aes(x = x2, xend = x2, y = y_pos * whisker_ratio, yend = y_pos),
        inherit.aes = FALSE,
        linewidth = 0.5
      ) +
      geom_text(
        data = sig_data,
        aes(x = (x1 + x2) / 2, y = y_pos * 1.01, label = significance),
        inherit.aes = FALSE,
        size = 6,
        vjust = 0
      )
  }
  
  if (return_stats) {
    return(list(plot = p, statistics = stats_results))
  } else {
    return(p)
  }
}

# RUN ---------------------------------------------------------------------


## CSF ---------------------------------------------------------------------

result_bw3_ab <- plot_biomarker_comparison(
  data = T1,
  group_var = "Apos_label",
  biomarkers = c("Ab42_40_csf_bridged", "Ab38_40_csf",
                 "pTau217_csf", "pTau181_csf",
                 "GFAP_csf_pgml", "NfL_csf_pgml"),
  title = "CSF Biomarkers",
  ncol = 3,
  show_violin = FALSE,
  trim_outliers = FALSE,
  color_palette = c("#EAAE80","#80B8D8","#D55E00","#0072B2"),
  pairwise_comparisons = list(
    c("ALS (Aβ-)", "DC (Aβ-)"),
    c("ALS (Aβ+)", "DC (Aβ+)"),
    c("ALS (Aβ-)", "ALS (Aβ+)"),
    c("DC (Aβ-)", "DC (Aβ+)")
  ),
  return_stats = TRUE,
  show_n = FALSE
)

# View plot
result_bw3_ab
result_bw3_ab$plot


## Blood -------------------------------------------------------------------

result_bw3_ab_blood <- plot_biomarker_comparison(
  data = T1,
  group_var = "Apos_label",
  biomarkers = c("Cr", "CK"
  ),
  title = "Blood Biomarkers", #NULL
  ncol = 1,
  show_violin =  FALSE,
  show_ylab =FALSE,
  trim_outliers = TRUE,
  y_log = TRUE,
  #color_palette = c("#F7B2AD", "#CDEAC0", "#8C1C13", "#1B5E20"),
  color_palette = c("#EAAE80","#80B8D8","#D55E00","#0072B2"),
  pairwise_comparisons = list(
    c("ALS (Aβ-)", "DC (Aβ-)"),
    c("ALS (Aβ+)", "DC (Aβ+)"),
    c("ALS (Aβ-)", "ALS (Aβ+)"),
    c("DC (Aβ-)", "DC (Aβ+)")
  ),
  show_n=FALSE,
  return_stats = TRUE)

result_bw3_ab_blood
result_bw3_ab_blood$plot


# 📊 Figure 1 -----------------------------------------------------------------
# FDR correction is applied across all pairwise tests shown in this figure.

Fig1 <- (result_bw3_ab$plot|result_bw3_ab_blood$plot) +
  plot_layout(widths = c(3, 1))

Fig1
