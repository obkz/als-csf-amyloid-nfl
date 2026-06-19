library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(tidyr)

# Configuration Section
# Specify the biomarker to use for tertile stratification
KM_FACET_BIOMARKER <- "tert_cNfL"
KM_FACET_BIOMARKER_LABEL <- "CSF NfL"

# Timepoint for survival probability extraction (months)
SURVIVAL_TIMEPOINT <- 24

# Data preparation
km_facet_survALS_ab <- survALS %>%
  mutate(
    facet_group = case_when(
      Ab_status == "negative" ~ "Aβ (–) ALS", #"ALS Aβ-",
      Ab_status == "positive" ~ "Aβ (+) ALS", #"ALS Aβ+",
      TRUE ~ NA_character_
    )
  )

km_facet_survALS_all <- survALS %>%
  mutate(facet_group = "Overall ALS")

km_facet_survALS_combined <- bind_rows(km_facet_survALS_all, km_facet_survALS_ab) %>%
  filter(!is.na(facet_group)) %>%
  mutate(
    surv_time_months = surv_time_days / 30.44,
    biomarker_tert_raw = case_when(
      .data[[KM_FACET_BIOMARKER]] == "low" ~ "low",
      .data[[KM_FACET_BIOMARKER]] == "mid" ~ "middle",
      .data[[KM_FACET_BIOMARKER]] == "middle" ~ "middle",
      .data[[KM_FACET_BIOMARKER]] == "high" ~ "high",
      TRUE ~ as.character(.data[[KM_FACET_BIOMARKER]])
    ),
    biomarker_tert_label = factor(
      biomarker_tert_raw,
      levels = c("low", "middle", "high"),
      labels = c("Low", "Middle", "High")
    )
  )

# Calculate tertile cutoff values for the selected biomarker
# Map tertile column names to their corresponding continuous biomarker columns
biomarker_mapping <- c(
  "tert_cNfL" = "NfL_csf_pgml",
  "tert_cGFAP" = "GFAP_csf_pgml",
  "tert_pTau217" = "pTau217_csf_pgml",
  "tert_pTau181" = "pTau181_csf_pgml",
  "tert_Ab38" = "Ab38_csf_pgml",
  "tert_Ab40" = "Ab40_csf_pgml",
  "tert_Ab42" = "Ab42_csf_pgml",
  "tert_Ab4240" = "Ab4240_ratio"
)

# Get the continuous biomarker column name
biomarker_col <- biomarker_mapping[KM_FACET_BIOMARKER]

# Calculate tertile cutoffs from the original continuous biomarker values
tertile_cutoffs <- survALS %>%
  summarise(
    low_max = quantile(.data[[biomarker_col]], 1/3, na.rm = TRUE),
    middle_max = quantile(.data[[biomarker_col]], 2/3, na.rm = TRUE)
  )

# Create tertile range labels with integer values and line breaks
tertile_range_low <- sprintf("(≤%.0f)", tertile_cutoffs$low_max)
tertile_range_middle <- sprintf("(%.0f-%.0f)", tertile_cutoffs$low_max, tertile_cutoffs$middle_max)
tertile_range_high <- sprintf("(>%.0f)", tertile_cutoffs$middle_max)
tertile_axis_labels <- c(
  "High"   = paste0("High\n", tertile_range_high),
  "Middle" = paste0("Middle\n", tertile_range_middle),
  "Low"    = paste0("Low\n", tertile_range_low)
)

# Function to extract survival probability at specific timepoint
extract_surv_prob <- function(data, group_var, time_point) {
  fit <- survfit(Surv(surv_time_months, surv_time_status == 1) ~ get(group_var), data = data)
  summary_fit <- summary(fit, times = time_point)
  
  # Extract survival probabilities and confidence intervals
  result <- data.frame(
    surv_prob = summary_fit$surv,
    lower_ci = summary_fit$lower,
    upper_ci = summary_fit$upper,
    strata = summary_fit$strata
  )
  
  # Clean strata names
  result$group <- gsub("get\\(group_var\\)=", "", result$strata)
  
  return(result)
}

# Extract 24-month survival probabilities for each facet group with sample sizes
surv_24m_results <- km_facet_survALS_combined %>%
  group_by(facet_group) %>%
  group_modify(~ {
    result_df <- extract_surv_prob(.x, "biomarker_tert_label", SURVIVAL_TIMEPOINT)
    
    # Calculate sample size for each tertile
    n_counts <- .x %>%
      group_by(biomarker_tert_label) %>%
      summarise(n = n(), .groups = "drop") %>%
      rename(group = biomarker_tert_label) %>%
      mutate(group = as.character(group))
    
    result_df <- result_df %>%
      left_join(n_counts, by = "group")
    
    return(result_df)
  }) %>%
  ungroup() %>%
  mutate(
    tertile = factor(group, levels = c("Low", "Middle", "High")),
    facet_group = factor(facet_group, levels = c("Overall ALS", "Aβ (–) ALS", "Aβ (+) ALS"))
  )


# Visualize ---------------------------------------------------------------

# Configuration for p3 heatmap display
# Set to TRUE to show confidence intervals, FALSE to show sample sizes
HEATMAP_SHOW_CI <- TRUE

# Visualization 3: Heatmap with sample sizes or confidence intervals
surv_24m_heatmap <- surv_24m_results %>%
  mutate(
    label_with_n = sprintf("%.2f\n(n=%d)", surv_prob, n),
    label_with_ci = sprintf("%.2f\n(%.2f-%.2f)", surv_prob, lower_ci, upper_ci),
    # Reverse the order: High -> Middle -> Low (top to bottom)
    tertile_reversed = factor(tertile, levels = c("High", "Middle", "Low"))
  )

if (HEATMAP_SHOW_CI) {
  surv_24m_heatmap <- surv_24m_heatmap %>%
    mutate(label_display = label_with_ci)
  subtitle_text <- "Values shown as: Probability (95% CI)"
} else {
  surv_24m_heatmap <- surv_24m_heatmap %>%
    mutate(label_display = label_with_n)
  subtitle_text <- "Values shown as: Probability (sample size)"
}

Fig5B_pub <- ggplot(surv_24m_heatmap, aes(x = facet_group, y = tertile_reversed, fill = surv_prob)) +
  geom_tile(color = "white", linewidth = 1.5) +
  geom_text(aes(label = label_display), 
            size = 5.25, fontface = "bold", color = "black", lineheight = 0.9) +
  scale_y_discrete(labels = tertile_axis_labels) +
  scale_fill_gradient2(
    low  = "#F0A65A",
    mid  = "#FAFAFA",
    high = "#88C4BB",
    midpoint = 0.5,
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1)
  )+
  labs(
    title = "24-Month Survival Probability", 
    #subtitle = subtitle_text,
    x = "Patient Subgroup (Aβ status)",
    y = paste(KM_FACET_BIOMARKER_LABEL, "Tertile"),
    fill = "Survival\nProbability"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = ggplot2::margin(t = 10)),
    plot.subtitle = element_text(size = 12, hjust = 0.5, margin = ggplot2::margin(b = 1)),
    plot.margin = ggplot2::margin(t = 0, r = 10, b = 10, l = 10),
    axis.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(margin = ggplot2::margin(t = 20)),
    axis.title.y = element_text(margin = ggplot2::margin(r = 20)),
    axis.text = element_text(size = 14, color = "black", margin = ggplot2::margin(t = 5, b = 15)),
    axis.text.x = element_text(size = 16, hjust = 0.5, angle = 0),
    axis.text.y = element_text(size = 16, hjust = 0.5, angle=90),
    legend.title = element_text(size = 14, margin = ggplot2::margin(b = 10)),
    legend.text = element_text(size = 12),
    panel.grid = element_blank()
  )

Fig5B_pub
