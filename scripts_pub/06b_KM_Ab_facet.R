library(survival)
library(survminer)
library(ggplot2)
library(dplyr)

# NfL tertiles were defined in the full ALS survival-analysis cohort.
# Kaplan–Meier curves were then plotted separately within each Aβ subgroup.

# Function to create individual KM plot with risk table ------------------
create_km_facet_plot <- function(data, title, show_event_info = TRUE) {
  
  # Calculate total sample size
  n_total <- nrow(data)
  # Calculate sample size for each tertile
  n_by_group <- data %>%
    group_by(biomarker_tert_label) %>%
    summarise(n = n(), .groups = "drop")
  legend_labels <- paste0(c("Low", "Middle", "High"), " (n=", n_by_group$n, ")")
  # N events
  n_event <- sum(data$surv_time_status == 1, na.rm = TRUE)
  
  km_facet_fit <- survfit(Surv(surv_time_months, surv_time_status == 1) ~ biomarker_tert_label, data = data)
  
  p <- ggsurvplot(
    km_facet_fit,
    data = data,
    pval = TRUE,
    pval.method = TRUE,
    pval.method.coord = c(KM_FACET_PVAL_METHOD_X, KM_FACET_PVAL_METHOD_Y),
    pval.coord = c(KM_FACET_PVAL_X, KM_FACET_PVAL_Y),
    xlab = "Time (months)",
    xlim = c(0, KM_FACET_X_MAX),
    break.x.by = KM_FACET_X_BREAK,
    ylab = "Survival Probability",
    ggtheme = theme_classic(),
    legend.title = paste(KM_FACET_BIOMARKER_LABEL, "Tertile"),
    legend.labs = c("low", "mid", "high"),
    title = title,
    palette = km_facet_nature_colors,
    legend = "bottom",
    risk.table = TRUE,
    risk.table.title = "Number at risk",
    risk.table.y.text = FALSE,
    risk.table.height = 0.25,
    tables.theme = theme_cleantable(),
    conf.int.alpha = CONF_INT_ALPHA,
    conf.int = FALSE  # 95%CI
  )
  
  # Add event information inside plot if requested
  if (show_event_info) {
    p$plot <- p$plot +
      annotate(
        "text",
        x = KM_FACET_X_MAX * 0.6,
        y = 0.91,
        label = paste0("Events / Total: ", n_event, " / ", n_total),
        size = 5,
        hjust = 0,
        color = "black"
      )
  }
  
  return(p)
}

# Customize each plot with left and bottom axes only
customize_km_facet_plot <- function(plot_obj) {
  plot_obj$plot <- plot_obj$plot +
    geom_vline(
      xintercept = 24, 
      linetype = "dashed", 
      color = "gray50", 
      linewidth = 0.6
    ) +
    theme(
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      axis.title.x = element_text(size = 16, face = "bold", margin = ggplot2::margin(t = 15, unit = "pt")),
      axis.title.y = element_text(size = 16, face = "bold", margin = ggplot2::margin(r = 10, unit = "pt")),
      axis.text = element_text(size = 14, color = "black"),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 1, face = "bold", margin = ggplot2::margin(t = 10)),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_line(color = "black", linewidth = 0.5),
      panel.border = element_blank(),
      axis.ticks = element_line(color = "black")
    )
  
  # Customize risk table appearance
  plot_obj$table <- plot_obj$table +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(size = 14, face = "bold"),
      plot.margin = margin(l = 20, r = 5, t = 5, b = 5, unit = "pt")
    )
  
  return(plot_obj)
}


# Configuration Section =====

# Parameters
CONF_INT_ALPHA <- 0.15

# Specify the biomarker to use for tertile stratification
KM_FACET_BIOMARKER <- "tert_cNfL"  

# Specify the display name for the biomarker
KM_FACET_BIOMARKER_LABEL <- "CSF NfL"  

# Whether to show event information inside the plots
SHOW_EVENT_INFO <- FALSE

# Color palette for the KM plots
km_facet_nature_colors <- c("#009E73","#56B4E9","#D55E00")

# P-value position configuration
KM_FACET_PVAL_METHOD_X <- 60
KM_FACET_PVAL_METHOD_Y <- 0.84
KM_FACET_PVAL_X <- 78
KM_FACET_PVAL_Y <- 0.84

# X-axis configuration
KM_FACET_X_MAX <- 100
KM_FACET_X_BREAK <- 12

# Subgroup variable for faceting
KM_FACET_SUBGROUP_VAR      <- "Ab_status"   
KM_FACET_SUBGROUP_POS      <- "positive"
KM_FACET_SUBGROUP_NEG      <- "negative"
KM_FACET_SUBGROUP_POS_DISP <- "Aβ+"
KM_FACET_SUBGROUP_NEG_DISP <- "Aβ−"


# Data preparation --------------------------------------------------------

km_facet_survALS_ab <- survALS %>%
  mutate(
    facet_group = case_when(
      .data[[KM_FACET_SUBGROUP_VAR]] == KM_FACET_SUBGROUP_NEG ~ paste("ALS", KM_FACET_SUBGROUP_NEG_DISP),
      .data[[KM_FACET_SUBGROUP_VAR]] == KM_FACET_SUBGROUP_POS ~ paste("ALS", KM_FACET_SUBGROUP_POS_DISP),
      TRUE ~ NA_character_
    )
  )


# Create a dataset for "ALS all" facet
km_facet_survALS_all <- survALS %>%
  mutate(facet_group = "ALS all")

# Combine the datasets and create better labels for tertiles
km_facet_survALS_combined <- bind_rows(km_facet_survALS_all, km_facet_survALS_ab) %>%
  filter(!is.na(facet_group)) %>%
  mutate(
    surv_time_months = surv_time_days / 30.44
  )

# Create a generic tertile label column based on the selected biomarker
km_facet_survALS_combined <- km_facet_survALS_combined %>%
  mutate(
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



# Create plots for each facet group -------------------------------------------------------

# km_facet_plot_all <- create_km_facet_plot(
#   km_facet_survALS_combined %>% filter(facet_group == "ALS all"),
#   show_event_info = SHOW_EVENT_INFO,
#   "All ALS Patients"
# )

km_facet_plot_ab_neg <- create_km_facet_plot(
  km_facet_survALS_combined %>% filter(facet_group == paste("ALS", KM_FACET_SUBGROUP_NEG_DISP)),
  show_event_info = SHOW_EVENT_INFO,
  title = paste0(KM_FACET_SUBGROUP_NEG_DISP, " ALS")
)

km_facet_plot_ab_pos <- create_km_facet_plot(
  km_facet_survALS_combined %>% filter(facet_group == paste("ALS", KM_FACET_SUBGROUP_POS_DISP)),
  show_event_info = SHOW_EVENT_INFO,
  title = paste0(KM_FACET_SUBGROUP_POS_DISP, " ALS")
)


# Apply customization to all plots ----------------------------------------

km_facet_plot_ab_neg <- customize_km_facet_plot(km_facet_plot_ab_neg)
km_facet_plot_ab_pos <- customize_km_facet_plot(km_facet_plot_ab_pos)

km_facet_plot_ab_neg$plot <- km_facet_plot_ab_neg$plot +
  labs(tag = "D") +
  theme(plot.tag = element_text(face = "bold", size = 22, hjust = 0, vjust = 0))

km_facet_plot_ab_pos$plot <- km_facet_plot_ab_pos$plot +
  labs(tag = "E") +
  theme(plot.tag = element_text(face = "bold", size = 22, hjust = 0, vjust = 0))


# 📊 Figure 5DE -----------------------------------------------------------

Fig5DE <- survminer::arrange_ggsurvplots(
  list(km_facet_plot_ab_neg, km_facet_plot_ab_pos),
  ncol = 2, nrow = 1
)

Fig5DE
