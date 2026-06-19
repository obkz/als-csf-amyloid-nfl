library(ggplotify)

# ---- Validation cohort: data preparation for KM ----
# KM tertiles are cohort-specific tertiles of NfL, calculated separately within each cohort.

ext_km <- ext |>
  filter(survival_days > 0) |>
  mutate(
    surv_time_months     = survival_days / 30.44,
    surv_time_status     = as.integer(survival_status == 1),
    ab_status_km         = if_else(ab42_40_csf_ratio < 0.068, "positive", "negative"),
    nfl_z_km             = as.numeric(scale(log(nfl_csf_pgml))),
    tert_nfl_km          = dplyr::ntile(nfl_z_km, 3),
    biomarker_tert_label = factor(
      tert_nfl_km,
      levels = c(1, 2, 3),
      labels = c("Low", "Middle", "High")
    )
  ) |>
  filter(!is.na(ab_status_km), !is.na(biomarker_tert_label))

cat("Bologna KM - N:", nrow(ext_km),
    "| Events:", sum(ext_km$surv_time_status),
    "| Aβ+:", sum(ext_km$ab_status_km == "positive"),
    "| Aβ-:", sum(ext_km$ab_status_km == "negative"), "\n")


# ---- Bologna KM plots using same function ----
km_ext_plot_neg <- create_km_facet_plot(
  ext_km |> filter(ab_status_km == "negative"),
  show_event_info = SHOW_EVENT_INFO,
  title = paste0(KM_FACET_SUBGROUP_NEG_DISP, " ALS (Bologna)")
)

km_ext_plot_pos <- create_km_facet_plot(
  ext_km |> filter(ab_status_km == "positive"),
  show_event_info = SHOW_EVENT_INFO,
  title = paste0(KM_FACET_SUBGROUP_POS_DISP, " ALS (Bologna)")
)

# ---- Apply customization ----
km_ext_plot_neg <- customize_km_facet_plot(km_ext_plot_neg)
km_ext_plot_pos <- customize_km_facet_plot(km_ext_plot_pos)

# ---- Panel labels ----
km_facet_plot_ab_neg$plot <- km_facet_plot_ab_neg$plot +
  labs(tag = "A",
       title = "Aβ-negative ALS (Nagoya)",
       subtitle = NULL) +
  theme(plot.tag = element_text(face = "bold", size = 22))

km_facet_plot_ab_pos$plot <- km_facet_plot_ab_pos$plot +
  labs(tag = "B",
       title = "Aβ-positive ALS (Nagoya)"
  ) +
  theme(plot.tag = element_text(face = "bold", size = 22))

km_ext_plot_neg$plot <- km_ext_plot_neg$plot +
  labs(tag = "C",
       title = "Aβ-negative ALS (Bologna)",
  ) +
  theme(plot.tag = element_text(face = "bold", size = 22))

km_ext_plot_pos$plot <- km_ext_plot_pos$plot +
  labs(tag = "D",
       title= "Aβ-positive ALS (Bologna)",
  ) +
  theme(plot.tag = element_text(face = "bold", size = 22))



grob_int_neg <- ggplotify::as.grob(function() print(km_facet_plot_ab_neg))
grob_int_pos <- ggplotify::as.grob(function() print(km_facet_plot_ab_pos))
grob_ext_neg <- ggplotify::as.grob(function() print(km_ext_plot_neg))
grob_ext_pos <- ggplotify::as.grob(function() print(km_ext_plot_pos))

SuppleFig12_added <- gridExtra::grid.arrange(
  grob_int_neg, grob_int_pos,
  grob_ext_neg, grob_ext_pos,
  ncol = 2, nrow = 2
)

cat("Median survival (months) Nagoya:", med_surv_int, "\n")
cat("Median survival (months) Bologna:", med_surv_ext, "\n")


# 🔍 SuppleFig12 ------------------------------------------------------------------

SuppleFig12 <- SuppleFig12_added

