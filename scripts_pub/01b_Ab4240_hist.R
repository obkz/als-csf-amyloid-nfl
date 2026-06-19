library(dplyr)
library(ggplot2)
library(patchwork)

AB_cutoff_value <- 0.077

Ab_df_hist <- csfdf %>%
  select(Ab42_40_csf_bridged, Ab_status, Apos_label, ALS_label) %>%
  filter(is.finite(Ab42_40_csf_bridged)) %>%
  mutate(
    ALS_label = recode_factor(ALS_label,
                              `1` = "ALS",
                              `0` = "DC"),
    Apos_label = recode_factor(Apos_label,
                               "DC_Apos"  = "DC (Aβ+)",
                               "DC_Aneg"  = "DC (Aβ-)",
                               "ALS_Apos" = "ALS (Aβ+)",
                               "ALS_Aneg" = "ALS (Aβ-)"),
    Ab_status = recode_factor(Ab_status,
                              "positive" = "Aβ (+)" ,
                              "negative" = "Aβ (-)")
  )

base_theme_his <- theme_bw(base_size = 16) +
  theme(
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    plot.margin = ggplot2::margin(t = 14, r = 14, b = 14, l = 14),
    plot.title    = element_text(size=20, face = "bold", 
                                 margin = ggplot2::margin(b = 10)),
    plot.subtitle = element_text(margin = ggplot2::margin(b = 8)),
    plot.caption  = element_text(margin = ggplot2::margin(t = 10)),
    axis.title.x = element_text(margin = ggplot2::margin(t = 14)),
    axis.title.y = element_text(margin = ggplot2::margin(r = 10)),
    axis.text.x  = element_text(margin = ggplot2::margin(t = 6)),
    axis.text.y  = element_text(margin = ggplot2::margin(r = 6)),
    strip.text = element_text(#face = "bold", 
                              margin = ggplot2::margin(t = 6, b = 6), size = 16),
    panel.spacing = unit(12, "pt"),
    legend.box.margin = ggplot2::margin(b = 6)
  )

# A) Distribution of CSF Aβ42/40 by Aβ status (overall)
SuppleFig_Ab_h2 <- ggplot(Ab_df_hist, aes(x = Ab42_40_csf_bridged, fill = Ab_status)) +
  geom_histogram(bins = 40, position = "identity", alpha = 0.5, color = "black") +
  geom_vline(xintercept = AB_cutoff_value, linetype = "dashed", color = "darkred", linewidth = 1) +
  scale_fill_manual(
    values = c(
      "Aβ (-)" = "gray90",
      "Aβ (+)" = "red4"
    )
  ) +
  labs(
    title = "(A) Distribution of CSF Aβ42/40 by Aβ status (overall)",
    x = "Aβ42/40 (CSF)",
    y = "Count",
    fill = "Ab status"
  ) + 
  ylim(0,50) +
  base_theme_his +  theme(aspect.ratio = 1) 

# B) Distribution of CSF Aβ42/40 by Aβ status × ALS diagnosis
SuppleFig_Ab_h4 <- ggplot(Ab_df_hist, aes(x = Ab42_40_csf_bridged, fill = Apos_label)) +
  geom_histogram(bins = 40, position = "identity", alpha = 0.5, color = "black") +
  geom_vline(xintercept = AB_cutoff_value, linetype = "dashed", color = "darkred", linewidth = 1) +
  facet_wrap(~ALS_label, scales = "fixed") +
  scale_fill_manual(
    values = c(
      "ALS (Aβ-)" = "#F7DDBF",  # very light orange
      "DC (Aβ-)"  = "#D6E8F5",  # very light blue
      "ALS (Aβ+)" = "#D55E00",  # dark orange-red
      "DC (Aβ+)"  = "#0072B2"   # dark blue
    )
  ) +
  ylim(0,50) +
  labs(
    title = "(B) Distribution of CSF Aβ42/40 by Aβ status × ALS diagnosis",
    x = "Aβ42/40 (CSF)",
    y = "Count",
    fill = "Group"
  ) + base_theme_his + theme(aspect.ratio = 1)


# Display
SuppleFig1A <- SuppleFig_Ab_h2
SuppleFig1B <- SuppleFig_Ab_h4


# 🔍 Supple Fig.1 -------------------------------------------------------------

SuppleFig1A
SuppleFig1B
