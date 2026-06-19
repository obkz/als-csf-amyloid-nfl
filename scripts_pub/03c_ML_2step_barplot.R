library(cowplot)
library(grid)
library(gridExtra)

source("scripts_pub/03b_ML_2step_oof.R")


# Function ----------------------------------------------------------------

# Modified plot function with right-aligned y-axis labels
plot_waterfall_barplot_fixed <- function(data, 
                                         als_color = "#E07A2D",
                                         dc_color  = "#4FA3C7",
                                         title = NULL,
                                         subtitle = NULL,
                                         x_max = 440,
                                         bar_height = 0.2) {
  data <- data %>% mutate(group = factor(group, levels = rev(group))) %>%
    mutate(label = if("label" %in% names(.)) label else sprintf("ALS: %d (%.0f%%), DC: %d (%.0f%%)",
                                                                n_als, (n_als/n_total)*100,
                                                                n_dc,  (n_dc/n_total)*100))
  data_long <- data %>%
    select(stage, group, n_als, n_dc, n_total, label) %>%
    pivot_longer(cols = c(n_als, n_dc),
                 names_to = "class", values_to = "n") %>%
    mutate(class = factor(class, levels = c("n_dc","n_als"),
                          labels = c("DC","ALS")))
  
  ggplot(data_long, aes(x = n, y = group, fill = class)) +
    geom_col(position = "stack", width = bar_height, color = "gray99", linewidth = 0.2) +
    geom_text(data = data,
              aes(x = 20, y = group, label = label),
              inherit.aes = FALSE, hjust = 0, vjust = 2.5,
              size = 4, color = "gray20") +
    scale_fill_manual(
      values = c("DC" = dc_color, "ALS" = als_color),
      breaks = c("ALS", "DC"),
      limits = c("ALS", "DC"),
      drop   = FALSE
    ) +
    scale_x_continuous(limits = c(0, x_max),
                       breaks = seq(0, x_max, 100),
                       expand = c(0, 0)) +
    scale_y_discrete() +
    coord_cartesian(clip = "off") +
    labs(title = title, 
         subtitle = subtitle,
         x = "Number of patients", y = NULL, fill = NULL) +
    theme_classic(base_size = 16) +
    theme(axis.text.y = element_text(hjust = 1, size = 13, lineheight = 1.025),  # Changed hjust = 0 to hjust = 1
          axis.text.x = element_text(size = 13),
          axis.title.x = element_text(size = 14, margin = ggplot2::margin(t = 10)),
          legend.position = "bottom",
          plot.subtitle = element_text(size = 14, hjust = 0, margin = ggplot2::margin(b = 5)),
          legend.text = element_text(size = 12),
          plot.margin = ggplot2::margin(10,10,0,10),
          plot.title = element_text(margin = ggplot2::margin(b = 5)))
}

# Modified plot_all_bar with right-aligned y-axis labels
plot_all_bar_fixed <- function(data,
                               als_color = "#E07A2D",
                               dc_color  = "#4FA3C7",
                               x_max = 440,
                               bar_height = 0.2) {
  data <- data %>% mutate(group = factor(group, levels = group)) %>%
    mutate(label = sprintf("ALS: %d (%.0f%%), DC: %d (%.0f%%)",
                           n_als, (n_als/n_total)*100,
                           n_dc,  (n_dc/n_total)*100))
  data_long <- data %>%
    select(group, n_als, n_dc, label) %>%
    pivot_longer(cols = c(n_als, n_dc),
                 names_to = "class", values_to = "n") %>%
    mutate(class = factor(class, levels = c("n_dc","n_als"),
                          labels = c("DC","ALS")))
  
  ggplot(data_long, aes(x = n, y = group, fill = class)) +
    geom_col(position = "stack", width = bar_height, color = "gray99", linewidth = 0.2) +
    geom_text(data = data,
              aes(x = 20, y = group, label = label),
              inherit.aes = FALSE, hjust = 0, vjust = 2.5,
              size = 4, color = "gray20") +
    scale_fill_manual(
      values = c("DC" = dc_color, "ALS" = als_color),
      breaks = c("ALS", "DC"),
      limits = c("ALS", "DC"),
      drop   = FALSE,
      guide  = "none"
    ) +
    scale_x_continuous(limits = c(0, x_max),
                       breaks = seq(0, x_max, 100),
                       expand = c(0, 0)) +
    labs(x = NULL, y = NULL, fill = NULL) +
    theme_classic(base_size = 16) +
    theme(axis.text.y = element_text(hjust = 1, size = 13, lineheight = 1.05),  
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          plot.margin = ggplot2::margin(5,10,0,10))
}


# Create all plots with adjusted bar heights --------------------------------

## Step 1 -------------------------------------------------

p_nfl_sp95_all_fixed <- plot_all_bar_fixed(nfl_sp95_data$all, x_max = 440, bar_height = 0.125)

p_nfl_sp95_class_fixed <- plot_waterfall_barplot_fixed(
  nfl_sp95_data$classification,
  title = "NfL cut-off at 95% Specificity",
  subtitle = "PPV 91%, NPV 44%",
  x_max = 440,
  bar_height = 0.2
) + theme(legend.position = "none")

p_nfl_se95_class_fixed <- plot_waterfall_barplot_fixed(
  nfl_se95_data$classification,
  title = "NfL cut-off at 95% Sensitivity",
  subtitle = "PPV 76%, NPV 86%",
  x_max = 440,
  bar_height = 0.2
) + theme(legend.position = "none")

p_step1_buckets_fixed <- plot_waterfall_barplot_fixed(
  step1_data$buckets,
  title = "Step 1: NfL Classification",
  x_max = 440,
  bar_height = 0.3
) + theme(legend.position = "none")

step2_plot_fixed <- plot_waterfall_barplot_fixed(
  step2_data,
  title = "Step 2: Gray-zone Resolution (RF)",
  x_max = 440,
  bar_height = 0.2
) + theme(legend.position = "none")

p_overall_class_fixed <- plot_waterfall_barplot_fixed(
  overall_data$classification,
  title = "Overall",
  x_max = 440,
  bar_height = 0.3
) + theme(legend.position = "none")

# Title plots with reduced margins
# Title plots with tags
p_onestep_title <- ggplot() + 
  labs(title = "One-step dichotomization using a single NfL cut-off",
       tag = "A") +
  theme_void() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0, 
                                  margin = ggplot2::margin(t = 15, b = 5, l = -2.5)),
        plot.tag = element_text(size = 22, face = "bold", hjust = 0, vjust = 1),
        plot.tag.position = c(-0.15, 0.9),
        plot.margin = ggplot2::margin(0, 0, 0, 0))

p_twostep_title <- ggplot() + 
  labs(title = "Two-step framework using two NfL cut-offs",
       tag = "B") +
  theme_void() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0, 
                                  margin = ggplot2::margin(t = 15, b = 5, l = -2.5)),
        plot.tag = element_text(size = 22, face = "bold", hjust = 0, vjust = 1),
        plot.tag.position = c(-0.15, 0.9),
        plot.margin = ggplot2::margin(0, 0, 0, 0))


## Step 2 -------------------------------------------------

# Create Step2 dataset with 4 rows

step2_merged <- bind_rows(
  # Row 1: Step 1 Rule-in (95% Sp. high NfL)
  step1_data$buckets %>% dplyr::slice(1) %>%
    mutate(group = sprintf("Step 1 ruled-in\n(n=%d)", n_total)),
  # Row 2: Additional tests (+) Rule-in from gray-zone
  step2_data %>% dplyr::slice(1),
  # Row 3: Additional tests (-) Still in gray-zone
  step2_data %>% dplyr::slice(2),
  # Row 4: Step 1 Rule-out (95% Se. low NfL)
  step1_data$buckets %>% dplyr::slice(3) %>%
    mutate(group = sprintf("Step 1 ruled-out\n(n=%d)", n_total))) %>%
  mutate(
    stage = factor("Step 2"),
    group = factor(group, levels = group)
  )

step2_merged


# Simply pre-process the data with custom labels, then use existing function
step2_merged_withlabel <- step2_merged %>%
  mutate(
    label = case_when(
      row_number() %in% c(1, 4) ~ "",  # No text for Step 1 rows
      TRUE ~ sprintf("ALS: %d (%.0f%%), DC: %d (%.0f%%)",
                     n_als, (n_als/n_total)*100,
                     n_dc,  (n_dc/n_total)*100)
    )
  )

# Use existing function - it will use the pre-calculated 'label' column
step2_plot_fixed <- plot_waterfall_barplot_fixed(
  step2_merged_withlabel,
  title = "Step 2: Gray-zone Resolution (RF)",
  x_max = 440,
  bar_height = 0.4
) + theme(legend.position = "none")

step2_plot_fixed


# Combine plots into main figure ------------------------------------------

col1 <- plot_grid(
  p_onestep_title, p_nfl_sp95_all_fixed, p_nfl_sp95_class_fixed,
  NULL,
  p_twostep_title, p_step1_buckets_fixed,
  ncol = 1, rel_heights = c(0.35, 1, 2, 0.15, 0.35, 2),
  align = "v", axis = "lr"
)

col2 <- plot_grid(
  NULL, NULL, p_nfl_se95_class_fixed,
  NULL,
  NULL, step2_plot_fixed,
  ncol = 1, rel_heights = c(0.35, 1, 2, 0.15, 0.35, 2),
  align = "v", axis = "lr"
)

col3 <- plot_grid(
  NULL, NULL, NULL,
  NULL,
  NULL, p_overall_class_fixed,
  ncol = 1, rel_heights = c(0.35, 1, 2, 0.15, 0.35, 2),
  align = "v", axis = "lr"
)

# Combine columns
main_plot <- plot_grid(
  col1, col2, col3,
  ncol = 3, rel_widths = c(0.9, 1, 1),
  align = "h", axis = "tb"
)




# add fake legend ---------------------------------------------------------

make_fake_legend <- function(
    als_color = "#E07A2D",
    dc_color  = "#4FA3C7",
    label_als = "ALS",
    label_dc  = "DC",
    box_w = 0.22,
    box_h = 0.1
){
  leg_df <- data.frame(
    class = factor(c(label_als, label_dc), levels = c(label_als, label_dc)),
    y = c(2, 1)
  )
  
  leg <- ggplot(leg_df, aes(x = 1, y = y)) +
    geom_rect(aes(xmin = 0.00, xmax = 0.18, ymin = y-0.2, ymax = y+0.2, fill = class),
              color = "gray20", linewidth = 0.2) +
    geom_text(aes(x = 0.22, label = class), hjust = 0, size = 4) +
    scale_fill_manual(values = setNames(c(als_color, dc_color), c(label_als, label_dc))) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.5, 2.5), clip = "off") +
    theme_void() +
    theme(legend.position = "none",
          plot.margin = ggplot2::margin(2, 2, 2, 2))
  
  list(plot = leg, w = box_w, h = box_h)
}



legend_bg <- ggplot() +
  geom_rect(aes(xmin=0, xmax=1, ymin=0, ymax=1),
            fill="white", color="gray70", linewidth=0.4) +
  theme_void() +
  theme(plot.margin = ggplot2::margin(0,0,0,0))


fake_legend <- make_fake_legend()
fake_legend$plot


# 📊 Figure 3AB --------------------------------------------------------------

Fig3AB <- cowplot::ggdraw(main_plot) +
  cowplot::draw_plot(legend_bg, x=0.83, y=0.8, width=0.12, height=0.105) +
  cowplot::draw_plot(fake_legend$plot, x = 0.85, y = 0.8, width = fake_legend$w, height = fake_legend$h)

Fig3AB
