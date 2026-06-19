library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(purrr)
library(pwr)
library(patchwork)
library(grid)
library(ggsci)


# ────────────────────────────────────────────────
# Load & prepare data -----------------------------
# ────────────────────────────────────────────────
# Assumes df_slope_model exists in the environment and contains:
#   - numeric 'slope_total' (ALSFRS-R slope per patient)
#   - numeric 'cNfL' (for quantile-based enrichment)


df0 <- df_slope_model %>%
  dplyr::filter(
    !is.na(slope_total),
    !is.na(NfL_csf_pgml),
    !R_ElEscorial %in% c("suspected")
  ) %>%
  dplyr::mutate(
    quartile_cNfL = dplyr::ntile(NfL_csf_pgml, 4),
    quartile_cNfL = factor(
      quartile_cNfL,
      levels = 1:4,
      labels = c("Q1 (lowest)", "Q2", "Q3", "Q4 (highest)")
    )
  )

cohorts <- list(
  "All-comers"    = df0,
  "NfL ≥ median" = df0 %>%
    dplyr::filter(quartile_cNfL %in% c("Q3", "Q4 (highest)")),
  "NfL Q4"       = df0 %>%
    dplyr::filter(quartile_cNfL == "Q4 (highest)")
)


# ────────────────────────────────────────────────
# Sample size calculation function ---------------
# ────────────────────────────────────────────────
# Power-based required total sample size vs. relative effect X%.
# Assumptions:
#   - Treatment reduces the mean slope by X% of |mu_control|.
#   - Two-sample t-test, two-sided alpha=0.05, power=0.80, equal allocation.
#   - Cohen's d = |delta| / SD_control, computed from the cohort's mean/SD.
#   - Optional small-trial bump: if n_per_arm < 25, add +2 to total (≈ +1 per arm).

calc_n_vs_X <- function(df, X_seq = seq(10, 60, by = 2),
                        alpha = 0.05, power = 0.80,
                        add_small_trial_bump = FALSE
                        ) {
  mu  <- mean(df$slope_total, na.rm = TRUE)
  sd0 <- sd(df$slope_total,   na.rm = TRUE)
  
  # Guard against degenerate inputs
  if (!is.finite(mu) || !is.finite(sd0) || sd0 <= 0) {
    stop("Mean/SD not finite or SD <= 0 in the selected cohort.")
  }
  
  tibble(X = X_seq) %>%
    mutate(
      # Absolute difference in mean slope we aim to detect
      delta_abs = (X/100) * abs(mu),
      # Convert to Cohen's d for power calculation
      d_cohen   = delta_abs / sd0,
      # Required n per arm from pwr.t.test
      n_per_arm = purrr::map_dbl(
        d_cohen,
        ~ pwr.t.test(d = .x, sig.level = alpha, power = power,
                     type = "two.sample", alternative = "two.sided")$n
      ),
      # Round up to the next integer
      n_per_arm = ceiling(n_per_arm),
      n_total   = 2 * n_per_arm
    ) %>%
    mutate(
      # Add a small bump for very small trials (optional)
      n_total = if (add_small_trial_bump) {
        ifelse(n_per_arm < 25, n_total + 2, n_total)
      } else n_total
    )
}

# ────────────────────────────────────────────────
# Define enrichment cohorts ----------------------
# ────────────────────────────────────────────────
# Two enrichment strategies (plus All-comers):
#   1) All-comers (no enrichment)
#   2) Top 50% (median split): Q3 + Q4
#   3) Top quartile: Q4 only

# Quick transparency table (sample size, mean, SD) for each cohort
cohort_params <- lapply(names(cohorts), function(nm){
  d <- cohorts[[nm]]
  tibble(
    Cohort      = nm,
    N_available = nrow(d),
    mean_slope  = mean(d$slope_total, na.rm = TRUE),
    sd_slope    = sd(d$slope_total,   na.rm = TRUE)
  )
}) %>% bind_rows()

print(cohort_params)



# Run calculations -------------------------------

X_grid <- seq(28, 60, by = 2)

res_all_ssim <- lapply(names(cohorts), function(nm){
  calc_n_vs_X(cohorts[[nm]], X_seq = X_grid) %>%
    mutate(Cohort = nm)
}) %>% bind_rows()

# Extract 40% rows for labeling on the line plot
lab40 <- res_all_ssim %>%
  filter(X == 40) %>%
  mutate(
    label = paste0(Cohort, ": n=", n_total),
    # Manual vjust tweaks to reduce overlap (adjust as needed)
    vjust = dplyr::case_when(
      Cohort == "All-comers"                 ~ 0,
      Cohort == "NfL ≥ median"   ~  0.75,
      Cohort == "NfL Q4"         ~ 1.55
    )
  )

lab40

# ────────────────────────────────────────────────
# Grouped bar plot at X = {30, 40, 50} -----------
# ────────────────────────────────────────────────
X_bars <- c(30, 40, 50)

res_bar <- res_all_ssim %>%
  filter(X %in% X_bars) %>%
  mutate(X_fac = factor(X, levels = X_bars, labels = paste0(X_bars, "%")))

# Compute percent reduction vs All-comers for persuasive labeling
ref_tbl <- res_bar %>%
  filter(Cohort == "All-comers") %>%
  select(X, n_total_ref = n_total)

res_bar <- res_bar %>%
  left_join(ref_tbl, by = "X") %>%
  mutate(
    pct_reduction_vs_all = 100 * (n_total_ref - n_total) / n_total_ref,
    pct_reduction_vs_all = ifelse(pct_reduction_vs_all < 0, 0, pct_reduction_vs_all),
    label_pct = ifelse(Cohort == "All-comers", "",
                       sprintf("−%.0f%%", pct_reduction_vs_all))
  )

# ────────────────────────────────────────────────────────────
# Line plot (X% vs required total sample size) ---------------
# ────────────────────────────────────────────────────────────

#lab40

p_line <- ggplot(res_all_ssim, aes(x = X, y = n_total, color = Cohort)) +
  geom_vline(xintercept = 40, linetype = "dotted", alpha=0.35) +
  geom_point(data = lab40, size = 3.5, show.legend = FALSE) +             # keep only labels at 30%
  geom_smooth(se = FALSE, method = "loess", span = 0.5, size = 1.2) +  # smooth curve
  geom_label(
    data = lab40,
    aes(label = label),
    show.legend = FALSE,
    label.size = 0.2,
    label.padding = unit(0.25, "lines"),
    vjust = lab40$vjust,
    hjust = 0,
    nudge_x = 0.5,
    nudge_y = 10,
    color = "gray25",
    fill="white"
  ) +
  coord_cartesian(xlim = c(28, 55)) +
  labs(
    title = "Sample size requirements under NfL-based enrichment in ALS trials",
    #subtitle = "Two-sample t-test, alpha=0.05 (two-sided), power=0.80; 1:1 allocation",
    x = NULL,
    y = "Required total sample size",
    color = "Enrollment strategy"
  ) +
  theme_classic(base_size = 18) +
  theme(legend.position = "none",
        plot.title = element_text(size = 20, face = "bold")) +
  scale_color_manual(values = c("gray70","#5FAE9C","#2F7F6F"))+
  scale_x_continuous(labels = function(x) paste0(x, "%"))

p_line

# ────────────────────────────────────────────────
# Grouped bar plot at X = {20, 30, 40} -----------
# ────────────────────────────────────────────────
p_bar <- ggplot(res_bar, aes(x = X_fac, y = n_total, fill = Cohort)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(
    aes(label = label_pct),
    position = position_dodge(width = 0.7),
    #size = 4, vjust = -1,
    size = 4, vjust = -1, hjust=0.125, angle=15
  ) +
  labs(
    #title = "Required sample size at common effect levels",
    #subtitle = "Numbers above bars show % reduction vs All-comers",
    x = NULL,
    y = NULL, #y = "Required total sample size",
    fill = "Enrichment strategy"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c("right", "top"),
    axis.title = element_text(size = 16),
    legend.title = element_text(size=15),
  ) +
  ylim(0, 400) +
  #scale_fill_manual(values = c("#00876C","#4C72B0","#D55E00"))
  scale_fill_manual(values = c("gray70","#5FAE9C","#2F7F6F"))
  #scale_fill_manual(values =c("#009E73","#56B4E9","#E69F00"))

p_bar <- p_bar +
  theme(
    legend.background = element_rect(
      fill = "white",
      colour = NA 
    ),
    legend.box.background = element_rect(
      fill = "white",
      colour = "gray50",
      linewidth = 0.6
    )
  )


# 📊 Fig 6A -------------------------------------------------

p_line|p_bar


Fig6A <- (p_line | p_bar) +
  plot_annotation(
    caption = "Expected relative treatment effect on ALSFRS-R slope (%)"
  ) &
  theme(
    plot.caption = element_text(hjust = 0.5, size = 18, face = "plain",
                                margin = ggplot2::margin(t = 10, b = 0)),
    plot.margin  =ggplot2:: margin(5, 5, 5, 5)
  )

Fig6A

