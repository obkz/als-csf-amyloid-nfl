# Build pre-baseline (+α) extension to onset (y=48) ----------
# Assumes visitdf has:
#   SSBS_ID, DaysFromFirstVisit_months, ALSFRSR_Total,
#   DiseaseDuration_M (months at baseline),
#   deltaFS (points/month; larger = faster decline)

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(purrr)

theme_pub <- function() {
  theme(
    plot.title   = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16, margin = ggplot2::margin(t = 15)),
    axis.title.y = element_text(size = 16, margin = ggplot2::margin(r = 10)),
    axis.text    = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 14),
    legend.position = "bottom"
  )
}

# 1) identify each subject’s baseline row (earliest visit)
baseline_df <- visitdf %>%
  group_by(SSBS_ID) %>%
  arrange(DaysFromFirstVisit_months, .by_group = TRUE) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  select(
    SSBS_ID, SampleID,
    ALSFRSR_baseline = ALSFRSR_Total,
    t0_months       = DaysFromFirstVisit_months,
    DiseaseDuration_M,
    deltaFS,
    tert_cNfL,
    Ab_status, VC_Percent
  )


# 2) label “fast vs slow” by median deltaFS among ALS only
med_delta <- visitdf %>%
  filter(ALS_label == 1) %>%
  group_by(SSBS_ID) %>%
  arrange(DaysFromFirstVisit_months, .by_group = TRUE) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  pull(deltaFS) %>%
  stats::median(na.rm = TRUE)


baseline_df <- baseline_df %>%
  mutate(deltaFS_group = case_when(
    is.na(deltaFS)            ~ NA_character_,
    deltaFS >= med_delta      ~ "≥ median ΔFS at baseline",
    TRUE                      ~ "< median ΔFS at baseline"),
    VC_group = case_when(
      is.na(VC_Percent)            ~ NA_character_,
      VC_Percent >= median(VC_Percent, na.rm = TRUE)   ~ "≥ median %VC at baseline",
      TRUE                         ~ "< median %VC at baseline"
    )
  )

# 3) create a pseudo onset point at x = -DiseaseDuration_M, y = 48
# Assumption for visualization: ALSFRS-R total score is set to 48 at symptom onset.
pseudo_onset <- baseline_df %>%
  filter(!is.na(DiseaseDuration_M), DiseaseDuration_M > 0) %>%
  transmute(
    SSBS_ID,
    SampleID,
    rel_months       = -DiseaseDuration_M,  # time since onset (negative)
    ALSFRSR_Total    = 48,
    deltaFS_group,
    tert_cNfL,
    Ab_status,
    VC_group,
    segment_type     = "pre"                   # tag for styling
  )

# 4) add an explicit baseline anchor at x=0 to ensure a straight segment joins onset→baseline
baseline_anchor <- baseline_df %>%
  transmute(
    SSBS_ID,
    SampleID,
    rel_months       = 0,
    ALSFRSR_Total    = ALSFRSR_baseline,
    deltaFS_group,
    tert_cNfL,
    Ab_status,
    VC_group,
    segment_type     = "baseline"
  )

# 5) real (post-baseline) visits with time re-zeroed at baseline
post_baseline <- visitdf %>%
  left_join(baseline_df %>% select(SSBS_ID, t0_months, deltaFS_group, VC_group),
            by = "SSBS_ID") %>%
  mutate(rel_months = DaysFromFirstVisit_months - t0_months) %>%
  select(SSBS_ID, SampleID, rel_months, ALSFRSR_Total, deltaFS_group, tert_cNfL, Ab_status, VC_group) %>%
  mutate(segment_type = "post")

# 6) assemble extended dataset
spaghetti_ext <- bind_rows(pseudo_onset, baseline_anchor, post_baseline) %>%
  arrange(SSBS_ID, rel_months)



# DeltaFS: dashed pre-baseline segments + solid post-baseline ----------------------------

Fig4A <- ggplot() +
  # pre-baseline dashed line from onset→baseline (per subject)
  geom_line(
    data = spaghetti_ext %>% filter(segment_type %in% c("pre", "baseline")),
    aes(x = rel_months, y = ALSFRSR_Total, group = SSBS_ID,
        color = deltaFS_group, linetype = deltaFS_group),
    alpha = 0.25
  ) +
  # post-baseline solid
  geom_line(
    data = spaghetti_ext %>% filter(segment_type == "post"),
    aes(x = rel_months, y = ALSFRSR_Total, group = SSBS_ID,
        color = deltaFS_group, linetype = deltaFS_group),
    alpha = 0.48
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3, alpha = 0.6) +
  scale_x_continuous(
    limits =c(-48,128), 
    breaks = seq(
      floor(min(spaghetti_ext$rel_months, na.rm = TRUE)/12)*12,
      ceiling(max(spaghetti_ext$rel_months, na.rm = TRUE)/12)*12,
      by = 24
    )
  ) +
  scale_y_continuous(limits = c(0, 48), breaks = seq(0, 48, by = 12)) +
  labs(
    title = "ALSFRS-R trajectories by baseline ΔFS",
    #subtitle = "Dashed = onset → baseline (assumes 48 at onset, linear decline at baseline ΔFS)",
    x = "Months from baseline",
    y = "ALSFRS-R",
    color = "Progression speed",
    linetype = "Progression speed"
  ) +
  theme_classic() +
  scale_color_manual(values = c(
    "< median ΔFS at baseline"  = "#65C2AE",
    "≥ median ΔFS at baseline" = "#E07B00"    
  ), na.translate = FALSE
  ) +
  scale_linetype_manual(values = c(
    "< median ΔFS at baseline"  = "F2",
    "≥ median ΔFS at baseline"  = "solid"
  ), na.value = "blank", na.translate=FALSE)+
  theme_pub()


## 📊 Fig.4A -------------------------------------------------------------------

Fig4A

# Tertiles of cNfL instead of deltaFS_group -------------------

Fig4B <- ggplot() +
  geom_line(
    data = spaghetti_ext %>% filter(segment_type %in% c("pre", "baseline")),
    aes(x = rel_months, y = ALSFRSR_Total, group = SSBS_ID, color = tert_cNfL,
        linetype = tert_cNfL
    ),
    alpha = 0.25
  ) +
  geom_line(
    data = spaghetti_ext %>% filter(segment_type == "post"),
    aes(x = rel_months, y = ALSFRSR_Total, group = SSBS_ID, color = tert_cNfL,
        linetype = tert_cNfL),
    alpha = 0.48
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3, alpha = 0.6) +
  scale_x_continuous(
    limits = c(-48,128),
    breaks = seq(
      floor(min(spaghetti_ext$rel_months, na.rm = TRUE)/12)*12,
      ceiling(max(spaghetti_ext$rel_months, na.rm = TRUE)/12)*12,
      by = 24
    )
  ) +
  scale_y_continuous(limits = c(0, 48), breaks = seq(0, 48, by = 12)) +
  labs(
    title = "ALSFRS-R trajectories by CSF NfL tertile",
    x = "Months from baseline",
    y = "ALSFRS-R",
    color = "CSF NfL tertile",
    linetype = "CSF NfL tertile"
  ) +
  theme_classic() +
  scale_linetype_manual(values = c(
    "low"  = "longdash",
    "mid"  = "F2",
    "high" = "solid"
  ), na.value = "blank")+
  scale_color_manual(values = c(
    "low"  = "#65C2AE", 
    "mid"  = "gray50",
    "high" = "#E07B00" 
  ),na.translate = FALSE)+
  theme_pub()

## 📊 Fig.4B -------------------------------------------------------------------

Fig4B


# %VC ------------------------------------------------------------------------

SuppleFig9 <-ggplot() +
  geom_line(
    data = spaghetti_ext %>% filter(segment_type %in% c("pre", "baseline")),
    aes(x = rel_months, y = ALSFRSR_Total, group = SSBS_ID, color = VC_group,
        linetype = VC_group
    ),
    #linetype = "dashed", 
    alpha = 0.25
  ) +
  geom_line(
    data = spaghetti_ext %>% filter(segment_type == "post"),
    aes(x = rel_months, y = ALSFRSR_Total, group = SSBS_ID, color = VC_group,
        linetype = VC_group),
    alpha = 0.5
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3, alpha = 0.6) +
  scale_x_continuous(
    limits = c(-48,128),
    breaks = seq(
      floor(min(spaghetti_ext$rel_months, na.rm = TRUE)/12)*12,
      ceiling(max(spaghetti_ext$rel_months, na.rm = TRUE)/12)*12,
      by = 24
    )
  ) +
  scale_y_continuous(limits = c(0, 48), breaks = seq(0, 48, by = 12)) +
  labs(
    title = "ALSFRS-R trajectories by baseline %VC",
    x = "Months from baseline",
    y = "ALSFRS-R",
    color = "Respiratory function",
    linetype = "Respiratory function"
  ) +
  theme_classic() +
  scale_color_manual(
    values=c(
      "≥ median %VC at baseline" = "#65C2AE", #"#7FCDBB","#4BB8A3",#"#7FCDBB",#, 
      "< median %VC at baseline" = "#E07B00" #"#D55E00" #"#E69F00" #
    ),
    na.translate = FALSE
  ) +
  scale_linetype_manual(values = c(
    "< median %VC at baseline"  = "solid",
    "≥ median %VC at baseline"  = "F2"
  ), na.value = "blank", na.translate=FALSE)+
  theme_pub()


# 🔍 Supple Fig.9 ------------------------------------------------------------

SuppleFig9
