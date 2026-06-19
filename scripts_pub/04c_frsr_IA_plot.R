# ============================================================================
# NfL Modifier Interaction Visualization: Example Usage
# ============================================================================

# Load functions
source("scripts_pub/function/lmer_IAvisualize_func.R")

# Configuration ------------------------------------------------------------

# Base covariates (for main effects)
base_covariates_viz <- c("Age_init", "VC_Percent", "deltaFS", "REEC_definite", "OnsetSite")

# *Time Interaction covariates
interaction_covariates_viz = c("Age_init","deltaFS","VC_Percent","REEC_definite","OnsetSite")

# 📊 Figure 4C --------------------------------------------------------

fm_Ab_status_core <- fit_nfl_modifier_interaction(
  df_vis %>%
    mutate(Ab_status_bl = factor(Ab_status_bl, 
                                 levels = c("negative", "positive"))),
  modifier_var = "Ab_status_bl",
  modifier_type = "factor",
  categorize_nfl = FALSE,
  base_covariates = base_covariates_viz,
  interaction_covariates = interaction_covariates_viz,
  cut_re = FALSE
)

Fig_Ab_status_core <- plot_nfl_modifier_interaction(
  fm_Ab_status_core,
  modifier_var  = "Ab_status_bl",
  modifier_type = "factor",
  level_mode    = "quantile",
  nfL_probs     = c(.25,.5,.75),
  months_range  = c(0,36),
  months_by     = 3,
  base_font_size= 18,
  T1label       = T1label
)


Fig4C <- Fig_Ab_status_core +
  facet_wrap(
    ~ facet,
    labeller = as_labeller(c(
      "negative" = "Aβ-negative ALS",
      "positive" = "Aβ-positive ALS"
    ))
  )

Fig4C
