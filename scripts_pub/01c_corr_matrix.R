library(reshape2)
library(gtsummary)
library(purrr)


# Function  ------------------------------------------------------

plot_corr_matrix_combined <- function(data_clinical, 
                                      data_biomarker = NULL,
                                      method = c("pearson", "spearman"), 
                                      title = NULL, 
                                      text_size = 4,
                                      sig_only = FALSE,
                                      label_map = NULL) {
  method <- match.arg(method)
  
  # Compute p-values for correlations
  cor_mtest <- function(mat1, mat2 = NULL, ...) {
    mat1 <- as.matrix(mat1)
    if (is.null(mat2)) {
      # Square matrix case
      n <- ncol(mat1)
      p.mat <- matrix(NA, n, n)
      diag(p.mat) <- 0
      for (i in 1:n) {
        for (j in 1:n) {
          if (i != j) {
            tmp <- cor.test(mat1[, i], mat1[, j], method = method, ...)
            p.mat[i, j] <- tmp$p.value
          }
        }
      }
      rownames(p.mat) <- colnames(p.mat) <- colnames(mat1)
    } else {
      # Rectangular matrix case
      mat2 <- as.matrix(mat2)
      n1 <- ncol(mat1)
      n2 <- ncol(mat2)
      p.mat <- matrix(NA, n1, n2)
      for (i in 1:n1) {
        for (j in 1:n2) {
          tmp <- cor.test(mat1[, i], mat2[, j], method = method, ...)
          p.mat[i, j] <- tmp$p.value
        }
      }
      rownames(p.mat) <- colnames(mat1)
      colnames(p.mat) <- colnames(mat2)
    }
    return(p.mat)
  }
  
  fill_label <- ifelse(method == "pearson", "Pearson Correlation", "Spearman Correlation")
  if (is.null(title)) title <- paste("Correlation Matrix with p-values -", method)
  
  # Check if biomarker data is provided
  if (is.null(data_biomarker) || ncol(data_biomarker) == 0) {
    # Clinical-only mode: square matrix with all cells
    corr_combined <- cor(data_clinical, use = "pairwise.complete.obs", method = method)
    p_combined <- cor_mtest(data_clinical)
    clinical_vars <- colnames(data_clinical)
  } else {
    # Combined mode: clinical vs clinical + clinical vs biomarker
    corr_clinical <- cor(data_clinical, use = "pairwise.complete.obs", method = method)
    p_clinical <- cor_mtest(data_clinical)
    
    corr_biomarker <- cor(data_clinical, data_biomarker, 
                          use = "pairwise.complete.obs", method = method)
    p_biomarker <- cor_mtest(data_clinical, data_biomarker)
    
    corr_combined <- cbind(corr_clinical, corr_biomarker)
    p_combined <- cbind(p_clinical, p_biomarker)
    clinical_vars <- colnames(data_clinical)
  }
  
  # Reshape to long format
  melted_corr <- melt(corr_combined, varnames = c("Var1", "Var2"), value.name = "corr")
  melted_p <- melt(p_combined, varnames = c("Var1", "Var2"), value.name = "pvalue")
  melted_all <- left_join(melted_corr, melted_p, by = c("Var1", "Var2")) %>%
    mutate(Var1 = factor(Var1, levels = rownames(corr_combined)),
           Var2 = factor(Var2, levels = colnames(corr_combined))) %>%
    mutate(
      # Identify diagonal cells
      is_diagonal = (as.character(Var1) == as.character(Var2)) & 
        (as.character(Var2) %in% clinical_vars),
      # Create label
      label = ifelse(is_diagonal, "", 
                     paste0(round(corr, 2), "\n(",
                            ifelse(pvalue < 0.001, "p<0.001", paste0("p=", round(pvalue, 3))), ")"))
    )
  
  # Significance-based alpha adjustment
  if (sig_only) {
    melted_all <- melted_all %>%
      mutate(alpha_val = ifelse(is_diagonal, 0, ifelse(pvalue < 0.05, 1, 0.5)))
  } else {
    melted_all <- melted_all %>%
      mutate(alpha_val = ifelse(is_diagonal, 0, 1))
  }
  
  # Apply human-readable labels
  if (!is.null(label_map)) {
    label_fun <- function(x) ifelse(x %in% names(label_map), label_map[[x]], x)
    x_labels <- sapply(levels(melted_all$Var2), label_fun)
    y_labels <- sapply(levels(melted_all$Var1), label_fun)
  } else {
    x_labels <- levels(melted_all$Var2)
    y_labels <- levels(melted_all$Var1)
  }
  
  # Plot with diagonal masking
  p <- ggplot(melted_all, aes(x = Var2, y = Var1)) +
    # Base tiles
    geom_tile(aes(fill = corr, alpha = alpha_val), color = "gray80") +
    # Add diagonal lines only
    geom_segment(data = melted_all %>% filter(is_diagonal),
                 aes(x = as.numeric(Var2) - 0.5, y = as.numeric(Var1) - 0.5,
                     xend = as.numeric(Var2) + 0.5, yend = as.numeric(Var1) + 0.5),
                 color = "grey80", linewidth = 0.3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limit = c(-1, 1), space = "Lab",
                         name = fill_label) +
    geom_text(aes(label = label), color = "black", size = text_size) +
    scale_alpha_identity() +
    theme_minimal() +
    labs(title = title, x = "", y = "") +
    theme(
      legend.position = "none",
      text = element_text(size = 16),
      axis.text = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_text(size = 16),
      strip.text = element_text(size = 16),
      plot.title = element_text(size = 20, hjust = 0, face = "bold", margin = ggplot2::margin(b = 10)),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_x_discrete(labels = x_labels) +
    scale_y_discrete(labels = y_labels)
  
  return(p)
}

# Pearson ----------------------------------------------------------------
# Samples excluded from analyses requiring temporal alignment:
# - `second`: repeated CSF samples from the same participant
# - `NotValidSample`: samples with substantial timing mismatch (>90 days) between CSF collection and baseline clinical/demographic assessment

# Pearson: Biomarker-only analyses (square, all cells, log-transformed)
df_biomarker_als_pearson <- csfdf %>% 
  dplyr::filter(!SampleID %in% second & !SampleID %in% NotValidSample) %>%
  dplyr::filter(ALS_label == "1") %>%
  dplyr::select(Age_init,
                log_NfL_csf, log_GFAP_csf, 
                log_Ab38_csf, log_Ab40_csf, log_Ab42_csf,
                log_Ab4240_csf, log_Ab3840_csf,
                log_pTau217_csf, log_pTau181_csf,
                log_Cr, log_Cr_CysC_ratio, log_CK, log_Alb)

SuppleFig_Pearson_Corr_ALS <- plot_corr_matrix_combined(
  data_clinical = df_biomarker_als_pearson,
  data_biomarker = NULL,
  method = "pearson",
  sig_only = TRUE,
  label_map = T1label,
  title = "(A) Pearson correlation matrix among fluid biomarkers in ALS"
)

df_biomarker_dc_pearson <- csfdf %>% 
  dplyr::filter(!SampleID %in% second & !SampleID %in% NotValidSample) %>%
  dplyr::filter(ALS_label == "0") %>%
  dplyr::select(Age_init,
                log_NfL_csf, log_GFAP_csf, 
                log_Ab38_csf, log_Ab40_csf, log_Ab42_csf,
                log_Ab4240_csf, log_Ab3840_csf,
                log_pTau217_csf, log_pTau181_csf,
                log_Cr, log_Cr_CysC_ratio, log_CK, log_Alb)

SuppleFig_Pearson_Corr_DC <- plot_corr_matrix_combined(
  data_clinical = df_biomarker_dc_pearson,
  data_biomarker = NULL,
  method = "pearson",
  sig_only = TRUE,
  label_map = T1label,
  title = "(B) Pearson correlation matrix among fluid biomarkers in disease controls"
)


SuppleFig2AB <- SuppleFig_Pearson_Corr_ALS/SuppleFig_Pearson_Corr_DC
SuppleFig2AB

# Spearman ----------------------------------------------------------------

# df_slope_model is expected to contain one row per SampleID.
stopifnot(!anyDuplicated(df_slope_model$SampleID))

# Spearman: Combined Clinical vs Clinical + Clinical vs Biomarker (single rectangular plot)
df_clinical_als <- csfdf %>% 
  dplyr::filter(!SampleID %in% second & !SampleID %in% NotValidSample) %>%
  left_join(df_slope_model %>% dplyr::select(SampleID, slope_total), by = "SampleID") %>%
  dplyr::filter(ALS_label == "1") %>%
  dplyr::select(Age_init, slope_total, deltaFS, VC_Percent, MMSE, FAB)

df_biomarker_als <- csfdf %>% 
  dplyr::filter(!SampleID %in% second & !SampleID %in% NotValidSample) %>%
  left_join(df_slope_model %>% dplyr::select(SampleID, slope_total), by = "SampleID") %>%
  dplyr::filter(ALS_label == "1") %>%
  dplyr::select(NfL_csf_pgml, GFAP_csf_pgml, 
                Ab38_csf, Ab40_csf, Ab42_csf,
                Ab42_40_csf_bridged, Ab38_40_csf,
                pTau217_csf, pTau181_csf,
                Cr, Cr_CysC_ratio, CK, Alb)


SuppleFig_Spearman_Combined_ALS <- plot_corr_matrix_combined(
  data_clinical = df_clinical_als,
  data_biomarker = df_biomarker_als,
  method = "spearman",
  sig_only = TRUE,
  label_map = T1label_bio_overwritten,
  title = "(C) Spearman correlation matrix among fluid biomarker and clinical characteristics in ALS"
)


# 🔍 Supple Fig 2 -------------------------------------------


SuppleFig_Spearman_Combined_ALS


SuppleFig2 <- SuppleFig_Pearson_Corr_ALS/
  SuppleFig_Pearson_Corr_DC/
  SuppleFig_Spearman_Combined_ALS+
  plot_layout(ncol = 1, heights = c(1,1,0.5))

SuppleFig2

# MMSE/FAB scatterplot vs Aβ42/40 ratio --------------------------------
# Scatter plots of cognitive scores vs CSF AD markers in ALS patients,
# faceted by cognitive measure (MMSE and FAB)

cognitive_levels <- c("MMSE", "FAB")

als_plot_data <- csfdf %>%
  dplyr::filter(!SampleID %in% second & !SampleID %in% NotValidSample) %>%
  dplyr::filter(ALS_label == "1") %>%
  dplyr::filter(!is.na(Ab42_40_csf_bridged), !is.na(Apos_label)) %>%
  dplyr::mutate(
    Group = "ALS"
  ) %>%
  tidyr::pivot_longer(
    cols = c(MMSE, FAB),
    names_to = "Cognitive_measure",
    values_to = "Cognitive_score"
  ) %>%
  dplyr::filter(!is.na(Cognitive_score)) %>%
  dplyr::mutate(
    Cognitive_measure = factor(Cognitive_measure, levels = cognitive_levels)
  )

# Use the same source data as the correlation matrix
cor_labels <- df_clinical_als |>
  dplyr::bind_cols(df_biomarker_als |> dplyr::select(Ab42_40_csf_bridged)) |>
  tidyr::pivot_longer(
    cols = c(MMSE, FAB),
    names_to = "Cognitive_measure",
    values_to = "Cognitive_score"
  ) |>
  dplyr::filter(!is.na(Ab42_40_csf_bridged), !is.na(Cognitive_score)) |>
  # Ensure factor levels match the main data
  dplyr::mutate(Cognitive_measure = factor(Cognitive_measure, levels = cognitive_levels)) |>
  dplyr::group_by(Cognitive_measure) |>
  dplyr::summarise(
    r = cor.test(Ab42_40_csf_bridged, Cognitive_score, method = "spearman")$estimate,
    p = cor.test(Ab42_40_csf_bridged, Cognitive_score, method = "spearman")$p.value,
    .groups = "drop"
  ) |>
  dplyr::mutate(
    label = paste0("\u03c1 = ", round(r, 2), ", p = ", signif(p, 2))
  )

cog_Ab4240_scatter <- ggplot(
  als_plot_data,
  aes(x = Ab42_40_csf_bridged, y = Cognitive_score)
) +
  geom_point(aes(color = Apos_label), size = 1.4, alpha = 0.5) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black",
    linewidth = 1,
    alpha = 0.2
  ) +
  geom_text(
    data = cor_labels,
    aes(label = label),
    x = Inf, y = -Inf,
    hjust = 1.05, vjust = -1.5,
    size = 5,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ Cognitive_measure, scales = "free_y") +
  labs(
    x = "CSF A\u03b242/40",
    y = "Cognitive score",
    color = "A\u03b2 status"
  ) +
  scale_color_manual(values = c(
    "ALS (A\u03b2-)" = "gray60",
    "ALS (A\u03b2+)" = "#8C1C13"
  )) +
  ggh4x::facetted_pos_scales(
    y = list(
      Cognitive_measure == "MMSE" ~ scale_y_continuous(breaks = seq(0, 30, 10), limits = c(0, 30)),
      Cognitive_measure == "FAB"  ~ scale_y_continuous(breaks = seq(0, 18, 6),  limits = c(0, 18))
    )
  ) +
  theme_classic(base_size = 20) + ggcommon_theme

cog_Ab4240_scatter


# 🔍 Supple Fig 3 --------------------------------------------------------------

SuppleFig3_added <- cog_Ab4240_scatter
SuppleFig3 <- SuppleFig3_added 

# Linear regression -------------------------------------------------------

# Create a modeling dataset with the scaled Aβ42/40 variable
# Multiplying by 100 means a 1-unit increase in the scaled variable 
# corresponds to a 0.01-unit increase in the original CSF Aβ42/40 ratio.
T1_model <- T1 |> 
  dplyr::filter(ALS_label == 1) |> 
  dplyr::mutate(Ab42_40_csf_bridged_100scaled = Ab42_40_csf_bridged * 100)

# MMSE ~ Aβ42/40 + age + sex + education year in ALS
lm_mmse_als <- lm(MMSE ~ Ab42_40_csf_bridged_100scaled + Age_init + Sex + Education_Years,
                  data = T1_model)

# FAB ~ Aβ42/40 + age + sex + education year in ALS
lm_fab_als  <- lm(FAB ~ Ab42_40_csf_bridged_100scaled + Age_init + Sex + Education_Years,
                  data = T1_model)

n_mmse <- nobs(lm_mmse_als)
n_fab  <- nobs(lm_fab_als)

# MMSE and FAB regression tables for ALS only, merged side by side
tbl_mmse_als <- tbl_regression(
  lm_mmse_als,
  label = list(
    Ab42_40_csf_bridged_100scaled ~ "CSF A\u03b242/40 (per 0.01 increase)",
    Age_init                   ~ "Age at baseline",
    Sex                        ~ "Sex (male vs female)",
    Education_Years            ~ "Education (years)"
  ),
  show_single_row = "Sex"  # Compress to single row showing Male vs Female
) 

tbl_fab_als <- tbl_regression(
  lm_fab_als,
  label = list(
    Ab42_40_csf_bridged_100scaled ~ "CSF A\u03b242/40 (per 0.01 increase)",
    Age_init                   ~ "Age at baseline",
    Sex                        ~ "Sex (male vs female)",
    Education_Years            ~ "Education (years)"
  ),
  show_single_row = "Sex"
)

# Merge the two tables side by side with spanner headers
tbl_cog_regression <- tbl_merge(
  tbls = list(tbl_mmse_als, tbl_fab_als),
  tab_spanner = c("**MMSE**", "**FAB**")
) |>
  modify_caption("Multivariable linear regression of cognitive scores in ALS")


# 🔍 Supple T7 ----------------------------------------------------------------

SuppleT7_rev_pub  <- tbl_cog_regression |>
  as_flex_table() |>
  flextable::add_footer_lines(
    paste0(
      "A\u03b2, amyloid-\u03b2; \u03b2, regression coefficient per unit increase in the predictor. ",
      "CI, confidence interval; CSF, cerebrospinal fluid; ",
      "MMSE, Mini-Mental State Examination (n = ", n_mmse, "); ",
      "FAB, Frontal Assessment Battery (n = ", n_fab, "). ",
      "Analyses were based on complete cases for all variables included in each model."
    )
  ) 

SuppleT7 <- SuppleT7_rev_pub 

SuppleT7