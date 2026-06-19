library(pROC)
library(ggplot2)
library(dplyr)
library(tibble)
library(RColorBrewer)
library(scales)
library(patchwork)

# Function --------------------------------------------------------------------

# Fixed component colors
component_colors <- c(
  "CK"                  = "#4E79A7",  # blue-cyan
  "pTau217_csf"         = "#59A14F",  # muted green
  "pTau181_csf"         = "#EDC948",  # yellow-gold (lighter)
  "Ab42_40_csf_bridged" = "#F28E8C",  # soft coral
  "Ab38_40_csf"         = "#E15759",  # red-orang
  "cGFAP"               = "#7F7F7F",  # neutral gray
  "Cr"                  = "#B07AA1",  # muted magenta
  "cNfL"                = "#1B9E77"   # KEEP as is
)

plot_multi_roc_gg <- function(
    data,
    predictors,
    outcome     = "ALS_label",
    boot.n      = 2000,
    conf.level  = 0.95,
    title       = NULL
) {
  # helper to count safely
  .count_key <- function(tbl, key) if (key %in% names(tbl)) as.integer(tbl[[key]]) else 0L
  
  # automatic title
  tbl_out <- table(data[[outcome]])
  n_ALS <- .count_key(tbl_out, "ALS"); if (n_ALS == 0L) n_ALS <- .count_key(tbl_out, "1")
  n_DC  <- .count_key(tbl_out, "DC");  if (n_DC  == 0L) n_DC  <- .count_key(tbl_out, "0") + .count_key(tbl_out, "nonALS")
  if (n_DC == 0L) n_DC <- sum(tbl_out) - n_ALS
  default_title <- sprintf("ROC: ALS (n=%d) vs DC (n=%d)", n_ALS, n_DC)
  final_title   <- if (is.null(title)) default_title else title
  
  # storage
  roc_dfs <- list()
  auc_tbl <- tibble(Predictor = character(), AUC = double(), CI_low = double(), CI_high = double())
  
  # fit models and build ROC
  for (pred in predictors) {
    fmla  <- stats::as.formula(paste(outcome, "~", pred))
    model <- stats::glm(fmla, data = data, family = stats::binomial())
    mf    <- stats::model.frame(model)
    probs <- stats::predict(model, type = "response")
    
    # robust levels for pROC
    y <- mf[[outcome]]
    y_chr <- as.character(y)
    levs <- if (all(c("DC","ALS") %in% unique(y_chr))) c("DC", "ALS") else sort(unique(y_chr))
    
    roc_obj <- pROC::roc(y, probs, levels = levs, direction = "auto")
    ci_obj  <- pROC::ci.auc(roc_obj, method = "bootstrap", boot.n = boot.n, conf.level = conf.level)
    auc_val <- as.numeric(pROC::auc(roc_obj))
    
    auc_tbl <- bind_rows(
      auc_tbl,
      tibble(Predictor = pred, AUC = auc_val, CI_low = ci_obj[1], CI_high = ci_obj[3])
    )
    
    df <- tibble(
      Predictor = pred,
      FPR = 1 - roc_obj$specificities,
      TPR = roc_obj$sensitivities
    ) |>
      bind_rows(
        tibble(Predictor = pred, FPR = 0, TPR = 0),
        tibble(Predictor = pred, FPR = 1, TPR = 1)
      ) |>
      group_by(Predictor, FPR) |>
      summarize(TPR = max(TPR), .groups = "drop") |>
      arrange(FPR)
    
    roc_dfs[[pred]] <- df
  }
  
  plot_df <- bind_rows(roc_dfs)
  
  # order by AUC desc
  auc_tbl <- auc_tbl |>
    arrange(desc(AUC)) |>
    mutate(Predictor = factor(Predictor, levels = Predictor))
  predictor_levels <- levels(auc_tbl$Predictor)
  plot_df <- plot_df |> mutate(Predictor = factor(Predictor, levels = predictor_levels))
  
  # pretty labels for legend only
  .pretty_map <- c(
    "cNfL"                = "NfL",
    "cGFAP"               = "GFAP",
    "pTau181_csf"         = "pTau181",
    "pTau217_csf"         = "pTau217",
    "Ab42_40_csf_bridged" = "Aβ42/40",
    "Ab38_40_csf"         = "Aβ38/40",
    "Ab38_42_csf"         = "Aβ38/42",
    "Cr"                  = "Cr",
    "CK"                  = "CK",
    "Ab_status"           = "Aβ status"
  )
  .escape_regex <- function(x) gsub("([\\^$.|?*+()\\[\\]{}])", "\\\\\\1", x, perl = TRUE)
  .pretty_predictor <- function(x, map = .pretty_map) {
    x <- gsub("`", "", x, fixed = TRUE)
    if (!is.null(map) && length(map)) {
      for (k in names(map)) {
        pat <- paste0("\\b", .escape_regex(k), "\\b")
        x   <- gsub(pat, map[[k]], x, perl = TRUE)
      }
    }
    x
  }
  pred_pretty <- .pretty_predictor(as.character(auc_tbl$Predictor))
  legend_labels <- sprintf("%s: AUC=%.2f (%.2f–%.2f)", pred_pretty, auc_tbl$AUC, auc_tbl$CI_low, auc_tbl$CI_high)
  names(legend_labels) <- as.character(auc_tbl$Predictor)
  
  # color channel resolver
  channel_key <- function(pred) {
    s <- gsub("`", "", pred, fixed = TRUE)
    s <- gsub("\\s+", "", s)
    # handle "cNfL+X" or "X+cNfL"
    if (grepl("^cNfL\\+", s)) {
      return(sub("^cNfL\\+", "", s))
    }
    if (grepl("\\+cNfL$", s)) {
      return(sub("\\+cNfL$", "", s))
    }
    # single marker
    return(s)
  }
  
  # stable color mapping per predictor level
  color_values <- setNames(
    vapply(predictor_levels, function(p) {
      key <- channel_key(as.character(p))
      if (key %in% names(component_colors)) component_colors[[key]] else "#555555"
    }, character(1)),
    predictor_levels
  )
  
  # plot
  p <- ggplot(plot_df, aes(x = FPR, y = TPR, color = Predictor)) +
    geom_step(direction = "vh", linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    scale_x_continuous(labels = number_format(accuracy = 0.1), breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
    scale_y_continuous(labels = number_format(accuracy = 0.1), breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
    scale_color_manual(
      values = color_values,
      breaks = predictor_levels,
      labels = legend_labels[predictor_levels],
      name   = "Predictor (AUC, 95% CI)"
    ) +
    labs(title = final_title, x = "1 − Specificity", y = "Sensitivity") +
    coord_equal() +
    theme_classic(base_size = 15) +
    theme(
      legend.position = "right",
      legend.text  = element_text(size = 11),
      legend.title = element_text(size = 12),
      axis.title.x = element_text(margin=ggplot2::margin(t=10)),
      axis.title.y = element_text(margin=ggplot2::margin(r=10))
    )
  
  print(p)
  invisible(list(roc_data = plot_df, auc_table = auc_tbl, plot = p, colors = color_values))
}


# Predictors --------------------------------------------------------------

ROCpredictors <- c(
  "Cr",
  "CK",
  "cNfL",
  "cGFAP",
  "pTau181_csf",
  "pTau217_csf",
  "Ab42_40_csf_bridged",
  "Ab38_40_csf"
)

ROCpredictorsInt <- c(
  "cNfL + cGFAP",
  "cNfL + Ab42_40_csf_bridged",
  "cNfL + Ab38_40_csf",
  "cNfL + pTau181_csf",
  "cNfL + pTau217_csf",
  "cNfL + Cr",
  "cNfL + CK"
)


# 📊 Figure 2AB ---------------------------------------------------------------------

ROC_df %>% filter(category == "ALS") %>% nrow()
ROC_df %>% filter(!category == "ALS") %>% nrow()

# Bootstrap n=2000 for final figure
multiROC <- plot_multi_roc_gg(
  ROC_df %>% rename(cNfL = NfL_csf_pgml, cGFAP=GFAP_csf_pgml), 
  predictors = ROCpredictors,
  title="ALS vs DC: single-biomarker performance")

multiROC2 <- plot_multi_roc_gg(
  ROC_df %>% rename(cNfL = NfL_csf_pgml, cGFAP=GFAP_csf_pgml), 
  ROCpredictorsInt, 
  title="Two-biomarker models including NfL")

multiROC$plot <- multiROC$plot + theme(
  plot.title = element_text(size=20, hjust=0, margin=ggplot2::margin(b=15), face="bold"))

multiROC2$plot <- multiROC2$plot + theme(
  plot.title = element_text(size=20, hjust=0, margin=ggplot2::margin(b=15), face="bold"))

# Fig2 AB
Fig2AB <- multiROC$plot|multiROC2$plot

Fig2AB


# ALS vs mimics dataset ----------------------------

ROC_mimics <- ROC_df %>% mutate(
  category2 = case_when(
    category %in% c("ALS mimics", "Myopathy") ~ "ALS mimics (broader)",
    TRUE ~ category  # keep original for others
  )) %>% 
  filter(category2 %in% c("ALS", "ALS mimics (broader)"))

ROC_mimics %>% filter(category2 == "ALS") %>% nrow()
ROC_mimics %>% filter(category2 == "ALS mimics (broader)") %>% nrow()

# Bootstrap n=2000 
mimics_multiROC <- plot_multi_roc_gg(
  ROC_mimics %>% rename(cNfL = NfL_csf_pgml, cGFAP=GFAP_csf_pgml), 
  predictors = ROCpredictors,
  title="Single-biomarker performance")

mimics_multiROC 

# Bootstrap n=2000 
mimics_multiROC2 <- plot_multi_roc_gg(
  ROC_mimics %>% rename(cNfL = NfL_csf_pgml, cGFAP=GFAP_csf_pgml), 
  ROCpredictorsInt, 
  title="Two-biomarker models including NfL")

mimics_multiROC2

# Merge for publication
ROC_mimics_pub <-(mimics_multiROC$plot | mimics_multiROC2$plot) +
  plot_annotation(title = "(A) ALS vs ALS mimics (n=264 vs 71)",
                  theme = theme(plot.title = element_text(size=20, hjust=0, face="bold",
                                                          margin=ggplot2::margin(l=45, b=5)))
  )

ROC_mimics_pub


# ALS vs other DC excluding mimics -----------------------------------

ROC_otherDC <- ROC_df %>% 
  filter(!category2 %in% c("ALS mimics (broader)"))

ROC_otherDC %>% filter(!category == "ALS") %>% nrow()

otherDC_multiROC <- plot_multi_roc_gg(
  ROC_otherDC %>% rename(cNfL = NfL_csf_pgml, cGFAP=GFAP_csf_pgml), 
  predictors = ROCpredictors,
  title="Single-biomarker performance")

otherDC_multiROC2 <- plot_multi_roc_gg(
  ROC_otherDC %>% rename(cNfL = NfL_csf_pgml, cGFAP=GFAP_csf_pgml), 
  ROCpredictorsInt, 
  title="Two-biomarker models including NfL")


ROC_otherDC_pub <- (otherDC_multiROC$plot | otherDC_multiROC2$plot) +
  plot_annotation(title = "(B) ALS vs other neurological diseases (n=264 vs 93)",
                  theme = theme(plot.title = element_text(size=20, hjust=0, face="bold",
                                                          margin=ggplot2::margin(l=45, b=5))))


# 🔎 Supple FigS4 --------------------------------------------------------------

SuppleFig4 <- wrap_elements(ROC_mimics_pub) / wrap_elements(ROC_otherDC_pub)

SuppleFig4
