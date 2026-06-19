suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(tibble); library(scales)
})

# Function ---------------------------------------------------------------

# representative repeat: pick the one closest to median performance
pick_representative_repeat <- function(rep_out, variant = "youden", anchor_method = "rf") {
  runs <- rep_out$runs
  stopifnot(length(runs) > 0, anchor_method %in% names(runs[[1]]$per_method))
  rows <- lapply(seq_along(runs), function(i) {
    pm <- runs[[i]]$per_method[[anchor_method]]
    vo <- pm$pooled$variants_overall
    go <- pm$pooled$gray_overall
    bc <- pm$pooled$bucket_counts
    vo <- vo[vo$variant == variant, , drop = FALSE]
    go <- go[go$variant == variant, , drop = FALSE]
    if (!"PPV_gray" %in% names(go) && "PPV" %in% names(go)) {
      names(go)[names(go)=="PPV"] <- "PPV_gray"; names(go)[names(go)=="NPV"] <- "NPV_gray"
    }
    cov_g <- bc$prop[bc$bucket=="gray"]
    data.frame(
      repeat_id = i,
      AUC = vo$AUC[1], AUC_gray = vo$AUC_gray[1],
      PPV_gray = go$PPV_gray[1], NPV_gray = go$NPV_gray[1],
      gray_cov = cov_g
    )
  })
  df <- dplyr::bind_rows(rows)
  med <- vapply(df[, c("AUC","AUC_gray","PPV_gray","NPV_gray","gray_cov")],
                median, numeric(1), na.rm = TRUE)
  dist <- apply(df[, names(med)], 1, function(x) sqrt(sum((x - med)^2, na.rm = TRUE)))
  df$repeat_id[which.min(dist)]
}

# rep_run → km like structure
as_km_like <- function(one_run) {
  per_method <- one_run$per_method
  overall_rows <- lapply(names(per_method), function(m) {
    ro <- per_method[[m]]$pooled$variants_overall; ro$method <- m
    ro[, c("method","variant","AUC","AUC_gray","Sens","Spec","PPV","NPV")]
  })
  overall_cmp <- dplyr::bind_rows(overall_rows)
  
  gray_rows <- lapply(names(per_method), function(m) {
    rg <- per_method[[m]]$pooled$gray_overall
    if (!"PPV_gray" %in% names(rg) && "PPV" %in% names(rg)) {
      names(rg)[names(rg)=="PPV"] <- "PPV_gray"; names(rg)[names(rg)=="NPV"] <- "NPV_gray"
    }
    rg$method <- m
    rg[, c("method","variant","gray_n","TP","TN","FP","FN",
           "Sens_gray","Spec_gray","PPV_gray","NPV_gray")]
  })
  gray_cmp <- dplyr::bind_rows(gray_rows)
  
  any_m <- names(per_method)[1]
  list(
    per_method = per_method,
    compare = list(
      overall   = overall_cmp,
      gray_only = gray_cmp,
      buckets   = per_method[[any_m]]$pooled$bucket_counts,
      nfl       = list(baselines = per_method[[any_m]]$pooled$nfl_baselines)
    ),
    meta = list(
      outcome = one_run$meta$outcome,
      methods_gray = names(per_method),
      data = per_method[[any_m]]$data
    )
  )
}

# use km to get NfL single-cutoff Sens/Spec (OOF reconstructed)
get_nfl_sens_spec_from_km <- function(km, variant = c("sp95","se95","youden"),
                                      outcome = "ALS_label",
                                      nfl_var = "NfL_csf_pgml") {
  variant <- match.arg(variant)
  any_m <- km$meta$methods_gray[1]
  df_clean <- km$per_method[[any_m]]$data$df_clean
  folds    <- km$per_method[[any_m]]$data$folds
  base_tbl <- km$compare$nfl$baselines
  row_idx <- switch(variant,
                    youden = which(base_tbl$baseline == "nfl_youden"),
                    sp95   = which(base_tbl$baseline == "nfl_sp95"),
                    se95   = which(base_tbl$baseline == "nfl_se95"))
  th_list <- base_tbl$th_train_all[[row_idx]]
  
  pos <- "ALS"; neg <- "nonALS"
  pred_oof <- rep(NA_character_, nrow(df_clean))
  for (i in seq_along(folds)) {
    idx <- folds[[i]]; th <- th_list[i]
    pred_oof[idx] <- ifelse(df_clean[[nfl_var]][idx] >= th, pos, neg)
  }
  pred_oof <- factor(pred_oof, levels = c(neg,pos))
  truth    <- factor(df_clean[[outcome]], levels = c(neg,pos))
  
  TP <- sum(truth==pos & pred_oof==pos, na.rm=TRUE)
  TN <- sum(truth==neg & pred_oof==neg, na.rm=TRUE)
  FP <- sum(truth==neg & pred_oof==pos, na.rm=TRUE)
  FN <- sum(truth==pos & pred_oof==neg, na.rm=TRUE)
  c(Sens = TP/(TP+FN), Spec = TN/(TN+FP))
}

# Simulate PPV/NPV curve across prevalence
.sim_curve <- function(sens, spec, prev) {
  tibble(
    prevalence = prev,
    PPV = (sens*prevalence) / (sens*prevalence + (1-spec)*(1-prevalence)),
    NPV = (spec*(1-prevalence)) / (spec*(1-prevalence) + (1-sens)*prevalence)
  )
}


# Main function to plot priority curve
plot_priority_curve <- function(rep_out,
                                priority = c("rule_in_sp95","rule_out_se95"),
                                method = "rf",
                                anchor_for_rep = "rf",
                                prev_range = seq(0.01, 0.90, by = 0.01),
                                outcome = "ALS_label",
                                nfl_var = "NfL_csf_pgml") {
  priority <- match.arg(priority)
  
  # representative repeat
  rep_id <- pick_representative_repeat(rep_out, variant = "youden", anchor_method = anchor_for_rep)
  km <- as_km_like(rep_out$runs[[rep_id]])
  
  # Two-step variant table
  vtab <- km$per_method[[method]]$pooled$variants_overall
  
  if (priority == "rule_in_sp95") {
    # Priority == "rule_in_sp95" → PPV 
    ss <- vtab %>% filter(variant=="sp95") %>% dplyr::slice(1) %>% transmute(Sens=Sens, Spec=Spec) %>% as.list()
    nfl <- get_nfl_sens_spec_from_km(km, "youden", outcome, nfl_var) #youden
    
    two <- .sim_curve(ss$Sens, ss$Spec, prev_range) %>%
      select(prevalence, PPV) %>% mutate(kind="Two-step")
    one <- .sim_curve(nfl["Sens"], nfl["Spec"], prev_range) %>%
      select(prevalence, PPV) %>% mutate(kind="NfL single cut-off")
    
    plot_df <- bind_rows(two, one) %>% rename(value = PPV)
    ylab <- "Predictive Values (PPV/NPV)"
    subtitle <- sprintf("Rule-in–prioritized two-step strategy vs single NfL cut-off. Method: %s", method)
    
  } else {
    # priority == "rule_out_se95" → NPV
    ss <- vtab %>% filter(variant=="se95") %>% dplyr::slice(1) %>% transmute(Sens=Sens, Spec=Spec) %>% as.list()
    nfl <- get_nfl_sens_spec_from_km(km, "youden", outcome, nfl_var) #youden
    
    two <- .sim_curve(ss$Sens, ss$Spec, prev_range) %>%
      select(prevalence, NPV) %>% mutate(kind="Two-step")
    one <- .sim_curve(nfl["Sens"], nfl["Spec"], prev_range) %>%
      select(prevalence, NPV) %>% mutate(kind="NfL single cut-off")
    
    plot_df <- bind_rows(two, one) %>% rename(value = NPV)
    ylab <- "Predictive values (PPV/NPV)"
    subtitle <- sprintf("Rule-out–prioritized two-step strategy vs single NfL cut-off.  Method: %s", method)
  }
  
  # add vertical line for study prevalence derived from df_clean (representative repeat)
  df_clean <- km$per_method[[km$meta$methods_gray[1]]]$data$df_clean
  prev0 <- mean(df_clean[[outcome]] == "ALS", na.rm = TRUE)
  
  p <- ggplot(plot_df, aes(x = prevalence, y = value, colour = kind)) +
    geom_line(linewidth = 1.2) +
    geom_vline(xintercept = prev0, linetype = "dotted", colour = "gray40") +
    annotate("text", x = prev0, y = 0.5,
             label = sprintf("study prevalence (%.1f%%)", 100*prev0),
             angle = 90, vjust = -0.5, size = 4, colour = "gray30") +
    scale_x_continuous(labels = percent_format()) +
    scale_y_continuous(labels = percent_format(), limits = c(0,1)) +
    scale_colour_manual(values = c("NfL single cut-off"="#0072B2","Two-step"="#D55E00"), name = "Strategy") +
    labs(
      title = if (priority=="rule_in_sp95")
        "Positive Predictive Values across Disease Prevalence"
      else
        "Negative Predictive Values across Disease Prevalence",
      subtitle = subtitle,
      x = "Disease prevalence",
      y = ylab
    ) +
    theme_classic(base_size = 16) +
    #height-width ratio 
    coord_fixed(ratio = 0.75) +
    theme(legend.position = "bottom", 
          panel.grid.minor = element_blank(),
          axis.title.x = element_text(margin = margin(t = 10)),
          axis.title.y = element_text(margin = margin(r = 10))
    )
  
  list(rep_id = rep_id,
       priority = priority,
       sens_spec = list(
         twostep = unlist(ss),
         nfl = nfl
       ),
       data = plot_df,
       plot = p)
}

# add NfL se95 NPV curve overlay to existing plot
add_se95_npv_overlay <- function(res_rulein,
                                 rep_out,
                                 anchor_for_rep = "rf",
                                 outcome = "ALS_label",
                                 nfl_var = "NfL_csf_pgml",
                                 prev_range = seq(0.01, 0.90, by = 0.01)) {
  stopifnot(is.list(res_rulein), !is.null(res_rulein$rep_id))
  
  # representative repeat
  km <- as_km_like(rep_out$runs[[res_rulein$rep_id]])
  
  # get NfL se95 Sens/Spec
  nfl_se <- get_nfl_sens_spec_from_km(km, variant = "se95",
                                      outcome = outcome, nfl_var = nfl_var)
  
  # generate NPV curve data
  ref_df <- .sim_curve(nfl_se["Sens"], nfl_se["Spec"], prev_range) |>
    dplyr::transmute(prevalence, value = NPV)
  
  # add to existing plot
  res_rulein$plot +
    ggplot2::geom_line(
      data = ref_df,
      ggplot2::aes(x = prevalence, y = value),
      inherit.aes = FALSE,
      linewidth = 0.7, linetype = "dashed", color = "gray40"
    ) +
    ggplot2::annotate(
      "text", x = 0.08, y = 0.9,
      label = "NPV reference (NfL Se95)",
      hjust = 0, size = 4, color = "gray30"
    )
}

# Run ---------------------------------------------------------------------
# The single-cutoff reference uses the Youden threshold.

res_rulein <- plot_priority_curve(rep_out, priority = "rule_in_sp95", method = "rf", anchor_for_rep = "rf")
res_rulein

# add NfL se95 NPV curve overlay to rule-in plot
p_rulein_with_npv_ref <- add_se95_npv_overlay(res_rulein, rep_out)


# 🔍 Supple Fig S8 -----------------------------------------------------------

SuppleFig8 <- p_rulein_with_npv_ref+
  labs(
    title = "Predictive values across disease prevalence",
    subtitle = "Rule-in–prioritized two-step strategy vs single cut-off approach")

SuppleFig8

