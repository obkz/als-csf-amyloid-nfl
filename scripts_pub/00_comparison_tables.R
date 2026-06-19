library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(broom)
library(tibble)
library(flextable)
library(officer)
library(forcats)
library(gtsummary)

# This script was used to generate descriptive and supplementary tables.
# Individual-level clinical data are not shared due to ethical and privacy restrictions.
# Therefore, the script is provided for transparency, not for full reproducibility.
# Running this script requires non-public preprocessed analysis objects, including `T1` and `csfdf`.

flextable::set_flextable_defaults(footnote_reference_symbol = letters)

# baseline characteristics ------------------------------------------------

Apos <- T1 %>% filter(Ab_status=="positive")
Aneg <- T1 %>% filter(Ab_status=="negative")

T1label <- list(Age_init="Age at baseline",
                Age_csf="Age at sample collection",
                ALS_label="ALS status",
                APOE = "APOE haplotype",
                APOE_group = "APOE haplotype†",
                APOE_e4_carrier = "APOE ε4 status (carrier vs non-carrier)",
                APOE_e2_carrier = "APOE ε2 status (carrier vs non-carrier)",
                APOE_e4_vs_e33 = "APOE ε4 carrier vs ε3 homozygote",
                APOE_e2_vs_e33 = "APOE ε2 carrier vs ε3 homozygote",
                slope_total="ALSFRS-R Slope (OLS)",
                genetic_ALS="ALS-related genetic status",
                genetic_label="ALS-related genetic mutation ",
                OnsetSite="Onset site ",
                DiseaseDur_init="Disease duration, months",
                R_ElEscorial ="REEC classification at baseline",
                REEC_definite="REEC ",
                REEC_binary="REEC classification (definite vs others)",
                ALSFRS_init ="ALSFRS-R at baseline",
                ALSFamHistory = "ALS family history",
                Education_Years ="Education years",
                MMSE= "MMSE at baseline",
                FAB="FAB at baseline",
                ACE_R="ACE-R at baseline",
                MOCA="MoCA-J at baseline",
                BMI_init="BMI at baseline",
                deltaFS="ΔFS at baseline",
                sNfL="serum NfL",
                surv_time_day="survival time, days",
                VC_Percent="%VC at baseline",
                CSF_TP="CSF TP (g/L)",
                CSF_Alb="CSF Alb (mg/dL)",
                Ab38_csf="Aβ38 (pg/mL)",
                Ab38="Aβ38",
                Ab40="Aβ40",
                Ab42="Aβ42",
                Ab40_csf="Aβ40 (pg/mL)",
                Ab42_csf="Aβ42 (pg/mL)",
                Ab42_40_csf_bridged="Aβ42/40 ratio",
                Ab38_40_csf="Aβ38/40 ratio",
                Ab3840="Aβ38/40 ratio",
                Ab38_42_csf="Aβ38/42 ratio",
                Ab42_38_csf="Aβ42/38 ratio",
                Ab4240_csf="Aβ42/40 ratio",
                Ab4240="Aβ42/40 ratio",
                Ab_status="Aβ status ",
                pTau181_csf="pTau181 (pg/mL)",
                pTau217_csf="pTau217 (pg/mL)",
                pTau181_ab40="pTau181/Aβ40 ratio",
                pTau217_ab40="pTau217/Aβ40 ratio",
                GFAP_csf_pgml="GFAP (pg/mL)",
                GFAP_csf="GFAP (pg/mL)",
                NfL_csf_pgml="NfL (pg/mL)",
                NfL_csf="NfL (pg/mL)",
                GFAP_csf_NTKbr="GFAP (NTK-br)",
                NfL_csf_NTKbr="NfL (NTK-br)",
                NfL_GFAP_csf="NfL/GFAP ratio",
                GFAP_ab40_csf="GFAP/Aβ40 ratio",
                Cr="Blood Cr (mg/dL)",
                CK="Blood CK (IU/L)",
                Alb="Blood Alb (g/dL)",
                cysC="Blood CysC (mg/L)",
                Cr_CysC_ratio="Cr/CysC ratio",
                CrCys_ratio = "Cr/CysC ratio",
                log_Alb="Blood Alb (ln)",
                log_Cr="Blood Cr (ln)",
                log_CK="Blood CK (ln)",
                log_Cr_CysC_ratio="Cr/CysC ratio (ln)",
                log_NfL_csf="NfL (ln)",
                log_GFAP_csf="GFAP (ln)",
                log_pTau181_csf="pTau181 (ln)",
                log_pTau217_csf="pTau217 (ln)",
                log_Ab3840_csf="Aβ38/40 ratio (ln)",
                log_Ab4240_csf="Aβ42/40 ratio (ln)",
                log_Ab42_csf="Aβ42 (ln)",
                log_Ab40_csf="Aβ40 (ln)",
                log_Ab38_csf="Aβ38 (ln)",
                surv_time_months="Follow-up duration, months",
                ever_Riluzole="Riluzole exposure (ever vs never)",
                ever_Edaravone="Edaravone exposure (ever vs never)"
                )

# Modify T1label for binary variables
T1label_mod <- T1label
T1label_mod$Sex <- "Sex (female)"
T1label_mod$Ab_status <- "Aβ positive"
T1label_mod$OnsetSite <- "Bulbar onset"

T1digits <- list(Age_init =c(0,0), 
                 MMSE = c(0,0),
                 FAB = c(0,0),
                 MOCA = c(0,0),
                 RCPM = c(0,0),
                 Education_Years=c(0,0),
                 ALSFRS_init = c(0,0),
                 BMI_init = c(1,1))

T1label_bio <- list(
  Alb = "Alb",
  Cr = "Cr",
  CK = "CK",
  Ab38_csf="Aβ38",
  Ab40_csf="Aβ40",
  Ab42_csf="Aβ42",
  Ab42_40_csf_bridged="Aβ42/40",
  Ab38_40_csf="Aβ38/40",
  pTau181_csf="pTau181",
  pTau217_csf="pTau217",
  GFAP_csf_pgml="GFAP",
  NfL_csf_pgml="NfL",
  log_NfL_csf="NfL (ln)",
  log_GFAP_csf="GFAP (ln)",
  log_pTau181_csf="pTau181 (ln)",
  log_pTau217_csf="pTau217 (ln)",
  log_Ab3840_csf="Aβ38/40 ratio (ln)",
  log_Ab4240_csf="Aβ42/40 ratio (ln)",
  log_Ab42_csf="Aβ42 (ln)",
  log_Ab40_csf="Aβ40 (ln)",
  log_Ab38_csf="Aβ38 (ln)",
  VC_Percent="%VC at baseline",
  deltaFS="ΔFS at baseline",
  OnsetSite.bulbar="Bulbar onset",
  REEC_definite.definite="REEC definite at baseline",
  SexM = "Sex (male)",
  Age_init = "Age at baseline"
)

T1label_bio_overwritten <- T1label
T1label_bio_overwritten[names(T1label_bio)] <- T1label_bio

variable_label_map <- unlist(T1label_bio)



## 📊 Table 1 ----------------------------------------------------------------

T1rough <- T1 %>%
  #ALS_label 0/1 to nonALS/ALS
  mutate(ALS_label = ifelse(ALS_label==1, "ALS", "DC")) %>%
  dplyr::select(ALS_label,Age_init,Sex,BMI_init,OnsetSite,ALSFamHistory,
                #genetic_ALS,
                genetic_label,
                ALSFRS_init, DiseaseDur_init,R_ElEscorial,VC_Percent,deltaFS,
                MMSE,FAB,
                ALSFamHistory,
                APOE_group,
                Ab_status
  ) %>% 
  dplyr::mutate(
    R_ElEscorial = forcats::fct_recode(
      R_ElEscorial,
      "probable, laboratory supported" = "probable_laboratory_supported"
    ),
    genetic_label = case_when(
      genetic_label == "positive" ~ "pathogenic variant*",
      genetic_label == "negative" ~ "no known pathogenic variant",
      is.na(genetic_label) ~ "test not consented"
    ),
    APOE_group = case_when(
      APOE_group == "e4_carrier" ~ "ε4 carrier",
      APOE_group == "e2_carrier" ~ "ε2 carrier",
      APOE_group == "e3e3" ~ "ε3 homozygote",
      is.na(APOE_group) ~ "genotyping not performed"
    ),
    dplyr::across(
      c(VC_Percent,deltaFS,ALSFRS_init,OnsetSite, DiseaseDur_init, R_ElEscorial, 
        ALSFamHistory, #genetic_ALS, 
        genetic_label,
        APOE_group
      ),
      ~ replace(.x, !ALS_label %in% c("ALS"), NA))
  ) %>% 
  tbl_summary(by = ALS_label,
              missing = "no",
              digits = T1digits,
              label = T1label_mod,
              value = list(
                Sex ~ "F",
                Ab_status ~ "positive",
                OnsetSite ~ "bulbar"
              )) %>%
  bold_labels() %>%
  modify_table_body(
    ~ .x %>%
      # NA to blank
      mutate(across(
        matches("^stat_"),
        ~ ifelse(. %in% c("0 (NA%)", "NA (NA, NA)"), "", .)
      ))
  )#%>% add_n() 

dc_breakdown <- T1 %>%
  filter(ALS_label == 0) %>%  # Select only DC group
  mutate(category2 = str_replace(category2, "ALS mimics \\(broader\\)", "ALS mimics")) %>%
  count(category2) %>%
  arrange(desc(n))

genetic_breakdown <- T1 %>%
  filter(ALS_label == "1") %>%  # Select only ALS group)
  count(genetic_ALS) %>%
  filter(!str_detect(genetic_ALS, "not performed")) %>%
  filter(!str_detect(genetic_ALS, "no known pathogenic variant")) %>%
  arrange(desc(n))

# Create legend text in format: category(n=x)
legend_text_dc_breakdown_raw <- paste0(
  dc_breakdown$category2, 
  " (n=", 
  dc_breakdown$n, 
  ")",
  collapse = ", "
)

legend_text_dc_breakdown <- paste0("DC (disease controls): ", 
                                   legend_text_dc_breakdown_raw,
                                   ". A detailed breakdown of ALS mimic diagnoses is shown in Supplementary Table S2.")

legend_text_genetic_breakdown <- paste0(
  genetic_breakdown$genetic_ALS, 
  " (n=", genetic_breakdown$n, ")", collapse = ", "
)
legend_text_genetic_breakdown <- paste0("*) ", 
                                        legend_text_genetic_breakdown,
                                        ".",
                                        "\n†) ε4 carriers include ε2/ε4 (n=4) and ε4/ε4 (n=2).")


T1rough %>% as_flex_table() %>% 
  flextable::add_footer_lines(values = c(legend_text_genetic_breakdown,
                                         legend_text_dc_breakdown))

Table1_pub <- T1rough %>% as_flex_table() %>% 
  flextable::add_footer_lines(values = c(legend_text_genetic_breakdown,
                                         legend_text_dc_breakdown)) %>% 
  set_caption(caption = flextable::as_paragraph(
    flextable::as_chunk("Table 1: Baseline Demographics of ALS patients and Disease Controls",
                        props = officer::fp_text(font.size   = 12, bold = TRUE, 
                                                 font.family = "Times New Roman"))),
    word_stylename = "Table Caption", 
    fp_p = officer::fp_par(text.align = "left", padding= 5)) %>% 
  flextable::font(fontname = "Times New Roman", part = "all") %>% 
  flextable::fontsize(size=9, part="all")


Table1_pub

# Demographic ---------------------------------------------------------

## 🔍Supple T1 ---------------------------------------------------------
nd_breakdown <- legend_text_dc_breakdown_raw %>%
  str_replace("ALS mimics \\(n=[0-9]+\\),? ?", "") %>%
  str_trim()

nd_caption <- paste0(
  "ND: neurological disease controls excluding ALS mimics, comprising ",
  nd_breakdown,
  "."
)

# Create table for ALS
T1rough_als <- T1 %>%
  dplyr::filter(category2 == "ALS") %>%
  dplyr::mutate(als_group = "ALS") %>%
  dplyr::select(
    als_group,
    Age_init, Sex, BMI_init,
    MMSE, FAB, Education_Years, Ab_status
  ) %>%
  tbl_summary(
    by = als_group,
    missing = "no",
    digits = T1digits,
    label = T1label_mod,
    type = list(Education_Years ~ "continuous"),
    value = list(Sex ~ "F", Ab_status ~ "positive")
  ) %>%
  modify_table_body(
    ~ .x %>%
      mutate(across(
        matches("^stat_"),
        ~ ifelse(. %in% c("0 (NA%)", "NA (NA, NA)"), "", .)
      ))
  )

# Create table for ALS mimics
T1rough_mimics <- T1 %>%
  dplyr::filter(category2 == "ALS mimics (broader)") %>%
  dplyr::mutate(mimics_group = "ALS mimics") %>%
  dplyr::select(
    mimics_group,
    Age_init, Sex, BMI_init,
    MMSE, FAB, Education_Years, Ab_status
  ) %>%
  tbl_summary(
    by = mimics_group,
    missing = "no",
    digits = T1digits,
    label = T1label_mod,
    type = list(Education_Years ~ "continuous"),
    value = list(Sex ~ "F", Ab_status ~ "positive")
  ) %>%
  modify_table_body(
    ~ .x %>%
      mutate(across(
        matches("^stat_"),
        ~ ifelse(. %in% c("0 (NA%)", "NA (NA, NA)"), "", .)
      ))
  )

# Create ND overall table
T1rough_nd_overall <- T1 %>%
  dplyr::filter(category2 %in% c("CIDP", "MSA", "PD", "iNPH")) %>%
  dplyr::mutate(nd_group = "ND (overall)") %>%
  dplyr::select(
    nd_group,
    Age_init, Sex, BMI_init,
    MMSE, FAB, Education_Years, Ab_status
  ) %>%
  tbl_summary(
    by = nd_group,
    missing = "no",
    digits = T1digits,
    label = T1label_mod,
    type = list(Education_Years ~ "continuous"),
    value = list(Sex ~ "F", Ab_status ~ "positive")
  ) %>%
  modify_table_body(
    ~ .x %>%
      mutate(across(
        matches("^stat_"),
        ~ ifelse(. %in% c("0 (NA%)", "NA (NA, NA)"), "", .)
      ))
  )

# Create ND subcategory table - reversed order: PD, MSA, iNPH, CIDP
T1rough_nd_sub <- T1 %>%
  dplyr::filter(category2 %in% c("CIDP", "MSA", "PD", "iNPH")) %>%
  dplyr::mutate(
    category2 = factor(category2, levels = c("PD", "MSA", "iNPH", "CIDP"))
  ) %>%
  dplyr::select(
    category2,
    Age_init, Sex, BMI_init,
    MMSE, FAB, Education_Years, Ab_status
  ) %>%
  tbl_summary(
    by = category2,
    missing = "no",
    digits = T1digits,
    label = T1label_mod,
    type = list(Education_Years ~ "continuous"),
    value = list(Sex ~ "F", Ab_status ~ "positive")
  ) %>%
  modify_table_body(
    ~ .x %>%
      mutate(across(
        matches("^stat_"),
        ~ ifelse(. %in% c("0 (NA%)", "NA (NA, NA)"), "", .)
      ))
  )

# Merge all tables horizontally
T1rough_combined <- tbl_merge(
  tbls = list(T1rough_als, T1rough_mimics, T1rough_nd_overall, T1rough_nd_sub),
  tab_spanner = c("", "", "", "")
)

# Convert to flextable
flex_table <- T1rough_combined %>% as_flex_table()

# Add custom header structure with uniform borders
SuppleT1_pub <- flex_table %>%
  flextable::delete_part(part = "header") %>%
  flextable::add_header_row(
    values = c("Characteristic", "ALS", "ALS mimics", "ND (overall)", "PD", "MSA", "iNPH", "CIDP"),
    colwidths = c(1, 1, 1, 1, 1, 1, 1, 1)
  ) %>%
  flextable::add_header_row(
    values = c("Characteristic", "ALS", "ALS mimics", "ND", "ND", "ND", "ND", "ND"),
    colwidths = c(1, 1, 1, 1, 1, 1, 1, 1),
    top = TRUE
  ) %>%
  flextable::merge_at(i = 1, j = 4:8, part = "header") %>%
  flextable::add_header_row(
    values = c("", paste0("N = ", table(T1$category2)["ALS"]), 
               paste0("N = ", table(T1$category2)["ALS mimics (broader)"]),
               paste0("N = ", sum(table(T1$category2)[c("CIDP", "iNPH", "MSA", "PD")])),
               paste0("N = ", table(T1$category2)["PD"]),
               paste0("N = ", table(T1$category2)["MSA"]),
               paste0("N = ", table(T1$category2)["iNPH"]),
               paste0("N = ", table(T1$category2)["CIDP"])),
    colwidths = c(1, 1, 1, 1, 1, 1, 1, 1),
    top = FALSE
  ) %>%
  flextable::hline_top(part = "header", border = officer::fp_border(width = 1)) %>%
  flextable::hline(i = 1, part = "header", border = officer::fp_border(width = 1)) %>%
  flextable::hline(i = 3, part = "header", border = officer::fp_border(width = 1)) %>%
  flextable::add_footer_lines(
    values = c(
      nd_caption,
      "Abbreviations: ALS, amyotrophic lateral sclerosis; BMI, body mass index; MMSE, Mini-Mental State Examination; FAB, Frontal Assessment Battery; Aβ, amyloid-β; PD, Parkinson’s disease; MSA, multiple system atrophy; iNPH, idiopathic normal pressure hydrocephalus; CIDP, chronic inflammatory demyelinating polyneuropathy."
    )
  ) %>%
  set_caption(
    caption = flextable::as_paragraph(
      flextable::as_chunk(
        "Supplementary Table S1: Baseline demographics by disease control categories",
        props = officer::fp_text(
          font.size = 12,
          bold = TRUE,
          font.family = "Times New Roman"
        )
      )
    ),
    word_stylename = "Table Caption",
    fp_p = officer::fp_par(text.align = "left", padding = 5)
  ) %>%
  flextable::font(fontname = "Times New Roman", part = "all") %>%
  flextable::fontsize(size = 8, part = "all") %>%
  flextable::align(align = "center", part = "header") %>%
  flextable::align(align = "left", part = "body", j = 1) %>%
  flextable::autofit()

SuppleT1_pub

## 🔍Supple T2: ALS mimics -----------------------------------------------------

tbl_mimics_summary <- flextable(mimics_summary_renamed) %>%
  autofit() %>%
  bold(part = "header") %>%
  bold(i = mimics_body_n, part = "body") %>%              
  border_remove() %>%
  hline_top(part = "header", border = officer::fp_border(color = "black", width = 1)) %>%
  hline_bottom(part = "header", border = officer::fp_border(color = "black", width = 1)) %>%
  hline(i = mimics_body_n - 1, part = "body", border = officer::fp_border(color = "black", width = 1)) %>%  
  hline_bottom(part = "body", border = officer::fp_border(color = "black", width = 1)) 

SuppleT2_pub <- tbl_mimics_summary %>% 
  set_caption(caption = flextable::as_paragraph(
    flextable::as_chunk("Supplementary Table S2: Diagnoses among ALS mimics",
                        props = officer::fp_text(font.size   = 12, bold = TRUE, 
                                                 font.family = "Times New Roman"))),
    word_stylename = "Table Caption", fp_p = officer::fp_par(text.align = "left", padding= 5)) %>% 
  flextable::font(fontname = "Times New Roman", part = "all") %>% 
  flextable::fontsize(size=8, part="all")

SuppleT2_pub

## 🔍Supple T3 -----------------------------------------------------------------

T1roughAB <- T1 %>%
  dplyr::select(Apos_label,Age_init,Sex,BMI_init,OnsetSite,
                ALSFRS_init, DiseaseDur_init,R_ElEscorial,VC_Percent,deltaFS,
                MMSE,FAB,Education_Years,
                #ALSFamHistory, 
                genetic_label, APOE_group,
  ) %>% 
  dplyr::mutate(
    R_ElEscorial = forcats::fct_recode(
      R_ElEscorial,
      "probable, laboratory supported" = "probable_laboratory_supported"
    ),
    genetic_label = case_when(
      genetic_label == "positive" ~ "pathogenic variant*",
      genetic_label == "negative" ~ "no known pathogenic variant",
      is.na(genetic_label) ~ "test not consented"
    ),
    APOE_group = case_when(
      APOE_group == "e4_carrier" ~ "ε4 carrier",
      APOE_group == "e2_carrier" ~ "ε2 carrier",
      APOE_group == "e3e3" ~ "ε3 homozygote",
      is.na(APOE_group) ~ "genotyping not performed"
    ),
    dplyr::across(
      c(VC_Percent,deltaFS,ALSFRS_init,OnsetSite, DiseaseDur_init, R_ElEscorial, #ALSFamHistory
        ,genetic_label,APOE_group
      ),
      ~ replace(.x, !Apos_label %in% c("ALS (Aβ+)", "ALS (Aβ-)"), NA)),
  ) %>% 
  tbl_summary(by = Apos_label,
              missing = "no",
              digits = T1digits,
              label = T1label_mod,
              value = list(
                Sex ~ "F",
                OnsetSite ~ "bulbar"
              )) %>% add_n() %>% 
  bold_labels() %>%
  modify_table_body(
    ~ .x %>%
      # NA -> blank
      mutate(across(
        matches("^stat_"),
        ~ ifelse(. %in% c("0 (NA%)", "NA (NA, NA)"), "", .)
      ))
  ) 

SuppleT3_pub <- T1roughAB %>% as_flex_table() %>%
  flextable::add_footer_lines(values = c(legend_text_genetic_breakdown,
                                         legend_text_dc_breakdown)) %>% 
  set_caption(caption = flextable::as_paragraph(
    flextable::as_chunk("Supplementary Table S3: Baseline demographics of patients by Aβ status",
                        props = officer::fp_text(font.size   = 12, bold = TRUE, 
                                                 font.family = "Times New Roman"))),
    word_stylename = "Table Caption", fp_p = officer::fp_par(text.align = "left", padding= 5)) %>% 
  flextable::font(fontname = "Times New Roman", part = "all") %>% 
  flextable::fontsize(size=8, part="all")

SuppleT3_pub

# Biomarker comparison tables ---------------------------------------------------------

## 🔍Supple T4 --------------------------------------------------------------------

T2rough <- T1 %>%
  #ALS_label 0/1 to nonALS/ALS
  mutate(ALS_label = ifelse(ALS_label==1, "ALS", "DC")) %>%
  dplyr::select(Age_init, 
                category3,
                Ab38_csf,Ab40_csf,Ab42_csf,
                Ab38_40_csf, 
                Ab42_40_csf_bridged, 
                GFAP_csf_pgml, 
                NfL_csf_pgml, 
                pTau181_csf,pTau217_csf,
                Cr, CK, Alb, Cr_CysC_ratio) %>% 
  tbl_summary(by=category3,
              missing="no",
              #missing_stat = "{N_miss}",
              digits=T1digits,
              label= T1label) %>% add_n() %>% 
  add_p(pvalue_fun = ~ style_pvalue(.x, digits = 2)) %>% add_q(method = "fdr") %>% 
  modify_header(q.value ~ "**p-value (FDR)**") #%>% bold_labels() 

T2rough %>% as_flex_table() %>% flextable::add_footer_lines(values = nd_caption)

SuppleT4_pub <- T2rough %>% as_flex_table() %>% 
  flextable::add_footer_lines(values = nd_caption) %>% 
  set_caption(caption = flextable::as_paragraph(
    flextable::as_chunk("Supplementary Table S4: Biomarker comparison by disease control categories",
                        props = officer::fp_text(font.size   = 12, bold = TRUE, 
                                                 font.family = "Times New Roman"))),
    word_stylename = "Table Caption", 
    fp_p = officer::fp_par(text.align = "left", padding= 5)) %>%
  flextable::font(fontname = "Times New Roman", part = "all") %>%
  flextable::fontsize(size=9, part="all")

SuppleT4_pub


## 🔍SuppleT5  -------------------------------------------------------------

T2ab_Aneg <- T1 %>% filter(Ab_status=="negative") %>% 
  dplyr::select(Apos_label, Age_init,
                Ab38_csf,Ab40_csf, Ab42_csf,
                Ab38_40_csf,
                Ab42_40_csf_bridged,
                GFAP_csf_pgml, NfL_csf_pgml, 
                pTau181_csf,pTau217_csf,
                Cr, CK, Alb, Cr_CysC_ratio) %>% 
  tbl_summary(by=Apos_label,
              missing="no",
              #missing_stat = "{N_miss}",
              digits=T1digits,
              label= T1label) %>% add_n() %>% add_p(pvalue_fun = ~ style_pvalue(.x, digits = 2)) %>% add_q(method="fdr") %>% #%>% bold_labels() #%>% add_n() %>% add_p(pvalue_fun = ~ style_pvalue(.x, digits = 2))#%>% add_overall()
  modify_header(q.value ~ "**p-value (FDR)**") %>% 
  modify_footnote(q.value ~ "FDR-adjusted p-values (Benjamini–Hochberg)")

T2ab_Aneg %>% as_flex_table()

SuppleT5_pub <- T2ab_Aneg %>% as_flex_table() %>% 
  set_caption(caption = flextable::as_paragraph(
    flextable::as_chunk("Supplementary Table S5: Biomarker comparison within Aβ-negative subgroup",
                        props = officer::fp_text(font.size   = 12, bold = TRUE, 
                                                 font.family = "Times New Roman"))),
    word_stylename = "Table Caption", 
    fp_p = officer::fp_par(text.align = "left", padding= 5)) %>% 
  flextable::font(fontname = "Times New Roman", part = "all") %>% 
  flextable::fontsize(size=8, part="all")

SuppleT5_pub


## 🔍SuppleT6 ----------------------------------------------------------------------

SuppleT_Aneg_bm_all <- T1 %>%
  mutate(category2 = dplyr::recode(category2,
                                   "ALS mimics (broader)"="ALS mimics")) %>%
  # Negative: category2 all ALS negative; put A(-)
  filter(Ab_status=="negative") %>% 
  mutate(category2 = paste0(category2, " (Aβ-)")) %>%
  dplyr::select(category2, 
                Age_csf,
                Ab38_csf,Ab40_csf,Ab42_csf,
                Ab38_40_csf,
                Ab42_40_csf_bridged, 
                GFAP_csf_pgml, 
                NfL_csf_pgml, 
                pTau181_csf,pTau217_csf,
                Cr, CK, Alb,Cr_CysC_ratio
  ) %>% 
  tbl_summary(by=category2,
              missing="no",
              #missing_stat = "{N_miss}",
              digits=T1digits,
              label= T1label) %>% add_n() %>% add_p(pvalue_fun = ~ style_pvalue(.x, digits = 2)) %>% add_q(method = "fdr") %>%
  modify_header(q.value ~ "**p-value (FDR)**") %>% 
  bold_labels() %>% 
  modify_table_body(
    ~ .x %>%
      # NA -> blank
      mutate(across(
        matches("^stat_"),
        ~ ifelse(. %in% c("0 (NA%)", "NA (NA, NA)"), "", .)
      ))
  )

SuppleT6_pub <- SuppleT_Aneg_bm_all %>% as_flex_table() %>% 
  flextable::width(j = 1:ncol_keys(.), width = 0.7) %>% 
  flextable::width(j = c(3,4,5,6,7,8), width = 1.0) %>%
  flextable::width(j = 1, width = 1.5) %>%
  flextable::padding(padding.left = 0, padding.right = 0, part = "all") %>% 
  set_caption(caption = flextable::as_paragraph(
    flextable::as_chunk("Supplementary Table S6: Biomarker comparison within Aβ-negative subgroup across all disease categories",
                        props = officer::fp_text(font.size   = 12, bold = TRUE, 
                                                 font.family = "Times New Roman"))),
    word_stylename = "Table Caption", 
    fp_p = officer::fp_par(text.align = "left", padding= 5)) %>% 
  flextable::font(fontname = "Times New Roman", part = "all") %>% 
  flextable::fontsize(size=9, part="all")

SuppleT6_pub


# APOE4 odds ratio ---------------------------------------------------

# logistic regression APOE e4 carrier vs Ab status in ALS
ApoEe4_Ab_ALS_logistic <-glm(Ab_status ~ APOE_e4_carrier + Age_init, 
                             data=csfdf %>% filter(ALS_label==1 & !SampleID %in% second) , family=binomial) 

ApoEe4_Ab_ALS_logistic %>% tbl_regression(exponentiate = TRUE,label = T1label) %>% 
  add_n() %>% modify_caption(caption = "CSF Aβ status in ALS ~ explanatory variables.") %>% as_flex_table()


# sensitivity analysis
ApoEe4_Ab_ALS_logistic_e33 <-glm(Ab_status ~ APOE_e4_vs_e33 + Age_init, 
                                 data=csfdf %>% filter(ALS_label==1 & !SampleID %in% second), family=binomial) 

ApoEe4_Ab_ALS_logistic_e33 %>% tbl_regression(exponentiate = TRUE,label = T1label) %>% 
  add_n() %>% modify_caption(caption = "CSF Aβ status in ALS ~ explanatory variables (logistic model).") %>% as_flex_table()



# ALS genetic +/- ---------------------------------------------------------
# Additional tables generated during peer review
# These tables were prepared in response to reviewer comments and are not part of the main analysis pipeline.

# ALS_genetic_demo <- T1 %>% filter(category=="ALS") %>% 
#   #genetic label rename
#   mutate(genetic_label = case_when(
#     genetic_label == "positive" ~ "pathogenic variant",
#     genetic_label == "negative" ~ "no known pathogenic variant",
#     is.na(genetic_label) ~ NA
#   )) %>%
#   dplyr::select(genetic_label,
#                 Age_init,Sex,BMI_init,OnsetSite,
#                 ALSFamHistory, #genetic_ALS,
#                 ALSFRS_init, DiseaseDur_init, R_ElEscorial, deltaFS,
#                 Education_Years, MMSE,FAB
#   ) %>% 
#   tbl_summary(by=genetic_label,
#               missing="no",
#               #missing_stat = "{N_miss}",
#               digits=T1digits,
#               label= T1label_mod,
#               value = list(
#                 Sex ~ "F",
#                 OnsetSite ~ "bulbar"
#               )
#   ) %>% bold_labels() %>% add_p(pvalue_fun = ~ style_pvalue(.x, digits = 2)) %>%   add_q(method = "fdr") %>%
#   modify_header(q.value ~ "**p-value (FDR)**") %>% 
#   modify_footnote(q.value ~ "FDR-adjusted p-values (Benjamini–Hochberg)") %>% 
#   modify_table_body(
#     ~ .x %>%
#       # NA -> blank
#       mutate(across(
#         matches("^stat_"),
#         ~ ifelse(. %in% c("0 (NA%)", "NA (NA, NA)"), "-", .)
#       ))
#   ) %>% add_n()
# 
# ALS_genetic_demo %>% as_flex_table()
# 
# 
# SuppleTable_rev2_5_2 <- ALS_genetic_demo %>% as_flex_table() %>% 
#   set_caption(caption = flextable::as_paragraph(
#     flextable::as_chunk("Demographic and clinical characteristics of ALS patients by genetic status",
#                         props = officer::fp_text(font.size   = 12, bold = TRUE, 
#                                                  font.family = "Times New Roman"))),
#     word_stylename = "Table Caption", fp_p = officer::fp_par(text.align = "left", padding= 5)) %>% 
#   flextable::font(fontname = "Times New Roman", part = "all") %>% 
#   flextable::fontsize(size=8, part="all")
# 
# SuppleTable_rev2_5_2
# 
# 
# ALS_genetic_bm <- T1 %>% filter(category=="ALS") %>% 
#   mutate(genetic_label = case_when(
#     genetic_label == "positive" ~ "pathogenic variant",
#     genetic_label == "negative" ~ "no known pathogenic variant",
#     is.na(genetic_label) ~ NA
#   )) %>%
#   dplyr::select(genetic_label,
#                 Age_init,
#                 Ab38_csf,Ab40_csf,Ab42_csf,
#                 Ab38_40_csf,Ab42_40_csf_bridged,Ab_status,
#                 GFAP_csf_pgml, NfL_csf_pgml, pTau181_csf,pTau217_csf,
#                 Cr, CK, Alb, Cr_CysC_ratio
#   ) %>% 
#   tbl_summary(by=genetic_label,
#               missing="no",
#               #missing_stat = "{N_miss}",
#               digits=T1digits,
#               label= T1label,
#   ) %>% bold_labels() %>% add_p(pvalue_fun = ~ style_pvalue(.x, digits = 2)) %>%   add_q(method = "fdr") %>%
#   modify_header(q.value ~ "**p-value (FDR)**") %>% 
#   modify_footnote(q.value ~ "FDR-adjusted p-values (Benjamini–Hochberg)") %>% 
#   modify_table_body(
#     ~ .x %>%
#       # NA -> blank
#       mutate(across(
#         matches("^stat_"),
#         ~ ifelse(. %in% c("0 (NA%)", "NA (NA, NA)"), "-", .)
#       ))
#   ) %>% add_n()
# 
# 
# SuppleTable_rev2_5_1 <- ALS_genetic_bm %>% as_flex_table() %>% 
#   set_caption(caption = flextable::as_paragraph(
#     flextable::as_chunk("Biomarker comparison of ALS patients by genetic status",
#                         props = officer::fp_text(font.size   = 12, bold = TRUE, 
#                                                  font.family = "Times New Roman"))),
#     word_stylename = "Table Caption", fp_p = officer::fp_par(text.align = "left", padding= 5)) %>% 
#   flextable::font(fontname = "Times New Roman", part = "all") %>% 
#   flextable::fontsize(size=8, part="all")
# 
# SuppleTable_rev2_5_1 
# 
# # Export
# # SuppleTable_rev2_5_1 %>% save_as_docx(path = "output/SuppleTable_rev2_5_1.docx")
