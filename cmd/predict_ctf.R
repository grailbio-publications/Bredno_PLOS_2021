#!/usr/bin/env Rscript
requireNamespace("dplyr", quietly = TRUE)
`%>%` <- magrittr::`%>%`
source("../R/biophysical_modeling.R")

user_input <- function(prompt) {
  if (interactive()) {
    entered <- readline(prompt)
  } else {
    cat(sprintf("%s\n", prompt))
    entered <- readLines("stdin", n = 1)
  }
}

main <- function() {
  # Get user selection of the cancer type. Read resepective models and
  # non-index lesion size imputation information.
  cancer_type <- user_input(
    "Cancer type\n[b]reast cancer / [l]ung cancer [c]olorectal cancer")
  stopifnot(cancer_type %in% c("b", "l", "c"))
  if (cancer_type == "b") {
    l_imputation <- readRDS("bc_impute_size.Rds")
    prediction_model <- readRDS("bc_prediction_model.Rds")
    quantitative_model <- readRDS("bc_quantitative_model.Rds")
    size_dim <- 3
  }
  if (cancer_type == "l") {
    l_imputation <- readRDS("luc_impute_size.Rds")
    prediction_model <- readRDS("luc_prediction_model.Rds")
    quantitative_model <- readRDS("luc_quantitative_model.Rds")
    size_dim <- 3
  }
  if (cancer_type == "c") {
    l_imputation <- readRDS("crc_impute_size.Rds")
    prediction_shallow <- readRDS("crc_shallow_prediction_model.Rds")
    prediction_deep <- readRDS("crc_deep_prediction_model.Rds")
    quantitative_shallow <- readRDS("crc_shallow_quantitative_model.Rds")
    quantitative_deep <- readRDS("crc_deep_quantitative_model.Rds")
    size_dim <- 2
  }

  # Collect user input for the size of all primary lesion foci and derive total
  # size of primary tumor volume (lung and breast cancer) or surface area
  # (colorectal cancer).
  n_primary_lesions <- as.integer(user_input(
    "Number of primary lesions: "
  ))
  stopifnot(!is.na(n_primary_lesions))
  stopifnot(n_primary_lesions > 0)
  stopifnot(n_primary_lesions <= 6)
  df_case <- data.frame(n_primary_lesions = n_primary_lesions)
  for (lesion_idx in 1:n_primary_lesions) {
    lesion_size <- as.numeric(user_input(
      sprintf("Max. size of lesion %d in mm (n for unknown): ", lesion_idx)))
    stopifnot(!is.na(lesion_size) || lesion_idx > 1)
    df_case[[sprintf("size_%d", lesion_idx)]] <- lesion_size
  }
  l_ret <- derive_primary_size(df_case,
                               l_imputation,
                               draw_or_fixed = "fixed",
                               dimension = size_dim)

  # Query cancer-tye specific clinical parameters for the models.
  df_case <- l_ret$df_model
  if (cancer_type == "b") {
    ki67_pos <- as.numeric(user_input("% Ki-67 pos: "))
    stopifnot(!is.na(ki67_pos))
    stopifnot(ki67_pos >= 0)
    stopifnot(ki67_pos <= 100)
    df_case$ki67_pos <- ki67_pos
    df_case$tmitv <- df_case$primary_volume * df_case$ki67_pos / 100
  }
  if (cancer_type == "l") {
    fdg_suv <- as.numeric(user_input("PET SUV: "))
    stopifnot(!is.na(fdg_suv))
    stopifnot(fdg_suv >= 1)
    df_case$fdg_suv <- fdg_suv
    df_case$elg <- df_case$primary_volume * (df_case$fdg_suv - 1.0)
  }
  if (cancer_type == "c") {
    lookup_microinvasion <- tibble::tribble(
      ~quant_invasion, ~display_microinvasion, ~score_invasion,
      0, "Other/missing", NA,
      1, "Intramucosal", 0,
      2, "Submucosa", 0,
      3, "Muscularis propria", 0,
      4, "Subserosa", 1,
      5, "Penetrates serosa", 1,
      6, "Invades adjacent", 1) %>%
      dplyr::mutate(display = sprintf("[%d] %s",
                                      quant_invasion, display_microinvasion))
    microinvasion <- as.integer(user_input(
      sprintf("Depth of microinvasion from path report:\n%s\n",
              paste(lookup_microinvasion$display, collapse = "\n"))))
      stopifnot(!is.na(microinvasion))
      stopifnot(microinvasion >= 0)
      stopifnot(microinvasion <= 6)
      if (microinvasion == 0) {
        t_stage <- as.integer(user_input(
          "T-Stage (1-4): "))
        stopifnot(!is.na(t_stage))
        stopifnot(t_stage > 0)
        stopifnot(t_stage <= 4)
        is_deep <- t_stage >= 3
      } else {
        is_deep <- ((lookup_microinvasion %>%
                       dplyr::filter(quant_invasion == microinvasion) %>%
                       dplyr::pull(score_invasion)) == 1)
      }
      df_case$area_shallow <- ifelse(is_deep, 0, df_case$primary_area)
      df_case$area_deep <- ifelse(is_deep, df_case$primary_area, 0)
      if (is_deep) {
        prediction_model <- prediction_deep
        quantitative_model <- quantitative_deep
      } else {
        prediction_model <- prediction_shallow
        quantitative_model <- quantitative_shallow
      }
  }

  # Model prediction and confidence intervals for estimated tumor fraction and
  # estimated number of GEs in patient circulation
  ctf_prediction <- predict(prediction_model,
                            df_case,
                            interval = "confidence")
  quantitative_prediction <- predict(quantitative_model,
                                     df_case,
                                     interval = "confidence")
  ctf_prediction <- pmin(1.0, ctf_prediction)
  cat(sprintf(
    "The model predicts a circulating tumor fraction of %f\n  (95%% CI [%f; %f])\n",
    ctf_prediction[1], ctf_prediction[2], ctf_prediction[3]))
  cat(sprintf(
    "Estimated number of tumor-derived GEs in circulation: %.1f\n  (95%% CI [%.1f; %.1f])\n",
    quantitative_prediction[1],
    quantitative_prediction[2],
    quantitative_prediction[3]))
}

if (sys.nframe() == 0) {
  main()
}
