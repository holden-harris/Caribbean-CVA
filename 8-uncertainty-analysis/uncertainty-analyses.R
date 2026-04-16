################################################################################
##------------------------------------------------------------------------------
## Caribbean CVA
## Script 1. Bootstrap uncertainty analysis and leave-one-out influence analysis
##
## Purpose:
## - Read compiled Caribbean CVA inputs
## - Reproduce baseline Sensitivity, Exposure, and Vulnerability scores
## - Run bootstrap uncertainty analysis for overall vulnerability
## - Run leave-one-out influence analysis for Sensitivity attributes and
##   Exposure factors
## - Write analysis-ready outputs for QA and figure generation
##
## NOAA FCVA workflow
##
## Step 1. Read and standardize compiled score inputs
## Step 2. Reproduce baseline component and vulnerability scores
## Step 3. Bootstrap uncertainty for Sensitivity and Vulnerability
## Step 4. Run leave-one-out influence analysis
## Step 5. Save final tables for QA and figures
##
## Notes:
## - Sensitivity comes from reviewer-entered qualitative attribute scores
## - Exposure comes from finalized quantitative exposure factor scores
## - Bootstrap uncertainty resamples individual reviewer Sensitivity scores only
## - Exposure remains fixed at the baseline finalized score during bootstrap
## - Leave-one-out analyses are deterministic, not bootstrap-based
## - This script stops if reproduced baseline scores do not match the
##   finalized Caribbean CVA outputs

##------------------------------------------------------------------------------
## Setup

rm(list = ls()); gc()

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(tibble)

##------------------------------------------------------------------------------
## Directories

proj_dir <- "."
in_dir   <- file.path(proj_dir, "outputs", "final-scores-compiled")
out_dir  <- file.path(in_dir, "uncertainty-loo")

intermediate_dir <- file.path(out_dir, "intermediate")
final_dir        <- file.path(out_dir, "final-tables")

dir.create(out_dir,          recursive = TRUE, showWarnings = FALSE)
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_dir,        recursive = TRUE, showWarnings = FALSE)

##------------------------------------------------------------------------------
## Input files

f_reviewer_scores <- file.path(in_dir,
                               "final-attribute-scores",
                               "table_final_attribute_scores_all.csv")

f_attr_means <- file.path(in_dir,
                          "overall-vulnerability-rankings",
                          "attribute_means_uscar.csv")

f_vuln <- file.path(in_dir,
                          "overall-vulnerability-rankings",
                          "overall_vulnerability_scores_uscar.csv")

##------------------------------------------------------------------------------
## Output files

f_precheck_log     <- file.path(intermediate_dir, "01_preanalysis_check_log.csv")
f_base_comp_repro  <- file.path(intermediate_dir, "02_baseline_reproduced_component_scores.csv")
f_base_vuln_repro  <- file.path(intermediate_dir, "02_baseline_reproduced_vulnerability_scores.csv")

f_boot_stock <- file.path(final_dir, "table_bootstrap_uncertainty_stock.csv")
f_boot_sens  <- file.path(final_dir, "table_bootstrap_uncertainty_sensitivity_component.csv")
f_boot_iter  <- file.path(intermediate_dir, "table_bootstrap_iterations_long.csv")

f_loo_sens_long <- file.path(final_dir, "table_leave_one_out_sensitivity_long.csv")
f_loo_exp_long  <- file.path(final_dir, "table_leave_one_out_exposure_long.csv")
f_loo_sens_sum  <- file.path(final_dir, "table_leave_one_out_sensitivity_summary.csv")
f_loo_exp_sum   <- file.path(final_dir, "table_leave_one_out_exposure_summary.csv")

f_plot_exp  <- file.path(final_dir, "table_exposure_factor_scores_for_plot.csv")
f_plot_sens <- file.path(final_dir, "table_sensitivity_attribute_scores_for_plot.csv")

##------------------------------------------------------------------------------
## Analysis settings

n_boot               <- 10000
save_iteration_table <- FALSE
borderline_threshold <- 0.25
bootstrap_seed       <- 99

################################################################################
##------------------------------------------------------------------------------
## Helper functions

## FCVA logic model
## Thresholds match the main scoring script exactly:
##   Very High = more than 3 attribute means >= 3.5  (n_ge_3_5 > 3)
##   High      = more than 2 attribute means >= 3.0  (n_ge_3_0 > 2)
##   Moderate  = more than 2 attribute means >= 2.5  (n_ge_2_5 > 2)
##   Low       = all other cases

fcva_logic_model <- function(mean_scores) {
  n_ge_35 <- sum(mean_scores >= 3.5, na.rm = TRUE)
  n_ge_30 <- sum(mean_scores >= 3.0, na.rm = TRUE)
  n_ge_25 <- sum(mean_scores >= 2.5, na.rm = TRUE)

  component_rank <- dplyr::case_when(
    n_ge_35 > 3 ~ "Very High",
    n_ge_30 > 2 ~ "High",
    n_ge_25 > 2 ~ "Moderate",
    TRUE        ~ "Low"
  )

  component_score_numeric <- dplyr::case_when(
    component_rank == "Low"       ~ 1L,
    component_rank == "Moderate"  ~ 2L,
    component_rank == "High"      ~ 3L,
    component_rank == "Very High" ~ 4L
  )

  tibble::tibble(component_rank, component_score_numeric)
}

## Vulnerability rank from numeric product score
assign_vulnerability_rank <- function(vulnerability_score_numeric) {
  dplyr::case_when(
    vulnerability_score_numeric <= 3                                           ~ "Low",
    vulnerability_score_numeric >= 4  & vulnerability_score_numeric <= 6      ~ "Moderate",
    vulnerability_score_numeric >= 8  & vulnerability_score_numeric <= 9      ~ "High",
    vulnerability_score_numeric >= 12 & vulnerability_score_numeric <= 16     ~ "Very High",
    TRUE                                                                       ~ NA_character_
  )
}

## Column existence check
check_required_columns <- function(dat, required_cols, object_name) {
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop(paste0("Missing required columns in ", object_name, ": ",
                paste(missing_cols, collapse = ", ")))
  }
}

## Bootstrap one attribute: resample reviewer scores with replacement
bootstrap_one_attribute <- function(scores_vec) {
  samp <- sample(scores_vec, size = length(scores_vec), replace = TRUE)
  mean(samp, na.rm = TRUE)
}

## Bootstrap one stock
bootstrap_one_stock <- function(stock_name_i,
                                 sens_scores_std,
                                 baseline_exp_score_num,
                                 n_boot               = 10000,
                                 borderline_threshold = 0.25) {

  stock_scores <- sens_scores_std %>%
    dplyr::filter(stock_name == stock_name_i)

  att_names <- sort(unique(stock_scores$attribute_name))

  ## Pre-split scores by attribute to avoid repeated filter calls in the loop
  scores_by_att <- split(stock_scores$score, stock_scores$attribute_name)

  iter_results <- purrr::map(seq_len(n_boot), function(i) {
    boot_means <- vapply(att_names, function(att) {
      v <- scores_by_att[[att]]
      if (length(v) == 0 || all(is.na(v))) return(NA_real_)
      bootstrap_one_attribute(v)
    }, numeric(1))

    logic    <- fcva_logic_model(boot_means)
    vuln_num <- logic$component_score_numeric * baseline_exp_score_num

    tibble::tibble(
      sens_score_numeric = logic$component_score_numeric,
      sens_rank          = logic$component_rank,
      vuln_score         = vuln_num,
      vuln_rank          = assign_vulnerability_rank(vuln_num)
    )
  })

  iter_df <- dplyr::bind_rows(iter_results) %>%
    dplyr::mutate(stock_name = stock_name_i)

  stock_summary <- iter_df %>%
    dplyr::count(vuln_rank) %>%
    dplyr::mutate(
      prop       = n / n_boot,
      stock_name = stock_name_i
    )

  top_rank_prop <- max(stock_summary$prop)
  is_borderline <- top_rank_prop < (1 - borderline_threshold)

  sens_summary <- iter_df %>%
    dplyr::count(sens_rank) %>%
    dplyr::mutate(
      prop       = n / n_boot,
      stock_name = stock_name_i
    )

  list(
    iter_df       = iter_df,
    stock_summary = stock_summary %>% dplyr::mutate(borderline = is_borderline),
    sens_summary  = sens_summary
  )
}

## LOO: sensitivity — omit one attribute at a time
run_loo_sensitivity_one_stock <- function(stock_name_i,
                                          sens_means_std,
                                          baseline_exp_score_num,
                                          baseline_sens_rank,
                                          baseline_vuln_rank,
                                          baseline_vuln_score_num) {

  stock_means <- sens_means_std %>%
    dplyr::filter(stock_name == stock_name_i)

  att_names <- unique(stock_means$attribute_name)

  purrr::map_dfr(att_names, function(att_omit) {
    reduced       <- stock_means %>% dplyr::filter(attribute_name != att_omit)
    logic         <- fcva_logic_model(reduced$mean_score)
    new_vuln_num  <- logic$component_score_numeric * baseline_exp_score_num
    new_vuln_rank <- assign_vulnerability_rank(new_vuln_num)

    tibble::tibble(
      stock_name             = stock_name_i,
      attribute_omitted      = att_omit,
      baseline_sens_rank     = baseline_sens_rank,
      new_sens_rank          = logic$component_rank,
      new_sens_score_numeric = logic$component_score_numeric,
      baseline_vuln_rank     = baseline_vuln_rank,
      new_vuln_rank          = new_vuln_rank,
      baseline_vuln_score    = baseline_vuln_score_num,
      new_vuln_score         = new_vuln_num,
      rank_changed           = (new_vuln_rank != baseline_vuln_rank),
      rank_change_direction  = dplyr::case_when(
        new_vuln_score > baseline_vuln_score ~ "Higher",
        new_vuln_score < baseline_vuln_score ~ "Lower",
        TRUE                                 ~ "No change"
      )
    )
  })
}

## LOO: exposure — omit one factor at a time
run_loo_exposure_one_stock <- function(stock_name_i,
                                       exp_means_std,
                                       baseline_sens_score_num,
                                       baseline_exp_rank,
                                       baseline_vuln_rank,
                                       baseline_vuln_score_num) {

  stock_means  <- exp_means_std %>%
    dplyr::filter(stock_name == stock_name_i)

  factor_names <- unique(stock_means$exposure_factor)

  purrr::map_dfr(factor_names, function(fac_omit) {
    reduced       <- stock_means %>% dplyr::filter(exposure_factor != fac_omit)
    logic         <- fcva_logic_model(reduced$mean_score)
    new_vuln_num  <- baseline_sens_score_num * logic$component_score_numeric
    new_vuln_rank <- assign_vulnerability_rank(new_vuln_num)

    tibble::tibble(
      stock_name            = stock_name_i,
      factor_omitted        = fac_omit,
      baseline_exp_rank     = baseline_exp_rank,
      new_exp_rank          = logic$component_rank,
      new_exp_score_numeric = logic$component_score_numeric,
      baseline_vuln_rank    = baseline_vuln_rank,
      new_vuln_rank         = new_vuln_rank,
      baseline_vuln_score   = baseline_vuln_score_num,
      new_vuln_score        = new_vuln_num,
      rank_changed          = (new_vuln_rank != baseline_vuln_rank),
      rank_change_direction = dplyr::case_when(
        new_vuln_score > baseline_vuln_score ~ "Higher",
        new_vuln_score < baseline_vuln_score ~ "Lower",
        TRUE                                 ~ "No change"
      )
    )
  })
}

################################################################################
##------------------------------------------------------------------------------
## Step 1 - Read and standardize compiled score inputs

## Individual reviewer Sensitivity scores (for bootstrap)
reviewer_scores_raw <- readr::read_csv(f_reviewer_scores, show_col_types = FALSE)

## Attribute means (Sensitivity + Exposure combined; split below)
attr_means_raw <- readr::read_csv(f_attr_means, show_col_types = FALSE)

## Baseline vulnerability scores
vuln_raw <- readr::read_csv(f_vuln, show_col_types = FALSE)

##------------------------------------------------------------------------------
## Step 1A - Standardize and split inputs

## Individual reviewer Sensitivity + Rigidity scores
## Both attribute types feed into the Sensitivity component in attribute_means_uscar
## (rows 21–28 = Sensitivity, rows 30–35 = Rigidity, both labeled "Sensitivity" there)
sens_scores_std <- reviewer_scores_raw %>%
  dplyr::filter(Attribute_type %in% c("Sensitivity", "Rigidity")) %>%
  dplyr::transmute(
    stock_name     = stock_name,
    reviewer_id    = Scorer,
    attribute_name = str_squish(Attribute_name),
    score          = as.numeric(Final_score)
  ) %>%
  dplyr::filter(!is.na(score)) %>%
  dplyr::arrange(stock_name, reviewer_id, attribute_name)

## Sensitivity attribute means (qualitative reviewer scores)
sens_means_std <- attr_means_raw %>%
  dplyr::filter(attribute_type == "Sensitivity", score_type == "Qualitative") %>%
  dplyr::transmute(
    stock_name     = stock_name,
    attribute_name = str_squish(attribute_name),
    mean_score     = as.numeric(attribute_mean)
  ) %>%
  dplyr::arrange(stock_name, attribute_name)

## Exposure factor means (calculated quantitative scores)
exp_means_std <- attr_means_raw %>%
  dplyr::filter(attribute_type == "Exposure", score_type == "Calculated") %>%
  dplyr::transmute(
    stock_name      = stock_name,
    exposure_factor = str_squish(attribute_name),
    mean_score      = as.numeric(attribute_mean)
  ) %>%
  dplyr::arrange(stock_name, exposure_factor)

## Baseline vulnerability table — rename abbreviated columns
vuln_std <- vuln_raw %>%
  dplyr::rename(
    exposure_score_numeric     = Exp_score,
    exposure_rank              = Exp_rank,
    sensitivity_score_numeric  = Sens_score,
    sensitivity_rank           = Sens_rank,
    vulnerability_score_numeric = Vuln_score,
    vulnerability_rank         = Vuln_rank
  )

##------------------------------------------------------------------------------
## Step 1B - Check required columns

check_required_columns(
  dat           = sens_scores_std,
  required_cols = c("stock_name", "reviewer_id", "attribute_name", "score"),
  object_name   = "sens_scores_std"
)

check_required_columns(
  dat           = sens_means_std,
  required_cols = c("stock_name", "attribute_name", "mean_score"),
  object_name   = "sens_means_std"
)

check_required_columns(
  dat           = exp_means_std,
  required_cols = c("stock_name", "exposure_factor", "mean_score"),
  object_name   = "exp_means_std"
)

check_required_columns(
  dat           = vuln_std,
  required_cols = c("stock_name",
                    "exposure_score_numeric",   "exposure_rank",
                    "sensitivity_score_numeric", "sensitivity_rank",
                    "vulnerability_score_numeric", "vulnerability_rank"),
  object_name   = "vuln_std"
)

################################################################################
##------------------------------------------------------------------------------
## Step 2 - Pre-analysis checks

precheck_log <- tibble::tibble(
  check_name = character(),
  status     = character(),
  detail     = character()
)

## Duplicate reviewer × attribute rows
dup_sens_scores <- sens_scores_std %>%
  dplyr::count(stock_name, reviewer_id, attribute_name) %>%
  dplyr::filter(n > 1)

## Impossible score values (outside 1–4)
bad_sens_scores <- sens_scores_std %>%
  dplyr::filter(score < 1 | score > 4)

bad_sens_means <- sens_means_std %>%
  dplyr::filter(mean_score < 1 | mean_score > 4)

bad_exp_means <- exp_means_std %>%
  dplyr::filter(mean_score < 1 | mean_score > 4)

## Stocks in reviewer score table that are missing from baseline vulnerability
stocks_missing_in_vuln <- setdiff(
  unique(sens_scores_std$stock_name),
  unique(vuln_std$stock_name)
)

precheck_log <- dplyr::bind_rows(
  tibble::tibble(check_name = "duplicate_reviewer_attribute_rows",
                 status     = ifelse(nrow(dup_sens_scores) == 0, "PASS", "FAIL"),
                 detail     = paste("n =", nrow(dup_sens_scores))),
  tibble::tibble(check_name = "reviewer_scores_outside_1_to_4",
                 status     = ifelse(nrow(bad_sens_scores) == 0, "PASS", "FAIL"),
                 detail     = paste("n =", nrow(bad_sens_scores))),
  tibble::tibble(check_name = "sensitivity_means_outside_1_to_4",
                 status     = ifelse(nrow(bad_sens_means) == 0, "PASS", "FAIL"),
                 detail     = paste("n =", nrow(bad_sens_means))),
  tibble::tibble(check_name = "exposure_means_outside_1_to_4",
                 status     = ifelse(nrow(bad_exp_means) == 0, "PASS", "FAIL"),
                 detail     = paste("n =", nrow(bad_exp_means))),
  tibble::tibble(check_name = "stocks_missing_in_baseline_vulnerability_table",
                 status     = ifelse(length(stocks_missing_in_vuln) == 0, "PASS", "FAIL"),
                 detail     = paste("n =", length(stocks_missing_in_vuln)))
)

readr::write_csv(precheck_log, f_precheck_log)
print(precheck_log)

if (any(precheck_log$status == "FAIL")) {
  stop("Pre-analysis checks failed. Review ", f_precheck_log)
}

################################################################################
##------------------------------------------------------------------------------
## Step 3 - Reproduce baseline scores and verify against official outputs
##
## The reproduced baseline MUST match the finalized Caribbean CVA outputs exactly.
## If it does not, the script stops — which indicates the fcva_logic_model()
## thresholds above do not match those used in the main scoring script.

base_sens_repro <- sens_means_std %>%
  dplyr::group_by(stock_name) %>%
  dplyr::summarise(
    sens_logic = list(fcva_logic_model(mean_scores = mean_score)),
    .groups    = "drop"
  ) %>%
  tidyr::unnest(cols = c(sens_logic)) %>%
  dplyr::rename(
    sensitivity_rank_repro          = component_rank,
    sensitivity_score_numeric_repro = component_score_numeric
  )

base_exp_repro <- exp_means_std %>%
  dplyr::group_by(stock_name) %>%
  dplyr::summarise(
    exp_logic = list(fcva_logic_model(mean_scores = mean_score)),
    .groups   = "drop"
  ) %>%
  tidyr::unnest(cols = c(exp_logic)) %>%
  dplyr::rename(
    exposure_rank_repro          = component_rank,
    exposure_score_numeric_repro = component_score_numeric
  )

baseline_reproduced <- base_sens_repro %>%
  dplyr::left_join(base_exp_repro, by = "stock_name") %>%
  dplyr::mutate(
    vulnerability_score_numeric_repro =
      sensitivity_score_numeric_repro * exposure_score_numeric_repro,
    vulnerability_rank_repro =
      assign_vulnerability_rank(vulnerability_score_numeric_repro)
  )

baseline_compare <- vuln_std %>%
  dplyr::left_join(baseline_reproduced, by = "stock_name") %>%
  dplyr::mutate(
    sens_rank_match  = sensitivity_rank  == sensitivity_rank_repro,
    sens_num_match   = sensitivity_score_numeric  == sensitivity_score_numeric_repro,
    exp_rank_match   = exposure_rank     == exposure_rank_repro,
    exp_num_match    = exposure_score_numeric     == exposure_score_numeric_repro,
    vuln_rank_match  = vulnerability_rank  == vulnerability_rank_repro,
    vuln_num_match   = vulnerability_score_numeric == vulnerability_score_numeric_repro
  )

readr::write_csv(baseline_reproduced, f_base_vuln_repro)
readr::write_csv(
  baseline_reproduced %>%
    dplyr::select(stock_name,
                  exposure_rank_repro,   exposure_score_numeric_repro,
                  sensitivity_rank_repro, sensitivity_score_numeric_repro),
  f_base_comp_repro
)

mismatches <- baseline_compare %>%
  dplyr::filter(
    is.na(sens_rank_match) | !sens_rank_match |
    is.na(sens_num_match)  | !sens_num_match  |
    is.na(exp_rank_match)  | !exp_rank_match  |
    is.na(exp_num_match)   | !exp_num_match   |
    is.na(vuln_rank_match) | !vuln_rank_match |
    is.na(vuln_num_match)  | !vuln_num_match
  )

if (nrow(mismatches) > 0) {
  message("Mismatched stocks:")
  print(mismatches %>%
          dplyr::select(stock_name, dplyr::ends_with("_match")),
        n = 50)
  stop("Baseline reproduction failed. Review reproduced score tables in ",
       intermediate_dir)
}

message("✓ Baseline reproduced successfully for all ", nrow(vuln_std), " stocks.")

################################################################################
##------------------------------------------------------------------------------
## Step 4 - Bootstrap uncertainty analysis
##
## For each stock:
##   - Pool individual reviewer Sensitivity scores within each attribute
##   - Resample with replacement (n = n_reviewers per attribute per iteration)
##   - Recalculate bootstrap attribute mean scores
##   - Apply FCVA logic model for bootstrap Sensitivity component score
##   - Hold baseline Exposure component score fixed
##   - Recalculate vulnerability score and rank
## Summarize proportions across n_boot iterations.

set.seed(bootstrap_seed)

stock_list <- sort(unique(vuln_std$stock_name))

message("Running bootstrap (n_boot = ", n_boot, ") for ",
        length(stock_list), " stocks ...")

boot_results <- purrr::map(stock_list, function(s) {
  message("  → ", s)
  baseline_exp <- vuln_std %>%
    dplyr::filter(stock_name == s) %>%
    dplyr::pull(exposure_score_numeric)
  bootstrap_one_stock(
    stock_name_i          = s,
    sens_scores_std       = sens_scores_std,
    baseline_exp_score_num = baseline_exp,
    n_boot                = n_boot,
    borderline_threshold  = borderline_threshold
  )
})

## Collect summary tables
boot_stock_summary <- dplyr::bind_rows(
  purrr::map(boot_results, "stock_summary")
) %>%
  dplyr::arrange(stock_name, vuln_rank)

boot_sens_summary <- dplyr::bind_rows(
  purrr::map(boot_results, "sens_summary")
) %>%
  dplyr::arrange(stock_name, sens_rank)

readr::write_csv(boot_stock_summary, f_boot_stock)
readr::write_csv(boot_sens_summary,  f_boot_sens)

if (save_iteration_table) {
  boot_iter_long <- dplyr::bind_rows(purrr::map(boot_results, "iter_df"))
  readr::write_csv(boot_iter_long, f_boot_iter)
}

message("✓ Bootstrap complete. Tables written to ", final_dir)

################################################################################
##------------------------------------------------------------------------------
## Step 5A - Leave-one-out: Sensitivity attributes
##
## Omit one Sensitivity attribute at a time, recompute Sensitivity component
## score, hold Exposure fixed, recompute vulnerability score and rank.

message("Running LOO — Sensitivity attributes ...")

loo_sens_long <- purrr::map_dfr(stock_list, function(s) {
  row <- vuln_std %>% dplyr::filter(stock_name == s)
  run_loo_sensitivity_one_stock(
    stock_name_i            = s,
    sens_means_std          = sens_means_std,
    baseline_exp_score_num  = row$exposure_score_numeric,
    baseline_sens_rank      = row$sensitivity_rank,
    baseline_vuln_rank      = row$vulnerability_rank,
    baseline_vuln_score_num = row$vulnerability_score_numeric
  )
})

loo_sens_summary <- loo_sens_long %>%
  dplyr::group_by(attribute_omitted) %>%
  dplyr::summarise(
    n_stocks_tested   = dplyr::n(),
    n_rank_changed    = sum(rank_changed, na.rm = TRUE),
    n_rank_lower      = sum(rank_change_direction == "Lower",    na.rm = TRUE),
    n_rank_higher     = sum(rank_change_direction == "Higher",   na.rm = TRUE),
    prop_rank_changed = round(mean(rank_changed, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(n_rank_changed), attribute_omitted)

readr::write_csv(loo_sens_long,    f_loo_sens_long)
readr::write_csv(loo_sens_summary, f_loo_sens_sum)

message("✓ Sensitivity LOO complete — ",
        nrow(loo_sens_long), " rows written.")

##------------------------------------------------------------------------------
## Step 5B - Leave-one-out: Exposure factors
##
## Omit one Exposure factor at a time, recompute Exposure component score,
## hold Sensitivity fixed, recompute vulnerability score and rank.

message("Running LOO — Exposure factors ...")

loo_exp_long <- purrr::map_dfr(stock_list, function(s) {
  row <- vuln_std %>% dplyr::filter(stock_name == s)
  run_loo_exposure_one_stock(
    stock_name_i            = s,
    exp_means_std           = exp_means_std,
    baseline_sens_score_num = row$sensitivity_score_numeric,
    baseline_exp_rank       = row$exposure_rank,
    baseline_vuln_rank      = row$vulnerability_rank,
    baseline_vuln_score_num = row$vulnerability_score_numeric
  )
})

loo_exp_summary <- loo_exp_long %>%
  dplyr::group_by(factor_omitted) %>%
  dplyr::summarise(
    n_stocks_tested   = dplyr::n(),
    n_rank_changed    = sum(rank_changed, na.rm = TRUE),
    n_rank_lower      = sum(rank_change_direction == "Lower",    na.rm = TRUE),
    n_rank_higher     = sum(rank_change_direction == "Higher",   na.rm = TRUE),
    prop_rank_changed = round(mean(rank_changed, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(n_rank_changed), factor_omitted)

readr::write_csv(loo_exp_long,    f_loo_exp_long)
readr::write_csv(loo_exp_summary, f_loo_exp_sum)

message("✓ Exposure LOO complete — ",
        nrow(loo_exp_long), " rows written.")

################################################################################
##------------------------------------------------------------------------------
## Step 6 - Save figure-ready descriptive score tables

plot_exp_scores <- exp_means_std %>%
  dplyr::left_join(
    vuln_std %>% dplyr::select(stock_name, exposure_rank),
    by = "stock_name"
  ) %>%
  dplyr::rename(component_rank_baseline = exposure_rank)

plot_sens_scores <- sens_means_std %>%
  dplyr::left_join(
    vuln_std %>% dplyr::select(stock_name, sensitivity_rank),
    by = "stock_name"
  ) %>%
  dplyr::rename(component_rank_baseline = sensitivity_rank)

readr::write_csv(plot_exp_scores,  f_plot_exp)
readr::write_csv(plot_sens_scores, f_plot_sens)

##------------------------------------------------------------------------------
## Done

message("Script completed successfully.")
message("Outputs written to: ", out_dir)
