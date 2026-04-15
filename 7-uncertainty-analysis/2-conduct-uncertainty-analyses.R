################################################################################
##------------------------------------------------------------------------------
## Caribbean CVA
## Script 1. Run bootstrap uncertainty analysis and leave-one-out influence
## analysis
##
## Purpose:
## - Read compiled Caribbean CVA inputs
## - Standardize and validate analysis tables
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
## - Sensitivity comes from reviewer tally data for qualitative attributes
## - Exposure comes from finalized quantitative exposure factor scores
## - Bootstrap uncertainty resamples Sensitivity tallies only
## - Exposure remains fixed at the baseline finalized score during bootstrap
## - Leave-one-out analyses are deterministic, not bootstrap-based
## - This script should stop if reproduced baseline scores do not match the
##   finalized Caribbean CVA outputs

##------------------------------------------------------------------------------
## User setup

rm(list = ls()); gc()

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(tibble)
library(ggplot2)

##------------------------------------------------------------------------------
## User setup - analysis settings

n_boot                <- 10000
save_iteration_table  <- FALSE
borderline_threshold  <- 0.25
bootstrap_seed        <- 99


##------------------------------------------------------------------------------
## Directories

proj_dir <- "."
in_dir   <- file.path(proj_dir, "outputs", "final-scores-compiled")
out_dir  <- file.path(proj_dir, "outputs", "analyses")

input_dir        <- file.path(out_dir, "1-inputs")
intermediate_dir <- file.path(out_dir, "2-intermediate")
final_dir        <- file.path(out_dir, "3-final-tables")

dir.create(out_dir,          recursive = TRUE, showWarnings = FALSE)
dir.create(input_dir,        recursive = TRUE, showWarnings = FALSE)
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_dir,        recursive = TRUE, showWarnings = FALSE)




##------------------------------------------------------------------------------
## User setup - input files
##
## Expected inputs:
## 1) Sensitivity tally table: one row per stock x reviewer x attribute
## 2) Baseline Sensitivity attribute means
## 3) Baseline Exposure factor means
## 4) Baseline component scores
## 5) Baseline overall vulnerability scores
## 6) Optional stock metadata table

f_sens_tallies <- file.path(input_dir,
                            "table_sensitivity_tallies_long.csv")

f_sens_means   <- file.path(in_dir,
                            "final-attribute-scores",
                            "attribute_means_uscar.csv")

f_exp_means    <- file.path(in_dir,
                            "final-attribute-scores",
                            "exposure_factor_means_uscar.csv")

f_components   <- file.path(in_dir,
                            "final-attribute-scores",
                            "component_scores_uscar.csv")

f_vuln         <- file.path(in_dir,
                            "final-attribute-scores",
                            "overall_vulnerability_scores_uscar.csv")

f_metadata     <- file.path(input_dir,
                            "table_stock_metadata.csv")


##------------------------------------------------------------------------------
## Step 1 - Read compiled score inputs

sens_tallies_raw <- read_csv(f_sens_tallies, show_col_types = FALSE)
sens_means_raw   <- read_csv(f_sens_means,   show_col_types = FALSE)
exp_means_raw    <- read_csv(f_exp_means,    show_col_types = FALSE)
components_raw   <- read_csv(f_components,   show_col_types = FALSE)
vuln_raw         <- read_csv(f_vuln,         show_col_types = FALSE)

metadata_raw <- if(file.exists(f_metadata)) {
  read_csv(f_metadata, show_col_types = FALSE)
} else {
  tibble()
}

##------------------------------------------------------------------------------
## Step 1A - Check required columns

check_required_columns(
  dat = sens_tallies_raw,
  required_cols = c("stock_name", "reviewer_id", "attribute_name",
                    "attribute_type", "tally_low", "tally_moderate",
                    "tally_high", "tally_very_high"),
  object_name = "sens_tallies_raw"
)

check_required_columns(
  dat = sens_means_raw,
  required_cols = c("stock_name", "attribute_name", "attribute_type",
                    "mean_score"),
  object_name = "sens_means_raw"
)

check_required_columns(
  dat = exp_means_raw,
  required_cols = c("stock_name", "exposure_factor", "mean_score"),
  object_name = "exp_means_raw"
)

check_required_columns(
  dat = vuln_raw,
  required_cols = c("stock_name",
                    "exposure_score_numeric", "exposure_rank",
                    "sensitivity_score_numeric", "sensitivity_rank",
                    "vulnerability_score_numeric", "vulnerability_rank"),
  object_name = "vuln_raw"
)






##------------------------------------------------------------------------------
## User setup - output files
##
## Intermediate standardized inputs and checks

f_std_sens_tallies <- file.path(intermediate_dir,
                                "01_inputs_standardized_sensitivity_tallies.csv")

f_std_sens_means   <- file.path(intermediate_dir,
                                "01_inputs_standardized_sensitivity_means.csv")

f_std_exp_means    <- file.path(intermediate_dir,
                                "01_inputs_standardized_exposure_means.csv")

f_precheck_log     <- file.path(intermediate_dir,
                                "02_preanalysis_check_log.csv")

f_base_comp_repro  <- file.path(intermediate_dir,
                                "03_baseline_reproduced_component_scores.csv")

f_base_vuln_repro  <- file.path(intermediate_dir,
                                "03_baseline_reproduced_vulnerability_scores.csv")

## Additional Step 2 diagnostic outputs

f_precheck_stock_summary <- file.path(intermediate_dir,
                                      "02A_preanalysis_stock_summary.csv")

f_precheck_tally_row_issues <- file.path(intermediate_dir,
                                         "02B_preanalysis_tally_row_issues.csv")

f_precheck_tally_attribute_summary <- file.path(intermediate_dir,
                                                "02C_preanalysis_tally_attribute_summary.csv")

f_precheck_score_issues <- file.path(intermediate_dir,
                                     "02D_preanalysis_score_issues.csv")

f_precheck_crosswalk <- file.path(intermediate_dir,
                                  "02E_preanalysis_crosswalk_checks.csv")

## Additional Step 3 diagnostic outputs

f_base_compare <- file.path(intermediate_dir,
                            "03A_baseline_compare_official_vs_reproduced.csv")

f_base_attr_counts <- file.path(intermediate_dir,
                                "03B_baseline_component_input_counts.csv")

## Final tables

f_boot_stock <- file.path(final_dir,
                          "table_bootstrap_uncertainty_stock.csv")

f_boot_sens  <- file.path(final_dir,
                          "table_bootstrap_uncertainty_sensitivity_component.csv")

f_boot_iter  <- file.path(intermediate_dir,
                          "table_bootstrap_iterations_long.csv")

f_loo_sens_long <- file.path(final_dir,
                             "table_leave_one_out_sensitivity_long.csv")

f_loo_exp_long  <- file.path(final_dir,
                             "table_leave_one_out_exposure_long.csv")

f_loo_sens_sum  <- file.path(final_dir,
                             "table_leave_one_out_sensitivity_summary.csv")

f_loo_exp_sum   <- file.path(final_dir,
                             "table_leave_one_out_exposure_summary.csv")

f_plot_exp <- file.path(final_dir,
                        "table_exposure_factor_scores_for_plot.csv")

f_plot_sens <- file.path(final_dir,
                         "table_sensitivity_attribute_scores_for_plot.csv")

##------------------------------------------------------------------------------
## Helper functions
##
## These helper functions support:
## - name standardization
## - FCVA logic model scoring
## - rank assignment and comparison
## - one-attribute bootstrap resampling
## - one-stock bootstrap uncertainty analysis
## - one-stock leave-one-out analysis for Sensitivity and Exposure

##------------------------------------------------------------------------------
## Helper functions
##
## Keep only small, reusable helper functions here.
## The main stock-by-stock analyses will be written sequentially below in
## Step 4 and Step 5 rather than wrapped in larger functions.

standardize_stock_names <- function(x) {
  ## Purpose:
  ## - Apply light-touch stock-name standardization
  ## - Preserve intended labels as much as possible
  ##
  ## Workflow:
  ## - Coerce to character
  ## - Normalize punctuation
  ## - Replace underscores with spaces
  ## - Trim and squish whitespace
  ##
  ## Notes:
  ## - Add project-specific recodes only if needed for known mismatches
  
  x_std <- x %>%
    as.character() %>%
    str_replace_all("[\u2018\u2019]", "'") %>%
    str_replace_all("[\u201C\u201D]", "\"") %>%
    str_replace_all("[\u2013\u2014]", "-") %>%
    str_replace_all("_", " ") %>%
    str_squish()
  
  ## Optional project-specific recodes
  ## x_std <- case_when(
  ##   x_std == "Old name" ~ "New name",
  ##   TRUE ~ x_std
  ## )
  
  x_std
}

standardize_attribute_names <- function(x) {
  ## Purpose:
  ## - Standardize Sensitivity attribute names
  
  x_std <- x %>%
    as.character() %>%
    str_replace_all("[\u2018\u2019]", "'") %>%
    str_replace_all("[\u201C\u201D]", "\"") %>%
    str_replace_all("[\u2013\u2014]", "-") %>%
    str_replace_all("_", " ") %>%
    str_squish()
  
  ## Optional project-specific recodes
  x_std <- case_when(
    x_std == "Stock Size/Status" ~ "Stock Size Status",
    TRUE ~ x_std
  )
  
  x_std
}

standardize_factor_names <- function(x) {
  ## Purpose:
  ## - Standardize Exposure factor names
  
  x_std <- x %>%
    as.character() %>%
    str_replace_all("[\u2018\u2019]", "'") %>%
    str_replace_all("[\u201C\u201D]", "\"") %>%
    str_replace_all("[\u2013\u2014]", "-") %>%
    str_replace_all("_", " ") %>%
    str_squish()
  
  x_std
}

component_numeric_to_rank <- function(score_numeric) {
  ## Purpose:
  ## - Convert component numeric score to NOAA FCVA rank label
  
  case_when(
    score_numeric == 1 ~ "Low",
    score_numeric == 2 ~ "Moderate",
    score_numeric == 3 ~ "High",
    score_numeric == 4 ~ "Very High",
    TRUE               ~ NA_character_
  )
}

rank_to_index <- function(rank_chr) {
  ## Purpose:
  ## - Convert ordered rank labels to a simple numeric index so that rank
  ##   changes can be compared directionally
  
  case_when(
    rank_chr == "Low"       ~ 1,
    rank_chr == "Moderate"  ~ 2,
    rank_chr == "High"      ~ 3,
    rank_chr == "Very High" ~ 4,
    TRUE                    ~ NA_real_
  )
}

compare_rank_direction <- function(baseline_rank, new_rank) {
  ## Purpose:
  ## - Compare a recalculated rank to the baseline rank
  ##
  ## Returns:
  ## - "No Change"
  ## - "Lower"
  ## - "Higher"
  
  base_i <- rank_to_index(baseline_rank)
  new_i  <- rank_to_index(new_rank)
  
  case_when(
    is.na(base_i) | is.na(new_i) ~ NA_character_,
    new_i == base_i              ~ "No Change",
    new_i < base_i               ~ "Lower",
    new_i > base_i               ~ "Higher",
    TRUE                         ~ NA_character_
  )
}

fcva_logic_model <- function(mean_scores) {
  ## Table-based rule used in prior NOAA FCVAs:
  ## - Very High = at least 3 means >= 3.5
  ## - High      = at least 2 means >= 3.0
  ## - Moderate  = at least 2 means >= 2.5
  ## - Low       = all other cases
  ##
  ## Returns:
  ## - component_rank
  ## - component_score_numeric
  ##
  ## Notes:
  ## - Input should be a numeric vector of mean scores for one component
  ## - Missing values are ignored
  ## - If all values are missing, return NA
  
  if(length(mean_scores) == 0 || all(is.na(mean_scores))) {
    return(
      tibble(
        component_rank = NA_character_,
        component_score_numeric = NA_real_
      )
    )
  }
  
  n_ge_35 <- sum(mean_scores >= 3.5, na.rm = TRUE)
  n_ge_30 <- sum(mean_scores >= 3.0, na.rm = TRUE)
  n_ge_25 <- sum(mean_scores >= 2.5, na.rm = TRUE)
  
  component_rank <- case_when(
    n_ge_35 >= 3 ~ "Very High",
    n_ge_30 >= 2 ~ "High",
    n_ge_25 >= 2 ~ "Moderate",
    TRUE         ~ "Low"
  )
  
  component_score_numeric <- case_when(
    component_rank == "Low"       ~ 1,
    component_rank == "Moderate"  ~ 2,
    component_rank == "High"      ~ 3,
    component_rank == "Very High" ~ 4,
    TRUE                          ~ NA_real_
  )
  
  tibble(component_rank, component_score_numeric)
}

assign_vulnerability_rank <- function(vulnerability_score_numeric) {
  ## Standard NOAA FCVA vulnerability bins
  ##
  ## Product of Exposure numeric score x Sensitivity numeric score:
  ## - 1 to 3   = Low
  ## - 4 to 6   = Moderate
  ## - 8 to 9   = High
  ## - 12 to 16 = Very High
  ##
  ## Notes:
  ## - Values such as 7, 10, or 11 are not valid under the standard NOAA FCVA
  ##   multiplication product and therefore return NA
  
  case_when(
    is.na(vulnerability_score_numeric)                          ~ NA_character_,
    vulnerability_score_numeric <= 3                            ~ "Low",
    vulnerability_score_numeric >= 4 &
      vulnerability_score_numeric <= 6                          ~ "Moderate",
    vulnerability_score_numeric >= 8 &
      vulnerability_score_numeric <= 9                          ~ "High",
    vulnerability_score_numeric >= 12 &
      vulnerability_score_numeric <= 16                         ~ "Very High",
    TRUE                                                        ~ NA_character_
  )
}

check_required_columns <- function(dat, required_cols, object_name) {
  ## Purpose:
  ## - Stop early if a required table is missing expected columns
  
  missing_cols <- setdiff(required_cols, names(dat))
  
  if(length(missing_cols) > 0) {
    stop(
      paste0(
        "Missing required columns in ", object_name, ": ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
}

calc_weighted_mean_from_tallies <- function(tally_low,
                                            tally_moderate,
                                            tally_high,
                                            tally_very_high) {
  ## Purpose:
  ## - Calculate a weighted mean score from four tally bins
  ##
  ## Scoring:
  ## - Low       = 1
  ## - Moderate  = 2
  ## - High      = 3
  ## - Very High = 4
  
  total_tallies <- tally_low + tally_moderate + tally_high + tally_very_high
  
  if(is.na(total_tallies) || total_tallies == 0) {
    return(NA_real_)
  }
  
  ((tally_low * 1) +
      (tally_moderate * 2) +
      (tally_high * 3) +
      (tally_very_high * 4)) / total_tallies
}

bootstrap_one_attribute <- function(dat_att) {
  ## Purpose:
  ## - Resample one stock x one Sensitivity attribute from pooled reviewer
  ##   tallies and return one bootstrap weighted mean score
  ##
  ## Workflow:
  ## - Determine observed reviewers for this stock x attribute
  ## - Sum pooled tallies across reviewers for bins 1:4
  ## - Create draw pile of 1, 2, 3, 4 values
  ## - Sample with replacement using observed reviewers x 5 tallies
  ## - Return the mean sampled score
  ##
  ## Notes:
  ## - Pre-analysis checks should already have verified 5 tallies per reviewer
  
  check_required_columns(
    dat = dat_att,
    required_cols = c("reviewer_id", "tally_low", "tally_moderate",
                      "tally_high", "tally_very_high"),
    object_name = "dat_att"
  )
  
  n_reviewers <- n_distinct(dat_att$reviewer_id)
  draw_size   <- n_reviewers * 5
  
  n1 <- sum(dat_att$tally_low,       na.rm = TRUE)
  n2 <- sum(dat_att$tally_moderate,  na.rm = TRUE)
  n3 <- sum(dat_att$tally_high,      na.rm = TRUE)
  n4 <- sum(dat_att$tally_very_high, na.rm = TRUE)
  
  total_pooled_tallies <- n1 + n2 + n3 + n4
  
  if(total_pooled_tallies == 0) {
    return(NA_real_)
  }
  
  if(total_pooled_tallies != draw_size) {
    stop(
      paste0(
        "Bootstrap tally mismatch detected for one stock x attribute. ",
        "Expected pooled tallies = ", draw_size,
        ", observed pooled tallies = ", total_pooled_tallies, "."
      )
    )
  }
  
  draw_pile <- c(rep(1, n1),
                 rep(2, n2),
                 rep(3, n3),
                 rep(4, n4))
  
  samp <- sample(x = draw_pile,
                 size = draw_size,
                 replace = TRUE)
  
  mean(samp, na.rm = TRUE)
}



##------------------------------------------------------------------------------
## Step 1B - Standardize names and score columns
##
## Workflow:
## - Harmonize stock, attribute, and factor names
## - Coerce score columns to numeric
## - Add tally totals
## - Arrange consistently for reproducibility

sens_tallies_std <- sens_tallies_raw %>%
  mutate(
    stock_name     = standardize_stock_names(stock_name),
    attribute_name = standardize_attribute_names(attribute_name),
    tally_low       = as.numeric(tally_low),
    tally_moderate  = as.numeric(tally_moderate),
    tally_high      = as.numeric(tally_high),
    tally_very_high = as.numeric(tally_very_high),
    n_tallies       = tally_low + tally_moderate + tally_high + tally_very_high
  ) %>%
  arrange(stock_name, reviewer_id, attribute_name)

sens_means_std <- sens_means_raw %>%
  mutate(
    stock_name     = standardize_stock_names(stock_name),
    attribute_name = standardize_attribute_names(attribute_name),
    mean_score     = as.numeric(mean_score)
  ) %>%
  arrange(stock_name, attribute_name)

exp_means_std <- exp_means_raw %>%
  mutate(
    stock_name       = standardize_stock_names(stock_name),
    exposure_factor  = standardize_factor_names(exposure_factor),
    mean_score       = as.numeric(mean_score)
  ) %>%
  arrange(stock_name, exposure_factor)

components_std <- components_raw %>%
  mutate(stock_name = standardize_stock_names(stock_name))

vuln_std <- vuln_raw %>%
  mutate(stock_name = standardize_stock_names(stock_name))

##------------------------------------------------------------------------------
## Step 1C - Write standardized inputs used by this script

write_csv(sens_tallies_std, f_std_sens_tallies)
write_csv(sens_means_std,   f_std_sens_means)
write_csv(exp_means_std,    f_std_exp_means)


##------------------------------------------------------------------------------
## Step 2 - Run pre-analysis checks
##
## Purpose:
## - Run sequential, transparent QA checks before any bootstrap or leave-one-out
##   analysis begins
## - Write diagnostic tables that make it easy to identify which stock,
##   reviewer, attribute, or factor caused a problem
##
## Workflow:
## Step 2A. Row-level checks on the Sensitivity tally table
## Step 2B. Stock-level summary diagnostics
## Step 2C. Stock x attribute pooled tally diagnostics
## Step 2D. Mean-score diagnostics for Sensitivity and Exposure
## Step 2E. Cross-table matching checks
## Step 2F. Build one summary pre-analysis log and stop on fatal failures
##
## Notes:
## - The row-level tally table should have one row per stock x reviewer x
##   attribute
## - Each reviewer x attribute row should sum to exactly 5 tallies
## - Pooled tallies for one stock x attribute should equal:
##   observed reviewers for that stock x 5
## - These checks are designed to be informative, not just pass/fail

##------------------------------------------------------------------------------
## Step 2A - Row-level checks on the Sensitivity tally table
##
## Workflow:
## - Check for duplicate stock x reviewer x attribute rows
## - Check for missing reviewer IDs or attribute names
## - Check for missing tally values
## - Check for negative tally values
## - Check whether each row sums to exactly 5 tallies

sens_tally_row_checks <- sens_tallies_std %>%
  mutate(
    missing_stock_name =
      is.na(stock_name) | stock_name == "",
    missing_reviewer_id =
      is.na(reviewer_id) | reviewer_id == "",
    missing_attribute_name =
      is.na(attribute_name) | attribute_name == "",
    missing_tally_value =
      if_any(c(tally_low, tally_moderate, tally_high, tally_very_high), is.na),
    negative_tally_value =
      if_any(c(tally_low, tally_moderate, tally_high, tally_very_high), ~ .x < 0),
    row_tally_sum =
      tally_low + tally_moderate + tally_high + tally_very_high,
    row_tally_sum_ok =
      row_tally_sum == 5
  )

sens_tally_duplicates <- sens_tallies_std %>%
  count(stock_name, reviewer_id, attribute_name, name = "n_rows") %>%
  mutate(duplicate_row = n_rows > 1)

sens_tally_row_issues <- sens_tally_row_checks %>%
  left_join(
    sens_tally_duplicates,
    by = c("stock_name", "reviewer_id", "attribute_name")
  ) %>%
  mutate(
    duplicate_row = ifelse(is.na(duplicate_row), FALSE, duplicate_row),
    any_row_issue =
      missing_stock_name |
      missing_reviewer_id |
      missing_attribute_name |
      missing_tally_value |
      negative_tally_value |
      !row_tally_sum_ok |
      duplicate_row
  ) %>%
  filter(any_row_issue) %>%
  arrange(stock_name, attribute_name, reviewer_id)

write_csv(sens_tally_row_issues, f_precheck_tally_row_issues)

##------------------------------------------------------------------------------
## Step 2B - Stock-level summary diagnostics
##
## Workflow:
## - Summarize reviewer counts and row counts by stock
## - Check whether every attribute within a stock has the same reviewer count
## - Summarize row-level issues by stock

stock_reviewer_summary <- sens_tallies_std %>%
  group_by(stock_name) %>%
  summarise(
    n_rows_tally_table = n(),
    n_reviewers_observed = n_distinct(reviewer_id),
    n_attributes_observed = n_distinct(attribute_name),
    .groups = "drop"
  )

stock_attribute_reviewer_counts <- sens_tallies_std %>%
  group_by(stock_name, attribute_name) %>%
  summarise(
    n_reviewers_for_attribute = n_distinct(reviewer_id),
    .groups = "drop"
  )

stock_reviewer_consistency <- stock_attribute_reviewer_counts %>%
  group_by(stock_name) %>%
  summarise(
    min_reviewers_across_attributes = min(n_reviewers_for_attribute, na.rm = TRUE),
    max_reviewers_across_attributes = max(n_reviewers_for_attribute, na.rm = TRUE),
    reviewer_count_consistent =
      min_reviewers_across_attributes == max_reviewers_across_attributes,
    .groups = "drop"
  )

stock_row_issue_summary <- sens_tally_row_issues %>%
  group_by(stock_name) %>%
  summarise(
    n_problem_rows = n(),
    n_problem_attributes = n_distinct(attribute_name),
    n_problem_reviewers = n_distinct(reviewer_id),
    .groups = "drop"
  )

precheck_stock_summary <- stock_reviewer_summary %>%
  left_join(stock_reviewer_consistency, by = "stock_name") %>%
  left_join(stock_row_issue_summary, by = "stock_name") %>%
  mutate(
    n_problem_rows = coalesce(n_problem_rows, 0L),
    n_problem_attributes = coalesce(n_problem_attributes, 0L),
    n_problem_reviewers = coalesce(n_problem_reviewers, 0L)
  ) %>%
  arrange(stock_name)

write_csv(precheck_stock_summary, f_precheck_stock_summary)

##------------------------------------------------------------------------------
## Step 2C - Stock x attribute pooled tally diagnostics
##
## Workflow:
## - Pool reviewer tallies within each stock x attribute
## - Calculate expected pooled tally total = reviewers observed for that
##   stock x attribute x 5
## - Compare pooled tally total to expectation
## - Calculate pooled weighted mean score for reference
##
## Notes:
## - This is the key diagnostic table for the bootstrap analysis because the
##   bootstrap resamples from these pooled tallies

precheck_tally_attribute_summary <- sens_tallies_std %>%
  group_by(stock_name, attribute_name) %>%
  summarise(
    n_reviewers_observed = n_distinct(reviewer_id),
    tally_low = sum(tally_low, na.rm = TRUE),
    tally_moderate = sum(tally_moderate, na.rm = TRUE),
    tally_high = sum(tally_high, na.rm = TRUE),
    tally_very_high = sum(tally_very_high, na.rm = TRUE),
    pooled_tally_sum =
      tally_low + tally_moderate + tally_high + tally_very_high,
    expected_pooled_tally_sum =
      n_reviewers_observed * 5,
    pooled_tally_sum_ok =
      pooled_tally_sum == expected_pooled_tally_sum,
    pooled_mean_score =
      calc_weighted_mean_from_tallies(
        tally_low       = tally_low,
        tally_moderate  = tally_moderate,
        tally_high      = tally_high,
        tally_very_high = tally_very_high
      ),
    .groups = "drop"
  ) %>%
  arrange(stock_name, attribute_name)

write_csv(precheck_tally_attribute_summary, f_precheck_tally_attribute_summary)

##------------------------------------------------------------------------------
## Step 2D - Mean-score diagnostics for Sensitivity and Exposure
##
## Workflow:
## - Check for duplicate stock x attribute rows in Sensitivity means
## - Check for duplicate stock x factor rows in Exposure means
## - Check for impossible mean scores outside 1 to 4
## - Write all problem rows to one combined score-issues table

sens_mean_duplicates <- sens_means_std %>%
  count(stock_name, attribute_name, name = "n_rows") %>%
  filter(n_rows > 1) %>%
  mutate(issue_type = "duplicate_sensitivity_mean_row")

sens_mean_range_issues <- sens_means_std %>%
  filter(is.na(mean_score) | mean_score < 1 | mean_score > 4) %>%
  mutate(
    n_rows = 1L,
    issue_type = "sensitivity_mean_score_outside_1_to_4"
  ) %>%
  select(stock_name, attribute_name, n_rows, issue_type, mean_score)

exp_mean_duplicates <- exp_means_std %>%
  count(stock_name, exposure_factor, name = "n_rows") %>%
  filter(n_rows > 1) %>%
  mutate(issue_type = "duplicate_exposure_mean_row")

exp_mean_range_issues <- exp_means_std %>%
  filter(is.na(mean_score) | mean_score < 1 | mean_score > 4) %>%
  mutate(
    n_rows = 1L,
    issue_type = "exposure_mean_score_outside_1_to_4"
  ) %>%
  select(stock_name, exposure_factor, n_rows, issue_type, mean_score)

precheck_score_issues <- bind_rows(
  sens_mean_duplicates %>%
    rename(item_name = attribute_name) %>%
    mutate(table_name = "sens_means_std"),
  sens_mean_range_issues %>%
    rename(item_name = attribute_name) %>%
    mutate(table_name = "sens_means_std"),
  exp_mean_duplicates %>%
    rename(item_name = exposure_factor) %>%
    mutate(table_name = "exp_means_std"),
  exp_mean_range_issues %>%
    rename(item_name = exposure_factor) %>%
    mutate(table_name = "exp_means_std")
) %>%
  arrange(table_name, stock_name, item_name)

write_csv(precheck_score_issues, f_precheck_score_issues)

##------------------------------------------------------------------------------
## Step 2E - Cross-table matching checks
##
## Workflow:
## - Compare stock coverage across:
##   - Sensitivity tally table
##   - Sensitivity mean table
##   - Exposure mean table
##   - Baseline vulnerability table
## - Check stock x attribute matching between tally and Sensitivity mean tables
##
## Notes:
## - This catches cases where a stock exists in one table but not another
## - It also catches stock x attribute combinations present in one table but not
##   in the other

stocks_tallies <- sens_tallies_std %>%
  distinct(stock_name) %>%
  mutate(in_tallies = TRUE)

stocks_sens_means <- sens_means_std %>%
  distinct(stock_name) %>%
  mutate(in_sens_means = TRUE)

stocks_exp_means <- exp_means_std %>%
  distinct(stock_name) %>%
  mutate(in_exp_means = TRUE)

stocks_vuln <- vuln_std %>%
  distinct(stock_name) %>%
  mutate(in_vuln = TRUE)

stock_crosswalk <- stocks_tallies %>%
  full_join(stocks_sens_means, by = "stock_name") %>%
  full_join(stocks_exp_means,  by = "stock_name") %>%
  full_join(stocks_vuln,       by = "stock_name") %>%
  mutate(
    across(starts_with("in_"), ~ coalesce(.x, FALSE)),
    in_all_tables =
      in_tallies & in_sens_means & in_exp_means & in_vuln
  ) %>%
  arrange(stock_name)

sens_pairs_tallies <- sens_tallies_std %>%
  distinct(stock_name, attribute_name) %>%
  mutate(in_tally_pairs = TRUE)

sens_pairs_means <- sens_means_std %>%
  distinct(stock_name, attribute_name) %>%
  mutate(in_mean_pairs = TRUE)

sens_pair_crosswalk <- sens_pairs_tallies %>%
  full_join(sens_pairs_means, by = c("stock_name", "attribute_name")) %>%
  mutate(
    in_tally_pairs = coalesce(in_tally_pairs, FALSE),
    in_mean_pairs  = coalesce(in_mean_pairs, FALSE),
    pair_in_both   = in_tally_pairs & in_mean_pairs
  ) %>%
  arrange(stock_name, attribute_name)

precheck_crosswalk_checks <- bind_rows(
  stock_crosswalk %>%
    transmute(
      check_level = "stock",
      stock_name = stock_name,
      item_name = NA_character_,
      issue_type = ifelse(in_all_tables,
                          "ok",
                          "stock_missing_from_one_or_more_tables"),
      in_tallies = in_tallies,
      in_sens_means = in_sens_means,
      in_exp_means = in_exp_means,
      in_vuln = in_vuln
    ),
  sens_pair_crosswalk %>%
    filter(!pair_in_both) %>%
    transmute(
      check_level = "stock_attribute",
      stock_name = stock_name,
      item_name = attribute_name,
      issue_type = "stock_attribute_pair_missing_from_tally_or_mean_table",
      in_tallies = in_tally_pairs,
      in_sens_means = in_mean_pairs,
      in_exp_means = NA,
      in_vuln = NA
    )
) %>%
  arrange(check_level, stock_name, item_name)

write_csv(precheck_crosswalk_checks, f_precheck_crosswalk)

##------------------------------------------------------------------------------
## Step 2F - Build one summary pre-analysis log
##
## Workflow:
## - Collapse the detailed diagnostics above into one high-level log
## - Mark clearly which checks are PASS / FAIL / WARN
## - Stop the script if any fatal checks fail
##
## Fatal failures:
## - duplicate tally rows
## - missing / negative / non-5 tally rows
## - pooled stock x attribute tally mismatch
## - duplicate mean rows
## - mean scores outside 1 to 4
## - stocks missing from one or more core tables
##
## Warning-level checks:
## - inconsistent reviewer counts across attributes within a stock
## - stock x attribute pair mismatch between tally and mean tables

n_dup_tally_rows <- sens_tally_duplicates %>%
  filter(duplicate_row) %>%
  nrow()

n_bad_row_tally_sum <- sens_tally_row_checks %>%
  filter(!row_tally_sum_ok) %>%
  nrow()

n_missing_tally_values <- sens_tally_row_checks %>%
  filter(missing_tally_value) %>%
  nrow()

n_negative_tally_values <- sens_tally_row_checks %>%
  filter(negative_tally_value) %>%
  nrow()

n_missing_reviewer_or_attribute <- sens_tally_row_checks %>%
  filter(missing_reviewer_id | missing_attribute_name | missing_stock_name) %>%
  nrow()

n_pooled_tally_mismatch <- precheck_tally_attribute_summary %>%
  filter(!pooled_tally_sum_ok) %>%
  nrow()

n_inconsistent_reviewer_counts <- precheck_stock_summary %>%
  filter(!reviewer_count_consistent) %>%
  nrow()

n_sens_mean_duplicates <- sens_mean_duplicates %>% nrow()
n_exp_mean_duplicates  <- exp_mean_duplicates %>% nrow()

n_sens_mean_range_issues <- sens_mean_range_issues %>% nrow()
n_exp_mean_range_issues  <- exp_mean_range_issues %>% nrow()

n_stock_missing_core_table <- stock_crosswalk %>%
  filter(!in_all_tables) %>%
  nrow()

n_stock_attribute_pair_mismatch <- sens_pair_crosswalk %>%
  filter(!pair_in_both) %>%
  nrow()

precheck_log <- bind_rows(
  tibble(
    check_name = "duplicate_stock_reviewer_attribute_rows",
    status = ifelse(n_dup_tally_rows == 0, "PASS", "FAIL"),
    level = "row",
    n_problem = n_dup_tally_rows,
    detail = "Sensitivity tally table should have one row per stock x reviewer x attribute"
  ),
  tibble(
    check_name = "row_tally_sum_not_equal_5",
    status = ifelse(n_bad_row_tally_sum == 0, "PASS", "FAIL"),
    level = "row",
    n_problem = n_bad_row_tally_sum,
    detail = "Each reviewer x attribute row should sum to exactly 5 tallies"
  ),
  tibble(
    check_name = "missing_tally_values",
    status = ifelse(n_missing_tally_values == 0, "PASS", "FAIL"),
    level = "row",
    n_problem = n_missing_tally_values,
    detail = "Sensitivity tally bins should not be missing"
  ),
  tibble(
    check_name = "negative_tally_values",
    status = ifelse(n_negative_tally_values == 0, "PASS", "FAIL"),
    level = "row",
    n_problem = n_negative_tally_values,
    detail = "Sensitivity tally bins should not be negative"
  ),
  tibble(
    check_name = "missing_stock_reviewer_or_attribute_labels",
    status = ifelse(n_missing_reviewer_or_attribute == 0, "PASS", "FAIL"),
    level = "row",
    n_problem = n_missing_reviewer_or_attribute,
    detail = "Stock name, reviewer ID, and attribute name should all be present"
  ),
  tibble(
    check_name = "pooled_stock_attribute_tally_mismatch",
    status = ifelse(n_pooled_tally_mismatch == 0, "PASS", "FAIL"),
    level = "stock_attribute",
    n_problem = n_pooled_tally_mismatch,
    detail = "Pooled stock x attribute tallies should equal observed reviewers x 5"
  ),
  tibble(
    check_name = "inconsistent_reviewer_counts_across_attributes_within_stock",
    status = ifelse(n_inconsistent_reviewer_counts == 0, "PASS", "WARN"),
    level = "stock",
    n_problem = n_inconsistent_reviewer_counts,
    detail = "Reviewer count differs among attributes within a stock"
  ),
  tibble(
    check_name = "duplicate_sensitivity_mean_rows",
    status = ifelse(n_sens_mean_duplicates == 0, "PASS", "FAIL"),
    level = "score_table",
    n_problem = n_sens_mean_duplicates,
    detail = "Sensitivity mean table should have one row per stock x attribute"
  ),
  tibble(
    check_name = "duplicate_exposure_mean_rows",
    status = ifelse(n_exp_mean_duplicates == 0, "PASS", "FAIL"),
    level = "score_table",
    n_problem = n_exp_mean_duplicates,
    detail = "Exposure mean table should have one row per stock x factor"
  ),
  tibble(
    check_name = "sensitivity_mean_scores_outside_1_to_4",
    status = ifelse(n_sens_mean_range_issues == 0, "PASS", "FAIL"),
    level = "score_table",
    n_problem = n_sens_mean_range_issues,
    detail = "Sensitivity mean scores should be between 1 and 4"
  ),
  tibble(
    check_name = "exposure_mean_scores_outside_1_to_4",
    status = ifelse(n_exp_mean_range_issues == 0, "PASS", "FAIL"),
    level = "score_table",
    n_problem = n_exp_mean_range_issues,
    detail = "Exposure mean scores should be between 1 and 4"
  ),
  tibble(
    check_name = "stocks_missing_from_one_or_more_core_tables",
    status = ifelse(n_stock_missing_core_table == 0, "PASS", "FAIL"),
    level = "stock",
    n_problem = n_stock_missing_core_table,
    detail = "Each stock should appear in the tally, Sensitivity mean, Exposure mean, and baseline vulnerability tables"
  ),
  tibble(
    check_name = "stock_attribute_pair_mismatch_between_tally_and_mean_tables",
    status = ifelse(n_stock_attribute_pair_mismatch == 0, "PASS", "WARN"),
    level = "stock_attribute",
    n_problem = n_stock_attribute_pair_mismatch,
    detail = "Each stock x attribute pair in the tally table should also exist in the Sensitivity mean table, and vice versa"
  )
) %>%
  arrange(desc(status), level, check_name)

write_csv(precheck_log, f_precheck_log)

##------------------------------------------------------------------------------
## Stop on fatal failures
##
## Notes:
## - WARN-level checks are written to file but do not stop the script
## - FAIL-level checks stop the script before analysis begins

if(any(precheck_log$status == "FAIL")) {
  stop(
    paste0(
      "Pre-analysis checks failed. Review:\n",
      " - ", basename(f_precheck_log), "\n",
      " - ", basename(f_precheck_tally_row_issues), "\n",
      " - ", basename(f_precheck_tally_attribute_summary), "\n",
      " - ", basename(f_precheck_score_issues), "\n",
      " - ", basename(f_precheck_crosswalk)
    )
  )
}

##------------------------------------------------------------------------------
## Step 3 - Reproduce baseline component and vulnerability scores
##
## Purpose:
## - Recalculate baseline Sensitivity and Exposure component scores directly
##   from the compiled mean-score tables
## - Recalculate baseline overall vulnerability score and rank
## - Compare reproduced values against the official finalized baseline outputs
## - Write a clear stock-by-stock comparison table for diagnostics
##
## Workflow:
## Step 3A. Summarize the number of component inputs per stock
## Step 3B. Reproduce baseline Sensitivity component scores
## Step 3C. Reproduce baseline Exposure component scores
## Step 3D. Reproduce baseline Vulnerability scores
## Step 3E. Compare official vs reproduced values stock by stock
## Step 3F. Stop the script if any reproduced values do not match
##
## Notes:
## - Sensitivity reproduction uses sens_means_std
## - Exposure reproduction uses exp_means_std
## - Official baseline values come from vuln_std
## - This step should be exact; any mismatch indicates a logic or input problem

##------------------------------------------------------------------------------
## Step 3A - Summarize the number of component inputs per stock
##
## Workflow:
## - Count the number of Sensitivity attributes used per stock
## - Count the number of Exposure factors used per stock
## - Save as a simple diagnostic table for traceability

base_attr_counts <- sens_means_std %>%
  group_by(stock_name) %>%
  summarise(
    n_sensitivity_attributes = n(),
    .groups = "drop"
  ) %>%
  full_join(
    exp_means_std %>%
      group_by(stock_name) %>%
      summarise(
        n_exposure_factors = n(),
        .groups = "drop"
      ),
    by = "stock_name"
  ) %>%
  arrange(stock_name)

write_csv(base_attr_counts, f_base_attr_counts)

##------------------------------------------------------------------------------
## Step 3B - Reproduce baseline Sensitivity component scores
##
## Workflow:
## - Group Sensitivity mean scores by stock
## - Apply the NOAA FCVA logic model to the vector of attribute means
## - Write one row per stock with reproduced Sensitivity score and rank

base_sens_repro <- sens_means_std %>%
  group_by(stock_name) %>%
  summarise(
    sensitivity_logic = list(fcva_logic_model(mean_scores = mean_score)),
    sensitivity_mean_of_means = mean(mean_score, na.rm = TRUE),
    n_sensitivity_attributes = n(),
    .groups = "drop"
  ) %>%
  unnest(cols = c(sensitivity_logic)) %>%
  rename(
    sensitivity_rank_repro = component_rank,
    sensitivity_score_numeric_repro = component_score_numeric
  ) %>%
  arrange(stock_name)

##------------------------------------------------------------------------------
## Step 3C - Reproduce baseline Exposure component scores
##
## Workflow:
## - Group Exposure mean scores by stock
## - Apply the NOAA FCVA logic model to the vector of factor means
## - Write one row per stock with reproduced Exposure score and rank

base_exp_repro <- exp_means_std %>%
  group_by(stock_name) %>%
  summarise(
    exposure_logic = list(fcva_logic_model(mean_scores = mean_score)),
    exposure_mean_of_means = mean(mean_score, na.rm = TRUE),
    n_exposure_factors = n(),
    .groups = "drop"
  ) %>%
  unnest(cols = c(exposure_logic)) %>%
  rename(
    exposure_rank_repro = component_rank,
    exposure_score_numeric_repro = component_score_numeric
  ) %>%
  arrange(stock_name)

##------------------------------------------------------------------------------
## Step 3D - Reproduce baseline Vulnerability scores
##
## Workflow:
## - Join reproduced Sensitivity and Exposure component scores by stock
## - Multiply numeric component scores to get reproduced vulnerability score
## - Convert reproduced vulnerability numeric score to reproduced rank

baseline_reproduced <- base_sens_repro %>%
  full_join(base_exp_repro, by = "stock_name") %>%
  mutate(
    vulnerability_score_numeric_repro =
      sensitivity_score_numeric_repro * exposure_score_numeric_repro,
    vulnerability_rank_repro =
      assign_vulnerability_rank(vulnerability_score_numeric_repro)
  ) %>%
  arrange(stock_name)

##------------------------------------------------------------------------------
## Step 3E - Compare official vs reproduced baseline values
##
## Workflow:
## - Join the official finalized baseline outputs to the reproduced outputs
## - Create explicit match flags for each component and the overall product
## - Create a single all_match flag for quick checking
##
## Notes:
## - This table is the main diagnostic product for baseline reproduction
## - It should make it immediately obvious which stock failed and why

baseline_compare <- vuln_std %>%
  select(
    stock_name,
    exposure_score_numeric,
    exposure_rank,
    sensitivity_score_numeric,
    sensitivity_rank,
    vulnerability_score_numeric,
    vulnerability_rank
  ) %>%
  full_join(baseline_reproduced, by = "stock_name") %>%
  mutate(
    stock_in_official = !is.na(exposure_score_numeric) |
      !is.na(sensitivity_score_numeric) |
      !is.na(vulnerability_score_numeric),
    stock_in_reproduced = !is.na(exposure_score_numeric_repro) |
      !is.na(sensitivity_score_numeric_repro) |
      !is.na(vulnerability_score_numeric_repro),
    
    exp_num_match =
      exposure_score_numeric == exposure_score_numeric_repro,
    exp_rank_match =
      exposure_rank == exposure_rank_repro,
    
    sens_num_match =
      sensitivity_score_numeric == sensitivity_score_numeric_repro,
    sens_rank_match =
      sensitivity_rank == sensitivity_rank_repro,
    
    vuln_num_match =
      vulnerability_score_numeric == vulnerability_score_numeric_repro,
    vuln_rank_match =
      vulnerability_rank == vulnerability_rank_repro,
    
    exp_num_diff =
      exposure_score_numeric_repro - exposure_score_numeric,
    sens_num_diff =
      sensitivity_score_numeric_repro - sensitivity_score_numeric,
    vuln_num_diff =
      vulnerability_score_numeric_repro - vulnerability_score_numeric,
    
    all_match =
      stock_in_official &
      stock_in_reproduced &
      exp_num_match &
      exp_rank_match &
      sens_num_match &
      sens_rank_match &
      vuln_num_match &
      vuln_rank_match
  ) %>%
  arrange(stock_name)

##------------------------------------------------------------------------------
## Step 3F - Write reproduced and comparison outputs
##
## Workflow:
## - Write the reproduced component table
## - Write the reproduced vulnerability table
## - Write the full comparison table

baseline_component_repro <- baseline_reproduced %>%
  select(
    stock_name,
    n_sensitivity_attributes,
    sensitivity_mean_of_means,
    sensitivity_rank_repro,
    sensitivity_score_numeric_repro,
    n_exposure_factors,
    exposure_mean_of_means,
    exposure_rank_repro,
    exposure_score_numeric_repro
  ) %>%
  arrange(stock_name)

write_csv(baseline_component_repro, f_base_comp_repro)
write_csv(baseline_reproduced,      f_base_vuln_repro)
write_csv(baseline_compare,         f_base_compare)

##------------------------------------------------------------------------------
## Step 3G - Stop on any mismatch
##
## Notes:
## - Baseline reproduction should be exact
## - Any mismatch means the script should stop before bootstrap or leave-one-out

n_missing_official   <- baseline_compare %>% filter(!stock_in_official) %>% nrow()
n_missing_reproduced <- baseline_compare %>% filter(!stock_in_reproduced) %>% nrow()
n_nonmatching_stocks <- baseline_compare %>% filter(!all_match) %>% nrow()

if(n_missing_official > 0 || n_missing_reproduced > 0 || n_nonmatching_stocks > 0) {
  stop(
    paste0(
      "Baseline reproduction failed. Review:\n",
      " - ", basename(f_base_comp_repro), "\n",
      " - ", basename(f_base_vuln_repro), "\n",
      " - ", basename(f_base_compare)
    )
  )
}

##------------------------------------------------------------------------------
## Step 4 - Bootstrap uncertainty analysis
##
## Workflow:
## - Loop through stocks sequentially
## - For each stock:
##   - pull stock-level tally data
##   - identify observed reviewers and attributes
##   - pull official baseline Exposure and Vulnerability values
##   - run 10,000 bootstrap iterations
##   - summarize bootstrap probabilities for Sensitivity and Vulnerability
##   - append one-row summaries to output tables
##   - optionally append iteration-level output
##
## Output tables:
## - table_bootstrap_uncertainty_stock.csv
## - table_bootstrap_uncertainty_sensitivity_component.csv
## - optional table_bootstrap_iterations_long.csv

set.seed(bootstrap_seed)

stock_list <- sort(unique(vuln_std$stock_name))

boot_stock_summary <- tibble()
boot_sens_summary  <- tibble()
boot_iter_long     <- tibble()

for(stock_i in stock_list) {
  
  message("Bootstrap uncertainty: ", stock_i)
  
  ##--------------------------------------------------------------------------
  ## Pull stock-specific Sensitivity tally data
  
  dat_stock <- sens_tallies_std %>%
    filter(stock_name == stock_i)
  
  if(nrow(dat_stock) == 0) {
    stop(paste0("No Sensitivity tally data found for stock: ", stock_i))
  }
  
  ##--------------------------------------------------------------------------
  ## Derive stock-specific metadata
  
  n_reviewers_observed <- n_distinct(dat_stock$reviewer_id)
  
  attribute_list <- dat_stock %>%
    distinct(attribute_name) %>%
    arrange(attribute_name) %>%
    pull(attribute_name)
  
  ##--------------------------------------------------------------------------
  ## Pull official baseline scores for this stock
  ##
  ## Use the finalized baseline outputs already carried in vuln_std. These are
  ## the official baseline values that the bootstrap certainty metrics should
  ## reference.
  
  base_row <- vuln_std %>%
    filter(stock_name == stock_i)
  
  if(nrow(base_row) != 1) {
    stop(paste0("Expected exactly one baseline row in vuln_std for stock: ", stock_i))
  }
  
  baseline_exp_score_num  <- base_row$exposure_score_numeric
  baseline_exp_rank       <- base_row$exposure_rank
  baseline_sens_score_num <- base_row$sensitivity_score_numeric
  baseline_sens_rank      <- base_row$sensitivity_rank
  baseline_vuln_score_num <- base_row$vulnerability_score_numeric
  baseline_vuln_rank      <- base_row$vulnerability_rank
  
  ##--------------------------------------------------------------------------
  ## Preallocate vectors for bootstrap results
  ##
  ## This is faster and cleaner than growing vectors inside the loop.
  
  sens_score_num_boot <- rep(NA_real_, n_boot)
  sens_rank_boot      <- rep(NA_character_, n_boot)
  vuln_score_num_boot <- rep(NA_real_, n_boot)
  vuln_rank_boot      <- rep(NA_character_, n_boot)
  
  ##--------------------------------------------------------------------------
  ## Bootstrap loop
  ##
  ## For each iteration:
  ## - resample every Sensitivity attribute
  ## - apply logic model to bootstrapped attribute means
  ## - hold Exposure fixed at the finalized baseline numeric score
  ## - calculate bootstrapped vulnerability score and rank
  
  for(i in seq_len(n_boot)) {
    
    boot_mean_scores_i <- rep(NA_real_, length(attribute_list))
    
    for(j in seq_along(attribute_list)) {
      
      att_j <- attribute_list[j]
      
      dat_att_j <- dat_stock %>%
        filter(attribute_name == att_j)
      
      boot_mean_scores_i[j] <- bootstrap_one_attribute(dat_att = dat_att_j)
    }
    
    sens_logic_i <- fcva_logic_model(mean_scores = boot_mean_scores_i)
    
    sens_score_num_boot[i] <- sens_logic_i$component_score_numeric
    sens_rank_boot[i]      <- sens_logic_i$component_rank
    
    vuln_score_num_boot[i] <- sens_score_num_boot[i] * baseline_exp_score_num
    vuln_rank_boot[i]      <- assign_vulnerability_rank(vuln_score_num_boot[i])
  }
  
  ##--------------------------------------------------------------------------
  ## Summarize bootstrap Vulnerability rank probabilities
  
  p_low       <- mean(vuln_rank_boot == "Low",       na.rm = TRUE)
  p_moderate  <- mean(vuln_rank_boot == "Moderate",  na.rm = TRUE)
  p_high      <- mean(vuln_rank_boot == "High",      na.rm = TRUE)
  p_very_high <- mean(vuln_rank_boot == "Very High", na.rm = TRUE)
  
  ##--------------------------------------------------------------------------
  ## Summarize bootstrap Sensitivity rank probabilities
  
  p_sens_low       <- mean(sens_rank_boot == "Low",       na.rm = TRUE)
  p_sens_moderate  <- mean(sens_rank_boot == "Moderate",  na.rm = TRUE)
  p_sens_high      <- mean(sens_rank_boot == "High",      na.rm = TRUE)
  p_sens_very_high <- mean(sens_rank_boot == "Very High", na.rm = TRUE)
  
  ##--------------------------------------------------------------------------
  ## Certainty metrics
  ##
  ## Certainty is defined as the proportion of bootstrap iterations that match
  ## the finalized baseline rank.
  
  certainty_baseline_rank <- case_when(
    baseline_vuln_rank == "Low"       ~ p_low,
    baseline_vuln_rank == "Moderate"  ~ p_moderate,
    baseline_vuln_rank == "High"      ~ p_high,
    baseline_vuln_rank == "Very High" ~ p_very_high,
    TRUE                              ~ NA_real_
  )
  
  certainty_baseline_sensitivity_rank <- case_when(
    baseline_sens_rank == "Low"       ~ p_sens_low,
    baseline_sens_rank == "Moderate"  ~ p_sens_moderate,
    baseline_sens_rank == "High"      ~ p_sens_high,
    baseline_sens_rank == "Very High" ~ p_sens_very_high,
    TRUE                              ~ NA_real_
  )
  
  ##--------------------------------------------------------------------------
  ## Modal and secondary rank metrics
  ##
  ## These help identify borderline stocks where the bootstrap distribution is
  ## split across adjacent ranks.
  
  rank_prob_tbl <- tibble(
    rank = c("Low", "Moderate", "High", "Very High"),
    prob = c(p_low, p_moderate, p_high, p_very_high)
  ) %>%
    arrange(desc(prob), rank)
  
  modal_boot_rank <- rank_prob_tbl$rank[1]
  secondary_rank  <- rank_prob_tbl$rank[2]
  secondary_rank_probability <- rank_prob_tbl$prob[2]
  
  baseline_rank_index  <- rank_to_index(baseline_vuln_rank)
  secondary_rank_index <- rank_to_index(secondary_rank)
  
  is_borderline <- !is.na(baseline_rank_index) &&
    !is.na(secondary_rank_index) &&
    abs(secondary_rank_index - baseline_rank_index) == 1 &&
    secondary_rank_probability >= borderline_threshold
  
  ##--------------------------------------------------------------------------
  ## Append stock-level Vulnerability bootstrap summary
  
  boot_stock_summary <- bind_rows(
    boot_stock_summary,
    tibble(
      stock_name = stock_i,
      n_boot = n_boot,
      n_reviewers_observed = n_reviewers_observed,
      baseline_exposure_rank = baseline_exp_rank,
      baseline_exposure_score_numeric = baseline_exp_score_num,
      baseline_sensitivity_rank = baseline_sens_rank,
      baseline_sensitivity_score_numeric = baseline_sens_score_num,
      baseline_vulnerability_rank = baseline_vuln_rank,
      baseline_vulnerability_score_numeric = baseline_vuln_score_num,
      p_low = p_low,
      p_moderate = p_moderate,
      p_high = p_high,
      p_very_high = p_very_high,
      modal_boot_rank = modal_boot_rank,
      certainty_baseline_rank = certainty_baseline_rank,
      secondary_rank = secondary_rank,
      secondary_rank_probability = secondary_rank_probability,
      is_borderline = is_borderline,
      borderline_definition =
        paste0("Adjacent non-baseline rank probability >= ",
               borderline_threshold)
    )
  )
  
  ##--------------------------------------------------------------------------
  ## Append stock-level Sensitivity bootstrap summary
  
  sens_rank_prob_tbl <- tibble(
    rank = c("Low", "Moderate", "High", "Very High"),
    prob = c(p_sens_low, p_sens_moderate, p_sens_high, p_sens_very_high)
  ) %>%
    arrange(desc(prob), rank)
  
  boot_sens_summary <- bind_rows(
    boot_sens_summary,
    tibble(
      stock_name = stock_i,
      n_boot = n_boot,
      baseline_sensitivity_rank = baseline_sens_rank,
      baseline_sensitivity_score_numeric = baseline_sens_score_num,
      p_sens_low = p_sens_low,
      p_sens_moderate = p_sens_moderate,
      p_sens_high = p_sens_high,
      p_sens_very_high = p_sens_very_high,
      modal_sensitivity_rank = sens_rank_prob_tbl$rank[1],
      certainty_baseline_sensitivity_rank =
        certainty_baseline_sensitivity_rank
    )
  )
  
  ##--------------------------------------------------------------------------
  ## Optionally append iteration-level output
  ##
  ## This can be large, so only save if requested.
  
  if(save_iteration_table) {
    boot_iter_long <- bind_rows(
      boot_iter_long,
      tibble(
        stock_name = stock_i,
        boot_iter = seq_len(n_boot),
        sensitivity_score_numeric_boot = sens_score_num_boot,
        sensitivity_rank_boot = sens_rank_boot,
        vulnerability_score_numeric_boot = vuln_score_num_boot,
        vulnerability_rank_boot = vuln_rank_boot
      )
    )
  }
}

##------------------------------------------------------------------------------
## Write bootstrap outputs to file

write_csv(boot_stock_summary, f_boot_stock)
write_csv(boot_sens_summary,  f_boot_sens)

if(save_iteration_table) {
  write_csv(boot_iter_long, f_boot_iter)
}


##------------------------------------------------------------------------------
## Step 5A - Sensitivity leave-one-out influence analysis
##
## Workflow:
## - Loop through stocks sequentially
## - For each stock:
##   - pull Sensitivity mean scores
##   - pull official baseline scores
##   - omit one Sensitivity attribute at a time
##   - recalculate Sensitivity component score
##   - hold Exposure fixed
##   - recalculate overall vulnerability score and rank
##   - append long-format results
##
## Output tables:
## - table_leave_one_out_sensitivity_long.csv
## - table_leave_one_out_sensitivity_summary.csv

loo_sens_long <- tibble()

for(stock_i in stock_list) {
  
  message("Sensitivity leave-one-out: ", stock_i)
  
  ##--------------------------------------------------------------------------
  ## Pull stock-specific Sensitivity mean-score data
  
  dat_stock <- sens_means_std %>%
    filter(stock_name == stock_i)
  
  if(nrow(dat_stock) == 0) {
    stop(paste0("No Sensitivity mean-score data found for stock: ", stock_i))
  }
  
  attribute_list <- dat_stock %>%
    distinct(attribute_name) %>%
    arrange(attribute_name) %>%
    pull(attribute_name)
  
  ##--------------------------------------------------------------------------
  ## Pull official baseline values for this stock
  
  base_row <- vuln_std %>%
    filter(stock_name == stock_i)
  
  if(nrow(base_row) != 1) {
    stop(paste0("Expected exactly one baseline row in vuln_std for stock: ", stock_i))
  }
  
  baseline_exp_score_num  <- base_row$exposure_score_numeric
  baseline_sens_score_num <- base_row$sensitivity_score_numeric
  baseline_sens_rank      <- base_row$sensitivity_rank
  baseline_vuln_score_num <- base_row$vulnerability_score_numeric
  baseline_vuln_rank      <- base_row$vulnerability_rank
  
  ##--------------------------------------------------------------------------
  ## Attribute omission loop
  
  for(att_i in attribute_list) {
    
    dat_omit_i <- dat_stock %>%
      filter(attribute_name != att_i)
    
    loo_sens_logic_i <- fcva_logic_model(mean_scores = dat_omit_i$mean_score)
    
    loo_sens_rank_i      <- loo_sens_logic_i$component_rank
    loo_sens_score_num_i <- loo_sens_logic_i$component_score_numeric
    
    loo_vuln_score_num_i <- loo_sens_score_num_i * baseline_exp_score_num
    loo_vuln_rank_i      <- assign_vulnerability_rank(loo_vuln_score_num_i)
    
    loo_sens_long <- bind_rows(
      loo_sens_long,
      tibble(
        stock_name = stock_i,
        attribute_name = att_i,
        baseline_sensitivity_rank = baseline_sens_rank,
        baseline_sensitivity_score_numeric = baseline_sens_score_num,
        baseline_exposure_score_numeric = baseline_exp_score_num,
        baseline_vulnerability_rank = baseline_vuln_rank,
        baseline_vulnerability_score_numeric = baseline_vuln_score_num,
        loo_sensitivity_rank = loo_sens_rank_i,
        loo_sensitivity_score_numeric = loo_sens_score_num_i,
        loo_vulnerability_rank = loo_vuln_rank_i,
        loo_vulnerability_score_numeric = loo_vuln_score_num_i,
        rank_changed = loo_vuln_rank_i != baseline_vuln_rank,
        rank_change_direction =
          compare_rank_direction(
            baseline_rank = baseline_vuln_rank,
            new_rank      = loo_vuln_rank_i
          ),
        delta_sensitivity_score_numeric =
          loo_sens_score_num_i - baseline_sens_score_num,
        delta_vulnerability_score_numeric =
          loo_vuln_score_num_i - baseline_vuln_score_num
      )
    )
  }
}

##------------------------------------------------------------------------------
## Summarize Sensitivity leave-one-out results for plotting and interpretation

loo_sens_summary <- loo_sens_long %>%
  group_by(attribute_name) %>%
  summarise(
    n_rank_changed = sum(rank_changed, na.rm = TRUE),
    n_rank_lower   = sum(rank_change_direction == "Lower", na.rm = TRUE),
    n_rank_higher  = sum(rank_change_direction == "Higher", na.rm = TRUE),
    prop_rank_changed = mean(rank_changed, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_rank_changed), attribute_name)

write_csv(loo_sens_long, f_loo_sens_long)
write_csv(loo_sens_summary, f_loo_sens_sum)



##------------------------------------------------------------------------------
## Step 5B - Exposure leave-one-out influence analysis
##
## Workflow:
## - Loop through stocks sequentially
## - For each stock:
##   - pull Exposure mean scores
##   - pull official baseline scores
##   - omit one Exposure factor at a time
##   - recalculate Exposure component score
##   - hold Sensitivity fixed
##   - recalculate overall vulnerability score and rank
##   - append long-format results
##
## Output tables:
## - table_leave_one_out_exposure_long.csv
## - table_leave_one_out_exposure_summary.csv

loo_exp_long <- tibble()

for(stock_i in stock_list) {
  
  message("Exposure leave-one-out: ", stock_i)
  
  ##--------------------------------------------------------------------------
  ## Pull stock-specific Exposure mean-score data
  
  dat_stock <- exp_means_std %>%
    filter(stock_name == stock_i)
  
  if(nrow(dat_stock) == 0) {
    stop(paste0("No Exposure mean-score data found for stock: ", stock_i))
  }
  
  factor_list <- dat_stock %>%
    distinct(exposure_factor) %>%
    arrange(exposure_factor) %>%
    pull(exposure_factor)
  
  ##--------------------------------------------------------------------------
  ## Pull official baseline values for this stock
  
  base_row <- vuln_std %>%
    filter(stock_name == stock_i)
  
  if(nrow(base_row) != 1) {
    stop(paste0("Expected exactly one baseline row in vuln_std for stock: ", stock_i))
  }
  
  baseline_exp_score_num  <- base_row$exposure_score_numeric
  baseline_exp_rank       <- base_row$exposure_rank
  baseline_sens_score_num <- base_row$sensitivity_score_numeric
  baseline_vuln_score_num <- base_row$vulnerability_score_numeric
  baseline_vuln_rank      <- base_row$vulnerability_rank
  
  ##--------------------------------------------------------------------------
  ## Factor omission loop
  
  for(fac_i in factor_list) {
    
    dat_omit_i <- dat_stock %>%
      filter(exposure_factor != fac_i)
    
    loo_exp_logic_i <- fcva_logic_model(mean_scores = dat_omit_i$mean_score)
    
    loo_exp_rank_i      <- loo_exp_logic_i$component_rank
    loo_exp_score_num_i <- loo_exp_logic_i$component_score_numeric
    
    loo_vuln_score_num_i <- baseline_sens_score_num * loo_exp_score_num_i
    loo_vuln_rank_i      <- assign_vulnerability_rank(loo_vuln_score_num_i)
    
    loo_exp_long <- bind_rows(
      loo_exp_long,
      tibble(
        stock_name = stock_i,
        exposure_factor = fac_i,
        baseline_exposure_rank = baseline_exp_rank,
        baseline_exposure_score_numeric = baseline_exp_score_num,
        baseline_sensitivity_score_numeric = baseline_sens_score_num,
        baseline_vulnerability_rank = baseline_vuln_rank,
        baseline_vulnerability_score_numeric = baseline_vuln_score_num,
        loo_exposure_rank = loo_exp_rank_i,
        loo_exposure_score_numeric = loo_exp_score_num_i,
        loo_vulnerability_rank = loo_vuln_rank_i,
        loo_vulnerability_score_numeric = loo_vuln_score_num_i,
        rank_changed = loo_vuln_rank_i != baseline_vuln_rank,
        rank_change_direction =
          compare_rank_direction(
            baseline_rank = baseline_vuln_rank,
            new_rank      = loo_vuln_rank_i
          ),
        delta_exposure_score_numeric =
          loo_exp_score_num_i - baseline_exp_score_num,
        delta_vulnerability_score_numeric =
          loo_vuln_score_num_i - baseline_vuln_score_num
      )
    )
  }
}

##------------------------------------------------------------------------------
## Summarize Exposure leave-one-out results for plotting and interpretation

loo_exp_summary <- loo_exp_long %>%
  group_by(exposure_factor) %>%
  summarise(
    n_rank_changed = sum(rank_changed, na.rm = TRUE),
    n_rank_lower   = sum(rank_change_direction == "Lower", na.rm = TRUE),
    n_rank_higher  = sum(rank_change_direction == "Higher", na.rm = TRUE),
    prop_rank_changed = mean(rank_changed, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_rank_changed), exposure_factor)

write_csv(loo_exp_long, f_loo_exp_long)
write_csv(loo_exp_summary, f_loo_exp_sum)
##------------------------------------------------------------------------------
## Step 6 - Save figure-ready descriptive score tables
##
## These tables are used by Script 2 for the boxplots of score distributions
## across stocks.

plot_exp_scores <- exp_means_std %>%
  left_join(
    vuln_std %>%
      select(stock_name, exposure_rank),
    by = "stock_name"
  ) %>%
  rename(component_rank_baseline = exposure_rank)

plot_sens_scores <- sens_means_std %>%
  left_join(
    vuln_std %>%
      select(stock_name, sensitivity_rank),
    by = "stock_name"
  ) %>%
  rename(component_rank_baseline = sensitivity_rank)

write_csv(plot_exp_scores,  f_plot_exp)
write_csv(plot_sens_scores, f_plot_sens)

##------------------------------------------------------------------------------
## Script completed

message("Script 1 completed successfully.")