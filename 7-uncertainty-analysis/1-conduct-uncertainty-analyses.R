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

standardize_stock_names <- function(x) {
  ## Purpose:
  ## - Apply light-touch name standardization without changing intended labels
  ## - Preserve capitalization and wording as much as possible
  ##
  ## Workflow:
  ## - Coerce to character
  ## - Replace curly quotes/apostrophes with plain text
  ## - Replace en/em dashes with hyphens
  ## - Replace underscores with spaces
  ## - Trim leading/trailing whitespace
  ## - Squish repeated internal whitespace
  ##
  ## Notes:
  ## - Add project-specific recodes below if needed
  ## - Keep this conservative unless a known stock-name mismatch exists
  
  x_std <- x %>%
    as.character() %>%
    str_replace_all("[\u2018\u2019]", "'") %>%
    str_replace_all("[\u201C\u201D]", "\"") %>%
    str_replace_all("[\u2013\u2014]", "-") %>%
    str_replace_all("_", " ") %>%
    str_squish()
  
  ## Optional project-specific recodes
  ## Example:
  ## x_std <- case_when(
  ##   x_std == "Old stock name" ~ "New stock name",
  ##   TRUE ~ x_std
  ## )
  
  x_std
}

standardize_attribute_names <- function(x) {
  ## Purpose:
  ## - Standardize Sensitivity attribute names
  ##
  ## Workflow:
  ## - Coerce to character
  ## - Normalize punctuation and spacing
  ## - Apply any known project-specific recodes
  ##
  ## Notes:
  ## - This is the place to resolve minor label inconsistencies such as
  ##   "Stock Size/Status" vs "Stock Size Status", if needed
  
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
  ##
  ## Workflow:
  ## - Coerce to character
  ## - Normalize punctuation and spacing
  ## - Apply any known project-specific recodes
  
  x_std <- x %>%
    as.character() %>%
    str_replace_all("[\u2018\u2019]", "'") %>%
    str_replace_all("[\u201C\u201D]", "\"") %>%
    str_replace_all("[\u2013\u2014]", "-") %>%
    str_replace_all("_", " ") %>%
    str_squish()
  
  ## Optional project-specific recodes
  x_std <- case_when(
    TRUE ~ x_std
  )
  
  x_std
}

component_numeric_to_rank <- function(score_numeric) {
  ## Purpose:
  ## - Convert NOAA FCVA component numeric scores to rank labels
  ##
  ## Mapping:
  ## - 1 = Low
  ## - 2 = Moderate
  ## - 3 = High
  ## - 4 = Very High
  
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
  ## - Convert ordered rank labels to a comparable numeric scale
  ##
  ## Mapping:
  ## - Low       = 1
  ## - Moderate  = 2
  ## - High      = 3
  ## - Very High = 4
  
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
  ## - Describe the direction of rank change relative to the baseline
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
  ## - If all values are missing, NA is returned
  
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
  ## - Values like 7, 10, or 11 are not legal under the standard 1:4 x 1:4
  ##   NOAA FCVA component-score product and therefore return NA
  
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
  ## - Stop the script early if a required input table is missing columns
  
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
  ## - Determine observed number of reviewers
  ## - Sum pooled tallies across reviewers for bins 1:4
  ## - Create a draw pile of scored tallies
  ## - Sample with replacement using observed reviewers x 5 tallies
  ## - Return the mean sampled score
  ##
  ## Notes:
  ## - This function assumes one row per reviewer for a single attribute
  ## - Pre-analysis checks should already have verified that each reviewer
  ##   contributed 5 tallies
  
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

bootstrap_one_stock <- function(stock_name_i,
                                sens_tallies_std,
                                baseline_exp_score_num,
                                n_boot = 10000,
                                borderline_threshold = 0.25) {
  ## Purpose:
  ## - Run the bootstrap uncertainty analysis for a single stock
  ##
  ## Workflow:
  ## - Filter one stock from the standardized Sensitivity tally table
  ## - Calculate baseline attribute means directly from pooled tallies
  ## - Apply the FCVA logic model to determine baseline Sensitivity score
  ## - Hold baseline Exposure numeric score fixed
  ## - For each bootstrap iteration:
  ##   - Resample each attribute
  ##   - Recalculate bootstrap Sensitivity score
  ##   - Recalculate overall vulnerability score and rank
  ## - Summarize proportions across vulnerability bins
  ## - Flag borderline cases when an adjacent non-baseline rank occurs with
  ##   probability >= borderline_threshold
  ##
  ## Returns:
  ## - list(stock_summary, sens_summary, iter_long)
  
  dat_stock <- sens_tallies_std %>%
    filter(stock_name == stock_name_i)
  
  if(nrow(dat_stock) == 0) {
    stop(paste0("No Sensitivity tally data found for stock: ", stock_name_i))
  }
  
  ##--------------------------------------------------------------------------
  ## Derive stock-specific metadata from the tally table
  
  n_reviewers_observed <- n_distinct(dat_stock$reviewer_id)
  
  attribute_list <- dat_stock %>%
    distinct(attribute_name) %>%
    arrange(attribute_name) %>%
    pull(attribute_name)
  
  ##--------------------------------------------------------------------------
  ## Calculate baseline attribute means from pooled tallies
  ##
  ## This provides a tally-based baseline that can be compared with the
  ## finalized baseline outputs elsewhere in the script.
  
  baseline_attribute_means <- dat_stock %>%
    group_by(attribute_name) %>%
    summarise(
      tally_low       = sum(tally_low, na.rm = TRUE),
      tally_moderate  = sum(tally_moderate, na.rm = TRUE),
      tally_high      = sum(tally_high, na.rm = TRUE),
      tally_very_high = sum(tally_very_high, na.rm = TRUE),
      baseline_mean_score = calc_weighted_mean_from_tallies(
        tally_low       = tally_low,
        tally_moderate  = tally_moderate,
        tally_high      = tally_high,
        tally_very_high = tally_very_high
      ),
      .groups = "drop"
    ) %>%
    arrange(attribute_name)
  
  baseline_sens_logic <- fcva_logic_model(
    mean_scores = baseline_attribute_means$baseline_mean_score
  )
  
  baseline_sens_rank      <- baseline_sens_logic$component_rank
  baseline_sens_score_num <- baseline_sens_logic$component_score_numeric
  
  baseline_exp_rank <- component_numeric_to_rank(baseline_exp_score_num)
  
  baseline_vuln_score_num <- baseline_sens_score_num * baseline_exp_score_num
  baseline_vuln_rank      <- assign_vulnerability_rank(baseline_vuln_score_num)
  
  ##--------------------------------------------------------------------------
  ## Run bootstrap iterations
  ##
  ## Preallocate vectors for speed and reproducibility
  
  sens_score_num_boot <- rep(NA_real_, n_boot)
  sens_rank_boot      <- rep(NA_character_, n_boot)
  vuln_score_num_boot <- rep(NA_real_, n_boot)
  vuln_rank_boot      <- rep(NA_character_, n_boot)
  
  for(i in seq_len(n_boot)) {
    
    boot_mean_scores_i <- map_dbl(
      attribute_list,
      function(att_i) {
        dat_att_i <- dat_stock %>%
          filter(attribute_name == att_i)
        
        bootstrap_one_attribute(dat_att = dat_att_i)
      }
    )
    
    sens_logic_i <- fcva_logic_model(mean_scores = boot_mean_scores_i)
    
    sens_score_num_boot[i] <- sens_logic_i$component_score_numeric
    sens_rank_boot[i]      <- sens_logic_i$component_rank
    
    vuln_score_num_boot[i] <- sens_score_num_boot[i] * baseline_exp_score_num
    vuln_rank_boot[i]      <- assign_vulnerability_rank(vuln_score_num_boot[i])
  }
  
  ##--------------------------------------------------------------------------
  ## Summarize bootstrap vulnerability rank probabilities
  
  vuln_rank_levels <- c("Low", "Moderate", "High", "Very High")
  
  vuln_rank_table <- tibble(vulnerability_rank_boot = vuln_rank_boot) %>%
    mutate(vulnerability_rank_boot = factor(vulnerability_rank_boot,
                                            levels = vuln_rank_levels)) %>%
    count(vulnerability_rank_boot, name = "n_iter") %>%
    complete(vulnerability_rank_boot = factor(vuln_rank_levels,
                                              levels = vuln_rank_levels),
             fill = list(n_iter = 0)) %>%
    mutate(prob = n_iter / n_boot)
  
  p_low       <- vuln_rank_table %>% filter(vulnerability_rank_boot == "Low")       %>% pull(prob)
  p_moderate  <- vuln_rank_table %>% filter(vulnerability_rank_boot == "Moderate")  %>% pull(prob)
  p_high      <- vuln_rank_table %>% filter(vulnerability_rank_boot == "High")      %>% pull(prob)
  p_very_high <- vuln_rank_table %>% filter(vulnerability_rank_boot == "Very High") %>% pull(prob)
  
  ##--------------------------------------------------------------------------
  ## Summarize bootstrap Sensitivity rank probabilities
  
  sens_rank_table <- tibble(sensitivity_rank_boot = sens_rank_boot) %>%
    mutate(sensitivity_rank_boot = factor(sensitivity_rank_boot,
                                          levels = vuln_rank_levels)) %>%
    count(sensitivity_rank_boot, name = "n_iter") %>%
    complete(sensitivity_rank_boot = factor(vuln_rank_levels,
                                            levels = vuln_rank_levels),
             fill = list(n_iter = 0)) %>%
    mutate(prob = n_iter / n_boot)
  
  p_sens_low       <- sens_rank_table %>% filter(sensitivity_rank_boot == "Low")       %>% pull(prob)
  p_sens_moderate  <- sens_rank_table %>% filter(sensitivity_rank_boot == "Moderate")  %>% pull(prob)
  p_sens_high      <- sens_rank_table %>% filter(sensitivity_rank_boot == "High")      %>% pull(prob)
  p_sens_very_high <- sens_rank_table %>% filter(sensitivity_rank_boot == "Very High") %>% pull(prob)
  
  ##--------------------------------------------------------------------------
  ## Derive certainty and borderline metrics
  ##
  ## Certainty is defined here as the proportion of bootstrap iterations that
  ## match the baseline rank for that stock.
  
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
  
  rank_prob_tbl <- tibble(
    rank = vuln_rank_levels,
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
  ## Create stock-level summary table
  
  stock_summary <- tibble(
    stock_name = stock_name_i,
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
  
  ##--------------------------------------------------------------------------
  ## Create Sensitivity-component summary table
  
  sens_summary <- tibble(
    stock_name = stock_name_i,
    n_boot = n_boot,
    baseline_sensitivity_rank = baseline_sens_rank,
    baseline_sensitivity_score_numeric = baseline_sens_score_num,
    p_sens_low = p_sens_low,
    p_sens_moderate = p_sens_moderate,
    p_sens_high = p_sens_high,
    p_sens_very_high = p_sens_very_high,
    modal_sensitivity_rank =
      sens_rank_table %>%
      arrange(desc(prob), sensitivity_rank_boot) %>%
      slice(1) %>%
      pull(sensitivity_rank_boot) %>%
      as.character(),
    certainty_baseline_sensitivity_rank =
      certainty_baseline_sensitivity_rank
  )
  
  ##--------------------------------------------------------------------------
  ## Create optional iteration-level output table
  
  iter_long <- tibble(
    stock_name = stock_name_i,
    boot_iter = seq_len(n_boot),
    sensitivity_score_numeric_boot = sens_score_num_boot,
    sensitivity_rank_boot = sens_rank_boot,
    vulnerability_score_numeric_boot = vuln_score_num_boot,
    vulnerability_rank_boot = vuln_rank_boot
  )
  
  list(
    stock_summary = stock_summary,
    sens_summary  = sens_summary,
    iter_long     = iter_long
  )
}

run_loo_sensitivity_one_stock <- function(stock_name_i,
                                          sens_means_std,
                                          baseline_exp_score_num,
                                          baseline_sens_rank,
                                          baseline_vuln_rank,
                                          baseline_vuln_score_num) {
  ## Purpose:
  ## - Run the deterministic Sensitivity leave-one-out analysis for one stock
  ##
  ## Workflow:
  ## - Filter one stock from the standardized Sensitivity mean-score table
  ## - Recalculate the baseline Sensitivity score from all attributes
  ## - Omit one Sensitivity attribute at a time
  ## - Recompute Sensitivity component score
  ## - Hold Exposure numeric score fixed
  ## - Recompute overall vulnerability score and rank
  ##
  ## Returns:
  ## - One long-format tibble with one row per omitted attribute
  
  dat_stock <- sens_means_std %>%
    filter(stock_name == stock_name_i)
  
  if(nrow(dat_stock) == 0) {
    stop(paste0("No Sensitivity mean-score data found for stock: ", stock_name_i))
  }
  
  attribute_list <- dat_stock %>%
    distinct(attribute_name) %>%
    arrange(attribute_name) %>%
    pull(attribute_name)
  
  baseline_sens_logic <- fcva_logic_model(dat_stock$mean_score)
  
  baseline_sens_score_num <- baseline_sens_logic$component_score_numeric
  
  loo_out <- map_dfr(
    attribute_list,
    function(att_i) {
      
      dat_omit_i <- dat_stock %>%
        filter(attribute_name != att_i)
      
      loo_sens_logic_i <- fcva_logic_model(dat_omit_i$mean_score)
      
      loo_sens_rank_i      <- loo_sens_logic_i$component_rank
      loo_sens_score_num_i <- loo_sens_logic_i$component_score_numeric
      
      loo_vuln_score_num_i <- loo_sens_score_num_i * baseline_exp_score_num
      loo_vuln_rank_i      <- assign_vulnerability_rank(loo_vuln_score_num_i)
      
      tibble(
        stock_name = stock_name_i,
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
    }
  )
  
  loo_out
}

run_loo_exposure_one_stock <- function(stock_name_i,
                                       exp_means_std,
                                       baseline_sens_score_num,
                                       baseline_exp_rank,
                                       baseline_vuln_rank,
                                       baseline_vuln_score_num) {
  ## Purpose:
  ## - Run the deterministic Exposure leave-one-out analysis for one stock
  ##
  ## Workflow:
  ## - Filter one stock from the standardized Exposure mean-score table
  ## - Recalculate the baseline Exposure score from all factors
  ## - Omit one Exposure factor at a time
  ## - Recompute Exposure component score
  ## - Hold Sensitivity numeric score fixed
  ## - Recompute overall vulnerability score and rank
  ##
  ## Returns:
  ## - One long-format tibble with one row per omitted factor
  
  dat_stock <- exp_means_std %>%
    filter(stock_name == stock_name_i)
  
  if(nrow(dat_stock) == 0) {
    stop(paste0("No Exposure mean-score data found for stock: ", stock_name_i))
  }
  
  factor_list <- dat_stock %>%
    distinct(exposure_factor) %>%
    arrange(exposure_factor) %>%
    pull(exposure_factor)
  
  baseline_exp_logic <- fcva_logic_model(dat_stock$mean_score)
  
  baseline_exp_score_num <- baseline_exp_logic$component_score_numeric
  
  loo_out <- map_dfr(
    factor_list,
    function(fac_i) {
      
      dat_omit_i <- dat_stock %>%
        filter(exposure_factor != fac_i)
      
      loo_exp_logic_i <- fcva_logic_model(dat_omit_i$mean_score)
      
      loo_exp_rank_i      <- loo_exp_logic_i$component_rank
      loo_exp_score_num_i <- loo_exp_logic_i$component_score_numeric
      
      loo_vuln_score_num_i <- baseline_sens_score_num * loo_exp_score_num_i
      loo_vuln_rank_i      <- assign_vulnerability_rank(loo_vuln_score_num_i)
      
      tibble(
        stock_name = stock_name_i,
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
    }
  )
  
  loo_out
}









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
## Workflow:
## - Check for duplicated stock x reviewer x attribute rows
## - Check for impossible tally totals
## - Check for impossible mean scores outside 1 to 4
## - Check that all stocks in tally tables exist in baseline vulnerability table
##
## Notes:
## - Expected tally total per reviewer x attribute is 5
## - Total draw size for bootstrap is observed reviewers x 5
## - Caribbean workflow should infer observed reviewer count by stock

precheck_log <- tibble(
  check_name = character(),
  status     = character(),
  detail     = character()
)

## Placeholder duplicated row check
dup_sens_tallies <- sens_tallies_std %>%
  count(stock_name, reviewer_id, attribute_name) %>%
  filter(n > 1)

## Placeholder invalid tally check
bad_tally_total <- sens_tallies_std %>%
  filter(n_tallies != 5)

## Placeholder impossible mean-score check
bad_sens_means <- sens_means_std %>%
  filter(mean_score < 1 | mean_score > 4)

bad_exp_means <- exp_means_std %>%
  filter(mean_score < 1 | mean_score > 4)

## Placeholder stock matching check
stocks_missing_in_vuln <- setdiff(unique(sens_tallies_std$stock_name),
                                  unique(vuln_std$stock_name))

## Assemble precheck log
precheck_log <- bind_rows(
  tibble(check_name = "duplicate_sensitivity_tally_rows",
         status     = ifelse(nrow(dup_sens_tallies) == 0, "PASS", "FAIL"),
         detail     = paste("n =", nrow(dup_sens_tallies))),
  tibble(check_name = "invalid_tally_total_not_equal_5",
         status     = ifelse(nrow(bad_tally_total) == 0, "PASS", "FAIL"),
         detail     = paste("n =", nrow(bad_tally_total))),
  tibble(check_name = "sensitivity_means_outside_1_to_4",
         status     = ifelse(nrow(bad_sens_means) == 0, "PASS", "FAIL"),
         detail     = paste("n =", nrow(bad_sens_means))),
  tibble(check_name = "exposure_means_outside_1_to_4",
         status     = ifelse(nrow(bad_exp_means) == 0, "PASS", "FAIL"),
         detail     = paste("n =", nrow(bad_exp_means))),
  tibble(check_name = "stocks_missing_in_baseline_vulnerability_table",
         status     = ifelse(length(stocks_missing_in_vuln) == 0, "PASS", "FAIL"),
         detail     = paste("n =", length(stocks_missing_in_vuln)))
)

write_csv(precheck_log, f_precheck_log)

if(any(precheck_log$status == "FAIL")) {
  stop("Pre-analysis checks failed. Review 02_preanalysis_check_log.csv")
}

##------------------------------------------------------------------------------
## Step 3 - Reproduce baseline component and vulnerability scores
##
## Workflow:
## - Apply FCVA logic model to baseline Sensitivity mean scores
## - Apply FCVA logic model to baseline Exposure mean scores
## - Multiply numeric component scores to recreate baseline vulnerability score
## - Assign vulnerability rank
##
## Notes:
## - Sensitivity should come from finalized qualitative mean scores
## - Exposure should come from finalized quantitative factor mean scores
## - This reproduced baseline must match the official final Caribbean outputs
##   exactly, or the script should stop

base_sens_repro <- sens_means_std %>%
  group_by(stock_name) %>%
  summarise(
    sens_logic = list(fcva_logic_model(mean_scores = mean_score)),
    .groups = "drop"
  ) %>%
  unnest(cols = c(sens_logic)) %>%
  rename(sensitivity_rank_repro = component_rank,
         sensitivity_score_numeric_repro = component_score_numeric)

base_exp_repro <- exp_means_std %>%
  group_by(stock_name) %>%
  summarise(
    exp_logic = list(fcva_logic_model(mean_scores = mean_score)),
    .groups = "drop"
  ) %>%
  unnest(cols = c(exp_logic)) %>%
  rename(exposure_rank_repro = component_rank,
         exposure_score_numeric_repro = component_score_numeric)

baseline_reproduced <- base_sens_repro %>%
  left_join(base_exp_repro, by = "stock_name") %>%
  mutate(
    vulnerability_score_numeric_repro =
      sensitivity_score_numeric_repro * exposure_score_numeric_repro,
    vulnerability_rank_repro =
      assign_vulnerability_rank(vulnerability_score_numeric_repro)
  )

baseline_compare <- vuln_std %>%
  left_join(baseline_reproduced, by = "stock_name") %>%
  mutate(
    sens_rank_match =
      sensitivity_rank == sensitivity_rank_repro,
    sens_num_match =
      sensitivity_score_numeric == sensitivity_score_numeric_repro,
    exp_rank_match =
      exposure_rank == exposure_rank_repro,
    exp_num_match =
      exposure_score_numeric == exposure_score_numeric_repro,
    vuln_rank_match =
      vulnerability_rank == vulnerability_rank_repro,
    vuln_num_match =
      vulnerability_score_numeric == vulnerability_score_numeric_repro
  )

write_csv(baseline_reproduced, f_base_vuln_repro)

baseline_component_repro <- baseline_reproduced %>%
  select(stock_name,
         exposure_rank_repro, exposure_score_numeric_repro,
         sensitivity_rank_repro, sensitivity_score_numeric_repro)

write_csv(baseline_component_repro, f_base_comp_repro)

if(any(!baseline_compare$sens_rank_match) |
   any(!baseline_compare$sens_num_match)  |
   any(!baseline_compare$exp_rank_match)  |
   any(!baseline_compare$exp_num_match)   |
   any(!baseline_compare$vuln_rank_match) |
   any(!baseline_compare$vuln_num_match)) {
  stop("Baseline reproduction failed. Review reproduced score tables.")
}

##------------------------------------------------------------------------------
## Step 4 - Bootstrap uncertainty analysis
##
## Workflow:
## - For each stock:
##   - Pool reviewer tallies within each Sensitivity attribute
##   - Resample with replacement using observed reviewers x 5 tallies
##   - Recalculate bootstrap attribute mean scores
##   - Apply logic model for bootstrap Sensitivity score
##   - Hold baseline Exposure numeric score fixed
##   - Recalculate vulnerability score and rank
## - Summarize proportions across 10,000 iterations
##
## Output tables:
## - table_bootstrap_uncertainty_stock.csv
## - table_bootstrap_uncertainty_sensitivity_component.csv
## - optional table_bootstrap_iterations_long.csv

set.seed(bootstrap_seed)

stock_list <- sort(unique(vuln_std$stock_name))

boot_results <- map(
  stock_list,
  ~ bootstrap_one_stock(
    stock_name_i           = .x,
    sens_tallies_std       = sens_tallies_std,
    baseline_exp_score_num =
      vuln_std %>%
      filter(stock_name == .x) %>%
      pull(exposure_score_numeric),
    n_boot                 = n_boot,
    borderline_threshold   = borderline_threshold
  )
)

boot_stock_summary <- bind_rows(map(boot_results, "stock_summary"))
boot_sens_summary  <- bind_rows(map(boot_results, "sens_summary"))
boot_iter_long     <- bind_rows(map(boot_results, "iter_long"))

write_csv(boot_stock_summary, f_boot_stock)
write_csv(boot_sens_summary,  f_boot_sens)

if(save_iteration_table) {
  write_csv(boot_iter_long, f_boot_iter)
}

###------------------------------------------------------------------------------
## Step 5A. Sensitivity leave-one-out
##
## Workflow:
## - Omit one Sensitivity attribute at a time
## - Recompute Sensitivity component score
## - Hold baseline Exposure fixed
## - Recompute vulnerability score and rank

loo_sens_results <- map(
  stock_list,
  ~ run_loo_sensitivity_one_stock(
    stock_name_i            = .x,
    sens_means_std          = sens_means_std,
    baseline_exp_score_num  =
      vuln_std %>%
      filter(stock_name == .x) %>%
      pull(exposure_score_numeric),
    baseline_sens_rank      =
      vuln_std %>%
      filter(stock_name == .x) %>%
      pull(sensitivity_rank),
    baseline_vuln_rank      =
      vuln_std %>%
      filter(stock_name == .x) %>%
      pull(vulnerability_rank),
    baseline_vuln_score_num =
      vuln_std %>%
      filter(stock_name == .x) %>%
      pull(vulnerability_score_numeric)
  )
)

loo_sens_long <- bind_rows(loo_sens_results)

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
## Step 5B. Exposure leave-one-out
##
## Workflow:
## - Omit one Exposure factor at a time
## - Recompute Exposure component score
## - Hold baseline Sensitivity fixed
## - Recompute vulnerability score and rank

loo_exp_results <- map(
  stock_list,
  ~ run_loo_exposure_one_stock(
    stock_name_i            = .x,
    exp_means_std           = exp_means_std,
    baseline_sens_score_num =
      vuln_std %>%
      filter(stock_name == .x) %>%
      pull(sensitivity_score_numeric),
    baseline_exp_rank       =
      vuln_std %>%
      filter(stock_name == .x) %>%
      pull(exposure_rank),
    baseline_vuln_rank      =
      vuln_std %>%
      filter(stock_name == .x) %>%
      pull(vulnerability_rank),
    baseline_vuln_score_num =
      vuln_std %>%
      filter(stock_name == .x) %>%
      pull(vulnerability_score_numeric)
  )
)

loo_exp_long <- bind_rows(loo_exp_results)

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