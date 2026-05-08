################################################################################
##------------------------------------------------------------------------------
## Caribbean CVA
## Module 9 — Potential for Distributional Change
## Script 2: Bootstrap DCP uncertainty
##
## Mirrors the draw-pile bootstrap in Module 8 (uncertainty-analyses.R):
##   1. Filter sensitivity_tallies_long.csv to the four DCP attributes.
##   2. For the three inverted movement attributes, swap tally counts
##      (L↔VH, M↔H) before building draw piles. This is equivalent to
##      inverting each individual tally vote (5 - vote) and is cheaper
##      than inverting at the draw step.
##   3. Aggregate swapped tallies across reviewers → draw piles of 20 votes.
##   4. Baseline reproduction gate: apply fcva_logic_model() to draw-pile
##      means and verify they reproduce the Script 1 DCP ranks. Halt if any
##      stock disagrees.
##   5. Run 10,000 bootstrap iterations per stock. Each iteration resamples
##      all four draw piles with replacement and applies fcva_logic_model().
##   6. Summarize rank-proportion distribution; flag borderline stocks
##      (dominant rank < 75% of iterations, matching Module 8 convention).
##
## Inputs:
##   outputs/final-tallies-long/sensitivity_tallies_long.csv
##   outputs/distributional_change_potential_uscar.csv   (Script 1 output)
##
## Outputs:
##   outputs/distributional_change_bootstrap_uscar.csv
##     Columns: stock_name, prop_L, prop_M, prop_H, prop_VH,
##              dominant_rank, dominant_prop, borderline
##   outputs/distributional_change_full_uscar.csv
##     Scripts 1 and 2 joined — input for Script 3 plotting
##
## Dependencies: dplyr, tidyr, readr, stringr

##------------------------------------------------------------------------------
## Setup

rm(list = ls()); gc()

library(dplyr)
library(tidyr)
library(readr)
library(stringr)

##------------------------------------------------------------------------------
## Configuration — must match Script 1 and Module 8

rank_threshold       <- 1       ## FCVA logic model threshold
bootstrap_seed       <- 99      ## matches Module 8
borderline_threshold <- 0.25    ## flag stocks where dominant prop < 0.75
n_boot               <- 10000

## Target attributes — same as Script 1
target_attributes <- c(
  "Adult mobility",
  "Habitat specificity",
  "Mobility and dispersal or early life stages",
  "Species range"
)

## Three movement attributes whose tally counts are swapped before draw-pile
## construction. Swapping L↔VH and M↔H is equivalent to the 5-mean inversion
## applied in Script 1: L(1)↔VH(4), M(2)↔H(3).
inverted_attrs <- c(
  "Adult mobility",
  "Habitat specificity",
  "Mobility and dispersal or early life stages"
)

##------------------------------------------------------------------------------
## File paths

proj_dir     <- "."
tallies_dir  <- file.path(proj_dir, "outputs", "final-tallies-long")
out_dir      <- file.path(proj_dir, "outputs", "distribution-change-potential")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

f_sens_tallies <- file.path(tallies_dir, "sensitivity_tallies_long.csv")
f_dcp_baseline <- file.path(out_dir,     "distributional_change_potential_uscar.csv")
f_boot_out     <- file.path(out_dir,     "distributional_change_bootstrap_uscar.csv")
f_full_out     <- file.path(out_dir,     "distributional_change_full_uscar.csv")

##------------------------------------------------------------------------------
## FCVA logic model (copied from Module 8 — uncertainty-analyses.R)

fcva_logic_model <- function(mean_scores, rank_threshold) {

  n_ge_35 <- sum(mean_scores >= 3.5, na.rm = TRUE)
  n_ge_30 <- sum(mean_scores >= 3.0, na.rm = TRUE)
  n_ge_25 <- sum(mean_scores >= 2.5, na.rm = TRUE)

  component_rank <- dplyr::case_when(
    n_ge_35 > rank_threshold + 1 ~ "Very High",
    n_ge_30 > rank_threshold     ~ "High",
    n_ge_25 > rank_threshold     ~ "Moderate",
    TRUE                         ~ "Low"
  )

  component_score_numeric <- dplyr::case_when(
    component_rank == "Low"       ~ 1L,
    component_rank == "Moderate"  ~ 2L,
    component_rank == "High"      ~ 3L,
    component_rank == "Very High" ~ 4L
  )

  tibble::tibble(component_rank, component_score_numeric)
}

##------------------------------------------------------------------------------
## Read inputs

dcp_baseline     <- readr::read_csv(f_dcp_baseline, show_col_types = FALSE)
sens_tallies_raw <- readr::read_csv(f_sens_tallies,  show_col_types = FALSE)

##------------------------------------------------------------------------------
## Stock name normalization (matches Module 8 — uncertainty-analyses.R)
##
## Tally files use Title Case workbook names; compiled output files use
## sentence case. str_to_sentence() handles most stocks; four require manual
## overrides because their common names differ beyond simple capitalization.

tally_name_overrides <- c(
  "Atlantic Herring" = "Atlantic thread herring",
  "Diadema"          = "Long-spined sea urchin",
  "Redhind"          = "Red hind",
  "Sea Cucumber"     = "Sea cucumbers"
)

##------------------------------------------------------------------------------
## Filter tallies to target attributes and normalize stock names

dcp_tallies_std <- sens_tallies_raw %>%
  dplyr::filter(
    attribute_type == "Sensitivity",
    attribute_name %in% target_attributes
  ) %>%
  dplyr::transmute(
    stock_name     = stringr::str_to_sentence(
                       dplyr::recode(stock_name, !!!tally_name_overrides)
                     ),
    reviewer_id    = reviewer_id,
    attribute_name = stringr::str_squish(attribute_name),
    tally_L  = tidyr::replace_na(as.integer(tally_L),  0L),
    tally_M  = tidyr::replace_na(as.integer(tally_M),  0L),
    tally_H  = tidyr::replace_na(as.integer(tally_H),  0L),
    tally_VH = tidyr::replace_na(as.integer(tally_VH), 0L),
    n_tallies = as.integer(n_tallies)
  )

## Halt if any target attribute is missing from the tally file
attrs_in_tallies <- unique(dcp_tallies_std$attribute_name)
missing_attrs <- setdiff(target_attributes, attrs_in_tallies)
if (length(missing_attrs) > 0) {
  stop(
    "Target attributes not found in sensitivity_tallies_long.csv:\n  ",
    paste(missing_attrs, collapse = "\n  ")
  )
}
message("✓ All four target attributes found in sensitivity_tallies_long.csv.")

##------------------------------------------------------------------------------
## Invert tally counts for movement attributes
##
## For each inverted attribute, the swap L↔VH, M↔H converts the raw tally
## counts into counts of inverted votes. A draw pile built from these swapped
## counts has a mean equal to 5 - (original draw pile mean), which is the same
## transformation applied in Script 1 to the attribute_mean.

dcp_tallies_inv <- dcp_tallies_std %>%
  dplyr::mutate(
    tally_L_new  = dplyr::if_else(attribute_name %in% inverted_attrs, tally_VH, tally_L),
    tally_M_new  = dplyr::if_else(attribute_name %in% inverted_attrs, tally_H,  tally_M),
    tally_H_new  = dplyr::if_else(attribute_name %in% inverted_attrs, tally_M,  tally_H),
    tally_VH_new = dplyr::if_else(attribute_name %in% inverted_attrs, tally_L,  tally_VH)
  ) %>%
  dplyr::select(-tally_L, -tally_M, -tally_H, -tally_VH) %>%
  dplyr::rename(
    tally_L  = tally_L_new,
    tally_M  = tally_M_new,
    tally_H  = tally_H_new,
    tally_VH = tally_VH_new
  )

##------------------------------------------------------------------------------
## Aggregate (inverted) tallies across reviewers → draw pile totals

dcp_tallies_agg <- dcp_tallies_inv %>%
  dplyr::group_by(stock_name, attribute_name) %>%
  dplyr::summarise(
    total_L        = sum(tally_L),
    total_M        = sum(tally_M),
    total_H        = sum(tally_H),
    total_VH       = sum(tally_VH),
    draw_pile_size = sum(tally_L + tally_M + tally_H + tally_VH),
    .groups = "drop"
  )

cat("Draw pile size summary (expected = 20 per stock × attribute):\n")
print(summary(dcp_tallies_agg$draw_pile_size))

##------------------------------------------------------------------------------
## Baseline reproduction gate
##
## Compute the mean of each inverted draw pile, apply fcva_logic_model(),
## and verify the resulting DCP ranks match Script 1 output for all 25 stocks.
## A mismatch means the tally inversion or rank_threshold is inconsistent
## with Script 1 — do not bootstrap until resolved.

baseline_repro <- dcp_tallies_agg %>%
  dplyr::mutate(
    draw_pile_mean = (total_L * 1L + total_M * 2L +
                      total_H * 3L + total_VH * 4L) / draw_pile_size
  ) %>%
  dplyr::group_by(stock_name) %>%
  dplyr::summarise(
    dcp_logic = list(fcva_logic_model(draw_pile_mean, rank_threshold)),
    .groups = "drop"
  ) %>%
  tidyr::unnest(cols = dcp_logic) %>%
  dplyr::rename(
    dcp_rank_repro    = component_rank,
    dcp_numeric_repro = component_score_numeric
  )

gate_check <- dcp_baseline %>%
  dplyr::select(stock_name, dcp_rank) %>%
  dplyr::left_join(baseline_repro, by = "stock_name") %>%
  dplyr::mutate(rank_match = dcp_rank == dcp_rank_repro)

mismatches <- dplyr::filter(gate_check, is.na(rank_match) | !rank_match)

if (nrow(mismatches) > 0) {
  message("Baseline reproduction mismatches:")
  print(dplyr::select(mismatches, stock_name, dcp_rank, dcp_rank_repro), n = 50)
  stop(
    "Baseline reproduction gate failed for ", nrow(mismatches), " stock(s). ",
    "Verify that rank_threshold matches Script 1 and that the tally inversion is correct. ",
    "Check reproduced vs baseline ranks above."
  )
}
message("✓ Baseline reproduction gate passed for all ", nrow(gate_check), " stocks.")

##------------------------------------------------------------------------------
## Bootstrap

all_ranks  <- c("Low", "Moderate", "High", "Very High")
stock_list <- sort(unique(dcp_baseline$stock_name))
n_stocks   <- length(stock_list)

set.seed(bootstrap_seed)

boot_rows <- vector("list", n_stocks)

message("Running bootstrap (n_boot = ", n_boot, ") for ", n_stocks, " stocks ...")

for (si in seq_along(stock_list)) {

  s <- stock_list[si]
  message("  → ", s, "  (", si, " / ", n_stocks, ")")

  ## Aggregated inverted tally data for this stock
  stock_agg <- dcp_tallies_agg[dcp_tallies_agg$stock_name == s, ]
  att_names <- sort(stock_agg$attribute_name)
  n_atts    <- length(att_names)

  ## Build draw piles once outside the bootstrap loop
  draw_piles <- vector("list", n_atts)
  names(draw_piles) <- att_names

  for (ai in seq_len(n_atts)) {
    att <- att_names[ai]
    r   <- stock_agg[stock_agg$attribute_name == att, ]
    draw_piles[[att]] <- c(
      rep(1L, r$total_L),
      rep(2L, r$total_M),
      rep(3L, r$total_H),
      rep(4L, r$total_VH)
    )
  }

  ## Pre-allocate rank record
  dcp_rank_rec <- character(n_boot)

  ## Bootstrap loop: resample draw piles, compute means, apply logic model
  for (bi in seq_len(n_boot)) {
    boot_att_means <- numeric(n_atts)
    for (ai in seq_len(n_atts)) {
      pile <- draw_piles[[att_names[ai]]]
      boot_att_means[ai] <- mean(sample(pile, size = length(pile), replace = TRUE))
    }
    logic_result     <- fcva_logic_model(boot_att_means, rank_threshold)
    dcp_rank_rec[bi] <- logic_result$component_rank
  }

  ## Summarize rank distribution; ensure all four ranks appear even if count = 0
  counts <- table(dcp_rank_rec)
  n_vec  <- setNames(rep(0L, 4L), all_ranks)
  n_vec[names(counts)] <- as.integer(counts)

  dom_rank <- all_ranks[which.max(n_vec)]
  dom_prop <- max(n_vec) / n_boot

  boot_rows[[si]] <- data.frame(
    stock_name    = s,
    prop_L        = n_vec["Low"]       / n_boot,
    prop_M        = n_vec["Moderate"]  / n_boot,
    prop_H        = n_vec["High"]      / n_boot,
    prop_VH       = n_vec["Very High"] / n_boot,
    dominant_rank = dom_rank,
    dominant_prop = dom_prop,
    borderline    = dom_prop < (1 - borderline_threshold),
    stringsAsFactors = FALSE
  )
}

boot_summary <- dplyr::bind_rows(boot_rows)

##------------------------------------------------------------------------------
## Validate and summarize bootstrap results

## Proportions must sum to 1.0 per stock within floating-point tolerance
prop_sums  <- boot_summary$prop_L + boot_summary$prop_M +
              boot_summary$prop_H + boot_summary$prop_VH
bad_sums   <- which(abs(prop_sums - 1.0) > 1e-9)
if (length(bad_sums) > 0) {
  warning("Bootstrap proportion sums deviate from 1.0 for: ",
          paste(boot_summary$stock_name[bad_sums], collapse = ", "))
} else {
  message("✓ Bootstrap proportions sum to 1.000 for all stocks.")
}

cat("\nBorderline stocks (dominant rank < 75% of iterations):\n")
borderline_stocks <- dplyr::filter(boot_summary, borderline)
if (nrow(borderline_stocks) == 0) {
  cat("  None\n")
} else {
  print(
    dplyr::select(borderline_stocks, stock_name, dominant_rank, dominant_prop),
    row.names = FALSE
  )
}

##------------------------------------------------------------------------------
## Write outputs

readr::write_csv(boot_summary, f_boot_out)
cat("\nBootstrap output written to:", f_boot_out, "\n")

## Join Scripts 1 + 2 for plotting
dcp_full <- dcp_baseline %>%
  dplyr::left_join(boot_summary, by = "stock_name")

readr::write_csv(dcp_full, f_full_out)
cat("Joined output written to:", f_full_out, "\n")
