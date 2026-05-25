################################################################################
##------------------------------------------------------------------------------
## Caribbean CVA
## Module 10 — Results tables
## Script 5: Produce directional effect and data quality results tables
##
## Produces two CSV tables saved to outputs/tables/:
##
##   table_directional_effect_results.csv
##     One row per stock. Columns: Stock, Vuln_rank, Dir_effect, Wt_mean,
##     N_negative, N_neutral, N_positive, N_tallies,
##     Boot_Negative, Boot_Neutral, Boot_Positive, Dominant_prop, Borderline.
##     Stocks ordered to match Results Table 1 (Vuln_rank descending, then
##     bootstrap certainty descending within each rank group).
##
##   table_data_quality_results.csv
##     One row per stock. Columns: Stock, Vuln_rank, Data_quality_rank,
##     Prop_ge_2, Mean_score, N_adequate, N_limited, N_expert, N_nodata.
##     Same stock ordering as the directional effect table.
##
## Inputs:
##   outputs/final-scores-compiled/overall-vulnerability-rankings/
##     overall_vulnerability_scores_uscar.csv
##   outputs/final-tallies-long/
##     directional_effect_tallies_by_stock.csv
##   outputs/analyses/uncertainty-loo/final-tables/
##     table_directional_effect_bootstrap.csv
##     table_bootstrap_uncertainty_stock.csv
##   outputs/final-scores-compiled/data-quality/
##     overall_data_quality_summary_by_stock.csv
##
## Run from the Caribbean-CVA RStudio project root (.Rproj file).
##------------------------------------------------------------------------------

rm(list = ls()); gc()
source("config.R")

library(dplyr)
library(tidyr)
library(readr)

##------------------------------------------------------------------------------
## File paths

proj_dir <- "."
out_dir  <- file.path(proj_dir, "outputs", run_label, "tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

f_vuln      <- file.path(proj_dir, "outputs", run_label,
                         "final-scores-compiled", "overall-vulnerability-rankings",
                         "overall_vulnerability_scores_uscar.csv")
f_dir_tally <- file.path(proj_dir, "outputs", "final-tallies-long",
                         "directional_effect_tallies_by_stock.csv")
f_dir_boot  <- file.path(proj_dir, "outputs", run_label, "analyses", "uncertainty-loo",
                         "final-tables", "table_directional_effect_bootstrap.csv")
f_vuln_boot <- file.path(proj_dir, "outputs", run_label, "analyses", "uncertainty-loo",
                         "final-tables", "table_bootstrap_uncertainty_stock.csv")
f_dq        <- file.path(proj_dir, "outputs", "final-scores-compiled",
                         "data-quality", "overall_data_quality_summary_by_stock.csv")

##------------------------------------------------------------------------------
## Read inputs and apply canonical stock names

vuln      <- read_csv(f_vuln,      show_col_types = FALSE) %>%
  mutate(stock_name = recode(stock_name, !!!stock_name_recode))

dir_tally <- read_csv(f_dir_tally, show_col_types = FALSE) %>%
  mutate(stock_name = recode(stock_name, !!!stock_name_recode))

dir_boot  <- read_csv(f_dir_boot,  show_col_types = FALSE) %>%
  mutate(stock_name = recode(stock_name, !!!stock_name_recode))

vuln_boot <- read_csv(f_vuln_boot, show_col_types = FALSE) %>%
  mutate(stock_name = recode(stock_name, !!!stock_name_recode))

dq        <- read_csv(f_dq,        show_col_types = FALSE) %>%
  mutate(stock_name = recode(stock_name, !!!stock_name_recode))

##------------------------------------------------------------------------------
## Sort key: dominant vulnerability bootstrap proportion per stock
## Stocks with higher certainty sort first within each Vuln_rank group,
## matching the ordering in Results Table 1 and Figure 7.

sort_key <- vuln_boot %>%
  group_by(stock_name) %>%
  slice_max(prop, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(stock_name, vuln_dom_prop = prop)

##------------------------------------------------------------------------------
## TABLE 1 — Directional effect results
##------------------------------------------------------------------------------

## Pivot bootstrap from long (one row per stock × dir_rank) to wide
dir_boot_wide <- dir_boot %>%
  select(stock_name, dir_rank, prop, baseline_dir_rank, baseline_w_mean) %>%
  pivot_wider(
    id_cols      = c(stock_name, baseline_dir_rank, baseline_w_mean),
    names_from   = dir_rank,
    values_from  = prop,
    names_prefix = "Boot_"
  ) %>%
  mutate(
    Boot_Negative = replace_na(Boot_Negative, 0),
    Boot_Neutral  = replace_na(Boot_Neutral,  0),
    Boot_Positive = replace_na(Boot_Positive, 0),
    Dominant_prop = pmax(Boot_Negative, Boot_Neutral, Boot_Positive),
    Borderline    = Dominant_prop < borderline_prop
  )

table_dir <- vuln %>%
  select(stock_name, Vuln_rank, Vuln_score) %>%
  left_join(
    dir_tally %>% select(stock_name,
                         N_negative = tally_neg,
                         N_neutral  = tally_neut,
                         N_positive = tally_pos,
                         N_tallies  = n_tallies),
    by = "stock_name"
  ) %>%
  left_join(dir_boot_wide, by = "stock_name") %>%
  left_join(sort_key,      by = "stock_name") %>%
  mutate(
    Vuln_rank     = factor(Vuln_rank, levels = rank_levels),
    Wt_mean       = round(baseline_w_mean, 3),
    Boot_Negative = round(Boot_Negative, 4),
    Boot_Neutral  = round(Boot_Neutral,  4),
    Boot_Positive = round(Boot_Positive, 4)
  ) %>%
  arrange(desc(Vuln_rank), desc(vuln_dom_prop), stock_name) %>%
  select(
    Stock         = stock_name,
    Vuln_rank,
    Dir_effect    = baseline_dir_rank,
    Wt_mean,
    N_negative, N_neutral, N_positive, N_tallies,
    Boot_Negative, Boot_Neutral, Boot_Positive,
    Dominant_prop, Borderline
  )

write_csv(table_dir, file.path(out_dir, "table_directional_effect_results.csv"))
cat("Table 1 written:", nrow(table_dir), "rows ->",
    file.path(out_dir, "table_directional_effect_results.csv"), "\n")

##------------------------------------------------------------------------------
## TABLE 2 — Data quality results
##------------------------------------------------------------------------------

table_dq <- vuln %>%
  select(stock_name, Vuln_rank, Vuln_score) %>%
  left_join(dq,       by = "stock_name") %>%
  left_join(sort_key, by = "stock_name") %>%
  mutate(Vuln_rank = factor(Vuln_rank, levels = rank_levels)) %>%
  arrange(desc(Vuln_rank), desc(vuln_dom_prop), stock_name) %>%
  select(
    Stock             = stock_name,
    Vuln_rank,
    Data_quality_rank = data_quality_rank,
    Prop_ge_2         = prop_ge_2,
    Mean_score        = mean_score,
    N_adequate        = n_3,
    N_limited         = n_2,
    N_expert          = n_1,
    N_nodata          = n_0
  )

write_csv(table_dq, file.path(out_dir, "table_data_quality_results.csv"))
cat("Table 2 written:", nrow(table_dq), "rows ->",
    file.path(out_dir, "table_data_quality_results.csv"), "\n")

cat("=== Results tables complete ===\n")
