
##------------------------------------------------------------------------------
## `2-finalize-tally-tables.R`
##
## Produces two sets of tables in ./outputs/final-tallies-long/:
##
## Set 1 — reviewer-level long tables (one row per scorer x stock x category):
##   - sensitivity_tallies_long.csv
##   - directional_effect_tallies_long.csv
##   - exposure_tallies_long.csv  (qual: reviewer-level; quant: stock x factor)
##
## Set 2 — final summary tables (one row per stock x grouping variable):
##   - sensitivity_final.csv        (stock x attribute; tallies summed, n_scorers)
##   - directional_effect_final.csv (stock; tally_neg/neut/pos, n_tallies, n_scorers)
##   - exposure_final.csv           (stock x attribute; qual grouped + quant)

## Set up ----------------------------------------------------------------------

rm(list = ls()); gc()
library(dplyr)

## Directories -----------------------------------------------------------------
dir_in    <- "./outputs/analyses/1-inputs"
dir_quant <- "./outputs/final-scores-compiled"
dir_out   <- "./outputs/final-tallies-long"
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

## Canonical stock name lookup -------------------------------------------------
stock_name_recode <- c(
  "Atlantic thread herring" = "Atlantic Herring",
  "Long-spined sea urchin"  = "Diadema",
  "Red hind"                = "Redhind",
  "Sea cucumbers"           = "Sea Cucumber",
  "Ballyhoo"                = "Ballyhoo",
  "Blue runner"             = "Blue Runner",
  "Dolphinfish"             = "Dolphinfish",
  "Gray angelfish"          = "Gray Angelfish",
  "Hogfish"                 = "Hogfish",
  "King mackerel"           = "King Mackerel",
  "Lane snapper"            = "Lane Snapper",
  "Misty grouper"           = "Misty Grouper",
  "Mutton snapper"          = "Mutton Snapper",
  "Nassau grouper"          = "Nassau Grouper",
  "Queen conch"             = "Queen Conch",
  "Queen triggerfish"       = "Queen Triggerfish",
  "Rainbow parrotfish"      = "Rainbow Parrotfish",
  "Red grouper"             = "Red Grouper",
  "Redhind"                 = "Red Hind",
  "Sea cucumber"            = "Sea Cucumber",
  "Silk snapper"            = "Silk Snapper",
  "Spiny lobster"           = "Spiny Lobster",
  "Stoplight parrotfish"    = "Stoplight Parrotfish",
  "White mullet"            = "White Mullet",
  "Yellowfin grouper"       = "Yellowfin Grouper",
  "Yellowtail snapper"      = "Yellowtail Snapper"
)

################################################################################
## Sensitivity

sensitivity_tallies_long <- read.csv(
  file.path(dir_in, "table_sensitivity_tallies_long.csv")) %>%
  mutate(stock_name = recode(stock_name, !!!stock_name_recode))

sensitivity_final <- sensitivity_tallies_long %>%
  group_by(stock_name, attribute_name) %>%
  summarise(
    tally_L   = sum(tally_L,   na.rm = TRUE),
    tally_M   = sum(tally_M,   na.rm = TRUE),
    tally_H   = sum(tally_H,   na.rm = TRUE),
    tally_VH  = sum(tally_VH,  na.rm = TRUE),
    n_scorers = n_distinct(reviewer_id),
    .groups   = "drop"
  ) %>%
  mutate(n_tallies = tally_L + tally_M + tally_H + tally_VH)

################################################################################
## Directional effect

directional_effect_tallies_long <- read.csv(
  file.path(dir_in, "table_directional_effect_tallies_long.csv")) %>%
  mutate(stock_name = recode(stock_name, !!!stock_name_recode)) %>%
  select(reviewer_id, stock_name, effect_category, tally)

directional_effect_final <- directional_effect_tallies_long %>%
  group_by(stock_name) %>%
  summarise(
    tally_neg  = sum(tally[effect_category == "Negative"], na.rm = TRUE),
    tally_neut = sum(tally[effect_category == "Neutral"],  na.rm = TRUE),
    tally_pos  = sum(tally[effect_category == "Positive"], na.rm = TRUE),
    n_scorers  = n_distinct(reviewer_id),
    .groups    = "drop"
  ) %>%
  mutate(n_tallies = tally_neg + tally_neut + tally_pos)

################################################################################
## Exposure

qualitative_exposure_tallies_long <- read.csv(
  file.path(dir_in, "table_qualitative_exposure_tallies_long.csv"))

quantitative_exposure_scores <- read.csv(
  file.path(dir_quant, "quantitative-exposure-attribute-scores-all.csv"))

## Quant conform (shared by both long and final exposure tables)
quant_exp_conform <- quantitative_exposure_scores %>%
  filter(spatial_extent == "U.S. Caribbean") %>%
  mutate(
    attribute_type = "Quantitative Exposure",
    attribute_name = full_names,
    n_tallies      = tally_L + tally_M + tally_H + tally_VH
  ) %>%
  select(stock_name, attribute_type, attribute_name,
         tally_L, tally_M, tally_H, tally_VH, n_tallies)

## Long: qualitative rows kept at reviewer level
qual_exp_long <- qualitative_exposure_tallies_long %>%
  filter(attribute_name != "Coral cover") %>%
  mutate(stock_name = recode(stock_name, !!!stock_name_recode)) %>%
  select(reviewer_id, stock_name, attribute_type, attribute_name,
         tally_L, tally_M, tally_H, tally_VH) %>%
  mutate(n_tallies = tally_L + tally_M + tally_H + tally_VH)

exposure_tallies_long <- bind_rows(qual_exp_long, quant_exp_conform) %>%
  arrange(stock_name, attribute_type)

## Final: qualitative rows grouped to stock x attribute
qual_exp_final <- qualitative_exposure_tallies_long %>%
  filter(attribute_name != "Coral cover") %>%
  mutate(stock_name = recode(stock_name, !!!stock_name_recode)) %>%
  group_by(stock_name, attribute_type, attribute_name) %>%
  summarise(
    tally_L   = sum(tally_L,   na.rm = TRUE),
    tally_M   = sum(tally_M,   na.rm = TRUE),
    tally_H   = sum(tally_H,   na.rm = TRUE),
    tally_VH  = sum(tally_VH,  na.rm = TRUE),
    n_scorers = n_distinct(reviewer_id),
    .groups   = "drop"
  ) %>%
  mutate(n_tallies = tally_L + tally_M + tally_H + tally_VH)

exposure_final <- bind_rows(qual_exp_final, quant_exp_conform) %>%
  arrange(stock_name, attribute_type)

################################################################################
## Write outputs

write.csv(sensitivity_tallies_long,
          file.path(dir_out, "sensitivity_tallies_long.csv"),       row.names = FALSE)
write.csv(sensitivity_final,
          file.path(dir_out, "sensitivity_tallies_by_stock.csv"),   row.names = FALSE)
write.csv(directional_effect_tallies_long,
          file.path(dir_out, "directional_effect_tallies_long.csv"), row.names = FALSE)
write.csv(directional_effect_final,
          file.path(dir_out, "directional_effect_tallies_by_stock.csv"), row.names = FALSE)
write.csv(exposure_tallies_long,
          file.path(dir_out, "exposure_tallies_long.csv"),           row.names = FALSE)
write.csv(exposure_final,
          file.path(dir_out, "exposure_tallies_by_stock.csv"),     row.names = FALSE)

################################################################################
## Row count summary

cat("--- Long tables (reviewer-level) ---\n")
cat("sensitivity_tallies_long:        ", nrow(sensitivity_tallies_long),        "rows\n")
cat("directional_effect_tallies_long: ", nrow(directional_effect_tallies_long), "rows\n")
cat("exposure_tallies_long:           ", nrow(exposure_tallies_long),           "rows\n")
cat("\n--- Final tables (grouped) ---\n")
cat("sensitivity_final:               ", nrow(sensitivity_final),               "rows\n")
cat("directional_effect_final:        ", nrow(directional_effect_final),        "rows\n")
cat("exposure_final:                  ", nrow(exposure_final),                  "rows\n")
cat("  distinct attribute names:      ", n_distinct(exposure_final$attribute_name), "\n")
cat("\nDone. Tables written to:", dir_out, "\n")
