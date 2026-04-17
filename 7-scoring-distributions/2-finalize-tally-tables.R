
##------------------------------------------------------------------------------
## `2-finalize-tally-tables.R`
##
## Produces three analysis-ready tally tables written to ./outputs/final-tallies-long/:
##   - sensitivity_tallies_long.csv        (reviewer-level; stock names recoded)
##   - directional_effect_tallies_long.csv (one row per stock; tally_neg/neut/pos, n_tallies, n_scorers)
##   - exposure_tallies_long.csv           (stock x attribute; qual + quant, U.S. Caribbean,
##                                          summed across reviewers)

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

## Read and recode sensitivity and directional-effect --------------------------
sensitivity_tallies_long <- read.csv(
  file.path(dir_in, "table_sensitivity_tallies_long.csv")) %>%
  mutate(stock_name = recode(stock_name, !!!stock_name_recode))

directional_effect_tallies_long <- read.csv(
  file.path(dir_in, "table_directional_effect_tallies_long.csv")) %>%
  mutate(stock_name = recode(stock_name, !!!stock_name_recode)) %>%
  group_by(stock_name) %>%
  summarise(
    tally_neg  = sum(tally[effect_category == "Negative"], na.rm = TRUE),
    tally_neut = sum(tally[effect_category == "Neutral"],  na.rm = TRUE),
    tally_pos  = sum(tally[effect_category == "Positive"], na.rm = TRUE),
    n_scorers  = n_distinct(reviewer_id),
    .groups    = "drop"
  ) %>%
  mutate(n_tallies = tally_neg + tally_neut + tally_pos)

## Conform and bind exposure tally tables --------------------------------------
qualitative_exposure_tallies_long <- read.csv(
  file.path(dir_in, "table_qualitative_exposure_tallies_long.csv"))

quantitative_exposure_scores <- read.csv(
  file.path(dir_quant, "quantitative-exposure-attribute-scores-all.csv"))

qual_exp_conform <- qualitative_exposure_tallies_long %>%
  filter(attribute_name != "Coral cover") %>%
  mutate(stock_name = recode(stock_name, !!!stock_name_recode)) %>%
  group_by(stock_name, attribute_type, attribute_name) %>%
  summarise(
    tally_L  = sum(tally_L,  na.rm = TRUE),
    tally_M  = sum(tally_M,  na.rm = TRUE),
    tally_H  = sum(tally_H,  na.rm = TRUE),
    tally_VH = sum(tally_VH, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  mutate(n_tallies = tally_L + tally_M + tally_H + tally_VH)

quant_exp_conform <- quantitative_exposure_scores %>%
  filter(spatial_extent == "U.S. Caribbean") %>%
  mutate(
    attribute_type = "Quantitative Exposure",
    attribute_name = full_names,
    n_tallies      = tally_L + tally_M + tally_H + tally_VH
  ) %>%
  select(stock_name, attribute_type, attribute_name,
         tally_L, tally_M, tally_H, tally_VH, n_tallies)

exposure_tallies_long <- bind_rows(qual_exp_conform, quant_exp_conform) %>%
  arrange(stock_name, attribute_type)

## Row count summary -----------------------------------------------------------
cat("sensitivity_tallies_long:        ", nrow(sensitivity_tallies_long),        "rows\n")
cat("directional_effect_tallies_long: ", nrow(directional_effect_tallies_long), "rows\n")
cat("exposure_tallies_long:           ", nrow(exposure_tallies_long),           "rows\n")
cat("  distinct attribute names:      ", n_distinct(exposure_tallies_long$attribute_name), "\n")

## Write outputs ---------------------------------------------------------------
write.csv(sensitivity_tallies_long,
          file.path(dir_out, "sensitivity_tallies_long.csv"),        row.names = FALSE)
write.csv(directional_effect_tallies_long,
          file.path(dir_out, "directional_effect_tallies_long.csv"), row.names = FALSE)
write.csv(exposure_tallies_long,
          file.path(dir_out, "exposure_tallies_long.csv"),           row.names = FALSE)
cat("\nDone. Tables written to:", dir_out, "\n")
