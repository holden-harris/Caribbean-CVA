rm(list = ls()); gc()

##------------------------------------------------------------------------------
## Setup

in_dir  <- "./outputs/final-scores-compiled/final-attribute-scores"
out_dir <- "./outputs/final-scores-compiled/overall-vulnerability-rankings/"

library(dplyr)
library(tidyr)

## Load compiled qualitative CVA scores
score_table <- read.csv(
  file.path(in_dir, "table_final_attribute_scores_all.csv"),
  stringsAsFactors = FALSE
)

## Load calculated quantitative exposure scores
exposure_scores <- read.csv(
  file.path(in_dir, "quantitative-exposure-attribute-scores-all.csv"),
  stringsAsFactors = FALSE
)

##------------------------------------------------------------------------------
## Harmonize quantitative exposure factor names

exposure_name_key <- c(
  "bs"     = "Bottom salinity",
  "bt"     = "Bottom temperature",
  "chl"    = "Chlorophyll-a concentration",
  "mld"    = "Mixed layer depth",
  "msstg"  = "Mean sea surface temperature gradient",
  "o200"   = "Oxygen at 200m",
  "ph"     = "Surface pH",
  "pp"     = "Primary production",
  "precip" = "Precipitation",
  "sso"    = "Sea surface oxygen",
  "sss"    = "Sea surface salinity",
  "sst"    = "Sea surface temperature",
  "swsm"   = "Surface wind speed magnitude"
)

##------------------------------------------------------------------------------
## Harmonize stock names to match the compiled reviewer score table

stock_name_key <- c(
  "Atlantic Herring"      = "Atlantic thread herring",
  "Blue Runner"           = "Blue runner",
  "Diadema"               = "Long-spined sea urchin",
  "Gray Angelfish"        = "Gray angelfish",
  "King Mackerel"         = "King mackerel",
  "Lane Snapper"          = "Lane snapper",
  "Misty Grouper"         = "Misty grouper",
  "Mutton Snapper"        = "Mutton snapper",
  "Nassau Grouper"        = "Nassau grouper",
  "Queen Conch"           = "Queen conch",
  "Queen Snapper"         = "Queen snapper",
  "Queen Triggerfish"     = "Queen triggerfish",
  "Rainbow Parrotfish"    = "Rainbow parrotfish",
  "Red Grouper"           = "Red grouper",
  "Redhind"               = "Red hind",
  "Sea Cucumber"          = "Sea cucumbers",
  "Silk Snapper"          = "Silk snapper",
  "Spiny Lobster"         = "Spiny lobster",
  "Stoplight Parrotfish"  = "Stoplight parrotfish",
  "White Mullet"          = "White mullet",
  "Yellowfin Grouper"     = "Yellowfin grouper",
  "Yellowtail Snapper"    = "Yellowtail snapper"
)

##------------------------------------------------------------------------------
## Standardize the qualitative reviewer score table
##
## region is not applicable to qualitative scores, so set to NA

score_table_qualitative <- score_table %>%
  mutate(
    attribute_type = case_when(
      Attribute_type == "Rigidity" ~ "Sensitivity",
      Attribute_type == "Sensitivity" ~ "Sensitivity",
      Attribute_type == "Qualitative Exposure Factors" ~ "Exposure",
      Attribute_type == "Exposure" ~ "Exposure",
      TRUE ~ Attribute_type
    ),
    score_type     = "Qualitative",
    region         = NA_character_,
    attribute_name = Attribute_name,
    scorer         = Scorer,
    score          = Final_score
  ) %>%
  filter(attribute_type %in% c("Exposure", "Sensitivity")) %>%
  select(
    stock_name,
    region,
    attribute_type,
    score_type,
    attribute_name,
    scorer,
    score
  )

##------------------------------------------------------------------------------
## Standardize the calculated exposure table
##
## region is based on spatial_extent:
## - U.S. Caribbean -> U.S. Caribbean
## - Caribbean Sea  -> Wider Caribbean
## - Western Atlantic -> Western Atlantic

if ("full_names" %in% names(exposure_scores)) {
  exposure_scores <- exposure_scores %>%
    mutate(attribute_name = full_names)
} else {
  exposure_scores <- exposure_scores %>%
    mutate(
      attribute_name = ifelse(
        quantitative_exposure_factor %in% names(exposure_name_key),
        exposure_name_key[quantitative_exposure_factor],
        quantitative_exposure_factor
      )
    )
}

score_table_calculated <- exposure_scores %>%
  mutate(
    stock_name = ifelse(
      stock_name %in% names(stock_name_key),
      stock_name_key[stock_name],
      stock_name
    ),
    region = case_when(
      spatial_extent == "U.S. Caribbean"   ~ "U.S. Caribbean",
      spatial_extent == "Caribbean Sea"    ~ "Wider Caribbean",
      spatial_extent == "Western Atlantic" ~ "Western Atlantic",
      TRUE ~ NA_character_
    ),
    attribute_type = "Exposure",
    score_type     = "Calculated",
    scorer         = "Calculated",
    score          = attribute_score
  ) %>%
  select(
    stock_name,
    region,
    attribute_type,
    score_type,
    attribute_name,
    scorer,
    score
  )

##------------------------------------------------------------------------------
## Append qualitative and calculated scores into one compiled table
score_table_all <- bind_rows(
  score_table_qualitative,
  score_table_calculated
) %>%
  arrange(stock_name, attribute_type, score_type, region, attribute_name, scorer)

## Make table for U.S. Caribbean Only ------------------------------------------
score_table_uscar <- score_table_all %>%
  filter(!region %in% c("Wider Caribbean", "Western Atlantic"))

head(score_table_uscar, n = 75)

##------------------------------------------------------------------------------
## QA checks
head(score_table_all, n = 100) ## Print out first stock

## Count rows by score type and region
score_table_all %>%
  count(score_type, region) ## Qualitative rows should all have NA region

## QA check using expected counts by score_type --------------------------------
## - Calculated should have 1 score
## - Qualitative should have 4 scores
qa_problem_scores <- score_table_uscar %>%
  group_by(stock_name, attribute_type, score_type, attribute_name) %>%
  summarise(
    n_scores = sum(!is.na(score)),
    scorers  = paste(sort(unique(scorer[!is.na(score)])), collapse = ", "),
    .groups = "drop"
  ) %>%
  filter(
    (score_type == "Calculated"  & n_scores != 1) |
      (score_type == "Qualitative" & n_scores != 4)
  ) %>%
  arrange(stock_name, attribute_type, attribute_name)
print(qa_problem_scores) ## Empty table means success :)


## QA check: count rows per stock x attribute x score_type
qa_attribute_counts <- score_table_uscar %>%
  group_by(stock_name, attribute_type, attribute_name, score_type) %>%
  summarise(
    n_rows = n(),
    n_scorers = n_distinct(scorer),
    .groups = "drop"
  ); print(qa_attribute_counts, n = 200)


##------------------------------------------------------------------------------
## Write final compiled table
write.csv(
  score_table_all,
  file.path(out_dir, "final_scores_compiled.csv"),
  row.names = FALSE
)

## Write final compiled table for U.S. Caribbean
write.csv(
  score_table_uscar,
  file.path(out_dir, "final_scores_uscar.csv"),
  row.names = FALSE
)