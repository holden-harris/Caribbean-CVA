rm(list = ls()); gc()

##------------------------------------------------------------------------------
## User setup

in_dir  <- "./outputs/final-scores-compiled/final-attribute-scores"
out_dir <- "./outputs/final-scores-compiled/overall-vulnerability-rankings/"

library(dplyr)
library(tidyr)

##------------------------------------------------------------------------------
## Load compiled qualitative CVA scores

score_table <- read.csv(
  file.path(in_dir, "table_final_attribute_scores_all.csv"),
  stringsAsFactors = FALSE
)

##------------------------------------------------------------------------------
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

##------------------------------------------------------------------------------
## QA checks
head(score_table_all, n = 110)

## Check region values
unique(score_table_all$region)

## Count rows by score type and region
score_table_all %>%
  count(score_type, region) ## Qualitative rows should all have NA region


##------------------------------------------------------------------------------
## Write final compiled table
write.csv(
  score_table_all,
  file.path(out_dir, "final_scores_compiled.csv"),
  row.names = FALSE
)

##------------------------------------------------------------------------------
## Make table for U.S. Caribbean Only
score_table_uscar <- score_table_all %>%
  filter(!region %in% c("Wider Caribbean", "Western Atlantic"))

head(score_table_uscar, n = 110)

##------------------------------------------------------------------------------
## Write final compiled table for U.S. Caribbean
write.csv(
  score_table_uscar,
  file.path(out_dir, "final_scores_uscar.csv"),
  row.names = FALSE
)


##------------------------------------------------------------------------------
## QA check: count rows per stock x attribute x score_type

qa_attribute_counts <- score_table_uscar %>%
  group_by(stock_name, attribute_type, attribute_name, score_type) %>%
  summarise(
    n_rows = n(),
    n_scorers = n_distinct(scorer),
    .groups = "drop"
  )

print(qa_attribute_counts, n = 200)


##------------------------------------------------------------------------------
## QA check using expected counts by score_type
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
qa_problem_scores

##------------------------------------------------------------------------------
## Pull the underlying rows from score_table_uscar for the problem cases

qa_problem_rows <- score_table_uscar %>%
  inner_join(
    qa_problem_scores %>%
      select(stock_name, attribute_type, score_type, attribute_name),
    by = c("stock_name", "attribute_type", "score_type", "attribute_name")
  ) %>%
  arrange(stock_name, attribute_type, attribute_name, scorer)

################################################################################
##------------------------------------------------------------------------------
## NOAA FCVA workflow
##
## Step 1. Calculate average score for each individual attribute/factor
## Step 2. Apply logic model to score overall Sensitivity and Exposure per species
## Step 3. Multiply component numeric scores to get overall vulnerability rank
##
## Notes:
## - Sensitivity comes from reviewer-entered final attribute scores
## - Exposure comes from the appended calculated quantitative exposure scores
## - The mean Exposure and mean Sensitivity columns below are descriptive only;
##   they are NOT used directly in the final vulnerability product


library(dplyr)
library(tidyr)

################################################################################
##------------------------------------------------------------------------------
## Step 1 - Calculate average score for each  attribute
##
## These are the cross-reviewer means for each stock x sensitivity attribute.
## This gives one average score per attribute per stock.
##
## Step 1B - Categorize each attribute mean into LMHV bins
## component_numeric_score:
##     4 = Very High
##     3 = High
##     2 = Moderate
##     1 = Low
##
## This is mainly for inspection / QA.
## The actual FCVA logic model below uses threshold counts directly.


attribute_means_uscar <- score_table_uscar %>%
  group_by(stock_name, attribute_type, score_type, attribute_name) %>%
  summarise(
    attribute_mean = mean(score, na.rm = TRUE),
    attribute_sd   = sd(score, na.rm = TRUE),
    n_scores       = sum(!is.na(score)),
    .groups = "drop"
  ) %>%
  mutate(
    attribute_mean_lmhv = case_when(
      is.na(attribute_mean) ~ NA_character_,
      attribute_mean < 2 ~ "Low",
      attribute_mean >= 2 & attribute_mean < 3.0 ~ "Moderate",
      attribute_mean >= 3.0 & attribute_mean < 3.5 ~ "High",
      attribute_mean >= 3.5 ~ "Very High"
    )
  )

attribute_means_uscar$attribute_mean <-  round(attribute_means_uscar$attribute_mean, 2)
attribute_means_uscar$attribute_sd <-  round(attribute_means_uscar$attribute_sd, 2)
print(attribute_means_uscar, n = 60)

##------------------------------------------------------------------------------
## Write attribute means table for U.S. Caribbean
write.csv(
  attribute_means_uscar,
  file.path(out_dir, "attribute_means_uscar.csv"),
  row.names = FALSE
)

##------------------------------------------------------------------------------
## Step 2 - Apply the FCVA logic model to determine overall component score
## per species for Sensitivity and Exposure
##
## Table-based rule used in prior FCVAs:
## - Very High = more than 3 attribute means >= 3.5
## - High      = more than 2 attribute means >= 3.0
## - Moderate  = more than 2 attribute means >= 2.5
## - Low       = all other cases
##
## Notes:
## - Sensitivity should come from qualitative scores
## - Exposure should come from calculated scores

component_scores <- attribute_means_uscar %>%
  filter(
    (attribute_type == "Sensitivity" & score_type == "Qualitative") |
      (attribute_type == "Exposure"    & score_type == "Calculated")
  ) %>%
  group_by(stock_name, attribute_type) %>%
  summarise(
    ## Descriptive summary of the attribute means within each component
    mean_attribute_score = mean(attribute_mean, na.rm = TRUE),
    sd_attribute_score   = sd(attribute_mean, na.rm = TRUE),
    n_attributes         = n(),
    
    ## Count how many attribute means meet each FCVA threshold
    n_ge_2_5 = sum(attribute_mean >= 2.5, na.rm = TRUE),
    n_ge_3_0 = sum(attribute_mean >= 3.0, na.rm = TRUE),
    n_ge_3_5 = sum(attribute_mean >= 3.5, na.rm = TRUE),
    
    ## Apply the FCVA logic model to assign a numeric component score
    component_score = case_when(
      n_ge_3_5 > 3 ~ 4,
      n_ge_3_0 > 2 ~ 3,
      n_ge_2_5 > 2 ~ 2,
      TRUE         ~ 1
    ),
    
    ## Assign the matching component rank
    component_rank = case_when(
      n_ge_3_5 > 3 ~ "Very High",
      n_ge_3_0 > 2 ~ "High",
      n_ge_2_5 > 2 ~ "Moderate",
      TRUE         ~ "Low"
    ),
    
    .groups = "drop"
  ); print(component_scores)

##------------------------------------------------------------------------------
## QA: each stock should have one Sensitivity row and one Exposure row

component_scores %>%
  count(stock_name, attribute_type) 

## Review the threshold counts that drove each component score
component_scores %>%
  arrange(attribute_type, desc(component_numeric_score), stock_name) %>%
  print(n = 50)

##------------------------------------------------------------------------------
## Write component scores table for U.S. Caribbean
write.csv(
  component_scores,
  file.path(out_dir, "component_scores_uscar.csv"),
  row.names = FALSE
)

################################################################################
##------------------------------------------------------------------------------
## Step 3 - Convert Sensitivity and Exposure component scores to one row per stock
## and calculate final overall vulnerability as the product of the two component
## numeric scores

stock_vulnerability <- component_scores %>%
  select(
    stock_name,
    attribute_type,
    component_score,
    component_rank
  ) %>%
  pivot_wider(
    names_from  = attribute_type,
    values_from = c(component_score, component_rank)
  ) %>%
  rename(
    Exp_score       = component_score_Exposure,
    Exp_rank    = component_rank_Exposure,
    Sens_score    = component_score_Sensitivity,
    Sens_rank = component_rank_Sensitivity
  ) %>%
  mutate(
    Vuln_score = Exp_score * Sens_score,
    Vuln_rank = case_when(
      Vuln_score <= 3 ~ "Low",
      Vuln_score >= 4  & Vuln_score <= 6  ~ "Moderate",
      Vuln_score >= 8  & Vuln_score <= 9  ~ "High",
      Vuln_score >= 12 & Vuln_score <= 16 ~ "Very High",
      TRUE ~ NA_character_
    )
  ); print(stock_vulnerability, n = 25)

##------------------------------------------------------------------------------
## Write component scores table for U.S. Caribbean
write.csv(
  stock_vulnerability,
  file.path(out_dir, "overall_vulnerability_scores_uscar.csv"),
  row.names = FALSE
)



