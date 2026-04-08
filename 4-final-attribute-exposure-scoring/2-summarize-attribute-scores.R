rm(list = ls()); gc()

##------------------------------------------------------------------------------
## User setup

in_dir  <- "./outputs/final-scores-compiled/final-attribute-scores"
out_dir <- "./outputs/final-scores-compiled/final-attribute-scores/overall-vulnerability-rankings/"

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
  score_table_all,
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

attribute_means <- score_table_uscar %>%
  group_by(stock_name, attribute_type, score_type, attribute_name) %>%
  summarise(
    attribute_mean = mean(score, na.rm = TRUE),
    attribute_sd   = sd(score, na.rm = TRUE),
    n_scores       = sum(!is.na(score)),
    .groups = "drop"
  ) 

print(attribute_means, n = 60)

## Categorize attribute score LMHV

################################################################################
##------------------------------------------------------------------------------
## Step 2 - Apply NOAA FCVA logic model 
##
## Input:
##   x = vector of weighted-average attribute/factor scores for one stock and one
##       component (Sensitivity or Exposure)
##
## Output:
##   component_numeric_score:
##     4 = Very High
##     3 = High
##     2 = Moderate
##     1 = Low
##
## Logic rule (following Morrison / HMS CVA):
##   Very High = 3 or more factors with mean >= 3.5
##   High      = 2 or more factors with mean >= 3.0
##   Moderate  = 2 or more factors with mean >= 2.5
##   Low       = all other cases

fcva_logic_model <- function(x) {
  
  x <- x[is.finite(x)]
  
  n_ge_2_5 <- sum(x >= 2.5, na.rm = TRUE)
  n_ge_3_0 <- sum(x >= 3.0, na.rm = TRUE)
  n_ge_3_5 <- sum(x >= 3.5, na.rm = TRUE)
  
  if (length(x) == 0) {
    return(
      tibble(
        component_numeric_score = NA_real_,
        component_rank = NA_character_,
        n_ge_2_5 = NA_integer_,
        n_ge_3_0 = NA_integer_,
        n_ge_3_5 = NA_integer_
      )
    )
  }
  
  if (n_ge_3_5 >= 3) {
    out_score <- 4
    out_rank  <- "Very High"
  } else if (n_ge_3_0 >= 2) {
    out_score <- 3
    out_rank  <- "High"
  } else if (n_ge_2_5 >= 2) {
    out_score <- 2
    out_rank  <- "Moderate"
  } else {
    out_score <- 1
    out_rank  <- "Low"
  }
  
  tibble(
    component_numeric_score = out_score,
    component_rank = out_rank,
    n_ge_2_5 = n_ge_2_5,
    n_ge_3_0 = n_ge_3_0,
    n_ge_3_5 = n_ge_3_5
  )
}

##------------------------------------------------------------------------------
## Helper function: overall vulnerability rank from product of component scores
##
## Product categories used in prior FCVAs:
##   1-3   = Low
##   4-6   = Moderate
##   8-9   = High
##   12-16 = Very High
fcva_overall_rank <- function(x) {
  case_when(
    is.na(x) ~ NA_character_,
    x <= 3 ~ "Low",
    x >= 4  & x <= 6  ~ "Moderate",
    x >= 8  & x <= 9  ~ "High",
    x >= 12 & x <= 16 ~ "Very High",
    TRUE ~ NA_character_
  )
}




################################################################################
##------------------------------------------------------------------------------
## Step 2. Apply the FCVA logic model separately to Sensitivity and Exposure
##
## For each stock and component:
## - count how many attribute/factor means are >= 2.5, 3.0, and 3.5
## - assign component numeric score and rank using the FCVA decision rules
##
## Also calculate the simple mean of all attribute means for reference only.
## This mean is useful to inspect, but it is not the final component score used
## in the NOAA FCVA product step.

component_scores <- component_attribute_means %>%
  group_by(stock_name, Component) %>%
  group_modify(~{
    
    logic_out <- fcva_logic_model(.x$attribute_mean)
    
    tibble(
      mean_attribute_score = mean(.x$attribute_mean, na.rm = TRUE),
      sd_attribute_score   = sd(.x$attribute_mean, na.rm = TRUE),
      n_attributes         = nrow(.x),
      component_numeric_score = logic_out$component_numeric_score,
      component_rank          = logic_out$component_rank,
      n_ge_2_5                = logic_out$n_ge_2_5,
      n_ge_3_0                = logic_out$n_ge_3_0,
      n_ge_3_5                = logic_out$n_ge_3_5
    )
  }) %>%
  ungroup()

print(component_scores, n = 50)

################################################################################
##------------------------------------------------------------------------------
## Step 3. Convert Sensitivity and Exposure component scores to one row per stock
## and calculate final overall vulnerability as the product of the two component
## numeric scores

attribute_scores_wide <- component_scores %>%
  mutate(Component = tolower(Component)) %>%
  select(
    stock_name,
    Component,
    mean_attribute_score,
    sd_attribute_score,
    n_attributes,
    component_numeric_score,
    component_rank,
    n_ge_2_5,
    n_ge_3_0,
    n_ge_3_5
  ) %>%
  pivot_wider(
    names_from = Component,
    values_from = c(
      mean_attribute_score,
      sd_attribute_score,
      n_attributes,
      component_numeric_score,
      component_rank,
      n_ge_2_5,
      n_ge_3_0,
      n_ge_3_5
    )
  ) %>%
  mutate(
    Vulnerability = case_when(
      is.na(sensitivity_component_numeric_score) |
        is.na(exposure_component_numeric_score) ~ NA_real_,
      TRUE ~ sensitivity_component_numeric_score * exposure_component_numeric_score
    ),
    Overall_rank = fcva_overall_rank(Vulnerability)
  ) %>%
  rename(
    mean_sensitivity_score = sensitivity_mean_attribute_score,
    mean_exposure_score    = exposure_mean_attribute_score,
    sd_sensitivity_score   = sensitivity_sd_attribute_score,
    sd_exposure_score      = exposure_sd_attribute_score,
    n_sensitivity_attributes = sensitivity_n_attributes,
    n_exposure_factors      = exposure_n_attributes,
    Sensitivity_score = sensitivity_component_numeric_score,
    Exposure_score    = exposure_component_numeric_score,
    Sensitivity_rank  = sensitivity_component_rank,
    Exposure_rank     = exposure_component_rank
  ) %>%
  arrange(desc(Vulnerability), stock_name)

print(attribute_scores_wide, n = 25)

################################################################################
##------------------------------------------------------------------------------
## Optional QA checks

## Check the number of attributes/factors contributing to each stock
attribute_scores_wide %>%
  select(
    stock_name,
    n_sensitivity_attributes,
    n_exposure_factors,
    mean_sensitivity_score,
    mean_exposure_score,
    Sensitivity_score,
    Exposure_score,
    Vulnerability,
    Overall_rank
  ) %>%
  print(n = 25)

## Inspect the FCVA threshold counts used to assign each component score
attribute_scores_wide %>%
  select(
    stock_name,
    starts_with("sensitivity_n_ge_"),
    starts_with("exposure_n_ge_"),
    Sensitivity_score,
    Sensitivity_rank,
    Exposure_score,
    Exposure_rank
  ) %>%
  print(n = 25)

## Count final overall ranks
attribute_scores_wide %>%
  count(Overall_rank, sort = TRUE)

################################################################################
##------------------------------------------------------------------------------
## Optional: write outputs

write.csv(
  component_attribute_means,
  file.path(out_dir, "table_component_attribute_means.csv"),
  row.names = FALSE
)

write.csv(
  component_scores,
  file.path(out_dir, "table_component_scores.csv"),
  row.names = FALSE
)

write.csv(
  attribute_scores_wide,
  file.path(out_dir, "table_stock_vulnerability_fcva_logic_model.csv"),
  row.names = FALSE
)