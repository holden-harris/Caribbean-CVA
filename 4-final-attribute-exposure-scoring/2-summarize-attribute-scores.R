##------------------------------------------------------------------------------
## User setup

in_dir      <- "./outputs/final-scores-compiled/final-attribute-scores"
out_dir     <- "./outputs/final-scores-compiled/final-attribute-scores"

## Breakpoints for assessing overall score
low_cut <- 3
moderate_cut <- 6
high_cut <- 9

## Load compiled CVA score table
score_table <- read.csv(file.path(in_dir, "table_final_attribute_scores_all.csv"), stringsAsFactors = FALSE)

## Set up score types: "Sensitivity" or "Exposure"

score_table$Attribute_type[score_table$Attribute_type == "Rigidity"] <- "Sensitivity" ## Update attribute names: ## All "Rigidity" scores are sensitivity
score_table$Attribute_type[score_table$Attribute_type == "Qualitative Exposure Factors"] <- "Exposure" ## Shorten name for "Exposure"

##------------------------------------------------------------------------------
## Join Exposure scores from the prior exposure analyses
##
## More info: https://github.com/holden-harris/Caribbean-CVA/tree/main/2-exposure-anomalies
## Exposure overlap for all species is available here: 
## https://github.com/holden-harris/Caribbean-CVA/tree/main/outputs/exposure-overlap-12panel

library(dplyr)
exposure_scores <- read.csv(file.path(in_dir, "quantitative-exposure-attribute-scores-all.csv"))

## -----------------------------------------------------------------------------
## Harmonize quantitative exposure factor names
## Only recode the short-code names; keep full names as-is

exposure_name_key <- c(
  "bs"    = "Bottom salinity",
  "bt"    = "Bottom temperature",
  "msstg" = "Mean sea surface temperature gradient",
  "precip"= "Precipitation",
  "sso"   = "Sea surface oxygen",
  "swsm"  = "Surface wind speed magnitude"
)

## -----------------------------------------------------------------------------
## Harmonize stock names to match the compiled score table

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

## -----------------------------------------------------------------------------
## Filter to U.S. Caribbean only and reshape to match score_table

exposure_scores_uscar <- exposure_scores %>%
  filter(spatial_extent == "U.S. Caribbean") %>%
  mutate(
    ## Match stock names to the existing compiled table
    stock_name = ifelse(
      stock_name %in% names(stock_name_key),
      stock_name_key[stock_name],
      stock_name
    ),
    
    ## Recode short quantitative factor names; keep full names unchanged
    Attribute_name = ifelse(
      quantitative_exposure_factor %in% names(exposure_name_key),
      exposure_name_key[quantitative_exposure_factor],
      quantitative_exposure_factor
    ),
    
    ## Match columns in score_table
    Final_score = attribute_score,
    SourceFile = "quantitative-exposure-attribute-scores.csv",
    Scorer = "Calculated",
    row_idx = NA_integer_,
    
    ## Because you already recoded this field above, use "Exposure"
    ## If you want the original long name instead, change to:
    ## Attribute_type = "Qualitative Exposure Factors"
    Attribute_type = "Exposure"
  ) %>%
  select(
    Attribute_name,
    Final_score,
    SourceFile,
    Scorer,
    stock_name,
    row_idx,
    Attribute_type
  )

## Quick check
print(unique(exposure_scores_uscar$Attribute_type))
print(unique(exposure_scores_uscar$Scorer))
print(table(exposure_scores_uscar$Attribute_name))
print(nrow(exposure_scores_uscar))   ## should be 325 if all 25 stocks x 13 factors

## -----------------------------------------------------------------------------
## Append calculated quantitative exposure scores to full score table
score_table_all <- bind_rows(score_table, exposure_scores_uscar)

## Inspect
print(dim(score_table))
print(dim(exposure_scores_uscar))
print(dim(score_table_all))

head(score_table_all)
tail(score_table_all)

## Write out combined table: "table_final_attribute_exposure_uscar_combined.csv"
write.csv(
  score_table_all,
  file = file.path(out_dir, "table_final_attribute_exposure_uscar_combined.csv"),
  row.names = FALSE
)

################################################################################
##------------------------------------------------------------------------------
## Calculate mean attribute scores by type

attribute_type_scores <- score_table_all %>%
  group_by(stock_name, Attribute_type) %>%
  summarise(
    mean_score   = mean(Final_score, na.rm = TRUE),
    sd_score     = sd(Final_score, na.rm = TRUE),
    n_attributes = n(),
    .groups = "drop"
  )

print(attribute_type_scores, n = 50)

##------------------------------------------------------------------------------
## Convert attribute means to wide format
## Calculate vulnerability score as the sum of Exposure + Sensitivity
## Note: "Calculated" rows only have Exposure, so Vulnerability is left as NA
## unless both Exposure and Sensitivity are present

attribute_scores_wide <- attribute_type_scores %>%
  select(stock_name, Attribute_type, mean_score) %>%
  tidyr::pivot_wider(
    names_from  = Attribute_type,
    values_from = mean_score
  ) %>%
  rowwise() %>%
  mutate(
    Vulnerability = case_when(
      is.na(Exposure) & is.na(Sensitivity) ~ NA_real_,
      is.na(Exposure) | is.na(Sensitivity) ~ NA_real_,
      TRUE ~ Exposure * Sensitivity
    )
  ) %>%
  ungroup()

##------------------------------------------------------------------------------
## Add overall rank

attribute_scores_wide <- attribute_scores_wide %>%
  mutate(
    Overall_rank = case_when(
      is.na(Vulnerability)                                     ~ NA_character_,
      Vulnerability < low_cut                                  ~ "Low",
      Vulnerability >= low_cut & Vulnerability < moderate_cut  ~ "Moderate",
      Vulnerability >= moderate_cut & Vulnerability < high_cut ~ "High",
      Vulnerability >= high_cut                                ~ "Very High"
    )
  )

print(attribute_scores_wide, n = 25)

##------------------------------------------------------------------------------
## Look at overall ranks
## This excludes rows with NA Vulnerability, such as exposure-only "Calculated" rows

attribute_scores_wide %>%
  filter(!is.na(Overall_rank)) %>%
  count(Overall_rank, sort = TRUE)