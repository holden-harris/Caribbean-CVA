##------------------------------------------------------------------------------
## User setup

in_dir      <- "./outputs/final-scores-compiled/"
out_dir     <- "./outputs/final-scores-compiled/"

## Breakpoints for assessing overall score
low_cut <- 3
moderate_cut <- 6
high_cut <- 9

## Load compiled CVA score table
score_table <- read.csv(file.path(in_dir, "table_final_attribute_scores_all.csv"), stringsAsFactors = FALSE)

## Update attribute names: ## All "Rigidity" scores are sensitivity
score_table$Attribute_type[score_table$Attribute_type == "Rigidity"] <- "Sensitivity" 
## Shorten name for "Exposure"
score_table$Attribute_type[score_table$Attribute_type == "Qualitative Exposure Factors"] <- "Exposure" 

################################################################################
##------------------------------------------------------------------------------
## Join Exposure scores

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
      TRUE ~ Exposure + Sensitivity
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

##------------------------------------------------------------------------------
## QA checks

attribute_scores_wide %>%
  summarise(
    n_rows = n(),
    n_scorers = n_distinct(Scorer),
    n_stocks = n_distinct(stock_name)
  )

## Rows missing either Exposure or Sensitivity
## Expected for Scorer == "Calculated" unless you also create calculated sensitivity
attribute_scores_wide %>%
  filter(
    is.na(Exposure) |
      is.na(Sensitivity)
  )

## Optional: specifically inspect calculated rows
attribute_scores_wide %>%
  filter(Scorer == "Calculated")

##------------------------------------------------------------------------------
## Aggregate again by species
## Only use rows with non-missing Vulnerability

stock_vulnerability_summary <- attribute_scores_wide %>%
  filter(!is.na(Vulnerability)) %>%
  group_by(stock_name) %>%
  summarise(
    mean_vulnerability = mean(Vulnerability, na.rm = TRUE),
    sd_vulnerability   = sd(Vulnerability, na.rm = TRUE),
    count_low          = sum(Overall_rank == "Low", na.rm = TRUE),
    count_moderate     = sum(Overall_rank == "Moderate", na.rm = TRUE),
    count_high         = sum(Overall_rank == "High", na.rm = TRUE),
    count_very_high    = sum(Overall_rank == "Very High", na.rm = TRUE),
    n_reviewers        = n_distinct(Scorer),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_vulnerability), stock_name)

##------------------------------------------------------------------------------
## Determine dominant overall rank

stock_vulnerability_summary <- stock_vulnerability_summary %>%
  rowwise() %>%
  mutate(
    max_val = max(c(count_low, count_moderate, count_high, count_very_high), na.rm = TRUE),
    overall_rank = paste(
      c(
        if (count_low == max_val) "Low",
        if (count_moderate == max_val) "Moderate",
        if (count_high == max_val) "High",
        if (count_very_high == max_val) "Very High"
      ),
      collapse = "-"
    )
  ) %>%
  ungroup() %>%
  select(-max_val)

print(stock_vulnerability_summary, n = 25)


















################################################################################
##------------------------------------------------------------------------------
## Calculate mean attribute scores by type

attribute_type_scores <- score_table_all %>%
  group_by(stock_name, Scorer, Attribute_type) %>%
  summarise(
    mean_score = mean(Final_score, na.rm = TRUE),
    sd_score = sd(Final_score, na.rm = TRUE),
    n_attributes = n(),
    .groups = "drop"
  ); print(attribute_type_scores, n = 35)

##------------------------------------------------------------------------------
## Convert attribute means to wide format
## Calculate vulnerability score as product of Exposure, Rigidity, and Sensitivity

attribute_scores_wide <- attribute_type_scores %>%
  select(stock_name, Scorer, Attribute_type, mean_score) %>%
  pivot_wider(
    names_from = Attribute_type,
    values_from = mean_score
  ) %>%
  rename(Exposure = `Qualitative Exposure Factors`) %>%
  rowwise() %>%
  mutate(
    Vulnerability = if(all(is.na(c(Exposure, Rigidity, Sensitivity)))) NA_real_
    else sum(c(Exposure, Rigidity, Sensitivity), na.rm = TRUE)
  ) %>%
  ungroup()

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

print(attribute_scores_wide, n = 30)

## Look at overall ranks
attribute_scores_wide %>%
  count(Overall_rank, sort = TRUE)

##------------------------------------------------------------------------------
## QA checks

attribute_scores_wide %>%
  summarise(
    n_rows = n(),
    n_scorers = n_distinct(Scorer),
    n_stocks = n_distinct(stock_name)
  )

attribute_scores_wide %>%
  filter(
    is.na(Sensitivity) |
      is.na(Rigidity) |
      is.na(`Qualitative Exposure Factors`)
  )

##------------------------------------------------------------------------------
## Aggregate again by species

stock_vulnerability_summary <- attribute_scores_wide %>%
  group_by(stock_name) %>%
  summarise(
    mean_vulnerability = mean(Vulnerability, na.rm = TRUE),
    sd_vulnerability   = sd(Vulnerability, na.rm = TRUE),
    count_low          = sum(Overall_rank == "Low", na.rm = TRUE),
    count_moderate     = sum(Overall_rank == "Moderate", na.rm = TRUE),
    count_high         = sum(Overall_rank == "High", na.rm = TRUE),
    count_very_high    = sum(Overall_rank == "Very High", na.rm = TRUE),
    n_reviewers        = n_distinct(Scorer),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_vulnerability), stock_name)

## Determine dominant overall rank
stock_vulnerability_summary <- stock_vulnerability_summary %>%
  rowwise() %>%
  mutate(
    max_val = max(c(count_low, count_moderate, count_high, count_very_high), na.rm = TRUE),
    
    overall_rank = paste(
      c(
        if (count_low == max_val) "Low",
        if (count_moderate == max_val) "Moderate",
        if (count_high == max_val) "High",
        if (count_very_high == max_val) "Very High"
      ),
      collapse = "-"
    )
  ) %>%
  ungroup() %>%
  select(-max_val)


print(stock_vulnerability_summary, n = 25)
