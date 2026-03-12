##------------------------------------------------------------------------------
## User setup

## Breakpoints for assessing overall score
low_cut <- 3
moderate_cut <- 6
high_cut <- 9

##------------------------------------------------------------------------------
## Load compiled CVA score table

score_table <- read.csv(
  "./data/final-scores/final_score_table_all.csv",
  stringsAsFactors = FALSE
)

##------------------------------------------------------------------------------
## Calculate mean attribute scores by type

attribute_type_scores <- score_table %>%
  group_by(stock_name, Scorer, Attribute_type) %>%
  summarise(
    mean_score = mean(Final_score, na.rm = TRUE),
    sd_score = sd(Final_score, na.rm = TRUE),
    n_attributes = n(),
    .groups = "drop"
  )
print(attribute_type_scores)

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
