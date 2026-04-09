data_quality_table <- 

##------------------------------------------------------------------------------
## Summaries

## By stock x attribute
data_quality_attribute_summary <- data_quality_table %>%
  group_by(stock_name, Attribute_type, row_idx, Attribute_name) %>%
  summarise(
    n_scores    = n(),
    n_adequate  = sum(Data_quality_score == 3, na.rm = TRUE),
    n_limited   = sum(Data_quality_score == 2, na.rm = TRUE),
    n_expert    = sum(Data_quality_score == 1, na.rm = TRUE),
    n_none      = sum(Data_quality_score == 0, na.rm = TRUE),
    n_ge_2      = sum(Data_quality_score >= 2, na.rm = TRUE),
    prop_ge_2   = n_ge_2 / n_scores,
    mean_score  = mean(Data_quality_score, na.rm = TRUE),
    sd_score    = sd(Data_quality_score, na.rm = TRUE),
    data_quality_rank = case_when(
      is.na(prop_ge_2)   ~ NA_character_,
      prop_ge_2 >= 0.80  ~ "High",
      prop_ge_2 >= 0.50  ~ "Moderate",
      prop_ge_2 <  0.50  ~ "Poor"
    ),
    .groups = "drop"
  ) %>%
  arrange(stock_name, row_idx); data_quality_attribute_summary

## Overall by stock
data_quality_stock_summary <- data_quality_table %>%
  group_by(stock_name) %>%
  summarise(
    n_scores    = n(),
    n_adequate  = sum(Data_quality_score == 3, na.rm = TRUE),
    n_limited   = sum(Data_quality_score == 2, na.rm = TRUE),
    n_expert    = sum(Data_quality_score == 1, na.rm = TRUE),
    n_none      = sum(Data_quality_score == 0, na.rm = TRUE),
    n_ge_2      = sum(Data_quality_score >= 2, na.rm = TRUE),
    prop_ge_2   = n_ge_2 / n_scores,
    mean_score  = mean(Data_quality_score, na.rm = TRUE),
    sd_score    = sd(Data_quality_score, na.rm = TRUE),
    data_quality_rank = case_when(
      is.na(prop_ge_2)   ~ NA_character_,
      prop_ge_2 >= 0.80  ~ "High",
      prop_ge_2 >= 0.50  ~ "Moderate",
      prop_ge_2 <  0.50  ~ "Poor"
    ),
    .groups = "drop"
  ) %>%
  arrange(desc(prop_ge_2), stock_name); data_quality_stock_summary