## -----------------------------------------------------------------------------
## Summarize Data Quality Scores
##
## Data quality was ranked using the proportion of scores above two. 
## High data quality was defined as 80% of data quality scores 2 or higher, 
## moderate data quality as 50–79% scores 2 or higher, 
## and poor data quality as 0–50% of scores 2 or higher

##------------------------------------------------------------------------------
## Set up
in_dir      <- "./outputs/final-scores-compiled/data-quality"
out_dir     <- "./outputs/final-scores-compiled/data-quality"

data_quality_table <- read.csv(file.path(in_dir, "table_data_quality_scores_extracted.csv"))

##------------------------------------------------------------------------------
## Summaries

## By stock x attribute
data_quality_attribute_summary <- data_quality_table %>%
  group_by(stock_name, Attribute_type, row_idx, Attribute_name) %>%
  summarise(
    n_scores      = n(),
    n_3           = sum(Data_quality_score == 3, na.rm = TRUE), 
    n_2           = sum(Data_quality_score == 2, na.rm = TRUE),
    n_1           = sum(Data_quality_score == 1, na.rm = TRUE),
    n_0           = sum(Data_quality_score == 0, na.rm = TRUE),
    n_ge_2        = sum(Data_quality_score >= 2, na.rm = TRUE),
    prop_ge_2     = n_ge_2 / n_scores,
    mean_score    = mean(Data_quality_score, na.rm = TRUE),
    sd_score      = sd(Data_quality_score, na.rm = TRUE),
    data_qual_rank = case_when(
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
#    n_scores    = n(),
    n_3           = sum(Data_quality_score == 3, na.rm = TRUE), 
    n_2           = sum(Data_quality_score == 2, na.rm = TRUE),
    n_1           = sum(Data_quality_score == 1, na.rm = TRUE),
    n_0           = sum(Data_quality_score == 0, na.rm = TRUE),
    n_ge_2        = sum(Data_quality_score >= 2, na.rm = TRUE),
    n_ge_2      = sum(Data_quality_score >= 2, na.rm = TRUE),
    prop_ge_2   = round(n_ge_2 / n(), 2),
    mean_score  = round(mean(Data_quality_score, na.rm = TRUE), 2),
    sd_score    = round(sd(Data_quality_score, na.rm = TRUE), 2),
    data_quality_rank = case_when(
      is.na(prop_ge_2)   ~ NA_character_,
      prop_ge_2 >= 0.80  ~ "High",
      prop_ge_2 >= 0.50  ~ "Moderate",
      prop_ge_2 <  0.50  ~ "Poor"
    ),
    .groups = "drop"
  ) %>%
  arrange(desc(prop_ge_2), stock_name); print(data_quality_stock_summary, n = 25)

write.csv(
  data_quality_stock_summary, 
  file = file.path(out_dir, "overall_data_quality_summary_by_stock.csv"),
  row.names = FALSE
)
