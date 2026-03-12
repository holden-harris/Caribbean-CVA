##------------------------------------------------------------------------------
## CONFIG

library(dplyr)

in_dir      <- "./data/final-scores"
out_dir     <- "./outputs/final-score-aggregates"
directional_effect_table <- read.csv(file.path(in_dir, "directional_effect_table_all.csv"))

## -----------------------------------------------------------------------------
## Aggregate to one overall directional effect per stock_name × scorer

directional_effect_summary <- directional_effect_table %>%
  dplyr::select(Scorer, stock_name, Directional_effect, Directional_score) %>%
  tidyr::pivot_wider(
    names_from  = Directional_effect,
    values_from = Directional_score
  ) %>%
  dplyr::mutate(
    Positive = dplyr::coalesce(Positive, 0),
    Neutral  = dplyr::coalesce(Neutral, 0),
    Negative = dplyr::coalesce(Negative, 0)
  ) %>%
  dplyr::mutate(
    overall_directional_effect = dplyr::case_when(
      Positive >= 2 &  Negative < 2                    ~ "positive",
      Negative >= 2 &  Positive < 2                    ~ "negative",
      Neutral  >= 3 &  Positive < 2  & Negative  < 2   ~ "neutral",
      Neutral  == 2 &  Positive == 1 & Negative == 1   ~ "neutral",
      Neutral  == 2 &  Positive == 2                   ~ "neutral"
    )
  )

directional_effect_summary

write.csv(
  directional_effect_summary,
  file = file.path(out_dir, "directional_effect_aggregated.csv"),
  row.names = FALSE
)

## Count directions per stock
stock_directional_summary <- directional_effect_summary %>%
  group_by(stock_name) %>%
  summarise(
    n_positive = sum(overall_directional_effect == "positive", na.rm = TRUE),
    n_neutral  = sum(overall_directional_effect == "neutral",  na.rm = TRUE),
    n_negative = sum(overall_directional_effect == "negative", na.rm = TRUE),
    n_reviews  = n(),
    .groups = "drop"
  )

print(stock_directional_summary)

## Determine dominent directional effect
stock_directional_summary <- stock_directional_summary %>%
  rowwise() %>%
  mutate(
    max_val = max(c(n_positive, n_neutral, n_negative)),
    
    overall_directional_effect = paste(
      c(
        if (n_positive == max_val) "positive",
        if (n_neutral  == max_val) "neutral",
        if (n_negative == max_val) "negative"
      ),
      collapse = "-"
    )
  ) %>%
  ungroup() %>%
  select(-max_val)

print(stock_directional_summary, n = 26)

write.csv(
  stock_directional_summary,
  file = file.path(out_dir, "directional_effect_summary_by-stock.csv"),
  row.names = FALSE
)
