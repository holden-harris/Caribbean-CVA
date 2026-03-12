##------------------------------------------------------------------------------
## CONFIG
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
      TRUE                                             ~ "mixed"
    )
  )

directional_effect_summary

write.csv(
  directional_effect_summary,
  file = file.path(out_dir, "directional_effect_summary.csv"),
  row.names = FALSE
)
