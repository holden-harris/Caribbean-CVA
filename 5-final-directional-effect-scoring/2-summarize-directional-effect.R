##------------------------------------------------------------------------------
## Caribbean CVA – Summarize reviewer Directional Effect scores 

##------------------------------------------------------------------------------
## CONFIG

library(dplyr)

in_dir      <- "./outputs/final-scores-compiled/"
out_dir     <- "./outputs/final-scores-compiled/"
directional_effect_table <- read.csv(file.path(in_dir, "table_directional_effect_scores.csv"))

## -----------------------------------------------------------------------------
## Aggregate to one overall directional effect per stock_name × scorer
## using CVA weighted-average method

directional_effect_summary <- directional_effect_table %>%
  dplyr::select(Scorer, stock_name, Directional_effect, Directional_score) %>%
  tidyr::pivot_wider(
    names_from  = Directional_effect,
    values_from = Directional_score
  ) %>%
  dplyr::mutate(
    Positive = dplyr::coalesce(Positive, 0),
    Neutral  = dplyr::coalesce(Neutral,  0),
    Negative = dplyr::coalesce(Negative, 0)
  ) %>%
  dplyr::mutate(
    n_tallies = Positive + Neutral + Negative,
    
    wt_avg = dplyr::if_else(
      n_tallies > 0,
      ((Negative * -1) + (Neutral * 0) + (Positive * 1)) / n_tallies,
      NA_real_
    ),
    overall = dplyr::case_when(
      is.na(wt_avg)            ~ NA_character_,
      wt_avg <= -0.333         ~ "negative",
      wt_avg >=  0.333         ~ "positive",
      wt_avg >  -0.333 &
        wt_avg <  0.333        ~ "neutral"
    )
  ); print(directional_effect_summary, n = 30)

## QA checks
directional_effect_summary %>%
  dplyr::count(overall)

directional_effect_summary %>%
  dplyr::select(
    stock_name, Scorer, Positive, Neutral, Negative,
    wt_avg, overall
  ) %>%
  dplyr::arrange(stock_name, Scorer)

## Write out pivot table
write.csv(
  directional_effect_summary,
  file = file.path(out_dir, "directional_effect_wide_all.csv"),
  row.names = FALSE
)

## -----------------------------------------------------------------------------
## Summarize by stock

## -----------------------------------------------------------------------------
## Summarize by stock using raw tally totals from all reviewers
## and HMS weighted-average directional effect method

stock_directional_summary <- directional_effect_table %>%
  group_by(stock_name, Directional_effect) %>%
  summarise(
    total_tallies = sum(Directional_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from  = Directional_effect,
    values_from = total_tallies
  ) %>%
  dplyr::mutate(
    Positive = dplyr::coalesce(Positive, 0),
    Neutral  = dplyr::coalesce(Neutral, 0),
    Negative = dplyr::coalesce(Negative, 0)
  ) %>%
  dplyr::mutate(
    n_tallies = Positive + Neutral + Negative,
    wt_avg = dplyr::if_else(
      n_tallies > 0,
      ((Negative * -1) + (Neutral * 0) + (Positive * 1)) / n_tallies,
      NA_real_
    ),
    overall = dplyr::case_when(
      is.na(wt_avg)            ~ NA_character_,
      wt_avg <= -0.333         ~ "negative",
      wt_avg >=  0.333         ~ "positive",
      wt_avg >  -0.333 &
        wt_avg <  0.333        ~ "neutral"
    )
  ) %>%
  dplyr::rename(
    n_positive = Positive,
    n_neutral  = Neutral,
    n_negative = Negative
  ) %>%
  dplyr::left_join(
    directional_effect_table %>%
      dplyr::group_by(stock_name) %>%
      dplyr::summarise(
        n_reviewers = dplyr::n_distinct(Scorer),
        .groups = "drop"
      ),
    by = "stock_name"
  ) %>%
  dplyr::select(
    stock_name,
    n_positive,
    n_neutral,
    n_negative,
    n_tallies,
    n_reviewers,
    wt_avg,
    overall
  ) %>%
  dplyr::arrange(stock_name); print(stock_directional_summary, n = 26)

## Write out species/stock summary
write.csv(
  stock_directional_summary,
  file = file.path(out_dir, "directional_effect_summary_by-stock.csv"),
  row.names = FALSE
)
