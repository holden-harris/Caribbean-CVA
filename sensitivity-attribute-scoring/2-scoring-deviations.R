##------------------------------------------------------------------------------
## Caribbean CVA – Preliminary Sensitivity Attribute SD plots (HMS-style)
## Inputs:
##   - score_table_all: tidy table built from the 16 reviewer workbooks with cols:
##       SourceFile, Scorer, stock_name, row_idx,
##       Attribute_name, Attribute_type,
##       Data_quality, Scoring_rank_1..4
## Output:
##   - ./outputs/prework/plots/<species>.pdf     (stacked LMHV bars by attribute)
##   - ./outputs/prework/scorer_colors.csv       (Scorer → color key; legend hidden in plots)
##   - ./outputs/prework/preliminary_sd.csv      (per-species mean SD across attributes)

rm(list = ls()); gc()
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(forcats)
  library(purrr)
})

##------------------------------------------------------------------------------
## 1) Read the compiled scoring table

in_csv  <- "./data/preliminary-scores/score_table_all.csv"   
out_dir <- "./outputs/prework"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE, recursive = TRUE)

score_table_all <- readr::read_csv(in_csv, show_col_types = FALSE)

## Sanity
stopifnot(all(c("Scorer","stock_name","Attribute_name",
                "Scoring_rank_1","Scoring_rank_2","Scoring_rank_3","Scoring_rank_4") %in% names(score_table_all)))

##------------------------------------------------------------------------------
## 2) Shorten/standardize attribute names (match HMS figure labels)

short_map <- c(
  "Adult mobility"                                   = "Adult Mobility",
  "Habitat specificity"                              = "Habitat Specificity",
  "Mobility and dispersal of early life stages"      = "Mobility_Dispersal of ELS",
  "Predation and competition dynamics"               = "Pred & Comp Dynamics",
  "Prey specificity"                                 = "Prey Specificity",
  "Complexity in reproductive strategy"              = "Repro Strat Sensitivity",
  "Spawning characteristics"                         = "Reproductive Cycle",
  "Specificity in early life history requirements"   = "Specificity EL Hist REQs",
  "Stock size/status"                                = "Stock Size_Status",
  "Tolerance to ocean acidification"                 = "Sensitivity to OA",
  "Population growth rate"                           = "Pop Growth Rate",
  "Species range"                                    = "Species Range",
  "Other stressors"                                  = "Other Stressors",
  "Genetic diversity"                                = "Genetic Diversity"
)

score_table_all <- score_table_all |>
  mutate(Attribute_name = str_squish(Attribute_name)) |>
  mutate(Attribute_short = dplyr::recode(Attribute_name, !!!short_map, .default = Attribute_name))

## Optional: enforce attribute display order (edit as desired)
att_order <- c(
  "Adult Mobility","Habitat Specificity","Mobility_Dispersal of ELS",
  "Other Stressors","Pop Growth Rate","Prey Specificity","Reproductive Cycle",
  "Repro Strat Sensitivity","Sensitivity to OA","Sensitivity to Temp",
  "Site Fidelity","Specificity EL Hist REQs","Stock Size_Status",
  "Species Range","Genetic Diversity","Pred & Comp Dynamics"
)
score_table_all <- score_table_all |>
  mutate(Attribute_short = factor(Attribute_short, levels = att_order))

##------------------------------------------------------------------------------
## 3) Build Scorer → color key (saved to CSV). Legend stays hidden in plots.

scorers <- score_table_all %>% distinct(Scorer) %>% arrange(Scorer) %>% pull(Scorer)

## Use a large qualitative palette; recycle if needed
pal <- scales::hue_pal()(max(length(scorers), 3))
scorer_colors <- tibble(Scorer = scorers, color = pal[seq_along(scorers)])
readr::write_csv(scorer_colors, file.path(out_dir, "scorer_colors.csv"))

##------------------------------------------------------------------------------
## 4) Prepare long format for plotting and HMS calculations
##   - one row per Scorer × Species × Attribute × Category (L/M/H/V)

## Tallies long-form: one row per Scorer × Species × Attribute × category, val ∈ [0..5]
lmhv_long <- score_table_all %>%
  dplyr::select(Scorer, stock_name, Attribute_name, Attribute_short,
                Scoring_rank_1, Scoring_rank_2, Scoring_rank_3, Scoring_rank_4) %>%
  tidyr::pivot_longer(cols = dplyr::starts_with("Scoring_rank_"),
                      names_to = "category", values_to = "val") %>%
  dplyr::mutate(
    category = dplyr::recode(category,
                             "Scoring_rank_1" = "L",
                             "Scoring_rank_2" = "M",
                             "Scoring_rank_3" = "H",
                             "Scoring_rank_4" = "V"),
    val = suppressWarnings(as.numeric(val))
  ) %>%
  tidyr::replace_na(list(val = 0)) %>%
  dplyr::mutate(category = factor(category, levels = c("L","M","H","V")))


##------------------------------------------------------------------------------
## 5) Compute per-species × attribute HMS mean and SD
##   - L,M,H,V = sums across scorers
##   - num_experts = number of scorers who picked any category for that attribute
##   - Weighted mean (HMS): ((L*1 + M*2 + H*3 + V*4) / (num_experts * 5))
##   - SD: sd(rep(1,L), rep(2,M), rep(3,H), rep(4,V))

## num_experts = distinct scorers who put ANY tallies (>0) on that attribute
experts_per_attr <- lmhv_long %>%
  dplyr::group_by(stock_name, Attribute_name, Scorer) %>%
  dplyr::summarise(tot = sum(val, na.rm = TRUE), .groups = "drop") %>%
  dplyr::group_by(stock_name, Attribute_name) %>%
  dplyr::summarise(num_experts = dplyr::n_distinct(Scorer[tot > 0]), .groups = "drop")

## Sum tallies across scorers by category
lmhv_sums <- lmhv_long %>%
  dplyr::group_by(stock_name, Attribute_name, Attribute_short, category) %>%
  dplyr::summarise(n = sum(val, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = category, values_from = n, values_fill = 0) %>%
  dplyr::left_join(experts_per_attr, by = c("stock_name","Attribute_name")) %>%
  dplyr::mutate(
    HMS_mean = dplyr::if_else(num_experts > 0,
                              ((L*1) + (M*2) + (H*3) + (V*4)) / (num_experts * 5),
                              NA_real_)
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    HMS_sd = {
      vec <- c(rep(1, L), rep(2, M), rep(3, H), rep(4, V))
      if (length(vec) >= 2) stats::sd(vec) else NA_real_
    }
  ) %>%
  dplyr::ungroup()


##------------------------------------------------------------------------------
## 6) Plot function (internal) to make per-species panel 
##   - stacked bars L/M/H/V by scorer (colors), one facet per attribute
##   - annotate facet labels with m= and sd= (rounded to 2)
##   - legend hidden; separate scorer_colors.csv is saved already

make_species_plot <- function(species_name) {
  
  ## counts only for selected tallies
  df_counts <- lmhv_long %>%
    filter(stock_name == stock_name, val == 1) %>%
    mutate(category = factor(category, levels = c("L","M","H","V")))
  
  if (nrow(df_counts) == 0) {
    message("No data for: ", species_name)
    return(NULL)
  }
  
  ## m/sd annotations per attribute
  df_annot <- lmhv_sums %>%
    filter(stock_name == stock_name) %>%
    mutate(m_txt  = sprintf("m= %.1f", round(HMS_mean, 1)),
           sd_txt = sprintf("sd= %.2f", round(HMS_sd, 2)),
           lab    = paste0(m_txt, "  ", sd_txt)) %>%
    select(Attribute_name, Attribute_short, lab)
  
  df_counts <- df_counts %>%
    left_join(df_annot, by = "Attribute_name") %>%
    mutate(Attribute_short = fct_inorder(Attribute_short))
  
  ## y pos for labels inside facets (max stacked + margin)
  lab_pos <- df_counts %>%
    count(Attribute_short, category, name = "n") %>%
    group_by(Attribute_short) %>%
    summarise(ymax = max(n, na.rm = TRUE), .groups = "drop") %>%
    left_join(distinct(df_counts, Attribute_short, lab), by = "Attribute_short") %>%
    mutate(x = 2.5, y = ymax + 0.5)
  
  ## colors per scorer
  cols <- scorer_colors$color
  names(cols) <- scorer_colors$Scorer
  
  ggplot(df_counts, aes(x = category, fill = Scorer)) +
    geom_bar(width = 0.85, position = "stack") +
    facet_wrap(~ Attribute_short, ncol = 3) +
    labs(title = species_name, x = NULL, y = NULL) +
    scale_x_discrete(limits = c("L","M","H","V")) +
    scale_fill_manual(values = cols, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    theme_bw(base_size = 11) +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9),
      plot.title  = element_text(size = 14, face = "bold")
    ) +
    geom_text(data = lab_pos, aes(x = x, y = y, label = lab),
              inherit.aes = FALSE, size = 3)
}

## Test 
first_species <- lmhv_long %>%
  distinct(stock_name) %>%
  arrange(stock_name) %>%
  slice(1) %>%
  pull(stock_name)

message("First species: ", first_species)
make_species_plot(first_species)




## ------------------------------------------------------------
## (2) Compile plots into a (multi-page) PDF

write_species_plots_pdf <- function(species_vec, outfile, width = 8.5, height = 11) {
  dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
  grDevices::pdf(outfile, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  
  for (sp in species_vec) {
    p <- make_species_plot(sp)
    if (!is.null(p)) print(p)
  }
  message("Wrote PDF: ", outfile)
}


























##------------------------------------------------------------------------------
## 8) Export per-species preliminary SD (note HMS ‘meansd’ = mean of attribute SDs)
prelim_sd <- lmhv_sums %>%
  group_by(stock_name) %>%
  summarise(
    prelim_sd = round(mean(HMS_sd, na.rm = TRUE), 2),
    n_attributes = sum(!is.na(HMS_sd)),
    .groups = "drop"
  ) %>%
  arrange(stock_name)

readr::write_csv(prelim_sd, file.path(out_dir, "preliminary_sd.csv"))
message("Wrote: ", file.path(out_dir, "preliminary_sd.csv"))
message("Wrote: ", file.path(out_dir, "scorer_colors.csv"))
