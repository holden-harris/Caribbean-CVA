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
##------------------------------------------------------------------------------
scorers <- score_table_all %>% distinct(Scorer) %>% arrange(Scorer) %>% pull(Scorer)
## Use a large qualitative palette; recycle if needed
pal <- scales::hue_pal()(max(length(scorers), 3))
scorer_colors <- tibble(Scorer = scorers, color = pal[seq_along(scorers)])
readr::write_csv(scorer_colors, file.path(out_dir, "scorer_colors.csv"))

##------------------------------------------------------------------------------
## 4) Prepare long format for plotting and HMS calculations
##   - one row per Scorer × Species × Attribute × Category (L/M/H/V)
##------------------------------------------------------------------------------
lmhv_long <- score_table_all |>
  select(Scorer, stock_name, Attribute_name, Attribute_short,
         Scoring_rank_1, Scoring_rank_2, Scoring_rank_3, Scoring_rank_4) |>
  pivot_longer(cols = starts_with("Scoring_rank_"),
               names_to = "category", values_to = "val") |>
  mutate(category = recode(category,
                           "Scoring_rank_1" = "L",
                           "Scoring_rank_2" = "M",
                           "Scoring_rank_3" = "H",
                           "Scoring_rank_4" = "V")) |>
  mutate(val = as.numeric(val)) |>
  replace_na(list(val = 0)) |>
  filter(val %in% c(0,1))  ## ensure binary tallies

##------------------------------------------------------------------------------
## 5) Compute per-species × attribute HMS mean and SD
##   - L,M,H,V = sums across scorers
##   - num_experts = number of scorers who picked any category for that attribute
##   - Weighted mean (HMS): ((L*1 + M*2 + H*3 + V*4) / (num_experts * 5))
##   - SD: sd(rep(1,L), rep(2,M), rep(3,H), rep(4,V))
##------------------------------------------------------------------------------
lmhv_sums <- lmhv_long %>%
  group_by(stock_name, Attribute_name, Attribute_short, category) %>%
  summarise(n = sum(val, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = category, values_from = n, values_fill = 0) %>%
  left_join(
    lmhv_long %>%
      group_by(stock_name, Attribute_name) %>%
      summarise(num_experts = n_distinct(Scorer[val == 1]),
                .groups = "drop"),
    by = c("stock_name","Attribute_name")
  ) %>%
  mutate(
    HMS_mean = ifelse(num_experts > 0,
                      ((L*1) + (M*2) + (H*3) + (V*4)) / (num_experts * 5),
                      NA_real_),
    ## build vector for SD
    sd_vec = pmap(list(L,M,H,V), ~c(rep(1, ..1), rep(2, ..2), rep(3, ..3), rep(4, ..4))),
    HMS_sd  = sapply(sd_vec, function(x) if (length(x) >= 2) sd(x) else NA_real_)
  ) %>%
  select(-sd_vec)

##------------------------------------------------------------------------------
## 6) Plot function (internal) to make per-species panel PDF
##   - stacked bars L/M/H/V by scorer (colors), one facet per attribute
##   - annotate facet labels with m= and sd= (rounded to 2)
##   - legend hidden; separate scorer_colors.csv is saved already
##------------------------------------------------------------------------------
plot_species <- function(species_name) {
  df_counts <- lmhv_long %>%
    filter(stock_name == species_name) %>%
    mutate(category = factor(category, levels = c("L","M","H","V")))
  
  if (nrow(df_counts) == 0) return(invisible(NULL))
  
  df_annot <- lmhv_sums %>%
    filter(stock_name == species_name) %>%
    mutate(m_txt  = sprintf("m= %.1f", round(HMS_mean, 1)),
           sd_txt = sprintf("sd= %.2f", round(HMS_sd, 2))) %>%
    mutate(lab = paste0(m_txt, "  ", sd_txt)) %>%
    select(stock_name, Attribute_name, Attribute_short, HMS_mean, HMS_sd, lab)
  
  ## Merge short labels into counts
  df_counts <- df_counts %>%
    left_join(df_annot %>% select(Attribute_name, Attribute_short, lab), by = "Attribute_name") %>%
    mutate(Attribute_short = fct_inorder(Attribute_short))
  
  ## Set colors
  cols <- scorer_colors$color
  names(cols) <- scorer_colors$Scorer
  
  p <- ggplot(df_counts, aes(x = category, y = val, fill = Scorer)) +
    geom_bar(stat = "count", position = "stack", width = 0.85) +
    facet_wrap(~ Attribute_short, ncol = 3, scales = "fixed") +
    labs(x = NULL, y = NULL) +
    scale_x_discrete(limits = c("L","M","H","V")) +
    scale_fill_manual(values = cols, guide = "none") +   ## hide legend
    theme_bw(base_size = 11) +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9),
      plot.title = element_text(size = 14, face = "bold")
    ) +
    ggtitle(species_name)
  
  ## Add m / sd text as facet labels (below strips) — simplest as an overlay
  ## Compute label positions per facet:
  lab_pos <- df_counts %>%
    group_by(Attribute_short) %>%
    summarise(ymax = max(..count.., na.rm = TRUE), .groups = "drop")  ## placeholder; will adjust
  
  ## Quick save (without extra label overlays). For brevity & robustness:
  out_pdf <- file.path(out_dir, "plots", paste0(gsub("[^A-Za-z0-9]+","_", species_name), ".pdf"))
  ggsave(out_pdf, p, width = 8.5, height = 11, units = "in")
  message("Wrote plot: ", out_pdf)
}

##------------------------------------------------------------------------------
## 7) Loop species and save PDFs
##------------------------------------------------------------------------------
species_vec <- lmhv_long %>% distinct(stock_name) %>% arrange(stock_name) %>% pull(stock_name)

walk(species_vec, plot_species)

##------------------------------------------------------------------------------
## 8) Export per-species preliminary SD (HMS ‘meansd’ = mean of attribute SDs)
##------------------------------------------------------------------------------
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
