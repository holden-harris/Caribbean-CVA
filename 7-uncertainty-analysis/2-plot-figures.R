
##------------------------------------------------------------------------------
## `2-plot-figures.R`
##
## Figure 1 - Reviewer x stock coverage heatmap
## Figure 2 - Sensitivity tally distributions by attribute (by stock)
## Figure 3 - Directional effect summary by stock
## Figure 4 - Sensitivity attribute score distributions

## Set up ----------------------------------------------------------------------

## Libraries -------------------------------------------------------------------
rm(list = ls()); gc()
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(forcats)
library(scales)

## Directories -----------------------------------------------------------------
dir_in  <- "./outputs/analyses/1-inputs"
dir_out <- file.path(dir_in, "figures")
dir_compiled <- "./outputs/final-scores-compiled/overall-vulnerability-rankings"

dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)


## Output file paths
f_fig_coverage  <- file.path(dir_out, "fig_reviewer_stock_coverage.png")
f_fig_tallies   <- file.path(dir_out, "fig_sensitivity_tally_distributions.png")
f_fig_tallies_stock   <- file.path(dir_out, "fig_sensitivity_tally_distributions_by_stock.png")
f_fig_dir       <- file.path(dir_out, "fig_directional_effect_summary.png")
f_fig_sens_box  <- file.path(dir_out, "fig_sensitivity_attribute_score_boxplot.png")

## Input data
all_qualitative_tallies_long    <- read.csv(file.path(dir_in, "table_all_qualitative_tallies_long.csv"))
sensitivity_tallies_long        <- read.csv(file.path(dir_in, "table_sensitivity_tallies_long.csv"))
directional_effect_tallies_long <- read.csv(file.path(dir_in, "table_directional_effect_tallies_long.csv"))
attr_means                      <- read.csv(file.path(dir_compiled, "attribute_means_uscar.csv"))

##------------------------------------------------------------------------------
## Data prep

## Reviewer x stock attribute counts (used by Figure 1)
reviewer_stock_coverage <- all_qualitative_tallies_long %>%
  group_by(reviewer_id, stock_name) %>%
  summarise(n_attributes_scored = n_distinct(attribute_name), .groups = "drop")

n_attributes_expected <- n_distinct(all_qualitative_tallies_long$attribute_name)

## Directional effect totals per stock (used by Figure 3)
qa_dir_summary <- directional_effect_tallies_long %>%
  group_by(stock_name, effect_category) %>%
  summarise(total_tally = sum(tally, na.rm = TRUE), .groups = "drop")

##------------------------------------------------------------------------------
## Figure 1 - Reviewer x stock coverage heatmap
##   - Each tile = number of attributes scored by that reviewer for that stock
##   - Red = missing or incomplete; blue = fully scored
##   - Quickly shows which reviewer x stock combinations have gaps

p_coverage <- ggplot(
  reviewer_stock_coverage,
  aes(
    x    = fct_rev(fct_inorder(reviewer_id)),
    y    = fct_rev(fct_inorder(stock_name)),
    fill = n_attributes_scored
  )
) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = n_attributes_scored), size = 2.8, color = "white") +
  scale_fill_gradient(
    low  = "#d9534f",
    high = "#2c7bb6",
    name = "Attributes\nscored"
  ) +
  labs(
    x        = "Reviewer",
    y        = "Stock",
    title    = "Reviewer x stock coverage",
    subtitle = paste0("Cell value = attributes scored (max = ",
                      n_attributes_expected, ")")
  ) +
  theme_bw(base_size = 10) +
  theme(
    text        = element_text(color = "black"),
    axis.text   = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid  = element_blank()
  )

ggsave(f_fig_coverage, p_coverage,
       width  = 10,
       height = max(5, 0.3 * n_distinct(reviewer_stock_coverage$stock_name)),
       dpi    = 300)

##------------------------------------------------------------------------------
## Figure 2 - Sensitivity tally distributions by attribute 
##   - Pooled tallies across all reviewers and stocks
##   - Horizontal stacked bar: proportion Low / Moderate / High / Very High
##   - Ordered by pooled mean score ascending (lowest at bottom)

sens_pooled <- sensitivity_tallies_long %>%
  group_by(attribute_name) %>%
  summarise(
    tally_L       = sum(tally_L,       na.rm = TRUE),
    tally_M  = sum(tally_M,  na.rm = TRUE),
    tally_H      = sum(tally_H,      na.rm = TRUE),
    tally_VH = sum(tally_VH, na.rm = TRUE),
    pooled_sum      = tally_L + tally_M + tally_H + tally_VH,
    .groups = "drop"
  ) %>%
  filter(pooled_sum > 0) %>%
  mutate(
    p_low       = tally_L       / pooled_sum,
    p_moderate  = tally_M  / pooled_sum,
    p_high      = tally_H      / pooled_sum,
    p_very_high = tally_VH / pooled_sum,
    mean_score  = (tally_L * 1 + tally_M * 2 +
                     tally_H * 3 + tally_VH * 4) / pooled_sum
  )

rank_levels <- c("Low", "Moderate", "High", "Very High")
rank_colors <- c(
  "Low"       = "green3",
  "Moderate"  = "yellow2",
  "High"      = "orange2",
  "Very High" = "red3"
)

## attr_order: attribute ordering by overall mean score (used across all Figure 2 plots)
attr_order <- sens_pooled %>%
  arrange(mean_score) %>%
  pull(attribute_name)

## Short display labels for y-axis (≤20 chars, single line)
attr_short_names <- c(
  "Adult mobility"                                 = "Adult mobility",
  "Complexity in reproductive strategy"            = "Reprod. complexity",
  "Genetic diversity"                              = "Genetic diversity",
  "Habitat specificity"                            = "Habitat specificity",
  "Mobility and dispersal or early life stages"    = "Early life dispersal",
  "Other stressors"                                = "Other stressors",
  "Population growth rate"                         = "Pop. growth rate",
  "Predation and competition dynamics"             = "Pred./competition",
  "Prey specificity"                               = "Prey specificity",
  "Spawning characteristics"                       = "Spawning charact.",
  "Species range"                                  = "Species range",
  "Specificity in early life history requirements" = "Early life req.",
  "Stock Size Status"                              = "Stock size/status",
  "Tolerance to ocean acidification"               = "OA tolerance"
)

attr_order_short <- ifelse(attr_order %in% names(attr_short_names),
                           attr_short_names[attr_order], attr_order)

## Per-stock tallies (one row per stock x attribute x rank)
sens_pooled_stock <- sensitivity_tallies_long %>%
  group_by(stock_name, attribute_name) %>%
  summarise(
    tally_L    = sum(tally_L,  na.rm = TRUE),
    tally_M    = sum(tally_M,  na.rm = TRUE),
    tally_H    = sum(tally_H,  na.rm = TRUE),
    tally_VH   = sum(tally_VH, na.rm = TRUE),
    pooled_sum = tally_L + tally_M + tally_H + tally_VH,
    .groups = "drop"
  ) %>%
  filter(pooled_sum > 0) %>%
  mutate(
    p_low       = tally_L  / pooled_sum,
    p_moderate  = tally_M  / pooled_sum,
    p_high      = tally_H  / pooled_sum,
    p_very_high = tally_VH / pooled_sum
  )

sens_pooled_stock_long <- sens_pooled_stock %>%
  select(stock_name, attribute_name,
         p_low, p_moderate, p_high, p_very_high) %>%
  pivot_longer(
    cols      = starts_with("p_"),
    names_to  = "rank",
    values_to = "proportion"
  ) %>%
  mutate(
    rank = case_when(
      rank == "p_low"       ~ "Low",
      rank == "p_moderate"  ~ "Moderate",
      rank == "p_high"      ~ "High",
      rank == "p_very_high" ~ "Very High"
    ),
    rank = factor(rank, levels = rank_levels),
    attribute_name_short = factor(
      ifelse(attribute_name %in% names(attr_short_names),
             attr_short_names[attribute_name], attribute_name),
      levels = attr_order_short
    )
  )

p_tallies <- ggplot(
  sens_pooled_stock_long,
  aes(
    x    = proportion,
    y    = attribute_name_short,
    fill = rank
  )
) +
  geom_col(width = 0.77, position = position_stack(reverse = TRUE),
           color = "black", linewidth = 0.2) +
  scale_fill_manual(values = rank_colors, name = "Vulnerability rank:",
                    breaks = rank_levels) +
  scale_x_continuous(breaks = c(0.25, 0.50, 0.75, 1.00),
                     labels = percent_format(accuracy = 1),
                     expand  = c(0, 0, 0, 0.02)) +
  facet_wrap(~ stock_name, ncol = 5, axes = "margins") +
  labs(
    x     = "Proportion of Tallies (Sensitivity Attributes)",
    y     = NULL,
#    title = "Sensitivity attribute tally distributions by stock"
  ) +
  theme(
    panel.background = element_rect(fill = "white"),
    strip.background = element_rect(fill = "grey80", color = "black"),
    strip.text       = element_text(size = 9,  color = "black"),
    text             = element_text(size = 10,  color = "black"),
    axis.text.y        = element_text(size = 8.5,  color = "black"),
    axis.text.x        = element_text(size = 10,  color = "black", angle = 0),
    axis.line        = element_line(color = "black"),
    axis.title.x      = element_text(size = 12,  color = "black"),
    legend.title      = element_text(size = 11,  color = "black"),
    legend.text      = element_text(size = 11,  color = "black"),
    legend.position  = "bottom",
    panel.spacing    = unit(0.2, "lines")
  ); p_tallies

ggsave(f_fig_tallies_stock, p_tallies,
       width  = 12,
       height = 12,
       dpi    = 1200)

##------------------------------------------------------------------------------
## Figure 3 - Directional effect summary by stock
##   - Horizontal stacked bar: proportion positive / neutral / negative
##   - Ordered by proportion positive ascending

dir_prop <- qa_dir_summary %>%
  filter(!is.na(effect_category), !is.na(stock_name)) %>%
  group_by(stock_name) %>%
  mutate(prop = total_tally / sum(total_tally)) %>%
  ungroup() %>%
  mutate(
    effect_category = factor(effect_category,
                             levels = c("Positive", "Neutral", "Negative"))
  )
print(dir_prop, n = 75)

dir_colors <- c(
  "Positive" = "#2c7bb6",
  "Neutral"  = "bisque",
  "Negative" = "indianred4"
)

stock_order_dir <- dir_prop %>%
  filter(effect_category == "Negative") %>%
  arrange(prop) %>%
  pull(stock_name)

p_dir <- ggplot(
  dir_prop,
  aes(
    x    = prop,
    y    = factor(stock_name, levels = stock_order_dir),  
    fill = effect_category
  )
) +
  geom_col(width = 0.75, col = 'black') +
  scale_fill_manual(values = dir_colors, name = "Directional effect: ",
                    breaks = c("Negative", "Neutral", "Positive")) +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     expand  = c(0, 0, 0, 0.02)) +
  labs(
    x        = "Proportion of reviewer tallies: positive / neutral / negative",
    y        = NULL,
    title    = "Directional effect by stock",
  ) +
  theme(
    text            = element_text(size = 12, color = "black"),
    axis.text       = element_text(size = 11, color = "black"),
    legend.text     = element_text(size = 11, color = "black"),
    legend.position = "bottom",
    axis.line       = element_line(color = "black")
  ); p_dir

ggsave(f_fig_dir, p_dir,
       width  = 7,
       height = max(4, 0.3 * n_distinct(dir_prop$stock_name)),
       dpi    = 1200)

##------------------------------------------------------------------------------
## Figure 4 - Sensitivity attribute score distributions
##   - Horizontal boxplot per attribute
##   - x = mean score (1-4) across reviewers for each stock
##   - y = sensitivity attribute, ordered by median ascending
##   - Median bar, IQR box, 1.5x IQR whiskers, outlier points

sens_scores <- attr_means %>%
  filter(attribute_type == "Sensitivity") %>%
  filter(!is.na(attribute_mean))

attr_box_order <- sens_scores %>%
  group_by(attribute_name) %>%
  summarise(med = median(attribute_mean, na.rm = TRUE), .groups = "drop") %>%
  arrange(med) %>%
  pull(attribute_name)

p_sens_box <- ggplot(
  sens_scores,
  aes(
    x = attribute_mean,
    y = factor(attribute_name, levels = attr_box_order)
  )
) +
  geom_boxplot(
    fill          = "gray80",
    color         = "black",
    outlier.shape = 1,
    outlier.size  = 2
  ) +
  scale_x_continuous(
    breaks = 1:4,
    labels = c("1\nLow", "2\nModerate", "3\nHigh", "4\nVery High"),
    limits = c(0.5, 4.5)
  ) +
  labs(
    x     = "Score",
    y     = NULL,
    title = "Sensitivity attribute score distributions"
  ) +
  theme_bw(base_size = 11) +
  theme(
    text      = element_text(color = "black"),
    axis.text = element_text(color = "black")
  ); plot(p_sens_box)

ggsave(f_fig_sens_box, p_sens_box,
       width  = 7,
       height = max(4, 0.35 * n_distinct(sens_scores$attribute_name)),
       dpi    = 1200)
message("Figures written to: ", dir_out)

