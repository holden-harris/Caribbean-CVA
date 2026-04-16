
##------------------------------------------------------------------------------
## `2-plot-figures.R`
##
## Figure 1 - Reviewer x stock coverage heatmap
## Figure 2 - Sensitivity tally distributions by attribute (by stock)
## Figure 3 - Directional effect summary by stock
## Figure 4 - Sensitivity attribute score distributions
## Figure 5 - Exposure tally distributions by attribute (by stock)
## Figure 6 - Exposure attribute score distributions

## Set up ----------------------------------------------------------------------

## Libraries -------------------------------------------------------------------
rm(list = ls()); gc()
windows()
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(forcats)
library(patchwork)
library(scales)

## Directories -----------------------------------------------------------------
dir_in  <- "./outputs/analyses/1-inputs"
dir_compiled <- "./outputs/final-scores-compiled/overall-vulnerability-rankings"
dir_out <- file.path("./figures")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

## Output file paths
f_fig_coverage  <- file.path(dir_out, "fig_reviewer_stock_coverage.png")
f_fig_tallies   <- file.path(dir_out, "fig_sensitivity_tally_distributions.png")
f_fig_tallies_stock   <- file.path(dir_out, "fig_sensitivity_tally_distributions_by_stock.png")
f_fig_dir       <- file.path(dir_out, "fig_directional_effect_summary.png")
f_fig_sens_box         <- file.path(dir_out, "fig_sensitivity_attribute_score_boxplot.png")
f_fig_exp_tallies_stock <- file.path(dir_out, "fig_exposure_tally_distributions_by_stock.png")
f_fig_exp_box          <- file.path(dir_out, "fig_exposure_attribute_score_boxplot.png")

## Input data
all_qualitative_tallies_long    <- read.csv(file.path(dir_in, "table_all_qualitative_tallies_long.csv"))
sensitivity_tallies_long        <- read.csv(file.path(dir_in, "table_sensitivity_tallies_long.csv"))
directional_effect_tallies_long <- read.csv(file.path(dir_in, "table_directional_effect_tallies_long.csv"))
attr_means                      <- read.csv(file.path(dir_compiled, "attribute_means_uscar.csv"))
qualitative_exposure_tallies_long <- read.csv(file.path(dir_in, "table_qualitative_exposure_tallies_long.csv"))
quantitative_exposure_scores      <- read.csv(
  "./outputs/final-scores-compiled/quantitative-exposure-attribute-scores-all.csv")

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
## Conform and row-bind exposure tally tables 
qual_exp_conform <- qualitative_exposure_tallies_long %>%
  mutate(stock_name = recode(stock_name,
    "Atlantic thread herring" = "Atlantic Herring",
    "Long-spined sea urchin"  = "Diadema",
    "Red hind"                = "Redhind",
    "Sea cucumbers"           = "Sea Cucumber",
    "Ballyhoo"                = "Ballyhoo",
    "Blue runner"             = "Blue Runner",
    "Dolphinfish"             = "Dolphinfish",
    "Gray angelfish"          = "Gray Angelfish",
    "Hogfish"                 = "Hogfish",
    "King mackerel"           = "King Mackerel",
    "Lane snapper"            = "Lane Snapper",
    "Misty grouper"           = "Misty Grouper",
    "Mutton snapper"          = "Mutton Snapper",
    "Nassau grouper"          = "Nassau Grouper",
    "Queen conch"             = "Queen Conch",
    "Queen triggerfish"       = "Queen Triggerfish",
    "Rainbow parrotfish"      = "Rainbow Parrotfish",
    "Red grouper"             = "Red Grouper",
    "Spiny lobster"           = "Spiny Lobster",
    "Stoplight parrotfish"    = "Stoplight Parrotfish",
    "White mullet"            = "White Mullet",
    "Yellowfin grouper"       = "Yellowfin Grouper",
    "Yellowtail snapper"      = "Yellowtail Snapper"
  )) %>%
  group_by(stock_name, attribute_type, attribute_name) %>%
  summarise(
    tally_L  = sum(tally_L,  na.rm = TRUE),
    tally_M  = sum(tally_M,  na.rm = TRUE),
    tally_H  = sum(tally_H,  na.rm = TRUE),
    tally_VH = sum(tally_VH, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  mutate(n_tallies = tally_L + tally_M + tally_H + tally_VH)

quant_exp_conform <- quantitative_exposure_scores %>%
  filter(spatial_extent == "U.S. Caribbean") %>%
  mutate(attribute_type = "Quantitative Exposure",
         attribute_name  = full_names,
         n_tallies = tally_L + tally_M + tally_H + tally_VH) %>%
  select(stock_name, attribute_type, attribute_name,
         tally_L, tally_M, tally_H, tally_VH, n_tallies) 

exposure_tallies_long <- bind_rows(qual_exp_conform, quant_exp_conform) %>%
  arrange(stock_name, attribute_type)

##------------------------------------------------------------------------------
## Short display labels for y-axis (used by boxplots and tally figures)

attr_short_names <- c(
  "Adult mobility"                                 = "Adult mobility",
  "Complexity in reproductive strategy"            = "Reprod. complex.",
  "Genetic diversity"                              = "Genetic divers.",
  "Habitat specificity"                            = "Habitat specif.",
  "Mobility and dispersal or early life stages"    = "Early life disp.",
  "Other stressors"                                = "Other stressors",
  "Population growth rate"                         = "Pop. growth rate",
  "Predation and competition dynamics"             = "Pred. and compet.",
  "Prey specificity"                               = "Prey specif.",
  "Spawning characteristics"                       = "Spawning charact.",
  "Species range"                                  = "Species range",
  "Specificity in early life history requirements" = "Early life req.",
  "Stock Size Status"                              = "Stock size",
  "Tolerance to ocean acidification"               = "OA tolerance"
)

exp_attr_short_names <- c(
  "Bottom salinity"                       = "Bott. sal.",
  "Bottom temperature"                    = "Bott. temp.",
  "Chlorophyll-a concentration"           = "Chl-a",
  "Mean sea surface temperature gradient" = "SST gradient",
  "Mixed layer depth"                     = "Mixed layer",
  "Oxygen at 200m"                        = "O2 at 200m",
  "Precipitation"                         = "Precip.",
  "Primary production"                    = "Primary prod.",
  "Sea surface oxygen"                    = "SS 02",
  "Sea surface salinity"                  = "SS sal.",
  "Sea surface temperature"               = "SS tempe.",
  "Surface pH"                            = "SS pH",
  "Surface wind speed magnitude"          = "Wind speed",
  "Coral cover"                           = "Coral cover",
  "Sargassum influx"                      = "Sargassum ",
  "Thermocline depth"                     = "Thermocl. depth"
)

##------------------------------------------------------------------------------
## Figure 1 - Score distributions
##   - Horizontal boxplot per attribute
##   - x = mean score (1-4) across reviewers for each stock
##   - y = sensitivity attribute, ordered by median ascending
##   - Median bar, IQR box, 1.5x IQR whiskers, outlier points

##------------------------------------------------------------------------------
## Figure 1A - Sensitivity attribute distributions

sens_scores <- attr_means %>%
  filter(attribute_type == "Sensitivity") %>%
  filter(!is.na(attribute_mean)) %>%
  group_by(attribute_name) %>%
  mutate(attr_median = median(attribute_mean, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(display_name = ifelse(attribute_name %in% names(attr_short_names),
                               attr_short_names[attribute_name], attribute_name))

attr_box_order <- sens_scores %>%
  group_by(attribute_name, display_name) %>%
  summarise(med = median(attribute_mean, na.rm = TRUE), .groups = "drop") %>%
  arrange(med) %>%
  pull(display_name)

p_sens_box <- ggplot(
  sens_scores,
  aes(
    x    = attribute_mean,
    y    = factor(display_name, levels = attr_box_order),
    fill = attr_median
  )
) +
  geom_boxplot(
    color         = "black",
    outlier.shape = 1,
    outlier.size  = 2
  ) +
  scale_fill_gradientn(
    colors = c("green3", "yellow2", "orange2", "red3"),
    values = scales::rescale(c(1, 2, 3, 4)),
    limits = c(1, 4),
    guide  = "none"
  ) +
  scale_x_continuous(
    breaks = 1:4,
    labels = c("1\nLow", "2\nModerate", "3\nHigh", "4\nVery High"),
    limits = c(0.5, 4.5)
  ) +
  labs(
    x     = NULL,
    y     = NULL,
    title = "Biological sensitivity attributes"
  ) +
  theme_bw(base_size = 11) +
  theme(
    text      = element_text(color = "black"),
    axis.text = element_text(color = "black")
  ); plot(p_sens_box)

## Write out sensitivity attribute distribution plot
ggsave(f_fig_sens_box, p_sens_box,
       width  = 7,
       height = max(4, 0.35 * n_distinct(sens_scores$attribute_name)),
       dpi    = 1200)

##------------------------------------------------------------------------------
## Figure 1B - Exposure attribute score distributions

exp_scores <- attr_means %>%
  filter(attribute_type == "Exposure") %>%
  filter(!is.na(attribute_mean)) %>%
  group_by(attribute_name) %>%
  mutate(attr_median = median(attribute_mean, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(display_name = ifelse(attribute_name %in% names(exp_attr_short_names),
                               exp_attr_short_names[attribute_name], attribute_name))

exp_box_order <- exp_scores %>%
  group_by(attribute_name, display_name) %>%
  summarise(med = median(attribute_mean, na.rm = TRUE), .groups = "drop") %>%
  arrange(med) %>%
  pull(display_name)

p_exp_box <- ggplot(
  exp_scores,
  aes(
    x    = attribute_mean,
    y    = factor(display_name, levels = exp_box_order),
    fill = attr_median
  )
) +
  geom_boxplot(
    color         = "black",
    outlier.shape = 1,
    outlier.size  = 2
  ) +
  scale_fill_gradientn(
    colors = c("green3", "yellow2", "orange2", "red3"),
    values = scales::rescale(c(1, 2, 3, 4)),
    limits = c(1, 4),
    guide  = "none"
  ) +
  scale_x_continuous(
    breaks = 1:4,
    labels = c("1\nLow", "2\nModerate", "3\nHigh", "4\nVery High"),
    limits = c(0.5, 4.5)
  ) +
  labs(
    x     = NULL,
    y     = NULL,
    title = "Exposure factors"
  ) +
  theme_bw(base_size = 11) +
  theme(
    text      = element_text(color = "black"),
    axis.text = element_text(color = "black")
  ); p_exp_box

## Write out exposure distribution plot
ggsave(f_fig_exp_box, p_exp_box,
       width  = 7,
       height = max(4, 0.35 * n_distinct(exp_scores$attribute_name)),
       dpi    = 1200)

##------------------------------------------------------------------------------
## Combined panel: sensitivity (A) + exposure (B) attribute score distributions

p_combined_box <- p_sens_box + p_exp_box +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 13, color = "black", face = "bold"))
plot(p_combined_box)

f_fig_combined_box <- file.path(dir_out, "fig_attribute_score_boxplot_combined.png")

ggsave(f_fig_combined_box, p_combined_box,
       width  = 14,
       height = max(4, 0.35 * max(n_distinct(sens_scores$attribute_name),
                                  n_distinct(exp_scores$attribute_name))),
       dpi    = 1200)

message("Figures written to: ", dir_out)



##------------------------------------------------------------------------------
## Figure 2 - Directional effect summary by stock
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
                     expand  = c(0, 0, 0, 0.03)) +
  labs(
    x        = "Proportion of reviewer tallies",
    y        = NULL,
    title    = "Directional effect",
  ) +
  theme(
    panel.background = element_rect(fill = "white"),
    text            = element_text(size = 12, color = "black"),
    axis.text       = element_text(size = 11, color = "black"),
    legend.text     = element_text(size = 11, color = "black"),
    legend.position = "bottom",
    axis.line       = element_line(color = "black")
  ); p_dir

## Write out directional effect plot -------------------------------------------
ggsave(f_fig_dir, p_dir,
       width  = 7,
       height = max(4, 0.3 * n_distinct(dir_prop$stock_name)),
       dpi    = 1200)

message("Figures written to: ", f_fig_dir)



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
## Figure 5 - Exposure tally distributions by attribute (by stock)
##   - Combined qualitative + quantitative exposure factors
##   - Faceted 5x5 by stock; horizontal stacked bar per attribute
##   - Ordered by pooled mean score ascending (lowest at bottom)

exp_pooled <- exposure_tallies_long %>%
  group_by(attribute_name) %>%
  summarise(
    tally_L    = sum(tally_L,  na.rm = TRUE),
    tally_M    = sum(tally_M,  na.rm = TRUE),
    tally_H    = sum(tally_H,  na.rm = TRUE),
    tally_VH   = sum(tally_VH, na.rm = TRUE),
    pooled_sum = tally_L + tally_M + tally_H + tally_VH,
    .groups = "drop"
  ) %>%
  filter(pooled_sum > 0) %>%
  mutate(mean_score = (tally_L * 1 + tally_M * 2 + tally_H * 3 + tally_VH * 4) / pooled_sum)

exp_attr_order <- exp_pooled %>% arrange(mean_score) %>% pull(attribute_name)

exp_attr_order_short <- ifelse(exp_attr_order %in% names(exp_attr_short_names),
                               exp_attr_short_names[exp_attr_order], exp_attr_order)

exp_pooled_stock <- exposure_tallies_long %>%
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

exp_pooled_stock_long <- exp_pooled_stock %>%
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
      ifelse(attribute_name %in% names(exp_attr_short_names),
             exp_attr_short_names[attribute_name], attribute_name),
      levels = exp_attr_order_short
    )
  )

p_exp_tallies <- ggplot(
  exp_pooled_stock_long,
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
    x = "Proportion of Tallies (Exposure Factors)",
    y = NULL
  ) +
  theme(
    panel.background  = element_rect(fill = "white"),
    strip.background  = element_rect(fill = "grey80", color = "black"),
    strip.text        = element_text(size = 9,   color = "black"),
    text              = element_text(size = 10,  color = "black"),
    axis.text.y       = element_text(size = 8.5, color = "black"),
    axis.text.x       = element_text(size = 10,  color = "black", angle = 0),
    axis.line         = element_line(color = "black"),
    axis.title.x      = element_text(size = 12,  color = "black"),
    legend.title      = element_text(size = 11,  color = "black"),
    legend.text       = element_text(size = 11,  color = "black"),
    legend.position   = "bottom",
    panel.spacing     = unit(0.2, "lines")
  ); p_exp_tallies

ggsave(f_fig_exp_tallies_stock, p_exp_tallies,
       width  = 12,
       height = 14,
       dpi    = 1200)






##------------------------------------------------------------------------------
## QA Figure - Reviewer x stock coverage heatmap
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