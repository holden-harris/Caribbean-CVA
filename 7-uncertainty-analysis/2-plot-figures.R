
##------------------------------------------------------------------------------
## Step 6 - Figures
##
## Figure 1. Reviewer x stock coverage heatmap
##   - Each tile = number of attributes scored by that reviewer for that stock
##   - Red = missing or incomplete; blue = fully scored
##   - Quickly shows which reviewer x stock combinations have gaps
##
## Figure 2. Sensitivity tally distributions by attribute
##   - Pooled tallies across all reviewers and stocks
##   - Horizontal stacked bar: proportion Low / Moderate / High / Very High
##   - Ordered by pooled mean score ascending (lowest at bottom)
##
## Figure 3. Directional effect summary by stock
##   - Horizontal stacked bar: proportion positive / neutral / negative
##   - Ordered by proportion positive ascending

##------------------------------------------------------------------------------
## Figure 1 - Reviewer x stock coverage heatmap

n_attributes_expected <- n_distinct(all_qualitative_tallies_long$attribute_name)

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
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid  = element_blank()
  )

ggsave(f_fig_coverage, p_coverage,
       width  = 10,
       height = max(5, 0.3 * n_distinct(reviewer_stock_coverage$stock_name)),
       dpi    = 300)

##------------------------------------------------------------------------------
## Figure 2 - Sensitivity tally distributions by attribute

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
  "Low"       = "#2c7bb6",
  "Moderate"  = "#abd9e9",
  "High"      = "#fdae61",
  "Very High" = "#d7191c"
)

sens_pooled_long <- sens_pooled %>%
  select(attribute_name, mean_score,
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
    rank = factor(rank, levels = rank_levels)
  )

attr_order <- sens_pooled %>%
  arrange(mean_score) %>%
  pull(attribute_name)

p_tallies <- ggplot(
  sens_pooled_long,
  aes(
    x    = proportion,
    y    = factor(attribute_name, levels = attr_order),
    fill = rank
  )
) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = rank_colors, name = "Rank") +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     expand  = c(0, 0)) +
  labs(
    x        = "Proportion of tallies",
    y        = NULL,
    title    = "Sensitivity attribute tally distributions",
    subtitle = "Pooled across all stocks and reviewers; ordered by mean score"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")

ggsave(f_fig_tallies, p_tallies,
       width  = 8,
       height = max(4, 0.35 * n_distinct(sensitivity_tallies_long$attribute_name)),
       dpi    = 300)

##------------------------------------------------------------------------------
## Figure 3 - Directional effect summary by stock

dir_prop <- qa_dir_summary %>%
  group_by(stock_name) %>%
  mutate(prop = total_tally / sum(total_tally)) %>%
  ungroup() %>%
  mutate(
    effect_category = factor(effect_category,
                             levels = c("Positive", "Neutral", "Negative"))
  )

dir_colors <- c(
  "Positive" = "#2c7bb6",
  "Neutral"  = "#ffffbf",
  "Negative" = "#d7191c"
)

stock_order_dir <- dir_prop %>%
  filter(effect_category == "Positive") %>%
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
  geom_col(width = 0.75) +
  scale_fill_manual(values = dir_colors, name = "Effect") +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     expand  = c(0, 0)) +
  labs(
    x        = "Proportion of tallies",
    y        = NULL,
    title    = "Directional effect by stock",
    subtitle = "Proportion of reviewer tallies: positive / neutral / negative"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")

ggsave(f_fig_dir, p_dir,
       width  = 7,
       height = max(4, 0.3 * n_distinct(dir_prop$stock_name)),
       dpi    = 300)

message("QA tables written to:  ", qa_dir)
message("Figures written to:    ", fig_dir)

