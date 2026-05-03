################################################################################
##------------------------------------------------------------------------------
## Caribbean CVA — Uncertainty analysis figures
##
## Figure 1 — Leave-one-out influence analysis bar plots
##   Panel A — Number of vulnerability rank changes when each EXPOSURE FACTOR
##             is omitted (tallied across all 25 stocks). Bars ordered by
##             influence (most influential at top).
##   Panel B — Same metric for each SENSITIVITY ATTRIBUTE.
##
## Figure 2 — Bootstrap resampling uncertainty
##   Panel A — Proportion of 10,000 bootstrap iterations falling in each
##             VULNERABILITY RANK (Low / Moderate / High / Very High) for
##             every stock, shown as a horizontal stacked bar.
##   Panel B — Proportion of iterations in each DIRECTIONAL EFFECT category
##             (Negative / Neutral / Positive) for every stock.
##
##   Stocks in Figure 2 are ordered by baseline vulnerability rank (Very High
##   at top, Moderate at bottom). Within each rank group, stocks are sorted by
##   their dominant bootstrap proportion ascending, so the most uncertain
##   (borderline) stocks cluster at each rank-group boundary — exactly where
##   the scientific uncertainty is most relevant to the reader.
##
## Outputs:
##   figures/fig_loo_bar_plots.png          (6.5 × 9 in,  300 dpi)
##   figures/fig_bootstrap_uncertainty.png  (7.5 × 12 in, 300 dpi)
##
## Run from the Caribbean-CVA RStudio project root (.Rproj file location).
##------------------------------------------------------------------------------

rm(list = ls()); gc()

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(forcats)
library(patchwork)

##------------------------------------------------------------------------------
## Directories

proj_dir <- "."
dir_out  <- file.path(proj_dir, "figures")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

## All analysis outputs live under this subdirectory
loo_dir <- file.path(proj_dir, "outputs", "analyses", "uncertainty-loo",
                     "final-tables")

##------------------------------------------------------------------------------
## Input file paths — Figure 1 (LOO)

f_loo_sens_sum <- file.path(loo_dir, "table_leave_one_out_sensitivity_summary.csv")
f_loo_exp_sum  <- file.path(loo_dir, "table_leave_one_out_exposure_summary.csv")

##------------------------------------------------------------------------------
## Input file paths — Figure 2 (Bootstrap)

f_boot_stock   <- file.path(loo_dir, "table_bootstrap_uncertainty_stock.csv")
f_dir_boot     <- file.path(loo_dir, "table_directional_effect_bootstrap.csv")

## Baseline vulnerability scores — used to establish the stock ordering for
## Figure 2 (primary sort key = baseline Vuln_rank)
f_vuln <- file.path(proj_dir, "outputs", "final-scores-compiled",
                    "overall-vulnerability-rankings",
                    "overall_vulnerability_scores_uscar.csv")

##------------------------------------------------------------------------------
## Output file paths

f_fig1 <- file.path(dir_out, "fig_loo_bar_plots.png")
f_fig2 <- file.path(dir_out, "fig_bootstrap_uncertainty.png")

##------------------------------------------------------------------------------
## Short display name lookup tables  (used only by Figure 1)
##
## Copied from 7-scoring-distributions/3-plot-figures.R (lines 104-137).
## Keys   = full attribute/factor names as stored in the LOO CSV files.
## Values = abbreviated labels for the y-axis of the bar charts.

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
  "Sea surface oxygen"                    = "Surf. O2",
  "Sea surface salinity"                  = "Surf. sal.",
  "Sea surface temperature"               = "Surf. temp.",
  "Surface pH"                            = "Surf. pH",
  "Surface wind speed magnitude"          = "Wind speed",
  "Sargassum influx"                      = "Sargassum",
  "Thermocline depth"                     = "Thermocl. depth"
)

##------------------------------------------------------------------------------
## Color palettes and factor level vectors  (shared across both figures)
##
## Vulnerability rank colors match 7-scoring-distributions/3-plot-figures.R
## for visual consistency across all CVA figures.

rank_levels <- c("Low", "Moderate", "High", "Very High")
rank_colors <- c(
  "Low"       = "green3",
  "Moderate"  = "yellow2",
  "High"      = "orange2",
  "Very High" = "red3"
)

## Directional effect colors use the Okabe-Ito colorblind-safe palette.
## Stacking order: Negative first (anchored at x = 0) so the dominant
## negative signal is directly readable as bar length.
dir_levels <- c("Negative", "Neutral", "Positive")
dir_colors <- c(
  "Negative" = "brown3",
  "Neutral"  = "bisque3",
  "Positive" = "turquoise4"
)

##------------------------------------------------------------------------------
## Shared ggplot theme  (used by all panels in both figures)
##
## theme_classic() gives a clean two-axis frame with no background grid,
## matching the HMS CVA reference style.

fig_theme <- ggplot2::theme_classic() +
  ggplot2::theme(
    axis.text.y       = ggplot2::element_text(size = 8.5, color = "black"),
    axis.text.x       = ggplot2::element_text(size = 8,   color = "black"),
    axis.title.x      = ggplot2::element_text(size = 9.5,
                                               margin = ggplot2::margin(t = 6)),
    ## No y-axis title — stock/attribute names are self-explanatory labels
    axis.title.y      = ggplot2::element_blank(),
    plot.tag          = ggplot2::element_text(face = "bold", size = 12),
    plot.tag.position = "topleft"
  )

################################################################################
##------------------------------------------------------------------------------
## Figure 1 — Leave-one-out influence analysis bar plots

##--- Read and prep: exposure LOO summary ---------------------------------------

exp_sum <- readr::read_csv(f_loo_exp_sum, show_col_types = FALSE)

exp_sum <- exp_sum %>%
  dplyr::mutate(
    short_name = dplyr::recode(factor_omitted, !!!exp_attr_short_names),
    ## fct_reorder: ascending n_rank_changed places the most influential
    ## factor at the top of the panel
    short_name = forcats::fct_reorder(short_name, n_rank_changed)
  )

##--- Read and prep: sensitivity LOO summary ------------------------------------
##
## "Stock size/status" in the CSV vs. "Stock Size Status" in the lookup key —
## recode before applying the short-name lookup.

sens_sum <- readr::read_csv(f_loo_sens_sum, show_col_types = FALSE)

sens_sum <- sens_sum %>%
  dplyr::mutate(
    attribute_omitted = dplyr::recode(attribute_omitted,
                                      "Stock size/status" = "Stock Size Status"),
    short_name = dplyr::recode(attribute_omitted, !!!attr_short_names),
    short_name = forcats::fct_reorder(short_name, n_rank_changed)
  )

loo_x_max_exp  <- ceiling(max(exp_sum$n_rank_changed,  na.rm = TRUE) / 4) * 4
loo_x_max_sens <- ceiling(max(sens_sum$n_rank_changed, na.rm = TRUE) / 4) * 4

##--- Panel A — Exposure factors ------------------------------------------------

p_exp <- ggplot2::ggplot(exp_sum,
                         ggplot2::aes(x = n_rank_changed, y = short_name)) +
  ggplot2::geom_col(fill = "black", width = 0.65) +
  ggplot2::scale_x_continuous(
    name   = "Number of Changes in Climate Vulnerability",
    limits = c(0, loo_x_max_exp),
    breaks = seq(0, loo_x_max_exp, by = 5),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::labs(tag = "A") +
  fig_theme

##--- Panel B — Sensitivity attributes -----------------------------------------

p_sens <- ggplot2::ggplot(sens_sum,
                          ggplot2::aes(x = n_rank_changed, y = short_name)) +
  ggplot2::geom_col(fill = "black", width = 0.65) +
  ggplot2::scale_x_continuous(
    name   = "Number of Changes in Climate Vulnerability",
    limits = c(0, loo_x_max_sens),
    breaks = seq(0, loo_x_max_sens, by = 1),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::labs(tag = "B") +
  fig_theme

##--- Combine and save: Figure 1 ------------------------------------------------

## Heights proportional to bar count so each bar is the same visual height:
## 13 exposure factors, 14 sensitivity attributes.
fig1 <- p_exp / p_sens +
  patchwork::plot_layout(heights = c(13, 14))

ggplot2::ggsave(
  filename = f_fig1,
  plot     = fig1,
  width    = 6.5, height = 9, units = "in", dpi = 300
)
cat("Figure 1 saved to:", f_fig1, "\n")

################################################################################
##------------------------------------------------------------------------------
## Figure 2 — Bootstrap resampling uncertainty

##--- Read inputs ---------------------------------------------------------------

boot_stock <- readr::read_csv(f_boot_stock, show_col_types = FALSE)
boot_dir   <- readr::read_csv(f_dir_boot,   show_col_types = FALSE)
vuln_base  <- readr::read_csv(f_vuln,        show_col_types = FALSE)

##--- Compute stock ordering ----------------------------------------------------
##
## Step 1: For each stock find the "dominant" bootstrap vulnerability rank
##         (the rank that captured the most iterations) and its proportion.
##         This proportion drives the secondary sort within each baseline rank
##         group, so borderline stocks (low dominant_prop) sit near the rank
##         boundaries where readers will look for uncertainty.

boot_dom <- boot_stock %>%
  dplyr::filter(n > 0) %>%             ## exclude zero-count rows before slice_max
  dplyr::group_by(stock_name) %>%
  ## with_ties = FALSE: defensive against exact ties (very unlikely with 10,000
  ## iterations, but ensures exactly one row per stock)
  dplyr::slice_max(prop, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(stock_name,
                dominant_rank = vuln_rank,
                dominant_prop = prop,
                borderline)

## Step 2: Build the ordered factor that all Figure 2 geom layers will use.

stock_order <- vuln_base %>%
  dplyr::select(stock_name, baseline_vuln_rank = Vuln_rank) %>%
  dplyr::left_join(boot_dom, by = "stock_name") %>%
  dplyr::mutate(
    ## Integer primary sort key: lower = higher vulnerability.
    ## Moderate gets 3 (→ bottom of chart), Very High gets 1 (→ top).
    rank_num = dplyr::case_when(
      baseline_vuln_rank == "Very High" ~ 1L,
      baseline_vuln_rank == "High"      ~ 2L,
      baseline_vuln_rank == "Moderate"  ~ 3L,
      baseline_vuln_rank == "Low"       ~ 4L
    ),
    ## Append asterisk to borderline stock labels; explained in figure caption
    label = dplyr::if_else(borderline,
                           paste0(stock_name, "*"),
                           stock_name)
  ) %>%
  ## desc(rank_num): Moderate first in vector (→ bottom of chart), VH last (→ top).
  ## dominant_prop ascending: least certain sits at the bottom of each group, so
  ## uncertainty is visible exactly at the Moderate/High and High/VH transitions.
  dplyr::arrange(desc(rank_num), dominant_prop) %>%
  dplyr::mutate(stock_fct = factor(label, levels = label))

## Named vector: raw stock_name → display label (with * for borderline stocks).
## Applied identically to both bootstrap data frames so their stock_fct factors
## use the same levels and the panels align row-for-row.
name_to_label <- stats::setNames(stock_order$label, stock_order$stock_name)

n_stocks <- nrow(stock_order)   ## used to position top-margin annotations

## Y-positions for horizontal rank-group separator lines.
## Layout: Moderate group (bottom) | High group (middle) | VH group (top).
n_mod  <- sum(stock_order$baseline_vuln_rank == "Moderate")
n_high <- sum(stock_order$baseline_vuln_rank == "High")

y_sep_mod_high <- n_mod            + 0.5   ## between Moderate and High
y_sep_high_vh  <- n_mod + n_high   + 0.5   ## between High and Very High

## X-position for the baseline rank indicator squares.
## Drawn just past x = 1.0; coord_cartesian(clip = "off") allows rendering
## in the right margin without distorting the axis scale.
strip_x <- 1.062

##--- Prep Panel A data: vulnerability bootstrap --------------------------------

boot_vuln_plot <- boot_stock %>%
  dplyr::mutate(
    stock_fct = factor(name_to_label[stock_name],
                       levels = levels(stock_order$stock_fct)),
    ## Factor order Low → VH controls stacking: Low leftmost, VH rightmost
    vuln_rank = factor(vuln_rank, levels = rank_levels)
  )

## One-row-per-stock baseline rank for the indicator squares.
## The squares use the same scale_fill_manual as the bars — no extra legend.
baseline_vuln_strip <- stock_order %>%
  dplyr::select(stock_fct, baseline_vuln_rank) %>%
  dplyr::mutate(baseline_vuln_rank = factor(baseline_vuln_rank, levels = rank_levels))

##--- Build Panel A: vulnerability bootstrap ------------------------------------
##
## Reading guide:
##   Each bar = one stock. Segment widths = proportion of 10,000 bootstrap
##   iterations in each vulnerability rank (Low/Moderate/High/Very High).
##   Solid single-color bar → robust classification.
##   Mixed-color bar straddling the 75% dashed line → borderline (*).
##   Small filled square on the right = finalized baseline rank for that stock.
##   Gray horizontal lines divide the three baseline rank groups.

p_boot_vuln <- ggplot2::ggplot(
    boot_vuln_plot,
    ggplot2::aes(x = prop, y = stock_fct, fill = vuln_rank)
  ) +

  ## Rank-group separator lines — behind bars (drawn first)
  ggplot2::geom_hline(
    yintercept = c(y_sep_mod_high, y_sep_high_vh),
    color = "gray55", linewidth = 0.40, linetype = "solid"
  ) +

  ## Stacked proportion bars.
  ## position_stack(reverse = FALSE) stacks in factor-level order: Low is placed
  ## leftmost because it is the first level, Very High rightmost as the last.
  ggplot2::geom_col(
    position = ggplot2::position_stack(reverse = FALSE),
    width    = 0.78
  ) +

  ## Baseline rank indicator squares (shape 22 = filled square).
  ## fill maps to baseline_vuln_rank, which shares rank_levels and rank_colors
  ## with the bar fill → same scale, same legend entry.
  ggplot2::geom_point(
    data        = baseline_vuln_strip,
    mapping     = ggplot2::aes(x    = strip_x,
                               y    = stock_fct,
                               fill = baseline_vuln_rank),
    shape = 22, size = 3.5, stroke = 0.25, color = "gray25",
    inherit.aes = FALSE
  ) +

  ## "Baseline" column header, rendered above the top bar via clip = "off"
  ggplot2::annotate(
    "text",
    x = strip_x, y = n_stocks + 0.85,
    label = "Baseline", size = 2.1, hjust = 0.5, color = "gray30"
  ) +

  ## Dashed vertical reference line at the 75% borderline threshold
  ggplot2::geom_vline(
    xintercept = 0.75, linetype = "dashed",
    color = "gray45", linewidth = 0.38
  ) +

  ## Label for the reference line, also in the top margin
  ggplot2::annotate(
    "text",
    x = 0.753, y = n_stocks + 0.85,
    label = "75%", size = 2.1, hjust = 0, color = "gray45"
  ) +

  ## Single fill scale: covers both bar segments and baseline squares.
  ## reverse = TRUE shows Very High at the top of the legend, matching the
  ## visual top-of-chart position of Very High stocks.
  ggplot2::scale_fill_manual(
    values = rank_colors, breaks = rank_levels,
    name   = "Vulnerability rank",
    guide  = ggplot2::guide_legend(reverse = TRUE)
  ) +

  ggplot2::scale_x_continuous(
    name   = "Proportion of bootstrap iterations",
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    labels = c("0%", "25%", "50%", "75%", "100%"),
    ## No padding left of zero so bars start flush with the y-axis
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +

  ## xlim constrains what is shown; clip = "off" allows geoms past x = 1.0
  ## (indicator squares and their header) to render in the right margin.
  ggplot2::coord_cartesian(xlim = c(0, 1.0), clip = "off") +

  ggplot2::labs(tag = "A") +

  fig_theme +
  ggplot2::theme(
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.title     = ggplot2::element_text(size = 8),
    legend.text      = ggplot2::element_text(size = 7.5),
    ## Extra right margin so indicator squares and "Baseline" text are not
    ## clipped by the outer figure boundary (they sit ~6% past x = 1.0)
    plot.margin      = ggplot2::margin(t = 8, r = 48, b = 4, l = 4, unit = "pt")
  ); p_boot_vuln

##--- Prep Panel B data: directional effect bootstrap --------------------------

## Select only what is needed for plotting; complete() guards against any
## missing stock × direction combinations (all three should exist per stock,
## but this ensures zero-prop rows are present for the stacked bars).
boot_dir_plot <- boot_dir %>%
  dplyr::mutate(
    stock_fct = factor(name_to_label[stock_name],
                       levels = levels(stock_order$stock_fct)),
    dir_rank  = factor(dir_rank, levels = dir_levels)
  ) %>%
  dplyr::select(stock_fct, dir_rank, n, prop) %>%
  tidyr::complete(stock_fct, dir_rank, fill = list(n = 0L, prop = 0))

## Baseline directional rank per stock — one row produced by distinct().
baseline_dir_strip <- boot_dir %>%
  dplyr::distinct(stock_name, baseline_dir_rank) %>%
  dplyr::mutate(
    stock_fct         = factor(name_to_label[stock_name],
                               levels = levels(stock_order$stock_fct)),
    baseline_dir_rank = factor(baseline_dir_rank, levels = dir_levels)
  )

##--- Build Panel B: directional effect bootstrap ------------------------------
##
## Reading guide:
##   Same stock ordering as Panel A for row-for-row comparison.
##   Stacking: Negative anchored at x = 0 (leftmost), Positive at the right.
##   A near-solid vermillion bar = strong negative directional consensus.
##   A mostly-gray bar = classifier is uncertain (Neutral dominant).
##   Baseline directional rank shown by indicator square, same as Panel A.

p_boot_dir <- ggplot2::ggplot(
    boot_dir_plot,
    ggplot2::aes(x = prop, y = stock_fct, fill = dir_rank)
  ) +

  ## Separator lines at the same y-positions as Panel A for visual alignment
  ggplot2::geom_hline(
    yintercept = c(y_sep_mod_high, y_sep_high_vh),
    color = "gray55", linewidth = 0.40, linetype = "solid"
  ) +

  ggplot2::geom_col(
    position = ggplot2::position_stack(reverse = FALSE),
    width    = 0.78
  ) +

  ggplot2::geom_point(
    data        = baseline_dir_strip,
    mapping     = ggplot2::aes(x    = strip_x,
                               y    = stock_fct,
                               fill = baseline_dir_rank),
    shape = 22, size = 3.5, stroke = 0.25, color = "gray25",
    inherit.aes = FALSE
  ) +

  ggplot2::annotate(
    "text",
    x = strip_x, y = n_stocks + 0.85,
    label = "Baseline", size = 2.1, hjust = 0.5, color = "gray30"
  ) +

  ggplot2::scale_fill_manual(
    values = dir_colors, breaks = dir_levels,
    name   = "Directional effect"
  ) +

  ggplot2::scale_x_continuous(
    name   = "Proportion of bootstrap iterations",
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    labels = c("0%", "25%", "50%", "75%", "100%"),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +

  ggplot2::coord_cartesian(xlim = c(0, 1.0), clip = "off") +

  ggplot2::labs(tag = "B") +

  fig_theme +
  ggplot2::theme(
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.title     = ggplot2::element_text(size = 8),
    legend.text      = ggplot2::element_text(size = 7.5),
    plot.margin      = ggplot2::margin(t = 8, r = 48, b = 4, l = 4, unit = "pt")
  )

##--- Combine panels and save: Figure 2 ----------------------------------------
##
## Equal heights (25 stocks per panel).
## plot_annotation adds a figure-level footnote defining the asterisk notation
## and the rank-group separator lines.

fig2 <- p_boot_vuln / p_boot_dir +
  patchwork::plot_layout(heights = c(1, 1)) +
  patchwork::plot_annotation(
    caption = paste0(
      "* Dominant bootstrap rank proportion < 75%: stock classification is",
      " sensitive to reviewer score uncertainty.\n",
      "Horizontal gray lines separate baseline vulnerability rank groups",
      " (Moderate, High, Very High, from bottom to top)."
    ),
    theme = ggplot2::theme(
      plot.caption = ggplot2::element_text(
        size   = 7,
        hjust  = 0,
        color  = "gray30",
        margin = ggplot2::margin(t = 6)
      )
    )
  )
fig2


ggplot2::ggsave(
  filename = f_fig2,
  plot     = fig2,
  width    = 7.5, height = 12, units = "in", dpi = 300
)
cat("Figure 2 saved to:", f_fig2, "\n")
