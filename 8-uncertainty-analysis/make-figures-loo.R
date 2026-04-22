################################################################################
##------------------------------------------------------------------------------
## Caribbean CVA — Leave-one-out bar plot figure
##
## Produces a two-panel horizontal bar chart:
##   Panel A — Exposure factors: number of stocks that change vulnerability rank
##             when each exposure factor is omitted from the LOO analysis.
##   Panel B — Sensitivity attributes: same metric for each sensitivity attribute.
##
## Bars are ordered by n_rank_changed (ascending so the most influential
## attribute sits at the top of each panel). Short display names are applied
## from the same lookup tables used by 7-scoring-distributions/3-plot-figures.R.
##
## Output: figures/fig_loo_bar_plots.png  (6.5 x 9 in, 300 dpi)
##
## Run from the Caribbean-CVA RStudio project root (.Rproj file location).
##------------------------------------------------------------------------------

rm(list = ls()); gc()

library(dplyr)
library(readr)
library(ggplot2)
library(forcats)
library(patchwork)

##------------------------------------------------------------------------------
## Directories

proj_dir <- "."
dir_out  <- file.path(proj_dir, "figures")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

##------------------------------------------------------------------------------
## Input file paths

f_loo_sens_sum <- file.path(proj_dir, "outputs", "analyses", "uncertainty-loo",
                            "final-tables",
                            "table_leave_one_out_sensitivity_summary.csv")

f_loo_exp_sum  <- file.path(proj_dir, "outputs", "analyses", "uncertainty-loo",
                            "final-tables",
                            "table_leave_one_out_exposure_summary.csv")

##------------------------------------------------------------------------------
## Output file path

f_out <- file.path(dir_out, "fig_loo_bar_plots.png")

##------------------------------------------------------------------------------
## Short display name lookup tables
##
## Copied from 7-scoring-distributions/3-plot-figures.R (lines 104-137).
## Keys are the full attribute/factor names stored in the LOO CSV files.
## Values are the abbreviated labels used on the y-axis.

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
## Read and prep exposure summary

exp_sum <- readr::read_csv(f_loo_exp_sum, show_col_types = FALSE)

## Apply short names via lookup; fct_reorder puts the highest bar at the top
## of the panel (ascending order so ggplot flips work correctly).
exp_sum <- exp_sum %>%
  dplyr::mutate(
    short_name = dplyr::recode(factor_omitted, !!!exp_attr_short_names),
    short_name = forcats::fct_reorder(short_name, n_rank_changed)
  )

##------------------------------------------------------------------------------
## Read and prep sensitivity summary
##
## The LOO CSV stores "Stock size/status" (sentence-case with slash) but the
## short-name lookup key uses "Stock Size Status" (Title Case, no slash) —
## recode before applying the lookup.

sens_sum <- readr::read_csv(f_loo_sens_sum, show_col_types = FALSE)

sens_sum <- sens_sum %>%
  dplyr::mutate(
    ## Align the one mismatched key to the lookup table
    attribute_omitted = dplyr::recode(attribute_omitted,
                                      "Stock size/status" = "Stock Size Status"),
    short_name = dplyr::recode(attribute_omitted, !!!attr_short_names),
    short_name = forcats::fct_reorder(short_name, n_rank_changed)
  )

##------------------------------------------------------------------------------
## Shared ggplot theme

loo_theme <- ggplot2::theme_classic() +
  ggplot2::theme(
    ## y-axis (factor/attribute names) — keep compact but readable
    axis.text.y  = ggplot2::element_text(size = 8.5, color = "black"),
    axis.text.x  = ggplot2::element_text(size = 8, color = "black"),
    axis.title.x = ggplot2::element_text(size = 9.5, margin = ggplot2::margin(t = 6)),
    ## No y-axis title — labels are self-explanatory
    axis.title.y = ggplot2::element_blank(),
    ## Bold panel tag (A / B)
    plot.tag     = ggplot2::element_text(face = "bold", size = 12),
    plot.tag.position = "topleft"
  )

##------------------------------------------------------------------------------
## Panel A — Exposure factors

p_exp <- ggplot2::ggplot(exp_sum,
                         ggplot2::aes(x = n_rank_changed, y = short_name)) +
  ggplot2::geom_col(fill = "black", width = 0.65) +
  ggplot2::scale_x_continuous(
    name   = "Number of Changes in Climate Vulnerability",
    limits = c(0, 25),
    breaks = seq(0, 25, by = 5),
    ## No gap between axis line and bars at zero
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::labs(tag = "A") +
  loo_theme

##------------------------------------------------------------------------------
## Panel B — Sensitivity attributes

p_sens <- ggplot2::ggplot(sens_sum,
                          ggplot2::aes(x = n_rank_changed, y = short_name)) +
  ggplot2::geom_col(fill = "black", width = 0.65) +
  ggplot2::scale_x_continuous(
    name   = "Number of Changes in Climate Vulnerability",
    limits = c(0, 25),
    breaks = seq(0, 25, by = 5),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::labs(tag = "B") +
  loo_theme

##------------------------------------------------------------------------------
## Combine panels and save

## Stack A over B; patchwork handles spacing automatically.
## Heights are proportional to number of bars: exp = 13 bars, sens = 14 bars.
fig_loo <- p_exp / p_sens +
  patchwork::plot_layout(heights = c(13, 14))
fig_loo

ggplot2::ggsave(
  filename = f_out,
  plot     = fig_loo,
  width    = 6.5, height = 9, units = "in", dpi = 300
)

cat("Figure saved to:", f_out, "\n")
