################################################################################
##------------------------------------------------------------------------------
## Caribbean CVA
## Module 9 — Potential for Distributional Change
## Script 3: Plot distributional change potential figures
##
## Produces two publication figures saved to figures/:
##
##   Figure A — fig_distributional_change_ranks.png
##     Column chart matching Craig et al. (2025) South Atlantic CVA Fig. 5:
##     one solid-colored column per DCP rank (Low / Moderate / High / Very High).
##     Column height = number of stocks. Stock names listed inside the column
##     from top (highest certainty) to bottom (lowest certainty); alphabetical
##     within each certainty tier.
##     Label format: "Stock name (V)" where V = overall vulnerability abbreviation
##       (L = Low, M = Moderate, H = High, VH = Very High).
##     Font and color encode bootstrap certainty:
##       > 0.95 dominant_prop → black bold    (Very High certainty)
##       0.90–0.95            → black italic  (High certainty)
##       0.67–0.89            → white bold    (Moderate certainty)
##       ≤ 0.66               → white italic  (Low certainty)
##
##   Figure B — fig_distributional_change_vs_vulnerability.png
##     4×4 cross-plot grid: x = DCP rank, y = Overall Vulnerability rank.
##     Cells colored by overall vulnerability rank. Stock names positioned
##     within cells with certainty-encoded font. The High/VH vulnerability +
##     Low/Moderate DCP quadrant ("stuck and exposed") is the primary
##     management-concern view.
##
## Inputs:
##   outputs/distributional_change_full_uscar.csv           (Scripts 1+2 joined)
##   outputs/final-scores-compiled/overall-vulnerability-rankings/
##     overall_vulnerability_scores_uscar.csv
##
## Outputs:
##   figures/fig_distributional_change_ranks.png
##   figures/fig_distributional_change_vs_vulnerability.png
##
## Dependencies: ggplot2, dplyr, readr, stringr

##------------------------------------------------------------------------------
## Setup

rm(list = ls()); gc()

library(ggplot2)
library(dplyr)
library(readr)
library(stringr)

##------------------------------------------------------------------------------
## File paths
## proj_dir = "." assumes the script is run from the RStudio project root
## (C:/Repos/Caribbean-CVA), which is the default when using the .Rproj file.

proj_dir     <- "."
compiled_dir <- file.path(proj_dir, "outputs", "final-scores-compiled",
                          "overall-vulnerability-rankings")
out_dir      <- file.path(proj_dir, "outputs", "distribution-change-potential")
fig_dir      <- file.path(proj_dir, "figures")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

f_dcp_full  <- file.path(out_dir,      "distributional_change_full_uscar.csv")
f_vuln      <- file.path(compiled_dir, "overall_vulnerability_scores_uscar.csv")
f_fig_ranks <- file.path(fig_dir,      "fig_distributional_change_ranks.png")
f_fig_cross <- file.path(fig_dir,      "fig_distributional_change_vs_vulnerability.png")

##------------------------------------------------------------------------------
## Shared constants

rank_order  <- c("Low", "Moderate", "High", "Very High")
rank_levels <- c("Low", "Moderate", "High", "Very High")

## Standardized CVA color palette (matches CVA project convention)
rank_colors <- c(
  "Low"       = "green3",
  "Moderate"  = "yellow4",
  "High"      = "orange2",
  "Very High" = "red3"
)

rank_num <- c("Low" = 1, "Moderate" = 2, "High" = 3, "Very High" = 4)

## Overall vulnerability abbreviations for stock labels
vuln_abbrev_map <- c(
  "Low"       = "L",
  "Moderate"  = "M",
  "High"      = "H",
  "Very High" = "VH"
)

## Canonical stock name lookup: CSV sentence-case names → display names.
## Applied to all tables on read so figure labels use project-standard names.
stock_name_recode <- c(
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
  "Queen snapper"           = "Queen Snapper",
  "Queen triggerfish"       = "Queen Triggerfish",
  "Rainbow parrotfish"      = "Rainbow Parrotfish",
  "Red grouper"             = "Red Grouper",
  "Redhind"                 = "Red Hind",
  "Sea cucumber"            = "Sea Cucumber",
  "Silk snapper"            = "Silk Snapper",
  "Spiny lobster"           = "Spiny Lobster",
  "Stoplight parrotfish"    = "Stoplight Parrotfish",
  "White mullet"            = "White Mullet",
  "Yellowfin grouper"       = "Yellowfin Grouper",
  "Yellowtail snapper"      = "Yellowtail Snapper"
)

##------------------------------------------------------------------------------
## Read data

dcp_full <- readr::read_csv(f_dcp_full, show_col_types = FALSE)
vuln     <- readr::read_csv(f_vuln,     show_col_types = FALSE)

## Join overall vulnerability rank to DCP data
dcp_full <- dcp_full %>%
  dplyr::left_join(
    dplyr::select(vuln, stock_name, Vuln_rank),
    by = "stock_name"
  ) %>%
  ## Apply display name standardization after the join so the key matches
  dplyr::mutate(stock_name = dplyr::recode(stock_name, !!!stock_name_recode))

missing_vuln <- dplyr::filter(dcp_full, is.na(Vuln_rank))
if (nrow(missing_vuln) > 0) {
  warning("Stocks not joined to overall vulnerability table: ",
          paste(missing_vuln$stock_name, collapse = ", "))
}

##------------------------------------------------------------------------------
## Derive certainty categories, text encoding, and stock labels
##
## Certainty thresholds match Module 4 / Module 8 / legend description:
##   > 0.95 → Very High certainty: black bold
##   0.90–0.95 → High certainty: black italic
##   0.67–0.89 → Moderate certainty: white bold
##   ≤ 0.66 → Low certainty: white italic
##
## Label format: "Stock name (V)" where V = overall vulnerability abbreviation
##   e.g., "Atlantic thread herring (M)"

dcp_plot <- dcp_full %>%
  dplyr::mutate(
    dcp_rank  = factor(dcp_rank,  levels = rank_order),
    Vuln_rank = factor(Vuln_rank, levels = rank_order),
    certainty = dplyr::case_when(
      dominant_prop >  0.95 ~ "Very High",
      dominant_prop >= 0.90 ~ "High",
      dominant_prop >= 0.67 ~ "Moderate",
      TRUE                  ~ "Low"
    ),
    ## Numeric order for within-column sorting: 1 = most certain (top of column)
    certainty_order = dplyr::case_when(
      certainty == "Very High" ~ 1L,
      certainty == "High"      ~ 2L,
      certainty == "Moderate"  ~ 3L,
      TRUE                     ~ 4L
    ),
    text_face  = dplyr::if_else(
      certainty %in% c("Very High", "Moderate"), "bold", "italic"
    ),
    text_color = dplyr::if_else(
      certainty %in% c("Very High", "High"), "black", "white"
    ),
    vuln_abbrev = vuln_abbrev_map[as.character(Vuln_rank)],
    label       = paste0(stock_name, " (", vuln_abbrev, ")")
  )

if (nrow(dcp_plot) != 25) {
  warning("Expected 25 stocks; got ", nrow(dcp_plot))
}

################################################################################
##------------------------------------------------------------------------------
## Figure A — DCP rank categories (column chart matching SA CVA Fig. 5)
##
## Column height = number of stocks in each DCP rank category.
## Within each column, stocks are ordered highest certainty at top,
## then alphabetically within each certainty tier.
## Text is centered in each 1-unit vertical slot (y = y_pos - 0.5).

## Column heights for geom_col background
col_heights <- dcp_plot %>%
  dplyr::count(dcp_rank, name = "n_stocks")

## Stock label positions within each column
## Sorted: most certain first (top of column), then alphabetical within tier
fig_a_data <- dcp_plot %>%
  dplyr::group_by(dcp_rank) %>%
  dplyr::arrange(certainty_order, stock_name, .by_group = TRUE) %>%
  dplyr::mutate(
    n_in_rank = dplyr::n(),
    ## row 1 (highest certainty) → top of column; row N → bottom
    y_pos = n_in_rank - dplyr::row_number() + 1L
  ) %>%
  dplyr::ungroup()

## Y-axis ceiling: round up to nearest 5 above max count
y_max <- ceiling((max(col_heights$n_stocks) + 0.5) / 3) * 3

p_ranks <- ggplot() +
  ## Solid colored column — height = number of stocks in that DCP rank
  geom_col(
    data  = col_heights,
    aes(x = dcp_rank, y = n_stocks, fill = dcp_rank),
    color = "black", linewidth = 0.35,
    width = 0.85
  ) +
  ## Stock name labels centered in each 1-unit vertical slot (y_pos - 0.5)
  geom_text(
    data  = fig_a_data,
    aes(x        = dcp_rank,
        y        = y_pos - 0.5,
        label    = label,
        color    = text_color,
        fontface = text_face),
    size = 3.2, hjust = 0.5
  ) +
  scale_fill_manual(values = rank_colors, breaks = rank_order, guide = "none") +
  scale_color_identity() +
  scale_x_discrete(limits = rank_order, expand = expansion(add = 0.7)) +
  scale_y_continuous(
    limits = c(0, y_max),
    breaks = seq(0, y_max, by = 3),
    expand = c(0, 0)
  ) +
  labs(
    x = "Species Distribution Change Potential",
    y = "Number of species"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_blank(),
    axis.text.x        = element_text(size = 12, color = "black"),
    axis.text.y        = element_text(size = 11, color = "black"),
    axis.title.x       = element_text(size = 12, face = "bold", margin = margin(t = 8)),
    axis.title.y       = element_text(size = 12, face = "bold", margin = margin(r = 8)),
    legend.position    = "none",
    plot.margin        = margin(10, 15, 10, 10)
  ); p_ranks

ggsave(
  filename = f_fig_ranks,
  plot     = p_ranks,
  width    = 10, height = 8, units = "in",
  dpi      = 1200, bg = "white"
)
cat("Figure A written to:", f_fig_ranks, "\n")

################################################################################
##------------------------------------------------------------------------------
## Figure B — DCP vs. Overall Vulnerability cross-plot
##
## 4×4 tile grid: x = DCP rank, y = Overall Vulnerability rank.
## Cells colored by vulnerability rank. Stock names positioned within cells
## with certainty-encoded font; y-offset separates stocks in the same cell.

## 4×4 background tile grid — all 16 cells, colored by vulnerability rank.
## Both axes use numeric positions (1–4) so geom_tile renders correctly;
## axis labels are added via scale breaks. Mixing a discrete x factor with
## a continuous y in geom_tile causes a silent blank-layer failure in ggplot2.
tile_grid <- data.frame(
  dcp_num  = rep(1:4, times = 4),
  vuln_num = rep(1:4, each  = 4)
) %>%
  dplyr::mutate(
    vuln_rank = factor(rank_order[vuln_num], levels = rank_order)
  )

## Stock positions in the grid with within-cell y offsets
cross_data <- dcp_plot %>%
  dplyr::mutate(
    dcp_num  = rank_num[as.character(dcp_rank)],
    vuln_num = rank_num[as.character(Vuln_rank)]
  ) %>%
  dplyr::group_by(dcp_num, vuln_num) %>%
  dplyr::arrange(certainty_order, stock_name, .by_group = TRUE) %>%
  dplyr::mutate(
    n_cell    = dplyr::n(),
    cell_rank = dplyr::row_number(),
    y_offset  = ((n_cell + 1) / 2 - cell_rank) * (0.8 / max(n_cell, 1)),
    y_pos     = vuln_num + y_offset
  ) %>%
  dplyr::ungroup()

p_cross <- ggplot() +
  geom_tile(
    data = tile_grid,
    aes(x = dcp_num, y = vuln_num, fill = vuln_rank),
    color = "white", linewidth = 0.6,
    width = 1, height = 1
  ) +
  geom_text(
    data = cross_data,
    aes(x        = dcp_num,
        y        = y_pos,
        label    = label,
        color    = text_color,
        fontface = text_face),
    size = 3.0, hjust = 0.5
  ) +
  scale_fill_manual(
    values  = rank_colors, breaks = rank_order,
    name    = "Overall climate\nvulnerability"
  ) +
  scale_color_identity() +
  scale_x_continuous(
    breaks = 1:4, labels = rank_order,
    limits = c(0.5, 4.5), expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = 1:4, labels = rank_order,
    limits = c(0.5, 4.5), expand = c(0, 0)
  ) +
  labs(
    x = "Potential for Distributional Change",
    y = "Overall Climate Vulnerability"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text      = element_text(size = 11, color = "black"),
    axis.title.x   = element_text(size = 12, face = "bold", margin = margin(t = 8)),
    axis.title.y   = element_text(size = 12, face = "bold", margin = margin(r = 8)),
    panel.grid     = element_blank(),
    axis.line      = element_line(color = "black", linewidth = 0.8),
    legend.position = "right",
    plot.margin    = margin(10, 10, 10, 10)
  ); p_cross

ggsave(
  filename = f_fig_cross,
  plot     = p_cross,
  width    = 9, height = 7, units = "in",
  dpi      = 1000, bg = "white"
)
cat("Figure B written to:", f_fig_cross, "\n")

cat("\nDone. Both figures written to:", fig_dir, "\n")
