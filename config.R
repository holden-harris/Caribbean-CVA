################################################################################
## config.R — Caribbean CVA shared constants
##
## Source this file at the top of every pipeline script (after rm(list = ls())).
## Changing a value here propagates to all analyses and figures automatically.
################################################################################

## --- FCVA Logic Model ---------------------------------------------------------
##
## rank_threshold sets the minimum number of attributes that must meet each score
## cutoff before a component rank of Moderate or higher is assigned. It is applied
## per stock × component (Sensitivity and Exposure) in Modules 04, 08, and 09.
##
## Logic model (from 04-final-attribute-exposure-scoring/3-calculate-overall-vulnerability-scores.R):
##   Very High : n attributes with mean ≥ 3.5  >  rank_threshold + 1
##   High      : n attributes with mean ≥ 3.0  >  rank_threshold
##   Moderate  : n attributes with mean ≥ 2.5  >  rank_threshold
##   Low       : otherwise (evaluated in order; first condition met wins)
##
## Current setting — rank_threshold = 1 (standard NOAA FCVA model):
##   Very High  →  ≥ 3 attributes with mean ≥ 3.5
##   High       →  ≥ 2 attributes with mean ≥ 3.0
##   Moderate   →  ≥ 2 attributes with mean ≥ 2.5
##   Low        →  0 or 1 attribute meeting any threshold above
##
## To apply a more restrictive model (fewer High / Very High outcomes):
##   rank_threshold <- 2L  →  Very High ≥ 4 attrs; High / Moderate ≥ 3 attrs
##
## To apply a more permissive model (more High / Very High outcomes):
##   rank_threshold <- 0L  →  Very High ≥ 2 attrs; High / Moderate ≥ 1 attr
##
## After changing rank_threshold, re-run the full pipeline (run-all.R) to
## propagate updated ranks through all analyses and figures.

rank_threshold    <- 1L
dir_eff_threshold <- 1/3  ## |weighted mean| cutoff for Neg/Pos directional classification
borderline_prop   <- 0.75 ## bootstrap dominant-rank prop below which a stock is borderline

## --- Bootstrap certainty bins -------------------------------------------------

cert_very_high <- 0.95
cert_high      <- 0.90
cert_moderate  <- 0.67

## --- Rank factor levels -------------------------------------------------------

rank_levels <- c("Low", "Moderate", "High", "Very High")

## Standard vulnerability rank colors.
## Scripts with intentionally different palettes define local overrides:
##   10-figures/1-plot-scoring-distributions.R  uses lighter tally-bar colors
##   10-figures/4-plot-overall-vulnerability.R  uses hex vuln_colors for the tile grid
##   These colors allow white and black writing to show. The local colors use brighter 
##   and oranges where the white doesn't show up. 

rank_colors <- c(
  "Low"       = "green3",
  "Moderate"  = "yellow2",
  "High"      = "orange2",
  "Very High" = "red3"
)

## --- Directional effect -------------------------------------------------------

dir_levels <- c("Negative", "Neutral", "Positive")
dir_colors <- c(
  "Negative" = "brown3",
  "Neutral"  = "bisque3",
  "Positive" = "turquoise4"
)

## --- Canonical stock name lookup ----------------------------------------------
## Applied to any table that may carry variant spellings from workbooks or CSVs.

stock_name_recode <- c(
  "Atlantic thread herring" = "Atlantic Herring",
  "Long-spined sea urchin"  = "Diadema",
  "Red hind"                = "Red Hind",
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

## --- Attribute display labels (used by figures) -------------------------------

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
