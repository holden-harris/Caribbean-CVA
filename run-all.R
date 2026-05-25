################################################################################
## Caribbean CVA — Master orchestration script
##
## Runs the full pipeline in two phases:
##   Phase 1: All analysis scripts in dependency order (Modules 4–9)
##   Phase 2: All figure scripts (Module 10)
##
## Run from the Caribbean-CVA RStudio project root (.Rproj file location).
## Each sourced script calls rm(list = ls()) internally, so environment state
## does not carry between scripts.
##
## TO SWITCH RUNS: edit active_run in config.R, then re-run this script.
##   active_run <- "broadened_distribution"   (rank_threshold = 2L, 14 attributes)
##   active_run <- "cross_region_comparable"  (rank_threshold = 1L, 12 attributes)
## Outputs for each run land in outputs/{run_label}/ and figures/{run_label}/.
################################################################################

## =============================================================================
## PHASE 1: Analyses
## =============================================================================

## Module 4 — Final attribute/exposure scoring and overall vulnerability
source("04-final-attribute-exposure-scoring/1-extract-final-scores-from-all-reviewers.R")
source("04-final-attribute-exposure-scoring/2-summarize-attribute-scores.R")
source("04-final-attribute-exposure-scoring/3-calculate-overall-vulnerability-scores.R")

## Module 5 — Directional effect scoring
source("05-final-directional-effect-scoring/1-extract-directional-effect.R")
source("05-final-directional-effect-scoring/2-summarize-directional-effect.R")

## Module 6 — Data quality scoring
source("06-final-data-quality-scoring/1-extract-data-quality-scores.R")
source("06-final-data-quality-scoring/2-summarize-data-quality-scores.R")

## Module 7 — Scoring distributions (tally extraction and finalization)
source("07-scoring-distributions/1-extract-tally-scores.R")
source("07-scoring-distributions/2-finalize-tally-tables.R")

## Module 8 — Uncertainty analysis (bootstrap + leave-one-out)
source("08-uncertainty-analysis/uncertainty-analyses.R")

## Module 9 — Distributional change potential
source("09-distributional-change-potential/1-calculate-distributional-change-potential.R")
source("09-distributional-change-potential/2-bootstrap-distributional-change.R")

## =============================================================================
## PHASE 2: Figures
## =============================================================================

## Module 10 — All publication figures
source("10-figures/1-plot-scoring-distributions.R")    ## score distributions + tally panels
source("10-figures/2-plot-uncertainty-figures.R")      ## LOO bar charts + bootstrap uncertainty
source("10-figures/3-plot-distributional-change.R")    ## DCP ranks + DCP vs. vulnerability
source("10-figures/4-plot-overall-vulnerability.R")    ## overall vulnerability grid + directional effect
source("10-figures/5-produce-results-tables.R")        ## directional effect + data quality result tables

cat("\n=== All analyses and figures complete ===\n")
