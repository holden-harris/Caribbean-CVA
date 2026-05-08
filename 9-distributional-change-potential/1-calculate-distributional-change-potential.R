################################################################################
##------------------------------------------------------------------------------
## Caribbean CVA
## Module 9 — Potential for Distributional Change
## Script 1: Calculate baseline distributional change potential (DCP) scores
##
## For each of the 25 assessed stocks, computes a DCP rank
## (Low / Moderate / High / Very High) using four sensitivity attributes.
## Three movement-related attributes are inverted (5 - mean) before applying
## the FCVA logic model so that a high transformed score consistently indicates
## "able to shift distribution":
##
##   Inverted (5 - mean):  Adult mobility
##                         Habitat specificity
##                         Mobility and dispersal or early life stages
##   Not inverted:         Species range (analog for temperature sensitivity;
##                         wider latitudinal range = broader thermal tolerance =
##                         lower temperature sensitivity = lower shift propensity)
##
## NOTE — open methodological decision: prior NOAA CVAs (HMS, SA, GoM) include
## "Sensitivity to Temperature" as the fourth attribute. The Caribbean CVA does
## not score this attribute directly. "Species range" serves as the analog here
## because both use latitudinal extent as a proxy for thermal tolerance.
## This substitution should be confirmed with the project lead before finalizing.
## To exclude the temperature analog, remove "Species range" from target_attributes
## below — the script requires no other changes.
##
## Inputs:
##   outputs/final-scores-compiled/overall-vulnerability-rankings/
##     attribute_means_uscar.csv
##
## Outputs:
##   outputs/distributional_change_potential_uscar.csv
##     Columns: stock_name, adult_mobility_inverted, habitat_specificity_inverted,
##              early_life_dispersal_inverted, species_range, dcp_numeric, dcp_rank
##
## Dependencies: dplyr, tidyr, readr, stringr

##------------------------------------------------------------------------------
## Setup

rm(list = ls()); gc()

library(dplyr)
library(tidyr)
library(readr)
library(stringr)

##------------------------------------------------------------------------------
## Configuration

## FCVA logic model threshold — must match Module 4 Script 3 and Module 8
rank_threshold <- 2

## Exact attribute names as they appear in attribute_means_uscar.csv.
## If these strings do not match the CSV, the script halts with a clear error.
target_attributes <- c(
  "Adult mobility",
  "Habitat specificity",
  "Mobility and dispersal or early life stages",
  "Species range"
)

## Three movement-related attributes: transformed as 5 - mean before the logic model.
## After inversion, a higher score means greater propensity to shift.
inverted_attrs <- c(
  "Adult mobility",
  "Habitat specificity",
  "Mobility and dispersal or early life stages"
)

##------------------------------------------------------------------------------
## File paths
## proj_dir = "." assumes the script is run from the RStudio project root
## (C:/Repos/Caribbean-CVA), which is the default when using the .Rproj file.

proj_dir     <- "."
compiled_dir <- file.path(proj_dir, "outputs", "final-scores-compiled",
                          "overall-vulnerability-rankings")
out_dir      <- file.path(proj_dir, "outputs", "distribution-change-potential")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

f_attr_means <- file.path(compiled_dir, "attribute_means_uscar.csv")
f_out        <- file.path(out_dir,      "distributional_change_potential_uscar.csv")

##------------------------------------------------------------------------------
## FCVA logic model (copied from Module 8 — uncertainty-analyses.R)
##
## Converts a vector of attribute mean scores (each 1–4) to a component rank
## and numeric score using the same cascading threshold rules as Modules 4 and 8.
##
## With rank_threshold = 2 (attr_means_current):
##   Very High: > 3 attribute means >= 3.5
##   High:      > 2 attribute means >= 3.0
##   Moderate:  > 2 attribute means >= 2.5
##   Low:       all other cases

fcva_logic_model <- function(mean_scores, rank_threshold) {

  n_ge_35 <- sum(mean_scores >= 3.5, na.rm = TRUE)
  n_ge_30 <- sum(mean_scores >= 3.0, na.rm = TRUE)
  n_ge_25 <- sum(mean_scores >= 2.5, na.rm = TRUE)

  component_rank <- dplyr::case_when(
    n_ge_35 > rank_threshold + 1 ~ "Very High",
    n_ge_30 > rank_threshold     ~ "High",
    n_ge_25 > rank_threshold     ~ "Moderate",
    TRUE                         ~ "Low"
  )

  component_score_numeric <- dplyr::case_when(
    component_rank == "Low"       ~ 1L,
    component_rank == "Moderate"  ~ 2L,
    component_rank == "High"      ~ 3L,
    component_rank == "Very High" ~ 4L
  )

  tibble::tibble(component_rank, component_score_numeric)
}

##------------------------------------------------------------------------------
## Read attribute means

attr_means_raw <- readr::read_csv(f_attr_means, show_col_types = FALSE)

##------------------------------------------------------------------------------
## Validate attribute names
##
## Print all Sensitivity attribute names found in the source file, then halt if
## any target attribute is missing. Catches capitalization mismatches early.

sens_names_available <- unique(
  attr_means_raw$attribute_name[attr_means_raw$attribute_type == "Sensitivity"]
)

cat("Sensitivity attributes in attribute_means_uscar.csv:\n")
cat(sort(sens_names_available), sep = "\n")
cat("\n")

missing_attrs <- setdiff(target_attributes, sens_names_available)
if (length(missing_attrs) > 0) {
  stop(
    "Target attributes not found in attribute_means_uscar.csv:\n  ",
    paste(missing_attrs, collapse = "\n  "),
    "\nAvailable Sensitivity attributes:\n  ",
    paste(sort(sens_names_available), collapse = "\n  ")
  )
}
message("✓ All four target attributes found in attribute_means_uscar.csv.")

##------------------------------------------------------------------------------
## Filter to target attributes and validate completeness

dcp_means <- attr_means_raw %>%
  dplyr::filter(
    attribute_type == "Sensitivity",
    attribute_name %in% target_attributes
  ) %>%
  dplyr::transmute(
    stock_name     = stock_name,
    attribute_name = stringr::str_squish(attribute_name),
    mean_score     = as.numeric(attribute_mean)
  )

## Halt if any NA means
na_rows <- dplyr::filter(dcp_means, is.na(mean_score))
if (nrow(na_rows) > 0) {
  stop(
    "NA attribute means found for:\n",
    paste(paste(" ", na_rows$stock_name, "/", na_rows$attribute_name), collapse = "\n")
  )
}

## Halt if any stock is missing one or more target attributes
attr_counts <- dcp_means %>%
  dplyr::count(stock_name, name = "n_attrs")
incomplete <- dplyr::filter(attr_counts, n_attrs != length(target_attributes))
if (nrow(incomplete) > 0) {
  stop(
    "Stocks missing one or more target attributes:\n",
    paste(paste(" ", incomplete$stock_name,
                "(", incomplete$n_attrs, "of", length(target_attributes), ")"),
          collapse = "\n")
  )
}
message("✓ All stocks have complete data for all four target attributes.")

##------------------------------------------------------------------------------
## Apply inversion to movement attributes

dcp_means <- dcp_means %>%
  dplyr::mutate(
    transformed_mean = dplyr::if_else(
      attribute_name %in% inverted_attrs,
      5 - mean_score,
      mean_score
    )
  )

##------------------------------------------------------------------------------
## Apply FCVA logic model per stock

dcp_scores <- dcp_means %>%
  dplyr::group_by(stock_name) %>%
  dplyr::summarise(
    adult_mobility_inverted       = transformed_mean[attribute_name == "Adult mobility"],
    habitat_specificity_inverted  = transformed_mean[attribute_name == "Habitat specificity"],
    early_life_dispersal_inverted = transformed_mean[attribute_name == "Mobility and dispersal or early life stages"],
    species_range                 = transformed_mean[attribute_name == "Species range"],
    dcp_logic = list(fcva_logic_model(transformed_mean, rank_threshold)),
    .groups = "drop"
  ) %>%
  tidyr::unnest(cols = dcp_logic) %>%
  dplyr::rename(
    dcp_rank    = component_rank,
    dcp_numeric = component_score_numeric
  )

##------------------------------------------------------------------------------
## Validate output and summarize

n_stocks_out <- nrow(dcp_scores)
if (n_stocks_out != 25) {
  warning("Expected 25 stocks in output; got ", n_stocks_out, ". Investigate before proceeding.")
} else {
  message("✓ All 25 stocks present in output.")
}

cat("\nDistributional Change Potential — rank distribution:\n")
print(table(dcp_scores$dcp_rank))

##------------------------------------------------------------------------------
## Write output

readr::write_csv(dcp_scores, f_out)
cat("\nOutput written to:", f_out, "\n")
