################################################################################
##------------------------------------------------------------------------------
## Caribbean CVA
## Extract FINAL SCORE tally tables from reviewer scoring workbooks
##
## Purpose:
## - Read reviewer scoring workbooks
## - Extract FINAL SCORE tally columns from stock sheets
## - Build long-format tally tables needed for bootstrap uncertainty analysis
##
## Output tables:
## - table_sensitivity_tallies_long.csv
## - table_qualitative_exposure_tallies_long.csv
## - table_directional_effect_tallies_long.csv
## - table_all_qualitative_tallies_long.csv
##
## Notes:
## - This script extracts tallies from the FINAL SCORE section only
## - Qualitative attributes use rank tallies in columns M:P
## - Directional effect uses tallies in column M for rows 38:40
## - Quantitative exposure factors are not extracted here because they are not
##   tally-based in the same way and are typically held fixed in NOAA-style
##   vulnerability uncertainty analyses

##------------------------------------------------------------------------------
## Libraries

rm(list = ls()); gc()

library(openxlsx)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(stringr)
library(tibble)
library(ggplot2)
library(forcats)
library(scales)

##------------------------------------------------------------------------------
## Directories

proj_dir <- "."
in_dir   <- file.path(proj_dir, "data", "final-scores")
out_dir  <- file.path(proj_dir, "outputs", "analyses", "1-inputs")
qa_dir   <- file.path(out_dir, "qa")
fig_dir  <- file.path(out_dir, "figures")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

##------------------------------------------------------------------------------
## Output files

f_sens_tallies <- file.path(out_dir,
                            "table_sensitivity_tallies_long.csv")

f_qexp_tallies <- file.path(out_dir,
                            "table_qualitative_exposure_tallies_long.csv")

f_dir_tallies  <- file.path(out_dir,
                            "table_directional_effect_tallies_long.csv")

f_all_tallies  <- file.path(out_dir,
                            "table_all_qualitative_tallies_long.csv")

## QA tables

f_qa_row_checks   <- file.path(qa_dir, "qa_tally_row_checks.csv")
f_qa_coverage     <- file.path(qa_dir, "qa_reviewer_stock_coverage.csv")
f_qa_attr_summary <- file.path(qa_dir, "qa_attribute_tally_summary.csv")
f_qa_dir_summary  <- file.path(qa_dir, "qa_directional_effect_summary.csv")

## Figures

f_fig_coverage <- file.path(fig_dir, "fig_reviewer_stock_coverage.png")
f_fig_tallies  <- file.path(fig_dir, "fig_sensitivity_tally_distributions.png")
f_fig_dir      <- file.path(fig_dir, "fig_directional_effect_summary.png")

##------------------------------------------------------------------------------
## Workbook files
##
## Assumption:
## - One workbook per scorer/reviewer
## - Reviewer identity can be inferred from the filename
##
## Update the file pattern if needed

workbook_files <- list.files(
  path = in_dir,
  pattern = "\\.xlsx$",
  full.names = TRUE
)

if(length(workbook_files) == 0) {
  stop("No reviewer workbook files found in data/final-scores/")
}

##------------------------------------------------------------------------------
## Functions

standardize_stock_names <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("[\u2018\u2019]", "'") %>%
    str_replace_all("[\u201C\u201D]", "\"") %>%
    str_replace_all("[\u2013\u2014]", "-") %>%
    str_replace_all("_", " ") %>%
    str_squish()
}

standardize_attribute_names <- function(x) {
  x_std <- x %>%
    as.character() %>%
    str_replace_all("[\u2018\u2019]", "'") %>%
    str_replace_all("[\u201C\u201D]", "\"") %>%
    str_replace_all("[\u2013\u2014]", "-") %>%
    str_replace_all("_", " ") %>%
    str_squish()

  dplyr::case_when(
    x_std == "Stock Size/Status" ~ "Stock Size Status",
    x_std == "Stock size/status" ~ "Stock Size Status",
    TRUE ~ x_std
  )
}

parse_reviewer_id <- function(file_path) {
  ## Purpose:
  ## - Create a reviewer ID from the workbook filename
  ##
  ## Notes:
  ## - Replace this if your file-naming convention supports a better parser

  basename(file_path) %>%
    tools::file_path_sans_ext() %>%
    str_squish()
}

get_stock_sheets <- function(wb) {
  ## Purpose:
  ## - Return stock sheets only, excluding instruction/example sheets
  
  setdiff(
    names(wb),
    c("Instructions", "Data Quality", "Example")
  )
}

safe_numeric <- function(x) {
  ## Purpose:
  ## - Convert blank strings and NULL-like values to NA, otherwise numeric

  if(is.null(x)) return(NA_real_)

  x_chr <- as.character(x)

  if(length(x_chr) == 0 || is.na(x_chr) || str_squish(x_chr) == "") {
    return(NA_real_)
  }

  suppressWarnings(as.numeric(x_chr))
}

safe_read_cell <- function(wb, sheet_name, row, col) {
  ## Purpose:
  ## - Read a single cell from a workbook and return its value as a character
  ## - Returns NA_character_ when readWorkbook finds no data in the range
  ##
  ## Notes:
  ## - readWorkbook returns NULL (not an empty data frame) for fully empty rows,
  ##   which causes as.matrix() to call array(data = NULL, ...) and crash.
  ##   This wrapper catches that case before it reaches as.matrix().

  result <- suppressWarnings(
    readWorkbook(wb, sheet = sheet_name, rows = row, cols = col,
                 colNames = FALSE)
  )

  if(is.null(result) || nrow(result) == 0 || ncol(result) == 0) {
    return(NA_character_)
  }

  val <- result[[1]][[1]]

  if(is.null(val) || length(val) == 0) return(NA_character_)

  as.character(val)
}


extract_rank_tally_block <- function(wb,
                                     sheet_name,
                                     rows,
                                     attribute_col,
                                     attribute_type,
                                     reviewer_id,
                                     source_file) {
  ## Purpose:
  ## - Extract one block of FINAL SCORE rank tallies where rows correspond to
  ##   attributes and columns M:P correspond to final rank tallies 1:4
  ##
  ## FINAL SCORE layout:
  ## - K = final attribute score
  ## - L = final data quality index
  ## - M = rank 1 tally
  ## - N = rank 2 tally
  ## - O = rank 3 tally
  ## - P = rank 4 tally
  
  out <- map_dfr(
    rows,
    function(r) {

      attribute_name     <- safe_read_cell(wb, sheet_name, r, attribute_col)
      final_score        <- safe_read_cell(wb, sheet_name, r, 11)
      data_quality_index <- safe_read_cell(wb, sheet_name, r, 12)

      ## Read each tally column individually to preserve absolute column
      ## alignment. readWorkbook silently drops empty columns when reading
      ## a range, shifting remaining values left. Reading one cell at a time
      ## pins each value to its correct absolute column position.
      tally_L  <- safe_numeric(safe_read_cell(wb, sheet_name, r, 13))
      tally_M  <- safe_numeric(safe_read_cell(wb, sheet_name, r, 14))
      tally_H  <- safe_numeric(safe_read_cell(wb, sheet_name, r, 15))
      tally_VH <- safe_numeric(safe_read_cell(wb, sheet_name, r, 16))
      
      tibble(
        source_file = basename(source_file),
        reviewer_id = reviewer_id,
        stock_name = standardize_stock_names(sheet_name),
        sheet_name = sheet_name,
        attribute_type = attribute_type,
        attribute_name = standardize_attribute_names(attribute_name),
        row_num = r,
        final_score = safe_numeric(final_score),
        data_quality_index = safe_numeric(data_quality_index),
        tally_L = tally_L,
        tally_M = tally_M,
        tally_H = tally_H,
        tally_VH = tally_VH,
        n_tallies = sum(c(tally_L,
                          tally_M,
                          tally_H,
                          tally_VH), na.rm = TRUE)
      )
    }
  )
  
  out
}

extract_directional_effect_block <- function(wb,
                                             sheet_name,
                                             reviewer_id,
                                             source_file) {
  ## Purpose:
  ## - Extract final directional effect tallies from rows 38:40
  ##
  ## FINAL SCORE layout:
  ## - effect labels in column B
  ## - final tallies in column M
  ##
  ## Rows:
  ## - 38 = Positive
  ## - 39 = Neutral
  ## - 40 = Negative
  
  rows <- 38:40
  
  map_dfr(
    rows,
    function(r) {
      
      effect_category <- safe_read_cell(wb, sheet_name, r, 2)
      tally_val       <- safe_read_cell(wb, sheet_name, r, 13)
      
      tibble(
        source_file = basename(source_file),
        reviewer_id = reviewer_id,
        stock_name = standardize_stock_names(sheet_name),
        sheet_name = sheet_name,
        effect_category = str_squish(as.character(effect_category)),
        row_num = r,
        tally = safe_numeric(tally_val)
      )
    }
  ) %>%
    mutate(
      n_tallies = sum(tally, na.rm = TRUE)
    )
}

##------------------------------------------------------------------------------
## Step 1 - Loop through reviewer workbooks and extract tally blocks

sensitivity_tallies_long <- tibble()
qualitative_exposure_tallies_long <- tibble()
directional_effect_tallies_long <- tibble()

for(file_i in workbook_files) {
  
  message("Reading workbook: ", basename(file_i))
  
  reviewer_id_i <- parse_reviewer_id(file_i)
  wb_i <- loadWorkbook(file_i)
  
  stock_sheets_i <- get_stock_sheets(wb_i)
  
  for(sheet_i in stock_sheets_i) {
    
    message("  Extracting sheet: ", sheet_i)
    
    ##--------------------------------------------------------------------------
    ## Qualitative Exposure Factors
    ##
    ## Rows 17:19
    ## Attribute names in column C
    ## Final rank tallies in M:P
    
    qexp_i <- extract_rank_tally_block(
      wb = wb_i,
      sheet_name = sheet_i,
      rows = 17:19,
      attribute_col = 3,
      attribute_type = "Qualitative Exposure",
      reviewer_id = reviewer_id_i,
      source_file = file_i
    )

    ##--------------------------------------------------------------------------
    ## Sensitivity Attributes
    ##
    ## Rows 21:28
    ## Attribute names in column B
    ## Final rank tallies in M:P

    sens_i <- extract_rank_tally_block(
      wb = wb_i,
      sheet_name = sheet_i,
      rows = 21:28,
      attribute_col = 2,
      attribute_type = "Sensitivity",
      reviewer_id = reviewer_id_i,
      source_file = file_i
    )

    ##--------------------------------------------------------------------------
    ## Rigidity Attributes
    ##
    ## Rows 30:35
    ## Attribute names in column B
    ## Final rank tallies in M:P
    ##
    ## Note:
    ## - To align with the current Caribbean workflow, we relabel these as
    ##   Sensitivity because they are part of the overall qualitative
    ##   sensitivity-side component used in vulnerability scoring

    rigidity_i <- extract_rank_tally_block(
      wb = wb_i,
      sheet_name = sheet_i,
      rows = 30:35,
      attribute_col = 2,
      attribute_type = "Sensitivity",
      reviewer_id = reviewer_id_i,
      source_file = file_i
    )

    ##--------------------------------------------------------------------------
    ## Directional Effect
    ##
    ## Rows 38:40
    ## Labels in column B
    ## Final tallies in column M

    directional_i <- extract_directional_effect_block(
      wb = wb_i,
      sheet_name = sheet_i,
      reviewer_id = reviewer_id_i,
      source_file = file_i
    )
    
    ##--------------------------------------------------------------------------
    ## Append to running long tables
    
    qualitative_exposure_tallies_long <- bind_rows(
      qualitative_exposure_tallies_long,
      qexp_i
    )
    
    sensitivity_tallies_long <- bind_rows(
      sensitivity_tallies_long,
      sens_i,
      rigidity_i
    )
    
    directional_effect_tallies_long <- bind_rows(
      directional_effect_tallies_long,
      directional_i
    )
  }
}

message("Extraction complete.")
message("Sensitivity tally rows: ", nrow(sensitivity_tallies_long))
message("Qualitative exposure tally rows: ", nrow(qualitative_exposure_tallies_long))
message("Directional effect tally rows: ", nrow(directional_effect_tallies_long))


##------------------------------------------------------------------------------
## Step 2 - Cleanup
##
## Workflow:
## - Remove rows with missing attribute names
## - Arrange consistently
## - Remove stock x reviewer combinations where all tallies are zero
##   (reviewer sheet was read but no scores were entered for that stock)
## - Create combined qualitative tally table

qualitative_exposure_tallies_long <- qualitative_exposure_tallies_long %>%
  filter(!is.na(attribute_name), attribute_name != "") %>%
  arrange(stock_name, reviewer_id, row_num)

sensitivity_tallies_long <- sensitivity_tallies_long %>%
  filter(!is.na(attribute_name), attribute_name != "") %>%
  arrange(stock_name, reviewer_id, row_num)

directional_effect_tallies_long <- directional_effect_tallies_long %>%
  filter(!is.na(effect_category), effect_category != "") %>%
  arrange(stock_name, reviewer_id, row_num)

##------------------------------------------------------------------------------
## Step 2B - Remove unscored stock x reviewer combinations
##
## The loop reads every sheet in every workbook, producing rows for reviewer x
## stock combinations where the reviewer did not enter any scores. These appear
## as rows where every tally bin is NA and n_tallies == 0. Identify those
## combinations and remove them from all tables before writing output.

unscored_sens <- sensitivity_tallies_long %>%
  group_by(stock_name, reviewer_id) %>%
  summarise(total_tallies = sum(n_tallies, na.rm = TRUE), .groups = "drop") %>%
  filter(total_tallies == 0)

unscored_qexp <- qualitative_exposure_tallies_long %>%
  group_by(stock_name, reviewer_id) %>%
  summarise(total_tallies = sum(n_tallies, na.rm = TRUE), .groups = "drop") %>%
  filter(total_tallies == 0)

unscored_dir <- directional_effect_tallies_long %>%
  group_by(stock_name, reviewer_id) %>%
  summarise(total_tallies = sum(tally, na.rm = TRUE), .groups = "drop") %>%
  filter(total_tallies == 0)

message("Removing ", nrow(unscored_sens),
        " unscored stock x reviewer combinations from sensitivity tallies")
message("Removing ", nrow(unscored_qexp),
        " unscored stock x reviewer combinations from qualitative exposure tallies")
message("Removing ", nrow(unscored_dir),
        " unscored stock x reviewer combinations from directional effect tallies")

sensitivity_tallies_long <- sensitivity_tallies_long %>%
  anti_join(unscored_sens, by = c("stock_name", "reviewer_id"))

qualitative_exposure_tallies_long <- qualitative_exposure_tallies_long %>%
  anti_join(unscored_qexp, by = c("stock_name", "reviewer_id"))

directional_effect_tallies_long <- directional_effect_tallies_long %>%
  anti_join(unscored_dir, by = c("stock_name", "reviewer_id"))

all_qualitative_tallies_long <- bind_rows(
  qualitative_exposure_tallies_long,
  sensitivity_tallies_long
) %>%
  arrange(stock_name, reviewer_id, attribute_type, row_num)

print(all_qualitative_tallies_long)

##------------------------------------------------------------------------------
## Step 3 - QA
##
## Workflow:
## Row-level tally checks
##   - Flag rows where n_tallies != 5 (each reviewer x attribute must use
##     exactly 5 tallies under the NOAA FCVA scoring rules)
##   - Flag rows with any NA tally bins
##   - Flag duplicate stock x reviewer x attribute rows
## Reviewer x stock coverage matrix
##   - Count distinct attributes scored per reviewer x stock combination
##   - A reviewer absent from a stock shows as 0 in the wide output
## Attribute-level pooled tally summary
##   - Pool tallies across all reviewers for each stock x attribute
##   - Calculate pooled weighted mean score
##   - Flag combinations where pooled sum != n_reviewers x 5
##   - This is the key diagnostic for the bootstrap analysis in Script 2
## Directional effect summary
##   - Summarize total tally counts per stock x effect category

##------------------------------------------------------------------------------
## Step 3A - Row-level tally checks

qa_row_checks <- all_qualitative_tallies_long %>%
  mutate(
    n_tallies_ok  = n_tallies == 5,
    has_na_tally  = is.na(tally_L) | is.na(tally_M) |
      is.na(tally_H) | is.na(tally_VH),
    any_row_issue = !n_tallies_ok | has_na_tally
  ); print(qa_row_checks)

qa_duplicates <- all_qualitative_tallies_long %>%
  count(stock_name, reviewer_id, attribute_name, name = "n_rows") %>%
  filter(n_rows > 1); print(qa_duplicates)

n_bad_tally_sum <- sum(!qa_row_checks$n_tallies_ok, na.rm = TRUE)
n_na_tally      <- sum(qa_row_checks$has_na_tally,  na.rm = TRUE)
n_duplicates    <- nrow(qa_duplicates)

message("Row-level QA: ", n_bad_tally_sum, " rows with n_tallies != 5")
message("Row-level QA: ", n_na_tally,      " rows with NA tally values")
message("Row-level QA: ", n_duplicates,    " duplicate stock x reviewer x attribute rows")

#write_csv(qa_row_checks %>% filter(any_row_issue), f_qa_row_checks)

##------------------------------------------------------------------------------
## Step 3B - Reviewer x stock coverage matrix

reviewer_stock_coverage <- all_qualitative_tallies_long %>%
  group_by(reviewer_id, stock_name) %>%
  summarise(
    n_attributes_scored = n_distinct(attribute_name),
    n_rows              = n(),
    .groups = "drop"
  ); print(reviewer_stock_coverage, n = 100)

reviewer_stock_wide <- reviewer_stock_coverage %>%
  pivot_wider(
    id_cols     = stock_name,
    names_from  = reviewer_id,
    values_from = n_attributes_scored,
    values_fill = 0L
  ) %>%
  arrange(stock_name); print(reviewer_stock_wide)

#write_csv(reviewer_stock_wide, f_qa_coverage)

##------------------------------------------------------------------------------
## Step 3C - Attribute-level pooled tally summary

qa_attr_summary <- all_qualitative_tallies_long %>%
  group_by(stock_name, attribute_type, attribute_name) %>%
  summarise(
    n_reviewers        = n_distinct(reviewer_id),
    tally_L            = sum(tally_L,       na.rm = TRUE),
    tally_M            = sum(tally_M,  na.rm = TRUE),
    tally_H            = sum(tally_H,      na.rm = TRUE),
    tally_VH            = sum(tally_VH, na.rm = TRUE),
    pooled_tally_sum   = tally_L + tally_M + tally_H + tally_VH,
    expected_tally_sum = n_reviewers * 5L,
    tally_sum_ok       = pooled_tally_sum == expected_tally_sum,
    pooled_mean_score  = (tally_L * 1 + tally_M * 2 +
                            tally_H * 3 + tally_VH * 4) / pooled_tally_sum,
    .groups = "drop"
  ) %>%
  arrange(stock_name, attribute_type, attribute_name)
print(qa_attr_summary, n = 50)

n_tally_sum_issues <- sum(!qa_attr_summary$tally_sum_ok, na.rm = TRUE)
message("Attribute QA: ", n_tally_sum_issues,
        " stock x attribute combinations with pooled tally sum != n_reviewers x 5")

write_csv(qa_attr_summary, f_qa_attr_summary)

##------------------------------------------------------------------------------
## Step 3D - Directional effect summary

qa_dir_summary <- directional_effect_tallies_long %>%
  group_by(stock_name, effect_category) %>%
  summarise(
    n_reviewers = n_distinct(reviewer_id),
    total_tally = sum(tally, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(stock_name, effect_category)
print(qa_dir_summary, n = 75)

write_csv(qa_dir_summary, f_qa_dir_summary)

##------------------------------------------------------------------------------
## Step 4 - Write output tables

write_csv(sensitivity_tallies_long,           f_sens_tallies)
write_csv(qualitative_exposure_tallies_long,  f_qexp_tallies)
write_csv(directional_effect_tallies_long,    f_dir_tallies)
write_csv(all_qualitative_tallies_long,       f_all_tallies)


