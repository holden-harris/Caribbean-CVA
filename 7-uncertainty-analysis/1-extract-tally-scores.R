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

##------------------------------------------------------------------------------
## Directories

proj_dir <- "."
in_dir   <- file.path(proj_dir, "data", "final-scores")
out_dir  <- file.path(proj_dir, "outputs", "analyses", "1-inputs")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

safe_read_tally_row <- function(wb, sheet_name, row, cols) {
  ## Purpose:
  ## - Read a row of tally cells and return a named list of numeric values
  ## - Returns NA for each column when readWorkbook finds no data
  ##
  ## Notes:
  ## - Same NULL-return problem as safe_read_cell applies here

  result <- suppressWarnings(
    readWorkbook(wb, sheet = sheet_name, rows = row, cols = cols,
                 colNames = FALSE)
  )

  n_cols <- length(cols)

  if(is.null(result) || nrow(result) == 0) {
    return(rep(list(NA_real_), n_cols))
  }

  lapply(seq_len(n_cols), function(i) {
    if(i > ncol(result)) return(NA_real_)
    safe_numeric(result[[i]])
  })
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

      tallies <- safe_read_tally_row(wb, sheet_name, r, 13:16)

      tally_low       <- tallies[[1]]
      tally_moderate  <- tallies[[2]]
      tally_high      <- tallies[[3]]
      tally_very_high <- tallies[[4]]
      
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
        tally_low = tally_low,
        tally_moderate = tally_moderate,
        tally_high = tally_high,
        tally_very_high = tally_very_high,
        n_tallies = sum(c(tally_low,
                          tally_moderate,
                          tally_high,
                          tally_very_high), na.rm = TRUE)
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

##------------------------------------------------------------------------------
## Step 2 - Final cleanup
##
## Workflow:
## - Remove rows with missing attribute names if any
## - Arrange consistently
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

all_qualitative_tallies_long <- bind_rows(
  qualitative_exposure_tallies_long,
  sensitivity_tallies_long
) %>%
  arrange(stock_name, reviewer_id, attribute_type, row_num)

##------------------------------------------------------------------------------
## Step 3 - Write output tables

write_csv(sensitivity_tallies_long,           f_sens_tallies)
write_csv(qualitative_exposure_tallies_long,  f_qexp_tallies)
write_csv(directional_effect_tallies_long,    f_dir_tallies)
write_csv(all_qualitative_tallies_long,       f_all_tallies)

##------------------------------------------------------------------------------
## Step 4 - Basic extraction QA summaries
##
## Notes:
## - These are quick console checks only
## - Full pre-analysis QA would happen later in Script 1

message("Extraction complete.")
message("Sensitivity tally rows: ", nrow(sensitivity_tallies_long))
message("Qualitative exposure tally rows: ", nrow(qualitative_exposure_tallies_long))
message("Directional effect tally rows: ", nrow(directional_effect_tallies_long))