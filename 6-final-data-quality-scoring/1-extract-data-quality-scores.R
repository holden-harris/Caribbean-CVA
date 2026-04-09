## -----------------------------------------------------------------------------
## Compile data quality scores from final scoring workbooks
##
## Target cells:
## - L17:L18
## - L21:L28
##
## Workflow:
## 1. Read each workbook
## 2. Extract the target rows from column L
## 3. Pull the row label from column B for each target row
## 4. Build a long table of individual data quality scores
## 5. Summarize by stock x attribute
## 6. Summarize overall data quality by stock
##
## Notes:
## - This assumes blanks mean no score / not available
## - Scores ranged from 0 to 3, 
##   where 3 = Adequate Data, 2 = Limited Data, 1 = Expert Judgment, and 0 = No Data. 
## - Following the NOAA CVA approach, overall data quality was summarized as 
##   the proportion of scores greater than or equal to 2, 
##   with High defined as ≥80%, Moderate as 50–79%, and Poor as <50%.
## -----------------------------------------------------------------------------

rm(list = ls()); gc()

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(tibble)
})

##------------------------------------------------------------------------------
## Set up

in_dir      <- "./data/final-scores"
out_dir     <- "./outputs/final-scores-compiled/data-quality"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ignore_tabs <- c("Instructions", "Data Quality", "Example")

dq_rows <- c(17:18, 21:28)

##------------------------------------------------------------------------------
## Helper functions

## Parse scorer initials from filename
parse_scorer <- function(path){
  fn <- basename(path)
  fn_noext <- sub("\\.[^.]+$", "", fn)
  parts <- strsplit(fn_noext, "_")[[1]]
  scorer <- parts[length(parts)]
  scorer <- stringr::str_squish(scorer)
  scorer <- substr(scorer, 1, 2)
  ifelse(nchar(scorer) > 0, scorer, "Unknown Scorer")
}

## Get visible sheet names without manual unzip
safe_visible_sheet_names <- function(xlsx_path) {
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    sn  <- openxlsx::getSheetNames(xlsx_path)
    vis <- tryCatch(
      openxlsx::getSheetVisibility(xlsx_path),
      error = function(e) rep("visible", length(sn))
    )
    keep <- is.na(vis) | tolower(vis) == "visible"
    return(trimws(sn[keep]))
  }
  readxl::excel_sheets(xlsx_path)
}

## Safely read a range; return empty tibble on failure
safe_read_range <- function(path, sheet, range){
  tryCatch({
    out <- readxl::read_excel(
      path      = path,
      sheet     = sheet,
      range     = range,
      col_names = FALSE,
      col_types = "text"
    )
    
    if (ncol(out) < 1) {
      return(tibble::tibble())
    }
    
    out
  }, error = function(e){
    tibble::tibble()
  })
}

##------------------------------------------------------------------------------
## Read data-quality score cells from column L
##
## Reads final data-quality scores already entered in each workbook:
##   L17:L18 and L21:L28
##
## This follows the same robust pattern as the directional-effect extractor:
## - safely copies workbook to temp
## - reads only visible sheets
## - skips empty/problem tabs
## - returns one long table with one row per stock x scorer x attribute score
##
## Notes:
## - Attribute labels are pulled from column B in the same rows
## - If a label cell is blank, a fallback label is used
## - Scores are expected to be 0, 1, 2, or 3


read_data_quality_scores <- function(xlsx_path, tab){
  
  x1 <- safe_read_range(
    path  = xlsx_path,
    sheet = tab,
    range = "L17:L18"
  )
  
  x2 <- safe_read_range(
    path  = xlsx_path,
    sheet = tab,
    range = "L21:L28"
  )
  
  vals1 <- if (ncol(x1) >= 1) suppressWarnings(as.numeric(x1[[1]])) else numeric(0)
  vals2 <- if (ncol(x2) >= 1) suppressWarnings(as.numeric(x2[[1]])) else numeric(0)
  
  if (length(vals1) < 2) vals1 <- c(vals1, rep(NA_real_, 2 - length(vals1)))
  if (length(vals2) < 8) vals2 <- c(vals2, rep(NA_real_, 8 - length(vals2)))
  
  vals <- c(vals1[1:2], vals2[1:8])
  vals
}

## Read matching row labels from column B
read_data_quality_labels <- function(xlsx_path, tab){
  
  x1 <- safe_read_range(
    path  = xlsx_path,
    sheet = tab,
    range = "C17:C18"
  )
  
  x2 <- safe_read_range(
    path  = xlsx_path,
    sheet = tab,
    range = "B21:B28"
  )
  
  lab1 <- if (ncol(x1) >= 1) as.character(x1[[1]]) else character(0)
  lab2 <- if (ncol(x2) >= 1) as.character(x2[[1]]) else character(0)
  
  if (length(lab1) < 2) lab1 <- c(lab1, rep(NA_character_, 2 - length(lab1)))
  if (length(lab2) < 8) lab2 <- c(lab2, rep(NA_character_, 8 - length(lab2)))
  
  labs <- c(lab1[1:2], lab2[1:8])
  labs <- stringr::str_squish(labs)
  
  fallback <- paste0("Row_", dq_rows)
  labs[is.na(labs) | labs == ""] <- fallback[is.na(labs) | labs == ""]
  
  labs
}

##------------------------------------------------------------------------------
## Extract one tab
extract_data_quality_rows <- function(xlsx_path, tab, Scorer){
  
  vals <- read_data_quality_scores(xlsx_path, tab)
  labs <- read_data_quality_labels(xlsx_path, tab)
  
  out <- tibble::tibble(
    SourceFile          = basename(xlsx_path),
    Scorer              = Scorer,
    stock_name          = tab,
    row_idx             = dq_rows,
    Attribute_type      = c(rep("Exposure", 2), rep("Sensitivity", 8)),
    Attribute_name      = labs,
    Data_quality_score  = vals
  ) |>
    dplyr::mutate(
      score_cell = paste0("L", row_idx),
      label_cell = paste0("B", row_idx),
      Data_quality_label = dplyr::case_when(
        Data_quality_score == 3 ~ "Adequate Data",
        Data_quality_score == 2 ~ "Limited Data",
        Data_quality_score == 1 ~ "Expert Judgment",
        Data_quality_score == 0 ~ "No Data",
        TRUE ~ NA_character_
      ),
      score_ge_2 = Data_quality_score >= 2
    ) |>
    dplyr::filter(!is.na(Data_quality_score))
  
  out
}

##------------------------------------------------------------------------------
## Helper: does sheet contain any data-quality scores?
sheet_has_data_quality_scores <- function(xlsx_path, tab){
  vals <- read_data_quality_scores(xlsx_path, tab)
  any(!is.na(vals))
}

##------------------------------------------------------------------------------
## Safe workbook loader
build_data_quality_table_for_workbook <- function(xlsx_path){
  tryCatch({
    tmp <- tempfile(fileext = ".xlsx")
    ok  <- file.copy(xlsx_path, tmp, overwrite = TRUE)
    
    if (!ok) {
      message("⚠️  Could not copy (locked?): ", basename(xlsx_path))
      return(tibble::tibble())
    }
    
    Scorer <- parse_scorer(xlsx_path)
    tabs   <- safe_visible_sheet_names(tmp)
    stock_tabs <- setdiff(tabs, ignore_tabs)
    
    if (!length(stock_tabs)) {
      message("ℹ️  No stock tabs in: ", basename(xlsx_path))
      return(tibble::tibble())
    }
    
    per_tab <- lapply(stock_tabs, function(tab){
      tryCatch({
        if (sheet_has_data_quality_scores(tmp, tab)) {
          extract_data_quality_rows(tmp, tab, Scorer)
        } else {
          message("  - Skipping unscored/empty tab: ", tab)
          tibble::tibble()
        }
      }, error = function(e){
        message("  ⚠️  Skipping tab due to read issue: ", tab,
                " | ", conditionMessage(e))
        tibble::tibble()
      })
    })
    
    dplyr::bind_rows(per_tab)
    
  }, error = function(e){
    message("❌ Error in ", basename(xlsx_path), ": ", conditionMessage(e))
    tibble::tibble()
  })
}

##------------------------------------------------------------------------------
## COMPILE

xlsx_files <- list.files(
  in_dir,
  pattern = "(?i)\\.xl[a-z]+$",
  full.names = TRUE
)

message("Found files: ", length(xlsx_files))

per_file <- lapply(xlsx_files, function(f){
  message("→ Reading: ", basename(f))
  out <- build_data_quality_table_for_workbook(f)
  message("  rows: ", nrow(out), "  (scorer: ", parse_scorer(f), ")")
  out
})

data_quality_table_all <- dplyr::bind_rows(per_file) |>
  dplyr::arrange(Scorer, stock_name, row_idx)

message("Total rows: ", nrow(data_quality_table_all),
        " | scorers: ", dplyr::n_distinct(data_quality_table_all$Scorer),
        " | stocks: ",  dplyr::n_distinct(data_quality_table_all$stock_name))

##------------------------------------------------------------------------------
## QA/QC

n_scorers <- data_quality_table_all %>% distinct(Scorer) %>% nrow()
n_stocks  <- data_quality_table_all %>% distinct(stock_name) %>% nrow()
n_rows    <- data_quality_table_all %>% distinct(row_idx) %>% nrow()

n_scorers; n_stocks; n_rows

## Per stock: how many scorers, how many scored rows
stock_qa <- data_quality_table_all %>%
  group_by(stock_name) %>%
  summarise(
    n_scorers = n_distinct(Scorer),
    n_rows    = n_distinct(row_idx),
    .groups = "drop"
  ) %>%
  arrange(stock_name) %>%
  as.data.frame(); stock_qa

## Scorers per stock
species_reviews <- data_quality_table_all %>%
  filter(!is.na(stock_name), !is.na(Scorer)) %>%
  distinct(stock_name, Scorer) %>%
  group_by(stock_name) %>%
  summarise(
    n_reviews = n_distinct(Scorer),
    reviewers = paste(sort(unique(Scorer)), collapse = ", "),
    .groups   = "drop"
  ) %>%
  arrange(desc(n_reviews), stock_name) %>%
  as.data.frame(); species_reviews

##------------------------------------------------------------------------------
## Write out
write.csv(
  species_reviews,
  file = file.path(out_dir, "qa_data_quality_by_scorers.csv"),
  row.names = FALSE
)

write.csv(
  stock_qa,
  file = file.path(out_dir, "qa_data_quality_by_stock.csv"),
  row.names = FALSE
)

write.csv(
  data_quality_table_all,
  file = file.path(out_dir, "table_data_quality_scores_extracted.csv"),
  row.names = FALSE
)

