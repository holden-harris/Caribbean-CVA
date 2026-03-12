##------------------------------------------------------------------------------
## Caribbean CVA – Compile reviewer Directional Effect scores into one table
## Builds: directional_effect_table_all_clean
##
## Reads final directional-effect scores already computed in each workbook:
##   M38:M40 = Positive / Neutral / Negative
##
## FIX:
## Some tabs return 0 columns for M38:M40. This version safely skips those tabs
## instead of crashing the entire workbook import.

rm(list = ls()); gc()
suppressPackageStartupMessages({
  library(readxl)
  library(xml2)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(tibble)
})

##------------------------------------------------------------------------------
## CONFIG
in_dir      <- "./data/final-scores"
ignore_tabs <- c("Instructions","Data Quality","Example")

##------------------------------------------------------------------------------
## HELPERS

## Parse Scorer from filename
parse_scorer <- function(path){
  fn <- basename(path)
  fn_noext <- sub("\\.[^.]+$", "", fn)
  parts <- strsplit(fn_noext, "_")[[1]]
  scorer <- parts[length(parts)]
  scorer <- stringr::str_squish(scorer)
  ifelse(nchar(scorer) > 0, scorer, "Unknown Scorer")
}

## Get visible sheet names WITHOUT unzipping
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

##------------------------------------------------------------------------------
## NEW: safely read a range and always return either a tibble or empty tibble
safe_read_range <- function(path, sheet, range){
  tryCatch({
    out <- readxl::read_excel(
      path      = path,
      sheet     = sheet,
      range     = range,
      col_names = FALSE,
      col_types = "text"
    )
    
    ## Some tabs may return 0 columns for an effectively empty range
    if (ncol(out) < 1) {
      return(tibble::tibble())
    }
    
    out
  }, error = function(e){
    tibble::tibble()
  })
}

##------------------------------------------------------------------------------
## NEW: helper to safely pull numeric values from M38:M40
read_directional_scores <- function(xlsx_path, tab){
  
  x <- safe_read_range(
    path  = xlsx_path,
    sheet = tab,
    range = "M38:M40"
  )
  
  if (ncol(x) < 1) {
    return(rep(NA_real_, 3))
  }
  
  vals <- suppressWarnings(as.numeric(x[[1]]))
  
  ## pad to length 3 if needed
  if (length(vals) < 3) {
    vals <- c(vals, rep(NA_real_, 3 - length(vals)))
  }
  
  vals[1:3]
}

##------------------------------------------------------------------------------
## Extract directional effect rows for one tab
extract_directional_effect_rows <- function(xlsx_path, tab, Scorer){
  
  vals <- read_directional_scores(xlsx_path, tab)
  
  out <- tibble::tibble(
    SourceFile         = basename(xlsx_path),
    Scorer             = Scorer,
    stock_name         = tab,
    row_idx            = 38:40,
    Directional_effect = c("Positive", "Neutral", "Negative"),
    Directional_score  = vals
  ) |>
    dplyr::filter(!is.na(Directional_score))
  
  out
}

##------------------------------------------------------------------------------
## Helper: determine whether a sheet has any directional-effect scores
sheet_has_directional_scores <- function(xlsx_path, tab){
  vals <- read_directional_scores(xlsx_path, tab)
  any(!is.na(vals))
}

##------------------------------------------------------------------------------
## Safe loader: copy to temp (breaks OneDrive/Excel locks), no unzip, resilient
## FIX: skip bad tabs rather than failing the whole workbook
build_directional_effect_table_for_workbook <- function(xlsx_path){
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
        if (sheet_has_directional_scores(tmp, tab)) {
          extract_directional_effect_rows(tmp, tab, Scorer)
        } else {
          ## optional message for diagnosis
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
  out <- build_directional_effect_table_for_workbook(f)
  message("  rows: ", nrow(out), "  (scorer: ", parse_scorer(f), ")")
  out
})

directional_effect_table_all <- dplyr::bind_rows(per_file) |>
  dplyr::arrange(Scorer, stock_name, row_idx)

message("Total rows: ", nrow(directional_effect_table_all),
        " | scorers: ", dplyr::n_distinct(directional_effect_table_all$Scorer),
        " | stocks: ",  dplyr::n_distinct(directional_effect_table_all$stock_name))

##------------------------------------------------------------------------------
## QA/QC

n_scorers <- directional_effect_table_all %>% distinct(Scorer) %>% nrow()
n_stocks  <- directional_effect_table_all %>% distinct(stock_name) %>% nrow()
n_levels  <- directional_effect_table_all %>% distinct(Directional_effect) %>% nrow()

n_scorers; n_stocks; n_levels

## Per stock: how many scorers, how many directional categories populated
stock_qa <- directional_effect_table_all %>%
  group_by(stock_name) %>%
  summarise(
    n_scorers = n_distinct(Scorer),
    n_levels  = n_distinct(Directional_effect),
    .groups = "drop"
  ) %>%
  arrange(stock_name) %>% as.data.frame(); stock_qa

## Scorers per stock
species_reviews <- directional_effect_table_all %>%
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

## Write out reviewers
write.csv(
  species_reviews,
  file = file.path(in_dir, "n_reviews_directional-effect.csv"),
  row.names = FALSE
)

##------------------------------------------------------------------------------
## Write to CSV for downstream steps
write.csv(
  directional_effect_table_all,
  file = file.path(in_dir, "directional_effect_table_all.csv"),
  row.names = FALSE
)