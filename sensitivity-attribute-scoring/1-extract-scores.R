##------------------------------------------------------------------------------
## Caribbean CVA – Compile reviewer scoring workbooks into one table
## Builds: score_table_all  (rows from rows 21–35, skipping 29; B/E/F–I)
##------------------------------------------------------------------------------

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
##------------------------------------------------------------------------------
in_dir      <- "./data/preliminary-scores"    ## folder containing the 16 Excel files
ignore_tabs <- c("Instructions","Data Quality","Example")

## Canonical attribute sets (for Attribute_type tagging)
sens_attrs <- c(
  "Habitat specificity",
  "Prey specificity",
  "Tolerance to ocean acidification",
  "Complexity in reproductive strategy",
  "Species range",
  "Specificity in early life history requirements",
  "Stock size/status",
  "Other stressors"
)

rigid_attrs <- c(
  "Population growth rate",
  "Mobility and dispersal or early life stages",
  "Adult mobility",
  "Spawning characteristics",
  "Predation and competition dynamics",
  "Genetic diversity"
)

##------------------------------------------------------------------------------
## HELPERS
##------------------------------------------------------------------------------

## Parse Scorer from filename: last '_' token before extension
parse_scorer <- function(path){
  fn <- basename(path)
  fn_noext <- sub("\\.[^.]+$", "", fn)
  parts <- strsplit(fn_noext, "_")[[1]]
  scorer <- parts[length(parts)]
  scorer <- str_squish(scorer)
  ifelse(nchar(scorer) > 0, scorer, "Unknown Scorer")
}

## Get visible (not hidden/veryHidden) sheet names by reading xl/workbook.xml
visible_sheet_names <- function(xlsx_path){
  tmp_dir <- tempfile("xlsx_unzip_"); dir.create(tmp_dir)
  unzip(xlsx_path, exdir = tmp_dir)
  wb_xml <- file.path(tmp_dir, "xl", "workbook.xml")
  if (!file.exists(wb_xml)) {
    ## Fall back to readxl if not a standard .xlsx (e.g., old .xls)
    return(readxl::excel_sheets(xlsx_path))
  }
  doc <- xml2::read_xml(wb_xml)
  ns  <- xml2::xml_ns(doc)
  nodes <- xml2::xml_find_all(doc, ".//d:sheets/d:sheet", ns = c(d = ns[[1]]))
  nm   <- xml2::xml_attr(nodes, "name")
  st   <- xml2::xml_attr(nodes, "state")     ## NA = visible; "hidden"/"veryHidden" otherwise
  vis  <- nm[is.na(st)]
  trimws(vis)
}

## Extract the 21–35 rows (skip 29) from B..I for a single stock tab
extract_tab_rows <- function(xlsx_path, tab, Scorer){
  ## Read B21:I35 as text
  rng <- readxl::read_excel(
    path      = xlsx_path,
    sheet     = tab,
    range     = "B21:I35",
    col_names = FALSE,
    col_types = "text"
  )
  
  if (ncol(rng) < 8) return(tibble())  ## must have B..I
  
  rng <- rng |>
    mutate(row_idx = 21:35) |>
    filter(row_idx != 29)
  
  out <- rng |>
    transmute(
      SourceFile      = basename(xlsx_path),
      Scorer          = Scorer,
      stock_name      = tab,
      row_idx         = row_idx,
      Attribute_name  = str_squish(as.character(`...1`)),
      Data_quality    = suppressWarnings(as.numeric(`...4`)),
      Scoring_rank_1  = suppressWarnings(as.numeric(`...5`)),
      Scoring_rank_2  = suppressWarnings(as.numeric(`...6`)),
      Scoring_rank_3  = suppressWarnings(as.numeric(`...7`)),
      Scoring_rank_4  = suppressWarnings(as.numeric(`...8`))
    ) |>
    ## Tag attribute type
    mutate(
      Attribute_type = case_when(
        Attribute_name %in% sens_attrs  ~ "Sensitivity",
        Attribute_name %in% rigid_attrs ~ "Rigidity",
        TRUE                            ~ NA_character_
      )
    ) |>
    ## drop rows that are completely empty
    filter(!(is.na(Attribute_name) & if_all(starts_with("Scoring_rank_"), ~is.na(.x))))
  
  out
}

## Process a single workbook → one tibble like your rows_21_35_all
build_score_table_for_workbook <- function(xlsx_path){
  Scorer <- parse_scorer(xlsx_path)
  
  ## visible sheets minus admin tabs
  all_tabs_vis <- visible_sheet_names(xlsx_path)
  all_tabs_rd  <- readxl::excel_sheets(xlsx_path)
  stock_tabs   <- setdiff(intersect(all_tabs_vis, all_tabs_rd), ignore_tabs)
  
  if (!length(stock_tabs)) return(tibble())
  
  bind_rows(lapply(stock_tabs, function(tab){
    extract_tab_rows(xlsx_path, tab, Scorer)
  }))
}

##------------------------------------------------------------------------------
## RUN: loop all files in folder and compile into score_table_all
##------------------------------------------------------------------------------
xlsx_files <- list.files(
  in_dir,
  pattern = "(?i)\\.xl[a-z]+$",     ## xl* (xlsx, xls, xlsm, xlx typos)
  full.names = TRUE
)

length(xlsx_files)  ## should be 16

score_table_all <- bind_rows(lapply(xlsx_files, build_score_table_for_workbook)) |>
  arrange(Scorer, stock_name, row_idx)

## Quick QA
nrow(score_table_all)
dplyr::count(score_table_all, Scorer, stock_name) |> print(n = 50)

## Flag unknown attribute names
unknown_attrs <- setdiff(unique(score_table_all$Attribute_name), c(sens_attrs, rigid_attrs))
if (length(unknown_attrs)) {
  message("Unrecognized Attribute_name values:\n- ", paste(unknown_attrs, collapse = "\n- "))
}

## Write to CSV for downstream steps
write.csv(score_table_all, file = file.path(in_dir, "score_table_all.csv"), row.names = FALSE)

