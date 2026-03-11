##------------------------------------------------------------------------------
## Caribbean CVA – TEST extract reviewer final attribute scores
## Reads ONLY the first Excel workbook and first non-empty scored sheet
## Extracts final reviewer scores already computed in the workbook:
##   C17:C18 with K17:K18
##   B21:B28 with K21:K28
##   B30:B35 with K30:K35

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
in_dir      <- "./data/final-scores"   ## CHANGED: now points to final score workbooks
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

qual_exp_attrs <- c(
  "Sargassum influx",
  "Thermocline depth"
)

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
    vis <- tryCatch(openxlsx::getSheetVisibility(xlsx_path),
                    error = function(e) rep("visible", length(sn)))
    keep <- is.na(vis) | tolower(vis) == "visible"
    return(trimws(sn[keep]))
  }
  readxl::excel_sheets(xlsx_path)
}

##------------------------------------------------------------------------------
## CHANGED: helper to read one attribute-name block + one score block
read_name_score_block <- function(xlsx_path, tab, name_range, score_range, name_col = "Attribute_name"){
  nm <- readxl::read_excel(
    path      = xlsx_path,
    sheet     = tab,
    range     = name_range,
    col_names = FALSE,
    col_types = "text"
  )
  
  sc <- readxl::read_excel(
    path      = xlsx_path,
    sheet     = tab,
    range     = score_range,
    col_names = FALSE,
    col_types = "text"
  )
  
  if (nrow(nm) != nrow(sc)) {
    stop("Name and score ranges have different numbers of rows in tab: ", tab)
  }
  
  tibble::tibble(
    Attribute_name = stringr::str_squish(as.character(nm[[1]])),
    Final_score    = suppressWarnings(as.numeric(sc[[1]]))
  )
}

##------------------------------------------------------------------------------
## CHANGED: new extractor using exact final-score cells
extract_tab_rows_final <- function(xlsx_path, tab, Scorer){
  
  ## Top two attributes: names in C17:C18, scores in K17:K18
  block1 <- read_name_score_block(
    xlsx_path   = xlsx_path,
    tab         = tab,
    name_range  = "C17:C18",
    score_range = "K17:K18"
  )
  
  ## Main sensitivity block: names in B21:B28, scores in K21:K28
  block2 <- read_name_score_block(
    xlsx_path   = xlsx_path,
    tab         = tab,
    name_range  = "B21:B28",
    score_range = "K21:K28"
  )
  
  ## Bottom rigidity block: names in B30:B35, scores in K30:K35
  block3 <- read_name_score_block(
    xlsx_path   = xlsx_path,
    tab         = tab,
    name_range  = "B30:B35",
    score_range = "K30:K35"
  )
  
  out <- dplyr::bind_rows(block1, block2, block3) |>
    dplyr::mutate(
      SourceFile = basename(xlsx_path),
      Scorer     = Scorer,
      stock_name = tab,
      row_idx    = c(17:18, 21:28, 30:35)   ## CHANGED: keep source row numbers for tracking
    ) |>
    dplyr::mutate(
      Attribute_type = dplyr::case_when(
        Attribute_name %in% sens_attrs  ~ "Sensitivity",
        Attribute_name %in% rigid_attrs ~ "Rigidity",
        Attribute_name %in% qual_exp_attrs  ~ "Qualitative Exposure Factors",
        TRUE                            ~ NA_character_
      )
    ) |>
    dplyr::filter(!(is.na(Attribute_name) & is.na(Final_score))) |>
    dplyr::filter(!is.na(Attribute_name) & Attribute_name != "")
  
  out
}

##------------------------------------------------------------------------------
## CHANGED: simpler sheet test
## Keep first visible stock sheet with at least one non-NA final score
sheet_has_final_scores <- function(xlsx_path, tab){
  
  k1 <- readxl::read_excel(
    path      = xlsx_path,
    sheet     = tab,
    range     = "K17:K18",
    col_names = FALSE,
    col_types = "text"
  )[[1]]
  
  k2 <- readxl::read_excel(
    path      = xlsx_path,
    sheet     = tab,
    range     = "K21:K28",
    col_names = FALSE,
    col_types = "text"
  )[[1]]
  
  k3 <- readxl::read_excel(
    path      = xlsx_path,
    sheet     = tab,
    range     = "K30:K35",
    col_names = FALSE,
    col_types = "text"
  )[[1]]
  
  vals <- suppressWarnings(as.numeric(c(k1, k2, k3)))
  any(!is.na(vals))
}

##------------------------------------------------------------------------------
## CHANGED: test loader for one workbook and first scored sheet
build_score_table_for_workbook_test <- function(xlsx_path){
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
    
    ## diagnostic
    diag_tabs <- tibble::tibble(
      tab = stock_tabs,
      has_final_scores = vapply(stock_tabs, function(tab){
        sheet_has_final_scores(tmp, tab)
      }, logical(1))
    )
    print(diag_tabs)
    
    scored_tabs <- diag_tabs$tab[diag_tabs$has_final_scores]
    
    if (!length(scored_tabs)) {
      message("ℹ️  No scored stock tabs found in: ", basename(xlsx_path))
      return(tibble::tibble())
    }
    
    first_tab <- scored_tabs[1]
    message("Testing first scored sheet only: ", first_tab)
    
    extract_tab_rows_final(tmp, first_tab, Scorer)
    
  }, error = function(e){
    message("❌ Error in ", basename(xlsx_path), ": ", conditionMessage(e))
    tibble::tibble()
  })
}

##------------------------------------------------------------------------------
## TEST: first workbook only

xlsx_files <- list.files(
  in_dir,
  pattern = "(?i)\\.xl[a-z]+$",
  full.names = TRUE
)

message("Found files: ", length(xlsx_files))

first_file <- xlsx_files[1]   ## CHANGED: only use first workbook for testing
message("→ Testing workbook: ", basename(first_file))

score_table_test <- build_score_table_for_workbook_test(first_file) |>
  dplyr::arrange(Scorer, stock_name, row_idx)

score_table_test