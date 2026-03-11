##------------------------------------------------------------------------------
## Caribbean CVA – Compile reviewer FINAL scoring workbooks into one table
## Builds: final_score_table_all_clean
##
## Reads final attribute scores already computed in each workbook:
##   C17:C18 with K17:K18
##   B21:B28 with K21:K28
##   B30:B35 with K30:K35
##
## !!!NOTE!!! Sheets being read in cannot be open in Excel or they may fail
## to read in the code loop

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

## CHANGED: added qualitative exposure factors
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
read_name_score_block <- function(xlsx_path, tab, name_range, score_range){
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
  
  ## Main block: names in B21:B28, scores in K21:K28
  block2 <- read_name_score_block(
    xlsx_path   = xlsx_path,
    tab         = tab,
    name_range  = "B21:B28",
    score_range = "K21:K28"
  )
  
  ## Bottom block: names in B30:B35, scores in K30:K35
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
      row_idx    = c(17:18, 21:28, 30:35)   ## CHANGED: keep source row numbers
    ) |>
    dplyr::mutate(
      Attribute_type = dplyr::case_when(
        Attribute_name %in% sens_attrs      ~ "Sensitivity",
        Attribute_name %in% rigid_attrs     ~ "Rigidity",
        Attribute_name %in% qual_exp_attrs  ~ "Qualitative Exposure Factors",
        TRUE                                ~ NA_character_
      )
    ) |>
    dplyr::filter(!is.na(Attribute_name), Attribute_name != "")
  
  out
}

##------------------------------------------------------------------------------
## CHANGED: helper to determine whether a sheet has any final scores
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
## CHANGED: Safe loader for final score workbooks
## Copies to temp, reads visible stock tabs, and keeps only tabs with real scores
build_score_table_for_workbook <- function(xlsx_path){
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
    
    ## CHANGED: only retain tabs that actually contain final scores
    scored_tabs <- stock_tabs[vapply(stock_tabs, function(tab){
      sheet_has_final_scores(tmp, tab)
    }, logical(1))]
    
    if (!length(scored_tabs)) {
      message("ℹ️  No scored stock tabs in: ", basename(xlsx_path))
      return(tibble::tibble())
    }
    
    dplyr::bind_rows(lapply(scored_tabs, function(tab){
      extract_tab_rows_final(tmp, tab, Scorer)
    }))
    
  }, error = function(e){
    message("❌ Error in ", basename(xlsx_path), ": ", conditionMessage(e))
    tibble::tibble()
  })
}

##------------------------------------------------------------------------------
## COMPILE (with simple logging)

xlsx_files <- list.files(
  in_dir,
  pattern = "(?i)\\.xl[a-z]+$",   ## xlsx/xls/xlsm…
  full.names = TRUE
)

message("Found files: ", length(xlsx_files))

per_file <- lapply(xlsx_files, function(f){
  message("→ Reading: ", basename(f))
  out <- build_score_table_for_workbook(f)
  message("  rows: ", nrow(out), "  (scorer: ", parse_scorer(f), ")")
  out
})

final_score_table_all <- dplyr::bind_rows(per_file) |>
  dplyr::arrange(Scorer, stock_name, row_idx)

message("Total rows: ", nrow(final_score_table_all),
        " | scorers: ", dplyr::n_distinct(final_score_table_all$Scorer),
        " | stocks: ",  dplyr::n_distinct(final_score_table_all$stock_name))

##------------------------------------------------------------------------------
## RUN: loop all files in folder and compile into final_score_table_all

xlsx_files <- list.files(
  in_dir,
  pattern = "(?i)\\.xl[a-z]+$",
  full.names = TRUE
); length(xlsx_files)  ## should be 16

## --- Build master table ---
final_score_table_all <- bind_rows(lapply(xlsx_files, build_score_table_for_workbook)) |>
  arrange(Scorer, stock_name, row_idx)

## Scorers × number of stocks they scored
scorer_stock_counts <- final_score_table_all %>%
  distinct(Scorer, stock_name) %>%
  count(Scorer, name = "n_stocks"); scorer_stock_counts

##------------------------------------------------------------------------------
## CHANGED: remove rows where Final_score is NA
final_score_table_all <- final_score_table_all %>%
  dplyr::filter(!is.na(Final_score))

## Check again
scorer_stock_counts <- final_score_table_all %>%
  distinct(Scorer, stock_name) %>%
  count(Scorer, name = "n_stocks"); scorer_stock_counts

##------------------------------------------------------------------------------
## QA/QC
n_scorers <- final_score_table_all %>% distinct(Scorer) %>% nrow()
n_stocks  <- final_score_table_all %>% distinct(stock_name) %>% nrow()
n_attrs   <- final_score_table_all %>% distinct(Attribute_name) %>% nrow()

n_scorers; n_stocks; n_attrs
## expected approx. 16 scorers, 25 stocks, 16 attributes

## Per stock: how many scorers, how many attributes populated
stock_qa <- final_score_table_all %>%
  group_by(stock_name) %>%
  summarise(
    n_scorers = n_distinct(Scorer),
    n_attrs   = n_distinct(Attribute_name),
    .groups = "drop"
  ) %>%
  arrange(stock_name) %>% as.data.frame(); stock_qa

## Scorers per stock
species_reviews <- final_score_table_all %>%
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
## Optional QA: show any attributes not categorized
uncategorized_attributes <- final_score_table_all %>%
  filter(is.na(Attribute_type)) %>%
  distinct(Attribute_name); uncategorized_attributes

##------------------------------------------------------------------------------
## Cleaning: remove extra scores
to_remove <- tibble::tribble(
  ~stock_name,          ~Scorer,
  "Blue runner",        "RGarciaSais",
  "King mackerel",      "AAcosta",
  "Yellowtail snapper", "REsteves"
)

## CHECK: show rows that will be removed
will_remove <- final_score_table_all %>%
  dplyr::semi_join(to_remove, by = c("stock_name","Scorer"))

message("Rows to remove: ", nrow(will_remove))
print(will_remove %>% dplyr::count(stock_name, Scorer, name = "n_rows"))

final_score_table_all_clean <- final_score_table_all %>%
  dplyr::anti_join(to_remove, by = c("stock_name","Scorer"))

message("Rows before: ", nrow(final_score_table_all),
        " | after: ", nrow(final_score_table_all_clean),
        " | removed: ", nrow(final_score_table_all) - nrow(final_score_table_all_clean))

## Final check: scorers per stock
species_reviews_clean <- final_score_table_all_clean %>%
  filter(!is.na(stock_name), !is.na(Scorer)) %>%
  distinct(stock_name, Scorer) %>%
  group_by(stock_name) %>%
  summarise(
    n_reviews = n_distinct(Scorer),
    reviewers = paste(sort(unique(Scorer)), collapse = ", "),
    .groups   = "drop"
  ) %>%
  arrange(desc(n_reviews), stock_name) %>%
  as.data.frame(); species_reviews_clean

##------------------------------------------------------------------------------
## Write to CSV for downstream steps
write.csv(
  final_score_table_all_clean,
  file = file.path(in_dir, "final_score_table_all_clean.csv"),
  row.names = FALSE
)