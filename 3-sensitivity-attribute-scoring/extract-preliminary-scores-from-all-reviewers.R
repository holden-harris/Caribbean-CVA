##------------------------------------------------------------------------------
## Caribbean CVA – Compile reviewer scoring workbooks into one table
## Builds: score_table_all  (rows from rows 21–35, skipping 29; B/E/F–I)
##
## !!!NOTE!!! Sheets being read in cannot be open in Excel or they will fail to read in the code loop

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

## Extract rows 21–35 (skip 29) from B..I for one tab
extract_tab_rows <- function(xlsx_path, tab, Scorer){
  rng <- readxl::read_excel(
    path      = xlsx_path,
    sheet     = tab,
    range     = "B21:I35",
    col_names = FALSE,
    col_types = "text"
  )
  if (ncol(rng) < 8) return(tibble::tibble())
  
  rng <- rng |>
    dplyr::mutate(row_idx = 21:35) |>
    dplyr::filter(row_idx != 29)
  
  out <- rng |>
    dplyr::transmute(
      SourceFile      = basename(xlsx_path),
      Scorer          = Scorer,
      stock_name      = tab,
      row_idx         = row_idx,
      Attribute_name  = stringr::str_squish(as.character(`...1`)),
      Data_quality    = suppressWarnings(as.numeric(`...4`)),
      Scoring_rank_1  = suppressWarnings(as.numeric(`...5`)),
      Scoring_rank_2  = suppressWarnings(as.numeric(`...6`)),
      Scoring_rank_3  = suppressWarnings(as.numeric(`...7`)),
      Scoring_rank_4  = suppressWarnings(as.numeric(`...8`))
    ) |>
    dplyr::mutate(
      Attribute_type = dplyr::case_when(
        Attribute_name %in% sens_attrs  ~ "Sensitivity",
        Attribute_name %in% rigid_attrs ~ "Rigidity",
        TRUE                            ~ NA_character_
      )
    ) |>
    dplyr::filter(!(is.na(Attribute_name) &
                      dplyr::if_all(dplyr::starts_with("Scoring_rank_"), ~is.na(.x))))
  out
}

## Safe loader: copy to temp (breaks OneDrive/Excel locks), no unzip, resilient
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
    
    dplyr::bind_rows(lapply(stock_tabs, function(tab){
      extract_tab_rows(tmp, tab, Scorer)
    }))
  }, error = function(e){
    message("❌ Error in ", basename(xlsx_path), ": ", conditionMessage(e))
    tibble::tibble()
  })
}

## ---------- COMPILE (with simple logging) ----------

xlsx_files <- list.files(
  in_dir,
  pattern = "(?i)\\.xl[a-z]+$",   # xlsx/xls/xlsm…
  full.names = TRUE
)

message("Found files: ", length(xlsx_files))
per_file <- lapply(xlsx_files, function(f){
  message("→ Reading: ", basename(f))
  out <- build_score_table_for_workbook(f)
  message("  rows: ", nrow(out), "  (scorer: ", parse_scorer(f), ")")
  out
})

score_table_all <- dplyr::bind_rows(per_file) |>
  dplyr::arrange(Scorer, stock_name, row_idx)

message("Total rows: ", nrow(score_table_all),
        " | scorers: ", dplyr::n_distinct(score_table_all$Scorer),
        " | stocks: ",  dplyr::n_distinct(score_table_all$stock_name))

##------------------------------------------------------------------------------
## RUN: loop all files in folder and compile into score_table_all

xlsx_files <- list.files(
  in_dir,
  pattern = "(?i)\\.xl[a-z]+$",     ## xl* (xlsx, xls, xlsm, xlx typos)
  full.names = TRUE
); length(xlsx_files)  ## should be 16

## --- Build master table ---
score_table_all <- bind_rows(lapply(xlsx_files, 
                                    build_score_table_for_workbook)) |>
  arrange(Scorer, stock_name, row_idx)

## Scorers × number of stocks they scored
scorer_stock_counts <- 
  score_table_all %>%
  distinct(Scorer, stock_name) %>%
  count(Scorer, name = "n_stocks"); scorer_stock_counts

## Note this includes hidden sheets, so includes stocks not scored by the reviewer
## Next, we'll remove these rows where all 5 fields (data_quality, scoring_rank1..4 are NA)

## Remove rows where ALL five fields are NA
score_table_all <- score_table_all %>%
  dplyr::filter(
    !dplyr::if_all(
      dplyr::all_of(c("Data_quality","Scoring_rank_1","Scoring_rank_2","Scoring_rank_3","Scoring_rank_4")),
      ~ is.na(.x)
    )
  )

## Check again
scorer_stock_counts <- 
  score_table_all %>%
  distinct(Scorer, stock_name) %>%
  count(Scorer, name = "n_stocks"); scorer_stock_counts

##------------------------------------------------------------------------------
## QA/QC
n_scorers <- score_table_all %>% distinct(Scorer) %>% nrow() 
n_stocks <- score_table_all %>% distinct(stock_name) %>% nrow() 
n_attrs  <- score_table_all %>% distinct(Attribute_name) %>% nrow() ## How many unique attributes?
n_scorers; n_stocks; n_attrs ## --> Should be 16 scorers, 25 stocks, 14 attributes

## Per stock: how many scorers, how many attributes populated
stock_qa <- score_table_all %>%
  group_by(stock_name) %>%
  summarise(
    n_scorers = n_distinct(Scorer),
    n_attrs   = n_distinct(Attribute_name),
    .groups = "drop"
  ) %>%
  arrange(stock_name) %>% as.data.frame(); stock_qa

## Scorers pers stock
species_reviews <- score_table_all %>%
  filter(!is.na(stock_name), !is.na(Scorer)) %>%
  distinct(stock_name, Scorer) %>%                 # one row per stock × reviewer
  group_by(stock_name) %>%
  summarise(
    n_reviews = n_distinct(Scorer),
    reviewers = paste(sort(unique(Scorer)), collapse = ", "),
    .groups   = "drop"
  ) %>%
  arrange(desc(n_reviews), stock_name) %>% 
  as.data.frame(); species_reviews

## -----------------------------------------------------------------------------
## Cleaning: remove extra scores
to_remove <- tibble::tribble(
  ~stock_name,     ~Scorer,
  "Blue runner",   "RGarciaSais",
  "King mackerel", "AAcosta"
)

##CHECK: show rows that will be removed (distinct by scorer/species)
will_remove <- score_table_all %>%
  dplyr::semi_join(to_remove, by = c("stock_name","Scorer"))

message("Rows to remove: ", nrow(will_remove)); print(will_remove %>% dplyr::count(stock_name, Scorer, name = "n_rows"))

score_table_all_clean <- score_table_all %>%
  dplyr::anti_join(to_remove, by = c("stock_name","Scorer"))

message("Rows before: ", nrow(score_table_all),
        " | after: ", nrow(score_table_all_clean),
        " | removed: ", nrow(score_table_all) - nrow(score_table_all_clean))

## Final check:ccorers pers sotck
species_reviews_clean <- score_table_all_clean %>%
  filter(!is.na(stock_name), !is.na(Scorer)) %>%
  distinct(stock_name, Scorer) %>%                 # one row per stock × reviewer
  group_by(stock_name) %>%
  summarise(
    n_reviews = n_distinct(Scorer),
    reviewers = paste(sort(unique(Scorer)), collapse = ", "),
    .groups   = "drop"
  ) %>%
  arrange(desc(n_reviews), stock_name) %>% 
  as.data.frame(); species_reviews_clean


## -----------------------------------------------------------------------------
## Write to CSV for downstream steps
write.csv(score_table_all_clean, file = file.path(in_dir, "score_table_all_clean.csv"), row.names = FALSE)


