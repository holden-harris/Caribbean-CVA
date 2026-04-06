##------------------------------------------------------------------------------
## Caribbean CVA – Preliminary Species SDs for Reviewers’ Pre-Workshop Scoring
## Mirrors HMS logic for attribute means/SDs and rolls up to a species-level SD

rm(list = ls()); gc()
library(readxl)
library(stringr)
library(dplyr)
library(tibble)

xlsx_path <- "./data/preliminary-scores/Caribbean_CVA_Scoring_Template_2025_AAcosta.xlsx"  ## as given


##------------------------------------------------------------------------------
## Parse Scorer from the filename
##   - Take the last underscore-delimited token before the extension
##   - "Caribbean CVA Scoring Template_2025_AAcosta.xlx" --> "AAcosta"

fn <- basename(xlsx_path)
fn_noext <- str_remove(fn, "\\.[Xx][Ll][SsXx]$")   ## strip .xls/.xlsx/.xlx (case-insensitive)
parts <- str_split(fn_noext, "_", simplify = TRUE)
Scorer <- parts[, ncol(parts)] |> str_squish()
Scorer

##------------------------------------------------------------------------------
## Get stock names from tabs
##   - Keep as `stock_name`

all_tabs <- readxl::excel_sheets(xlsx_path)       ## <-- use xlsx_path
ignore_tabs <- c("Instructions", "Data Quality", "Example")
stock_tabs <- setdiff(all_tabs, ignore_tabs)

## Quick peek
all_tabs
stock_tabs

## Inspect sheet visibility
props <- readxl::sheet_properties(xlsx_path)
props[, c("name","hidden")]
# hidden is logical: FALSE (visible), TRUE (hidden), or "veryHidden" in some cases

## Keep only visible (not hidden) + drop admin tabs
ignore_tabs <- c("Instructions","Data Quality","Example")

visible_tabs <- props |>
  dplyr::filter(!hidden %in% c(TRUE, "veryHidden")) |>
  dplyr::pull(name)

## (Optional) normalize whitespace just in case
visible_tabs <- trimws(visible_tabs)

## Final stock tabs list
stock_tabs <- setdiff(visible_tabs, ignore_tabs)

## Quick peek
props[, c("name","hidden")]
stock_tabs


doc <- xml2::read_xml(wb_xml)

## workbook.xml uses namespaces; grab them and set a short alias
ns <- xml2::xml_ns(doc)
 names(ns)  # optional

 ## Select all <sheet> nodes under <sheets>
 sheets_nodes <- xml2::xml_find_all(doc, ".//d:sheets/d:sheet", ns = c(d = ns[[1]]))
 
 sheet_names <- xml2::xml_attr(sheets_nodes, "name")
 sheet_state <- xml2::xml_attr(sheets_nodes, "state")  # NA means visible; can be "hidden" or "veryHidden"
 
 props <- data.frame(
   name   = sheet_names,
   hidden = ifelse(is.na(sheet_state), FALSE, sheet_state),  # FALSE, "hidden", or "veryHidden"
   stringsAsFactors = FALSE
 )
 
 props

 ##------------------------------------------------------------------------------
 ## Step D. Keep only **visible** sheets, then drop admin tabs
 ##------------------------------------------------------------------------------
 ignore_tabs  <- c("Instructions","Data Quality","Example")
 
 visible_tabs <- props$name[ props$hidden == FALSE ]           # visible only
 visible_tabs <- trimws(visible_tabs)
 
 ## For safety, keep only names that readxl also reports (ordering/encoding guard)
 all_tabs <- readxl::excel_sheets(xlsx_path)
 visible_tabs <- intersect(all_tabs, visible_tabs)
 
 stock_tabs <- setdiff(visible_tabs, ignore_tabs)
 
 ## Quick peek
 props[, c("name","hidden")]
 stock_tabs
 
 
 ##------------------------------------------------------------------------------
 ## Step 13. For each visible stock tab, read rows 21–35 (skip row 29)
 ##          Columns needed:
 ##          - B: Attribute_name
 ##          - E: Data_quality
 ##          - F–I: Scoring_rank_1..4
 ##------------------------------------------------------------------------------
 
 library(readxl)
 library(dplyr)
 library(stringr)
 library(tidyr)
 library(tibble)
 
 ## Container for all tabs
 rows_21_35_all <- NULL
 
 for (tab in stock_tabs) {
   
   ## Read exactly the range B21:I35 (8 columns wide, 15 rows tall)
   rng <- readxl::read_excel(
     path      = xlsx_path,
     sheet     = tab,
     range     = "B21:I35",
     col_names = FALSE,
     col_types = "text"
   )
   
   ## Add absolute row index (21..35) and drop row 29
   rng <- rng |>
     mutate(row_idx = 21:35) |>
     filter(row_idx != 29)
   
   ## Rename columns we need:
   ## B=...1, C=...2, D=...3, E=...4, F=...5, G=...6, H=...7, I=...8
   sub <- rng |>
     transmute(
       Scorer          = Scorer,
       stock_name      = tab,
       row_idx         = row_idx,
       Attribute_name  = as.character(`...1`),
       Data_quality    = suppressWarnings(as.numeric(`...4`)),
       Scoring_rank_1  = suppressWarnings(as.numeric(`...5`)),
       Scoring_rank_2  = suppressWarnings(as.numeric(`...6`)),
       Scoring_rank_3  = suppressWarnings(as.numeric(`...7`)),
       Scoring_rank_4  = suppressWarnings(as.numeric(`...8`))
     ) |>
     ## trim whitespace on the attribute text
     mutate(Attribute_name = str_squish(Attribute_name)) |>
     ## drop completely empty rows (no attribute and no scores)
     filter(!(is.na(Attribute_name) & if_all(starts_with("Scoring_rank_"), ~is.na(.x))))
   
   rows_21_35_all <- bind_rows(rows_21_35_all, sub)
 }
 
 ## Quick peek
 nrow(rows_21_35_all)
 head(rows_21_35_all, 12)
 
 ## Optional: enforce canonical order of attributes (for QA)
 canonical_order <- c(
   # Sensitivity attributes
   "Habitat specificity",
   "Prey specificity",
   "Tolerance to ocean acidification",
   "Complexity in reproductive strategy",
   "Species range",
   "Specificity in early life history requirements",
   "Stock size/status",
   "Other stressors",
   # Rigidity attributes
   "Population growth rate",
   "Mobility and dispersal of early life stages",
   "Adult mobility",
   "Spawning characteristics",
   "Predation and competition dynamics",
   "Genetic diversity"
 )
 
 rows_21_35_all <- rows_21_35_all |>
   mutate(Attribute_name = ifelse(Attribute_name %in% canonical_order, Attribute_name, Attribute_name)) |>
   arrange(stock_name, row_idx)
 
 ##------------------------------------------------------------------------------
 ## Add Attribute_type = "Sensitivity" or "Rigidity"
 
 ## Canonical sets
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
 
 ## Tag the type
 rows_21_35_all <- rows_21_35_all |>
   dplyr::mutate(
     Attribute_type = dplyr::case_when(
       Attribute_name %in% sens_attrs  ~ "Sensitivity",
       Attribute_name %in% rigid_attrs ~ "Rigidity",
       TRUE                            ~ NA_character_
     )
   )
 
 ## Optional: sanity-check for any unexpected attribute names
 unknown_attrs <- setdiff(unique(rows_21_35_all$Attribute_name), c(sens_attrs, rigid_attrs))
 if (length(unknown_attrs)) {
   message("Unrecognized Attribute_name values:\n- ", paste(unknown_attrs, collapse = "\n- "))
 }
 
 ## Quick peek
 head(rows_21_35_all, 12)
 
