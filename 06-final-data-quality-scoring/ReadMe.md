## Purpose

`1-extract-data-quality-scores.R` compiles reviewer-assigned **data quality scores** from the final Caribbean CVA scoring workbooks into a single long-format table for downstream quality assessment and analysis. `2-summarize-data-quality.R` summarizes into a final table. 

## Background

In this workflow, reviewers assign a **data quality score** to each qualitative exposure factor and sensitivity attribute. Scores range from `0` to `3`:

| Score | Label           | Description                           |
|-------|-----------------| --------------------------------------|
| 3     | Adequate Data   | The score is based on data which have been observed, modeled or empirically measured for the species in question and comes from a reputable source. |
| 2     | Limited Data    | The score is based on data which has a higher degree of uncertainty. The data used to score the attribute may be based on related or similar species, come from outside the study area, or the reliability of the source may be limited. |
| 1     | Expert Judgment | The attribute score reflects the expert judgment of the reviewer and is based on their general knowledge of the species, or other related species, and their relative role in the ecosystem. |
| 0     | No Data         | No information to base an attribute score on. Very little is known about the species or related species and there is no basis for forming an expert opinion. |

Following the NOAA CVA-style approach used in the script, scores are also classified as:

- **data-based** if `score >= 2`
- **low-information** if `score <= 1`

## File structure

##### Inputs
```
    data/
    └── final-scores/
        ├── CVA_Scores_Final_AB.xlsx
        ├── CVA_Scores_Final_CD.xlsx
        ├── CVA_Scores_Final_EF.xlsx
        └── ...
```

##### Outputs
```
    outputs/
    └── final-scores-compiled/
        └── data-quality/
            ├── table_data_quality_scores_extracted.csv
            ├── qa_data_quality_by_stock.csv
            └── qa_data_quality_by_scorers.csv
```

# 1-extract-data-quality-scores.R

This script is the **first step** in the data quality workflow. Its main job is to:

1. read data quality scores from each reviewer workbook,
2. assemble them into one standardized table across scorers and stocks,
3. attach row labels and metadata, and
4. write out QA/QC and summary files for later scripts.

The script is designed to be robust to:
- empty or partially completed tabs,
- workbook lock issues,
- hidden sheets, and
- sheet-level read problems.

### Inputs

**Filename convention:** the scorer's initials are parsed from the last underscore-delimited segment of the filename (e.g., `CVA_Scores_Final_AB.xlsx` → scorer `AB`).

**Workbook structure:**

- One Excel tab per stock (species). Tabs named `Instructions`, `Data Quality`, or `Example` are ignored.
- Target cells: `L17:L18` (Exposure attributes) and `L21:L28` (Sensitivity attributes).
- Attribute labels are read from columns B/C in the corresponding rows.

### Output

All outputs are written to `./outputs/final-scores-compiled/data-quality/`:

| File | Description |
|------|-------------|
| `table_data_quality_scores_extracted.csv` | Long table — one row per stock × scorer × attribute |
| `qa_data_quality_by_stock.csv` | Count of scorers and scored rows per stock |
| `qa_data_quality_by_scorers.csv` | Reviewer list and review count per stock |

### Workflow

The script follows a workbook-by-workbook, tab-by-tab extraction process.

#### 1. Set directories and target rows

The script defines the input and output directories, the stock-tab rows to extract, and the set of tabs to ignore.

- Input directory: `./data/final-scores`
- Output directory: `./outputs/final-scores-compiled/data-quality`
- Ignored tabs: `Instructions`, `Data Quality`, `Example`
- Target rows: `17:18` and `21:28`

These target rows are stored in:

- `dq_rows <- c(17:18, 21:28)`

#### 2. Load all Excel workbooks

The script searches the input directory for Excel files using a case-insensitive pattern that matches standard workbook extensions.

Each workbook is then processed one at a time.

#### 3. Parse scorer initials from the filename

For each workbook, the script extracts reviewer initials using `parse_scorer()`.

This function:
- removes the file extension,
- splits the filename on underscores,
- takes the last segment,
- trims whitespace, and
- keeps the first two characters.

If scorer initials cannot be parsed, the script assigns `Unknown Scorer`.

#### 4. Copy each workbook to a temporary file

Before reading any sheets, the script copies the workbook to a temporary `.xlsx` file.

This helps avoid failures caused by:
- Excel file locks,
- OneDrive or cloud-sync conflicts,
- partially open workbooks.

If a workbook cannot be copied, it is skipped and the script moves on.

#### 5. Identify visible stock tabs

The script uses `safe_visible_sheet_names()` to retrieve visible sheet names.

- If `openxlsx` is available, it checks sheet visibility directly.
- Otherwise, it falls back to `readxl::excel_sheets()`.

After that, the script removes ignored tabs and treats the remaining visible tabs as stock tabs.

#### 6. Safely read score and label ranges

All cell-range reads are handled through `safe_read_range()`.

This helper function:
- attempts to read the requested sheet/range,
- returns an empty tibble if the range is unreadable or malformed,
- prevents one failed tab from crashing the full script.

This makes the extraction process much more robust across inconsistent workbooks.

#### 7. Read data quality scores from column L

For each stock tab, the script reads the final entered data-quality scores from:

- `L17:L18`
- `L21:L28`
- `L30:L35`

These are combined into one 10-element score vector by `read_data_quality_scores()`.

If fewer than 2 or 8 values are found in the two ranges, the script pads the result with `NA` so the row alignment remains stable.

#### 8. Read row labels from columns C and B

For the same stock tab, the script reads attribute labels from:

- `C17:C18`
- `B21:B28`
- `B30:B35`

These are combined into one 10-element label vector by `read_data_quality_labels()`.

Blank or missing labels are replaced with fallback values like:
- `Row_17`
- `Row_18`
- `Row_21`
- ...
- `Row_28`

This guarantees that each extracted score has an associated row label.

#### 9. Build one long-format table per stock tab

The function `extract_data_quality_rows()` combines the extracted values into a tidy table with one row per:

- stock × scorer × attribute row

It also attaches the following metadata: `SourceFile`, `Scorer`, `stock_name`, `row_idx`, `Attribute_type`, `Attribute_name`, `Data_quality_score`, `score_cell`, `label_cell`, `Data_quality_label`, `score_ge_2`

Attribute type is assigned as: `Exposure` for rows `17:18` and `Sensitivity` for rows `21:28` and `30:35` (Note these are labeled "Rigidity" in the scorer workbooks. 

Rows with missing scores are dropped at this stage.

#### 10. Skip empty or unscored tabs

Before extraction, the script checks whether a stock tab contains any usable data quality scores using `sheet_has_data_quality_scores()`. If all extracted score values are missing, the tab is skipped. This prevents empty tabs from contributing empty rows to the final output.

#### 11. Combine all extracted rows across all workbooks

After all files are processed, the script binds all per-file results together into one master object: `data_quality_table_all`

The final extracted table is sorted by: `Scorer`, `stock_name`, `row_idx`. 

This becomes the main long-format dataset for downstream analyses.

#### 12. Run QA/QC checks

The script performs several quick QA/QC checks:

- number of distinct scorers
- number of distinct stocks
- number of distinct row indices extracted

It also builds a stock-level QA table showing:

- `n_scorers`: number of unique scorers per stock
- `n_rows`: number of unique scored rows extracted per stock

This helps verify that the expected rows were found and that reviewer coverage is reasonable.

##### Summarize reviewer coverage by stock

A second QA table is created to show reviewer coverage for each stock.

For each stock, the script reports:

- `n_reviews`: number of distinct reviewers
- `reviewers`: comma-separated scorer initials

This is useful for checking participation and confirming which reviewers contributed to each stock.

#### 14. Write output files

Finally, the script writes three CSV files to the output directory:

1. `table_data_quality_scores_extracted.csv`  
   Main extracted long-format data table.

2. `qa_data_quality_by_stock.csv`  
   QA summary of scorers and extracted row counts by stock.

3. `qa_data_quality_by_scorers.csv`  
   Reviewer coverage table by stock.

### Main extracted table structure

The primary output table contains one row per extracted data-quality score and includes the following fields:

| Column | Description |
|--------|-------------|
| `SourceFile` | Workbook filename |
| `Scorer` | Reviewer initials parsed from the filename |
| `stock_name` | Stock/species name, based on the tab name |
| `row_idx` | Original Excel row number |
| `Attribute_type` | `Exposure` or `Sensitivity` |
| `Attribute_name` | Label read from the workbook |
| `Data_quality_score` | Numeric score from 0 to 3 |
| `score_cell` | Excel cell address of the score |
| `label_cell` | Excel cell address of the label |
| `Data_quality_label` | Interpreted score label |
| `score_ge_2` | Logical indicator for `Data_quality_score >= 2` |

### Notes and assumptions
- The script assumes that final data quality scores have already been entered into the target cells.
- Scores are expected to be numeric values in the set `0, 1, 2, 3`.
- Blank score cells are treated as missing and are excluded from the final extracted table.
- Exposure row labels are read from column `C` for rows `17:18`.
- Sensitivity row labels are read from column `B` for rows `21:28` and `30:35`.
- Stock names are taken directly from the worksheet names.
- Reviewer initials depend on consistent file naming.

# Summarize Data Quality Scores: `2-summarize-data-quality-scores.R`

Summarizes extracted data quality scores by stock, classifying overall data quality using the NOAA CVA proportion-based ranking system.

## Overview

This script reads the long-format table of individual data quality scores produced by the extraction step (`1-compile-data-quality-scores.R`) and computes per-stock summaries including score distributions, the proportion of scores ≥ 2, and an overall quality rank.

## Data Quality Ranking

Overall data quality per stock is determined by the proportion of scores ≥ 2:

| Rank     | Proportion of scores ≥ 2 |
|----------|--------------------------|
| High     | ≥ 80%                    |
| Moderate | 50–79%                   |
| Poor     | < 50%                    |

Scores range from 0 (No Data) to 3 (Adequate Data).

## Input

| File | Location |
|------|----------|
| `table_data_quality_scores_extracted.csv` | `./outputs/final-scores-compiled/data-quality/` |

This file is produced by the upstream extraction script.

## Output

| File | Description |
|------|-------------|
| `overall_data_quality_summary_by_stock.csv` | One row per stock with score counts, proportion ≥ 2, mean, SD, and quality rank |

Written to `./outputs/final-scores-compiled/data-quality/`.

## Output Columns

- `stock_name` — Species or stock identifier
- `n_3`, `n_2`, `n_1`, `n_0` — Count of scores at each level
- `n_ge_2` — Count of scores ≥ 2
- `prop_ge_2` — Proportion of scores ≥ 2
- `mean_score` — Mean data quality score
- `sd_score` — Standard deviation of scores
- `data_quality_rank` — High, Moderate, or Poor
