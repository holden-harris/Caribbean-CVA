# 4-final-directional-effect-scoring

## Overview

This folder compiles **directional effect scores** (Positive, Neutral, Negative) from individually scored expert review workbooks as part of the Caribbean Climate Vulnerability Assessment (CVA). Each reviewer independently scored fish stocks using a standardized Excel template, and this code extracts and consolidates those scores into a single, analysis-ready table.

Directional effect scores characterize whether climate change is expected to have a positive, neutral, or negative impact on a given fish stock. Each score is read from a fixed cell range (M38:M40) within species-specific tabs of each reviewer's workbook.

## Final Products

| File | Description |
|------|-------------|
| `table_directional_effect_scores.csv` | Long-format table of all directional effect scores by stock, reviewer, and effect direction |
| `n_reviews_directional-effect.csv` | Summary of how many reviewers scored each stock, with reviewer names listed |

## Setup

### Dependencies

- R (≥ 4.0)
- R packages: `readxl`, `xml2`, `dplyr`, `stringr`, `tidyr`, `tibble`, `openxlsx`

### Input Directory Structure
4-final-directional-effect-scoring/
├── data/
│   └── final-scores/          # Reviewer workbooks (*.xlsx)
│                                # Filename convention: ..._ScorerName.xlsx
├── outputs/
│   └── final-scores-compiled/  # Created automatically by the script
├── 1-extract-directional-effect.R
└── README.md

Each input workbook contains one tab per fish stock. Tabs named `Instructions`, `Data Quality`, and `Example` are excluded automatically. The scorer's name is parsed from the last underscore-delimited segment of the filename (e.g., `Caribbean_CVA_Jones.xlsx` → Scorer = `Jones`).

### Output Directory Structure
outputs/
└── final-scores-compiled/
├── table_directional_effect_scores.csv
└── n_reviews_directional-effect.csv

## Workflow and Code Description

### `1-extract-directional-effect.R`

1. **Discovery** — Scans `./data/final-scores/` for all Excel workbooks and identifies visible, non-template sheet tabs within each file.

2. **Extraction** — For each species tab in each workbook, reads cells M38:M40, which correspond to the Positive, Neutral, and Negative directional effect scores. A safe-read wrapper handles edge cases such as empty tabs, locked files, or malformed ranges by returning `NA` rather than halting execution.

3. **Compilation** — Row-binds all per-tab results into a single long-format table (`directional_effect_table_all`) with columns for source file, scorer, stock name, row index, directional effect category, and score.

4. **QA/QC** — Reports the number of unique scorers, stocks, and directional effect levels. Produces a per-stock summary of reviewer coverage to flag stocks with missing or incomplete reviews.

5. **Export** — Writes the full scores table and the reviewer summary table to `./outputs/final-scores-compiled/`.

Key design choices:
- Workbooks are copied to a temporary file before reading to avoid conflicts with OneDrive or Excel file locks.
- Tabs that fail to read are skipped individually so that one problematic tab does not prevent the rest of a workbook from being processed.
