# 4-final-directional-effect-scoring

## Overview

This folder compiles **directional effect scores** (Positive, Neutral, Negative) from individually scored expert review workbooks as part of the Caribbean Climate Vulnerability Assessment (CVA). Each reviewer independently scored fish stocks using a standardized Excel template, and this code extracts and consolidates those scores into a single, analysis-ready table.

Directional effect scores characterize whether climate change is expected to have a positive, neutral, or negative impact on a given fish stock. Each score is read from a fixed cell range (M38:M40) within species-specific tabs of each reviewer's workbook.

## Final Products

| File | Description |
|------|-------------|
| `table_directional_effect.csv` | Long-format table of all directional effect scores by stock, reviewer, and effect direction |
| `n_reviews_directional-effect.csv` | QA summary of how many reviewers scored each stock|
| `directional_effect_wide_all.csv` | Reviewer-level wide summary with Positive, Neutral, and Negative tallies, total tallies, weighted average, and overall directional effect |
|**`directional_effect_summary_by-stock.csv`** | **Stock-level summary combining raw tallies across all reviewers and assigning one final directional effect per stock** |


## Setup

### Input Directory Structure
```text
project/
├─ data/
│  └─ final-scores/             ## Final Caribbean CVA reviewer Excel workbooks
│     ├─ [workbook_1].xlsx
│     ├─ [workbook_2].xlsx
│     └─ ...
├─ 4-final-directional-effect-scoring/
│  ├─ 1-extract-directional-effect.R
│  └─ 2-summarize-directional-effect.R
└─ outputs/
   └─ final-scores-compiled/           
      └─ directional-effect/ ## Compiled tables written by this workflow
```

Each input workbook contains one tab per fish stock. Tabs named `Instructions`, `Data Quality`, and `Example` are excluded automatically. The scorer's name is parsed from the last underscore-delimited segment of the filename (e.g., `Caribbean_CVA_Jones.xlsx` → Scorer = `Jones`).

### Output Directory Structure
```
project/
└─ outputs/
   └─ final-scores-compiled/
      ├─ n_reviews_directional-effect.csv
      ├─ table_directional_effect_scores.csv
      ├─ directional_effect_wide_all.csv
      └─ directional_effect_summary_by-stock.csv
```

# Workflow and Code Description

## `1-extract-directional-effect.R`

1. **Discovery** — Scans `./data/final-scores/` for all Excel workbooks and identifies visible, non-template sheet tabs within each file.

2. **Extraction** — For each species tab in each workbook, reads cells M38:M40, which correspond to the Positive, Neutral, and Negative directional effect scores. A safe-read wrapper handles edge cases such as empty tabs, locked files, or malformed ranges by returning `NA` rather than halting execution.

3. **Compilation** — Row-binds all per-tab results into a single long-format table (`directional_effect_table_all`) with columns for source file, scorer, stock name, row index, directional effect category, and score.

4. **QA/QC** — Reports the number of unique scorers, stocks, and directional effect levels. Produces a per-stock summary of reviewer coverage to flag stocks with missing or incomplete reviews.

5. **Export** — Writes the full scores table and the reviewer summary table to `./outputs/final-scores-compiled/`.

Key design choices:
- Workbooks are copied to a temporary file before reading to avoid conflicts with OneDrive or Excel file locks.
- Tabs that fail to read are skipped individually so that one problematic tab does not prevent the rest of a workbook from being processed.

## `2-summarize-directional-effect.R`

This script takes the compiled long-format directional effect table produced by `1-extract-directional-effect.R` and summarizes it at two levels:

1. **Reviewer × stock level** — one overall directional effect per reviewer for each stock  
2. **Stock level** — one overall directional effect per stock after combining raw tallies across all reviewers  

The script reads:

`./outputs/final-scores-compiled/table_directional_effect_scores.csv`

and writes two additional output files back to the same directory.

### Reviewer-level summarization

1. **Reshape to wide format** — The script selects `Scorer`, `stock_name`, `Directional_effect`, and `Directional_score`, then pivots the long table wider so that each `stock_name × Scorer` combination has separate columns for:

   - `Positive`
   - `Neutral`
   - `Negative`

2. **Fill missing tallies** — Missing directional effect categories are replaced with zero using `dplyr::coalesce()`. This ensures that a reviewer who only assigned one or two nonzero categories still has a complete row for weighted averaging.

3. **Calculate total tallies and weighted average** — For each reviewer-stock combination, the script calculates:

   - `n_tallies` = `Positive + Neutral + Negative`
   - `wt_avg` using the Caribbean CVA weighted-average method:

     `wt_avg = ((Negative * -1) + (Neutral * 0) + (Positive * 1)) / n_tallies`

   This converts the three-category directional effect tally into a continuous index ranging from -1 to 1, where:
   - more Negative tallies pull the score downward
   - more Positive tallies pull the score upward
   - Neutral tallies contribute zero

4. **Assign overall reviewer-level directional effect** — The weighted average is converted to one categorical result:

   - `negative` if `wt_avg <= -0.333`
   - `positive` if `wt_avg >= 0.333`
   - `neutral` if `-0.333 < wt_avg < 0.333`
   - `NA` if no tallies are available

5. **QA checks** — The script prints:
   - the full reviewer-level summary table
   - a count of reviewer-level outcomes by category
   - a sorted subset showing `stock_name`, `Scorer`, raw tallies, weighted average, and final category

6. **Export reviewer-level summary** — The reviewer-level wide table is written to:

   `./outputs/final-scores-compiled/directional_effect_wide_all.csv`

### Stock-level summarization

1. **Aggregate raw tallies across reviewers** — The script groups the long-format table by `stock_name` and `Directional_effect`, then sums `Directional_score` across all reviewers to generate total Positive, Neutral, and Negative tallies for each stock.

2. **Reshape to wide format** — The stock-level grouped table is pivoted wider so that each stock has separate columns for:

   - `Positive`
   - `Neutral`
   - `Negative`

   Missing categories are replaced with zero.

3. **Calculate stock-level totals and weighted average** — For each stock, the script calculates:

   - `n_tallies` = total Positive + Neutral + Negative tallies across all reviewers
   - `wt_avg` using the same weighted-average method:

     `wt_avg = ((Negative * -1) + (Neutral * 0) + (Positive * 1)) / n_tallies`

4. **Assign overall stock-level directional effect** — The weighted average is translated into one final stock-level directional effect category using the same thresholds:

   - `negative` if `wt_avg <= -0.333`
   - `positive` if `wt_avg >= 0.333`
   - `neutral` if `-0.333 < wt_avg < 0.333`
   - `NA` if no tallies are available

**Export stock-level summary** — The final stock-level summary is written to:

   `./outputs/final-scores-compiled/directional_effect_summary_by-stock.csv`

### Interpretation of the weighted-average method

This script uses a simple directional scoring framework in which:

- `Negative = -1`
- `Neutral = 0`
- `Positive = +1`

This approach converts reviewer tallies into a single weighted-average score and then classifies the result into a final directional-effect category. It preserves more information than using only the dominant raw tally because it accounts for mixed reviewer tallies and reflects the balance of positive versus negative scoring.

For example:

- all Negative tallies produce `wt_avg = -1`
- all Neutral tallies produce `wt_avg = 0`
- all Positive tallies produce `wt_avg = 1`
- mixed tallies produce intermediate values

### Notes and assumptions

- The reviewer-level summary treats missing directional effect categories as zero rather than missing.
- The stock-level summary is based on **raw tally totals across all reviewers**, not a vote of reviewer-level final categories.
- The same threshold values are used at both reviewer and stock levels to maintain consistency.
- This workflow assumes that `table_directional_effect_scores.csv` has already been created by `1-extract-directional-effect.R`.

### Full output directory structure after both scripts are run

outputs/
└── final-scores-compiled/
    ├── table_directional_effect_scores.csv
    ├── n_reviews_directional-effect.csv
    ├── directional_effect_wide_all.csv
    └── directional_effect_summary_by-stock.csv
