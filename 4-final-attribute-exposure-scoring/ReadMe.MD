# 4-final-attribute-exposure-scoring

## Overview

This folder compiles and summarizes final attribute scores from individually scored expert review workbooks for the Caribbean Climate Vulnerability Assessment (CVA). The workflow first extracts final reviewer scores into a single compiled table, then summarizes those scores across attribute types and joins them with calculated quantitative exposure scores for downstream vulnerability scoring.

The first script in this folder focuses on extracting and consolidating final reviewer-entered attribute scores from the completed scoring workbooks.

## Final Products

| File | Description |
|------|-------------|
| `table_final_attribute_scores_all.csv` | Long-format compiled table of final attribute scores from all reviewer workbooks |
| `quantitative-exposure-attribute-scores-all.csv` | Quantitative exposure attribute scores used downstream in the summary step |
| Additional summary outputs | Created by `2-summarize-attribute-scores.R` in the same output directory |

## Setup

### Input Directory Structure

```text
4-final-attribute-exposure-scoring/
├── data/
│   └── final-scores/                     # Reviewer workbooks (*.xlsx)
├── outputs/
│   └── final-scores-compiled/           # Compiled and summary tables
├── 1-extract-final-scores-from-all-reviewers.R
├── 2-summarize-attribute-scores.R
└── README.md
```

### Output Directory Structure
```
outputs/
└── final-scores-compiled/
    ├── table_final_attribute_scores_all.csv
    ├── quantitative-exposure-attribute-scores-all.csv
    └── [additional summary files written by downstream scripts]
```

## Workflow and Code Description

### `1-extract-final-scores-from-all-reviewers.R`

This script compiles final reviewer-entered CVA attribute scores from all completed scoring workbooks into a single standardized, long-format table for downstream analysis.

#### Purpose

The script reads the final reviewer workbooks stored in:

`./data/final-scores/`

and extracts the final attribute scores recorded for each stock and attribute. These outputs are combined into a single analysis-ready table used by the downstream summary script.

#### Main workflow

1. **Workbook discovery**  
   The script scans the input directory for all reviewer scoring workbooks and processes them one at a time.

2. **Sheet identification**  
   Within each workbook, the script identifies stock-specific worksheets and excludes non-scoring template tabs such as instructions or example sheets.

3. **Final score extraction**  
   For each stock tab, the script reads the final reviewer-entered attribute scores from the standardized workbook template and converts them into a consistent tabular format.

4. **Standardization**  
   Extracted rows are standardized to include the following fields:
   - source workbook
   - reviewer name
   - stock name
   - attribute name
   - attribute type
   - final score
   - any row or worksheet identifiers needed for QA/QC

5. **Compilation**  
   All reviewer-stock tables are row-bound into a single long-format dataset containing the full set of final attribute scores across all reviewers.

6. **Quality control**  
   The script performs basic checks to confirm that workbooks, tabs, and extracted scores were read successfully and that the compiled table contains the expected reviewer and stock coverage.

7. **Export**  
   The compiled output is written to:

   `./outputs/final-scores-compiled/table_final_attribute_scores_all.csv`

#### Output table structure

The compiled long-format table is designed to support downstream summarization and contains one row per:

`reviewer × stock × attribute`

with columns describing the reviewer, stock, attribute, attribute type, and final score.

#### Workflow notes

- This script should be run before `2-summarize-attribute-scores.R`.
- Reviewer workbooks should remain unchanged in `./data/final-scores/`.
- The script assumes a consistent workbook structure across reviewers.
- If the workbook template changes, the extraction logic should be reviewed and updated.
- Any stock tabs with missing or malformed score entries should be flagged during QA/QC before downstream summaries are interpreted.

To maintain transparency and reproducibility:

- preserve raw reviewer workbooks in `./data/final-scores/`
- avoid editing compiled CSV outputs by hand
- rerun the extraction step whenever reviewer workbooks are added or revised
