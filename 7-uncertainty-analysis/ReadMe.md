# 1-extract-tally-scores.R

## Purpose                                                                                                                         
Reads all reviewer scoring workbooks for the Caribbean CVA and extracts the **FINAL SCORE tally columns** from every stock sheet. The output is a set of long-format tally tables that serve as the primary inputs for the bootstrap uncertainty analysis and leave-one-out influence analysis in subsequent scripts. Quantitative exposure factor scores are not extracted here because they are not tally-based and are held fixed during uncertainty analyses under the NOAA FCVA workflow.

## Directory structure
  ```
  project root/
  ├── data/
  │   └── final-scores/          # Input: one .xlsx workbook per reviewer
  └── outputs/
      └── analyses/
          └── 1-inputs/          # Output: tally tables (CSV)
              ├── qa/            # Output: QA diagnostic tables (CSV)
              └── figures/       # Output: exploratory figures (PNG)
```

## Input files

  | File | Description |
  |------|-------------|
  | `data/final-scores/*.xlsx` | One workbook per reviewer. Each workbook contains one sheet per stock, plus non-stock
  
  Sheets without data (`Instructions`, `Data Quality`, `Example`) are excluded by the code. 

  ### Workbook layout

  | Row(s) | Column(s) | Content |
  |--------|-----------|---------|
  | 17–19  | C         | Qualitative Exposure attribute names |
  | 21–28  | B         | Sensitivity attribute names |
  | 30–35  | B         | Rigidity attribute names (relabelled as Sensitivity) |
  | 38–40  | B         | Directional effect labels (Positive / Neutral / Negative) |
  | 17–35  | K         | Final attribute score |
  | 17–35  | L         | Data quality index |
  | 17–35  | M–P       | Rank tallies: Low (1), Moderate (2), High (3), Very High (4) |
  | 38–40  | M         | Directional effect tally |

  ---

  ## Workflow

  ### Step 1 — Extract tally blocks from all workbooks

  Loops over every reviewer workbook and every stock sheet. For each sheet, calls
  two extraction functions:

  - **`extract_rank_tally_block()`** — reads one block of FINAL SCORE rank tally
    rows. Each tally column (M–P) is read individually via `safe_read_cell()` to
    preserve absolute column alignment; reading a range at once causes
    `readWorkbook()` to silently drop empty columns and shift values left.
    Called three times per sheet: Qualitative Exposure (rows 17–19), Sensitivity
    (rows 21–28), and Rigidity (rows 30–35, relabelled as Sensitivity).

  - **`extract_directional_effect_block()`** — reads the three directional effect
    tally rows (38–40) from column M.

  Results are appended into three running long-format tables:
  `sensitivity_tallies_long`, `qualitative_exposure_tallies_long`, and
  `directional_effect_tallies_long`.

  ### Step 2 — Cleanup

  Applied to all three tables in sequence:

  1. **Drop rows with missing attribute or effect-category names** — removes blank
     rows returned for empty worksheet regions.
  2. **Drop fully empty rows** — removes rows where `final_score`,
     `data_quality_index`, and all four tally bins are all `NA` and `n_tallies`
     equals zero. For directional effect, drops rows where `tally` is `NA`.
  3. **Remove unscored stock × reviewer combinations** — the loop reads every
     sheet in every workbook, including sheets the reviewer did not score.
     Combinations where the sum of `n_tallies` across all rows is zero are
     identified and removed via `anti_join()`.
  4. **Drop redundant columns** — `sheet_name` (identical to `stock_name`) and
     `source_file` (redundant with `reviewer_id`) are removed.
  5. **Shorten `reviewer_id` to initials** — the full workbook filename stem (e.g.
     `Caribbean CVA Scoring Template_2025_AAcosta`) is reduced to the leading
     uppercase letters of the name segment after the final underscore (e.g. `AA`).

  A combined table `all_qualitative_tallies_long` is then built by row-binding the cleaned qualitative exposure and sensitivity tables.

  ### Step 3 — QA checks

  Four diagnostic tables are produced and written to `outputs/analyses/1-inputs/qa/`:

  - **3A Row-level tally checks** — flags rows where `n_tallies ≠ 5`, any tally
    bin is `NA`, or duplicate stock × reviewer × attribute combinations exist.
  - **3B Reviewer × stock coverage matrix** — counts distinct attributes scored
    per reviewer × stock combination; pivoted wide for easy inspection.
  - **3C Attribute-level pooled tally summary** — pools tallies across all
    reviewers for each stock × attribute, calculates a pooled weighted mean score,
    and flags combinations where the pooled tally sum does not equal
    `n_reviewers × 5`.
  - **3D Directional effect summary** — totals directional effect tallies per
    stock across all reviewers.

  ### Step 4 — Write output tables

  Writes the four final tally tables to `outputs/analyses/1-inputs/`.

  ---

  ## QA checks

  | Check | Table | Flag condition |
  |-------|-------|---------------|
  | Row tally sum | `qa_tally_row_checks.csv` | `n_tallies ≠ 5` |
  | NA tally bins | `qa_tally_row_checks.csv` | Any of `tally_L`, `tally_M`, `tally_H`, `tally_VH` is `NA` |
  | Duplicate rows | `qa_tally_row_checks.csv` | More than one row per stock × reviewer × attribute |
  | Pooled tally mismatch | `qa_attribute_tally_summary.csv` | Pooled sum ≠ `n_reviewers × 5` |

  Console messages report counts for all flag conditions after Step 2 cleanup and
  again after Step 3A–3C.

  ---

  ## Output files

  ### Tally tables — `outputs/analyses/1-inputs/`

  | File | Rows | Key columns |
  |------|------|-------------|
  | `table_sensitivity_tallies_long.csv` | One row per reviewer × stock × sensitivity attribute | `stock_name`, `reviewer_id`, `attribute_type`, `attribute_name`, `row_num`, `final_score`, `data_quality_index`, `tally_L`, `tally_M`, `tally_H`, `tally_VH`, `n_tallies` |
  | `table_qualitative_exposure_tallies_long.csv` | One row per reviewer × stock × qualitative exposure attribute | Same columns as above |
  | `table_all_qualitative_tallies_long.csv` | Row-bind of sensitivity and qualitative exposure tables | Same columns as above |
  | `table_directional_effect_tallies_long.csv` | One row per reviewer × stock × effect category | `stock_name`, `reviewer_id`, `effect_category`, `row_num`, `tally`, `n_tallies` |  


  ### QA tables — `outputs/analyses/1-inputs/qa/`

  | File | Contents |
  |------|----------|
  | `qa_tally_row_checks.csv` | Problem rows only: `n_tallies ≠ 5`, NA tally bins |
  | `qa_reviewer_stock_coverage.csv` | Wide matrix: stocks × reviewers, cells = attributes scored |
  | `qa_attribute_tally_summary.csv` | Pooled tallies and weighted mean score per stock × attribute |
  | `qa_directional_effect_summary.csv` | Total tally per stock across all reviewers |
