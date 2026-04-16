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

---

# 2-plot-figures.R

## Purpose
Produces all primary uncertainty and distribution figures for the Caribbean CVA. The script visualizes the spread of reviewer-assigned vulnerability scores and LMHV tally distributions across the 25 assessed stocks, for both biological sensitivity attributes and exposure factors (qualitative + quantitative). A combined score-distribution panel figure (Figure 1) is the primary publication output; the tally distribution figures (Figure 3) provide per-stock diagnostic detail.

## Directory structure
```
project root/
├── outputs/
│   ├── analyses/
│   │   └── 1-inputs/              # Input: tally tables from 1-extract-tally-scores.R
│   └── final-scores-compiled/
│       ├── overall-vulnerability-rankings/
│       │   └── attribute_means_uscar.csv   # Input: per-stock attribute mean scores
│       └── quantitative-exposure-attribute-scores-all.csv  # Input: CMIP6 factor scores
└── figures/                       # Output: all PNG figures
```

## Inputs

| File | Description |
|------|-------------|
| `outputs/analyses/1-inputs/table_sensitivity_tallies_long.csv` | Reviewer-level tallies for 14 biological sensitivity attributes |
| `outputs/analyses/1-inputs/table_qualitative_exposure_tallies_long.csv` | Reviewer-level tallies for 2 qualitative exposure attributes (Coral cover excluded) |
| `outputs/analyses/1-inputs/table_all_qualitative_tallies_long.csv` | Row-bind of sensitivity and qualitative exposure tally tables |
| `outputs/analyses/1-inputs/table_directional_effect_tallies_long.csv` | Reviewer directional effect tallies (Positive / Neutral / Negative) per stock |
| `outputs/final-scores-compiled/overall-vulnerability-rankings/attribute_means_uscar.csv` | Per-stock × attribute mean scores (U.S. Caribbean spatial extent); `attribute_type` ∈ `"Sensitivity"`, `"Exposure"` |
| `outputs/final-scores-compiled/quantitative-exposure-attribute-scores-all.csv` | LMHV grid-cell tally counts for 13 CMIP6 exposure factors × 25 stocks × 3 spatial extents |

## Workflow

### Data carpentry

**Stock name standardization.** A canonical named vector `stock_name_recode` maps legacy lowercase stock names (e.g., `"Atlantic thread herring"`, `"Long-spined sea urchin"`) to title-case canonical names (e.g., `"Atlantic Herring"`, `"Diadema"`). This recode is applied via `recode(stock_name, !!!stock_name_recode)` to all five input tables on read, ensuring consistent joins and facet labels throughout.

**Coral cover exclusion.** Rows with `attribute_name == "Coral cover"` are dropped from `all_qualitative_tallies_long` and from the qualitative exposure conform step. This attribute was assessed by only a small subset of reviewers and is excluded from all tally distribution figures.

**Exposure tally table.** A combined `exposure_tallies_long` table is built by conforming and row-binding two sources:
- *Qualitative exposure* (`qual_exp_conform`): tallies are summed across reviewers within each stock × attribute, producing one row per stock × attribute with pooled `tally_L / tally_M / tally_H / tally_VH` and `n_tallies`.
- *Quantitative exposure* (`quant_exp_conform`): filtered to `spatial_extent == "U.S. Caribbean"`, with `full_names` renamed to `attribute_name` and `attribute_type` set to `"Quantitative Exposure"`.

The resulting table has one row per stock × exposure factor (2 qualitative + 13 quantitative = 15 attributes × 25 stocks).

**Short display labels.** Two named lookup vectors (`attr_short_names`, `exp_attr_short_names`) map full attribute names to abbreviated display labels (≤ 16 characters) used on all figure y-axes. These are defined once in the data prep section and referenced by all figures.

---

### Figure 1 — Score distributions (combined panel)

Two horizontal boxplot figures are produced and then combined into a single two-panel figure using `patchwork`.

**Figure 1A — Biological sensitivity attributes** (`fig_sensitivity_attribute_score_boxplot.png`):
- One box per sensitivity attribute (14 attributes), ordered by median score ascending (bottom to top).
- x-axis: mean score across stocks (1 = Low, 4 = Very High); data from `attr_means` filtered to `attribute_type == "Sensitivity"`.
- Fill color: each box is filled by its attribute median score, mapped through the LMHV gradient (`green3` → `yellow2` → `orange2` → `red3`, anchored at 1–4).

**Figure 1B — Exposure factors** (`fig_exposure_attribute_score_boxplot.png`):
- Same structure as Figure 1A, filtered to `attribute_type == "Exposure"` (15 attributes).
- Attributes cover both qualitative (expert-scored) and quantitative (CMIP6-derived) factors; both use the 1–4 LMHV scale.

**Combined panel** (`fig_attribute_score_boxplot_combined.png`):
- Figure 1A (panel A) and Figure 1B (panel B) placed side by side at 14" × height.
- Panel tags added via `plot_annotation(tag_levels = "A")`.

<img src="https://raw.githubusercontent.com/holden-harris/Caribbean-CVA/main/7-uncertainty-analysis/figures/fig_attribute_score_boxplot_combined.png" width="800"/>

---

### Figure 2 — Directional effect summary by stock

(`fig_directional_effect_summary.png`)

Horizontal stacked bar chart showing the proportion of reviewer tallies classified as Positive, Neutral, or Negative for each stock. Stocks are ordered by proportion Negative (ascending). Colors: Positive = `#2c7bb6`, Neutral = `bisque`, Negative = `indianred4`. Data source: `directional_effect_tallies_long`, pooled across all reviewers per stock.

<img src="https://raw.githubusercontent.com/holden-harris/Caribbean-CVA/main/7-uncertainty-analysis/figures/fig_directional_effect_summary.png" width="600"/>

---

### Figure 3A — Sensitivity tally distributions by stock

(`fig_sensitivity_tally_distributions_by_stock.png`)

Faceted 5 × 5 panel figure (one panel per stock). Each panel shows one horizontal stacked bar per sensitivity attribute, ordered by pooled mean score ascending. Bar segments show the proportion of reviewer tallies in each LMHV category (Low / Moderate / High / Very High). Y-axis labels use abbreviated names from `attr_short_names`. Colors: Low = `green3`, Moderate = `yellow2`, High = `orange2`, Very High = `red3`. Saved at 12" × 12", 1200 dpi.

<img src="https://raw.githubusercontent.com/holden-harris/Caribbean-CVA/main/7-uncertainty-analysis/figures/fig_sensitivity_tally_distributions_by_stock.png" width="800"/>

---

### Figure 3B — Exposure tally distributions by stock

(`fig_exposure_tally_distributions_by_stock.png`)

Same structure as Figure 3A, using `exposure_tallies_long` (15 exposure attributes per stock). Tally units differ by attribute type — qualitative tallies count reviewer votes; quantitative tallies count LMHV grid cells — but both are expressed as proportions and are visually comparable. Y-axis labels use abbreviated names from `exp_attr_short_names`. Saved at 12" × 14", 1200 dpi.

<img src="https://raw.githubusercontent.com/holden-harris/Caribbean-CVA/main/7-uncertainty-analysis/figures/fig_exposure_tally_distributions_by_stock.png" width="800"/>

---

### Figure QA — Reviewer × stock coverage heatmap

(`fig_reviewer_stock_coverage.png`)

Tile heatmap with reviewers on x and stocks on y. Each cell is colored by the count of attributes scored (red = few, blue = full coverage) with the count printed in white. Used to identify reviewer × stock combinations with missing or incomplete assessments.

<img src="https://raw.githubusercontent.com/holden-harris/Caribbean-CVA/main/7-uncertainty-analysis/figures/fig_reviewer_stock_coverage.png" width="600"/>

---

## Outputs

### Figures — `figures/`

| File | Figure | Description |
|------|--------|-------------|
| `fig_sensitivity_attribute_score_boxplot.png` | 1A | Sensitivity attribute score distributions (boxplot) |
| `fig_exposure_attribute_score_boxplot.png` | 1B | Exposure factor score distributions (boxplot) |
| `fig_attribute_score_boxplot_combined.png` | 1 (combined) | Panels 1A and 1B side by side |
| `fig_directional_effect_summary.png` | 2 | Directional effect proportions by stock |
| `fig_sensitivity_tally_distributions_by_stock.png` | 3A | Per-stock sensitivity LMHV tally distributions |
| `fig_exposure_tally_distributions_by_stock.png` | 3B | Per-stock exposure LMHV tally distributions |
| `fig_reviewer_stock_coverage.png` | QA | Reviewer × stock attribute coverage heatmap |
