# Module 10 — Final Figures and Results Tables

## Overview

This module is the final output stage of the Caribbean CVA pipeline. The five R scripts here consume pre-computed CSV files produced by Modules 04–09 and generate all publication figures and summary tables. No FCVA scoring logic is applied. All ranks and scores are read directly from existing CSVs. Run all scripts from the RStudio project root (the `.Rproj` file location), which sets the working directory to `C:/Repos/Caribbean-CVA/`.

The shared configuration file `config.R` at the project root is sourced at the top of every script. Changing any threshold in `config.R` (e.g., `rank_threshold`, `borderline_prop`) requires re-running Modules 04–09 before Module 10 figures will reflect the updates. Figure scripts read only pre-computed CSVs and will not recalculate ranks on their own.

---

## Color Scheme

Two types of color assignments are used across the five scripts:

- **Config-based colors** — referenced by variable name from `config.R` (e.g., `rank_colors`, `dir_colors`). Changing the value in `config.R` automatically propagates to all scripts that use the variable.
- **Local overrides** — a named vector or hex string defined directly inside the script, not referencing the `config.R` variable. These exist where the standard palette does not provide sufficient contrast for overlaid stock-name text, or where a visual distinction between figure types is needed.

### Colors defined in `config.R`

| Variable | Purpose | Values |
|----------|---------|--------|
| `rank_colors` | Standard vulnerability rank fill | Low = `green3`, Moderate = `yellow2`, High = `orange2`, Very High = `red3` |
| `dir_colors` | Directional effect fill | Negative = `brown3`, Neutral = `bisque3`, Positive = `turquoise4` |
| `cert_very_high` / `cert_high` / `cert_moderate` | Bootstrap certainty thresholds for text encoding | 0.95 / 0.90 / 0.67 |

### Local color overrides (defined inside individual scripts)

These colors appear in the script as a literal named vector or hex string — without a `rank_colors` or `dir_colors` variable reference — so they are not affected by changes to `config.R`.

| Script | Variable | Figures affected | Values | Reason |
|--------|----------|-----------------|--------|--------|
| Script 1 | `rank_colors` (local, overwrites config value) | Tally bar charts | Low = `green2`, Moderate = `yellow2`, High = `orange1`, Very High = `red2` | Lighter shades for stacked bars; boxplot fill uses a separate continuous gradient matching the config palette |
| Script 1 | Inline gradient in `scale_fill_gradientn` | Boxplots | Continuous `green3` → `yellow2` → `orange2` → `red3` | Gradient scale keyed to attribute median; same hues as config but applied as a continuous scale, not a discrete vector |
| Script 1 | Inline gradient in `scale_fill_gradient` | Reviewer coverage heatmap (QA) | Low = `#d9534f` → High = `#2c7bb6` | Red-to-blue diverging scale to highlight incomplete reviewer × stock combinations |
| Script 3 | `rank_colors` (local, overwrites config value) | DCP column chart and cross-plot | Low = `green3`, Moderate = **`yellow4`**, High = `orange2`, Very High = `red3` | `yellow4` (darker than config's `yellow2`) visually distinguishes DCP figures from the bootstrap uncertainty figures in Script 2 |
| Script 4 | `vuln_colors` (new variable, does not overwrite) | Vulnerability tile grid (Panel A) | Low = `#2d9a27`, Moderate = `yellow3`, High = `#e87722`, Very High = `#cc2222` | Hex values provide sufficient contrast for both black and white overlaid text; config `rank_colors` are not used in Panel A |

---

## Script 1: 1-plot-scoring-distributions.R

**Purpose:** Generate score distribution figures: attribute-level boxplots, per-stock tally bar charts (sensitivity and exposure), a directional effect summary, and a reviewer coverage QA heatmap.

### Inputs

| File | Contents |
|------|---------|
| `outputs/analyses/1-inputs/table_all_qualitative_tallies_long.csv` | Raw per-reviewer qualitative tallies in long format; "Coral cover" rows are filtered out before use |
| `outputs/final-scores-compiled/overall-vulnerability-rankings/attribute_means_uscar.csv` | Mean reviewer-assigned score per stock × attribute (sensitivity and exposure) |
| `outputs/final-tallies-long/sensitivity_tallies_long.csv` | Pooled sensitivity tally counts per stock × attribute (`tally_L`, `tally_M`, `tally_H`, `tally_VH`) |
| `outputs/final-tallies-long/directional_effect_tallies_long.csv` | Pooled directional effect tallies per stock |
| `outputs/final-tallies-long/exposure_tallies_long.csv` | Pooled exposure tally counts per stock × factor |
| `data/exposure-factor-filter-long.csv` | Expert-approved exposure factor × stock pairs; only rows with `include == TRUE` are retained |

### Outputs

**`figures/fig_attribute_score_boxplot_combined_vertical.png`** — Exposure and Sensitivity boxplot panels stacked vertically. Saved at 1200 dpi.

<p align="center">
  <img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/figures/fig_attribute_score_boxplot_combined_vertical.png?raw=true"
       alt="Combined attribute score boxplots (vertical)"
       width="700"/>
</p>

**`figures/fig_exposure_tally_distributions_by_stock.png`** — Faceted horizontal stacked bars (one facet per stock, 5 columns). Within each facet, bars show the proportion of tallies in each vulnerability rank category for each expert-approved exposure factor × stock pair. Factors ordered on the y-axis by ascending pooled mean score. Uses local lighter tally-bar palette. Saved at 8.5 × 11 in, 900 dpi.

<p align="center">
  <img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/figures/fig_exposure_tally_distributions_by_stock.png?raw=true"
       alt="Exposure tally distributions by stock"
       width="700"/>
</p>

**`figures/fig_sensitivity_tally_distributions_by_stock.png`** — Same layout for biological sensitivity attributes. Up to 20 pooled tallies (4 reviewers × 5 tallies) per attribute per stock. Saved at 8.5 × 11 in, 900 dpi.

<p align="center">
  <img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/figures/fig_sensitivity_tally_distributions_by_stock.png?raw=true"
       alt="Sensitivity tally distributions by stock"
       width="700"/>
</p>

**`figures/fig_directional_effect_summary.png`** — Horizontal stacked bar per stock showing proportions of reviewer tallies assigned Negative, Neutral, or Positive directional effect. 16 pooled tallies per stock (4 reviewers × 4 tallies). Stocks ordered by ascending proportion negative (stock with lowest negative proportion at bottom). Uses `dir_colors` from `config.R`. Saved at 1200 dpi.

<p align="center">
  <img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/figures/fig_directional_effect_summary.png?raw=true"
       alt="Directional effect summary by stock"
       width="700"/>
</p>

**`figures/fig_reviewer_stock_coverage.png`** — QA heatmap of reviewer × stock showing the count of attributes scored per combination. Red (#d9534f) indicates incomplete scoring; blue (#2c7bb6) indicates a fully scored combination. Used to verify data collection completeness. Saved at 300 dpi.

<p align="center">
  <img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/figures/fig_reviewer_stock_coverage.png?raw=true"
       alt="Reviewer × stock coverage heatmap"
       width="700"/>
</p>

### Workflow
1. Source `config.R` for rank/directional effect levels, `rank_colors`, `dir_colors`, `attr_short_names`, `exp_attr_short_names`, and `stock_name_recode`.
2. Filter raw tallies to remove "Coral cover" rows; apply canonical stock name recodes.
3. Filter attribute means and exposure tallies to expert-approved factor × stock pairs only.
4. Build sensitivity and exposure boxplots using `scale_fill_gradientn` with an inline green-to-red gradient; attributes ordered by ascending median.
5. Locally override `rank_colors` with a lighter palette for the tally stacked bar charts (green2, yellow2, orange1, red2).
6. Pool tallies across reviewers; compute per-stock proportions; generate faceted tally bar charts for sensitivity and exposure.
7. Generate directional effect summary bar chart (uses `dir_colors` from config).
8. Generate reviewer × stock coverage heatmap.

---

## Script 2: 2-plot-uncertainty-figures.R

**Purpose:** Visualize the leave-one-out (LOO) influence of each attribute and factor on vulnerability rank, and the bootstrap-resampled rank distributions for vulnerability and directional effect.

### Inputs

| File | Contents |
|------|---------|
| `outputs/analyses/uncertainty-loo/final-tables/table_leave_one_out_sensitivity_summary.csv` | Per-attribute count of vulnerability rank changes across all stocks when that attribute is omitted |
| `outputs/analyses/uncertainty-loo/final-tables/table_leave_one_out_exposure_summary.csv` | Same metric for each exposure factor |
| `outputs/analyses/uncertainty-loo/final-tables/table_bootstrap_uncertainty_stock.csv` | Proportion of 10,000 bootstrap iterations in each vulnerability rank per stock |
| `outputs/analyses/uncertainty-loo/final-tables/table_directional_effect_bootstrap.csv` | Proportion of 10,000 bootstrap iterations in each directional effect category per stock |
| `outputs/final-scores-compiled/overall-vulnerability-rankings/overall_vulnerability_scores_uscar.csv` | Baseline vulnerability ranks — used to set stock display order |

### Outputs

**`figures/fig_loo_bar_plots.png`** — **Panel A:** Horizontal bar chart showing the number of stocks (out of 25) whose vulnerability rank changed when each exposure factor was omitted, ordered by descending influence (most influential at top). **Panel B:** Same for each biological sensitivity attribute. Bars are solid black. X-axes scale independently. Panels combined via patchwork with heights proportional to bar count (13 exposure factors, 14 sensitivity attributes). Saved at 6.5 × 9 in, 300 dpi.

<p align="center">
  <img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/figures/fig_loo_bar_plots.png?raw=true"
       alt="Leave-one-out influence analysis"
       width="700"/>
</p>

**`figures/fig_bootstrap_uncertainty.png`** — **Panel A:** Horizontal stacked bars showing proportion of 10,000 bootstrap iterations in each vulnerability rank (Low/Moderate/High/Very High) per stock. Dashed vertical line at 75% borderline threshold. Borderline stocks (dominant proportion < 75%) are marked with an asterisk (*). Right-margin filled square shows finalized baseline rank. Gray horizontal lines separate baseline vulnerability rank groups (Moderate/High/Very High, bottom to top). **Panel B:** Same layout for directional effect bootstrap (Negative/Neutral/Positive). Stocks ordered top-to-bottom by descending baseline vulnerability rank; within each rank group, stocks sorted by ascending dominant bootstrap proportion so the most uncertain stocks cluster at rank-group boundaries. Panels use equal heights (25 stocks each). Saved at 7.5 × 12 in, 300 dpi.

<p align="center">
  <img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/figures/fig_bootstrap_uncertainty.png?raw=true"
       alt="Bootstrap resampling uncertainty"
       width="700"/>
</p>

### Colors
- LOO bars: solid `black` fill — no rank-coded palette
- Bootstrap vulnerability bars (Panel A): `rank_colors` from `config.R`
- Bootstrap directional effect bars (Panel B): `dir_colors` from `config.R`

### Workflow
1. Source `config.R` for `rank_levels`, `rank_colors`, `dir_levels`, `dir_colors`, `attr_short_names`, `exp_attr_short_names`.
2. Recode "Stock size/status" → "Stock Size Status" before applying `attr_short_names` lookup (CSV column name differs from config key).
3. Build LOO figure: one panel each for exposure and sensitivity, ordered by descending influence.
4. Compute stock ordering for bootstrap figure: primary sort by baseline vulnerability rank (Very High → Low, top to bottom); secondary sort by dominant bootstrap proportion ascending within each rank group.
5. Build bootstrap vulnerability panel with baseline rank indicator squares and 75% reference line.
6. Build bootstrap directional effect panel using the same stock order as the vulnerability panel for row-for-row comparison.

---

## Script 3: 3-plot-distributional-change.R

**Purpose:** Produce two figures showing the potential for distributional change (DCP) and comparing DCP against overall vulnerability.

### Inputs

| File | Contents |
|------|---------|
| `outputs/distribution-change-potential/distributional_change_full_uscar.csv` | DCP rank and bootstrap dominant proportion per stock (output from Module 09) |
| `outputs/final-scores-compiled/overall-vulnerability-rankings/overall_vulnerability_scores_uscar.csv` | Baseline overall vulnerability rank per stock — joined to add `Vuln_rank` to DCP data |

### Outputs

**`figures/fig_distributional_change_ranks.png`** — One solid-colored column per DCP rank (Low/Moderate/High/Very High). Column height = number of stocks in that category. Stock names stacked inside each column, ordered by bootstrap certainty (highest certainty at top), then alphabetically within each certainty tier. Label format: "Stock name (V)" where V = vulnerability abbreviation (L/M/H/VH). Text color/face encodes bootstrap certainty: very high (>95%) = black bold, high (90–95%) = black italic, moderate (67–89%) = white bold, low (<67%) = white italic. Y-axis ceiling rounded up to the nearest multiple of 3. Saved at 10 × 8 in, 1200 dpi.

<p align="center">
  <img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/figures/fig_distributional_change_ranks.png?raw=true"
       alt="Distributional change potential by rank"
       width="700"/>
</p>

**`figures/fig_distributional_change_vs_vulnerability.png`** — 4×4 tile grid with DCP rank on the x-axis and overall vulnerability rank on the y-axis. Cells are filled by overall vulnerability rank. Stock names are placed within each cell with y-offsets to separate multiple stocks in the same cell. Same certainty-encoded text color and face as the column chart. Saved at 9 × 7 in, 1000 dpi.

<p align="center">
  <img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/figures/fig_distributional_change_vs_vulnerability.png?raw=true"
       alt="Distributional change potential vs. overall vulnerability"
       width="700"/>
</p>

### Colors
- Column and tile fill: local `rank_colors` override — Low = `green3`, Moderate = **`yellow4`**, High = `orange2`, Very High = `red3`. Only Moderate differs from `config.R` (yellow4 is darker than yellow2), deliberately chosen so the DCP figures are visually distinct from the bootstrap uncertainty bars in Script 2.
- Text color: black (very high / high certainty) or white (moderate / low certainty), derived from `dominant_prop` thresholds in `config.R`.

### Workflow
1. Source `config.R` for `rank_levels`, certainty thresholds (`cert_very_high`, `cert_high`, `cert_moderate`), and `stock_name_recode`.
2. Read DCP and vulnerability CSVs; join on `stock_name`; apply canonical stock name recode after the join.
3. Derive certainty category from `dominant_prop` and assign `text_color` / `text_face` accordingly.
4. Build column chart: compute y-positions for text labels within each column (most certain stock at top row).
5. Build cross-plot: assign numeric grid positions (1–4) for both axes; compute within-cell y-offsets for stocks sharing a cell.

---

## Script 4: 4-plot-overall-vulnerability.R

**Purpose:** Produce the main results figure showing overall climate vulnerability on a 2-D Exposure × Sensitivity grid and the dominant directional effect for all 25 stocks.

### Inputs

| File | Contents |
|------|---------|
| `outputs/final-scores-compiled/overall-vulnerability-rankings/overall_vulnerability_scores_uscar.csv` | Exposure, sensitivity, and vulnerability ranks and scores per stock |
| `outputs/final-tallies-long/directional_effect_tallies_by_stock.csv` | Directional effect tally counts and dominant direction per stock |
| `outputs/analyses/uncertainty-loo/final-tables/table_bootstrap_uncertainty_stock.csv` | Vulnerability bootstrap proportions — used to derive certainty for Panel A text encoding |
| `outputs/analyses/uncertainty-loo/final-tables/table_directional_effect_bootstrap.csv` | Directional effect bootstrap proportions — used to derive certainty for Panel B text encoding |

### Output

**`figures/fig_overall_vulnerability.png`** — **Panel A:** 4×4 tile grid (x = Climate Exposure rank, y = Biological Sensitivity rank). Background tile fill reflects the expected FCVA vulnerability rank for each Exposure × Sensitivity combination, calculated as the product of numeric rank values (score = exp_num × sens_num: ≤3 = Low, ≤6 = Moderate, ≤9 = High, >9 = Very High). The Very High sensitivity row is rendered as a thin tile (height 0.3, centered at y = 3.65) because no stocks fall there; all other rows use standard height 1. Stock names placed in cells with y-offsets for multiple stocks. **Panel B:** Three columns (Negative, Neutral, Positive) with stock names stacked in each directional effect category, ordered by descending bootstrap certainty. Both panels encode certainty via text color and face: very high (>95%) = black bold, high (90–95%) = black italic, moderate (67–89%) = white bold, low (<67%) = white italic. Panel A uses vulnerability bootstrap certainty; Panel B uses directional effect bootstrap certainty. Panels combined via patchwork (Panel A height 2×, Panel B height 1×). Saved at 8 × 11.5 in, 1000 dpi.

<p align="center">
  <img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/figures/fig_overall_vulnerability.png?raw=true"
       alt="Overall climate vulnerability and directional effect"
       width="700"/>
</p>

### Colors
- Panel A tile fill: local `vuln_colors` — Low = `#2d9a27`, Moderate = `yellow3`, High = `#e87722`, Very High = `#cc2222`. These values (not `rank_colors` from `config.R`) were selected to provide adequate contrast for overlaid white and black stock name text against the tile background.
- Panel B tile fill: `dir_colors` from `config.R` (Negative = `brown3`, Neutral = `bisque3`, Positive = `turquoise4`).
- Text color: black or white, determined by certainty tier; thresholds (`cert_very_high`, `cert_high`, `cert_moderate`) sourced from `config.R`.

### Workflow
1. Source `config.R` for rank levels, certainty thresholds, `dir_colors`, `dir_levels`, `stock_name_recode`.
2. Join vulnerability scores with bootstrap data; extract each stock's dominant vulnerability bootstrap proportion (rows where `vuln_rank == Vuln_rank`); derive certainty and text encoding.
3. Compute within-cell y-offsets for stocks sharing a grid cell (sorted by descending bootstrap proportion).
4. Build 4×4 background tile grid using FCVA product-score logic; render Very High sensitivity row as a thin tile.
5. Build Panel A (vulnerability grid) with `vuln_colors` tile fill and certainty-encoded stock name text.
6. Join directional effect tally data with directional bootstrap; derive directional certainty and text encoding.
7. Build Panel B (directional effect) with `dir_colors` tile fill and directional certainty-encoded text.
8. Combine panels with patchwork and save.

---

## Script 5: 5-produce-results-tables.R

**Purpose:** Generate two CSV results tables (directional effect and data quality) ready for inclusion in the manuscript.

### Inputs

| File | Contents |
|------|---------|
| `outputs/final-scores-compiled/overall-vulnerability-rankings/overall_vulnerability_scores_uscar.csv` | Vulnerability ranks and scores per stock |
| `outputs/final-tallies-long/directional_effect_tallies_by_stock.csv` | Raw tally counts and dominant directional effect per stock |
| `outputs/analyses/uncertainty-loo/final-tables/table_directional_effect_bootstrap.csv` | Bootstrap proportions (Negative/Neutral/Positive) per stock, long format |
| `outputs/analyses/uncertainty-loo/final-tables/table_bootstrap_uncertainty_stock.csv` | Vulnerability bootstrap proportions per stock — used to derive the stock sort order |
| `outputs/final-scores-compiled/data-quality/overall_data_quality_summary_by_stock.csv` | Data quality rank and score distribution per stock |

### Outputs

| File | Columns |
|------|---------|
| `outputs/tables/table_directional_effect_results.csv` | Stock, Vuln_rank, Dir_effect, Wt_mean, N_negative, N_neutral, N_positive, N_tallies, Boot_Negative, Boot_Neutral, Boot_Positive, Dominant_prop, Borderline |
| `outputs/tables/table_data_quality_results.csv` | Stock, Vuln_rank, Data_quality_rank, Prop_ge_2, Mean_score, N_adequate, N_limited, N_expert, N_nodata |

Both tables are sorted by descending vulnerability rank, then descending dominant vulnerability bootstrap proportion within each rank group, then alphabetical by stock name.

### Workflow
1. Source `config.R` for `rank_levels`, `borderline_prop`, and `stock_name_recode`.
2. Read all five input CSVs; apply canonical stock name recode to each.
3. Derive sort key: for each stock, find the dominant vulnerability bootstrap proportion (maximum `prop` across all rank bins in `table_bootstrap_uncertainty_stock.csv`).
4. Build directional effect table: pivot bootstrap from long to wide (Boot_Negative, Boot_Neutral, Boot_Positive columns); compute `Dominant_prop` and `Borderline` flag; join vulnerability scores and tally counts; apply sort order; rename columns.
5. Build data quality table: join data quality summary with vulnerability scores and sort key; apply sort order; rename columns.
6. Write both CSVs to `outputs/tables/`.
