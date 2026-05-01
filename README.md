# Caribbean-CVA

> **Summary:** This repository contains the complete analytical pipeline for a NOAA Fisheries Climate Vulnerability Assessment (CVA) of 25 fish and invertebrate stocks managed in the U.S. Caribbean. Stocks are ranked by overall climate vulnerability — the product of independently scored **sensitivity** (biological susceptibility) and **exposure** (projected habitat change) components — following the FCVA framework of Morrison et al. (2015). The pipeline runs from raw CMIP6 oceanographic projections and IUCN species ranges (Module 2) through expert reviewer workbook extraction (Modules 3–7) to bootstrap and leave-one-out uncertainty analyses (Module 8). All analyses are in R; final figures and ranked scores are written to `figures/` and `outputs/final-scores-compiled/`. Conducted by Harris Analytics & Research LLC in support of [Isla Mar 501c3](https://www.islamar.org/).

---

## Project Overview

This repository contains the full analytical pipeline for a **Climate Vulnerability Assessment (CVA)** of 25 fish and invertebrate stocks managed in the U.S. Caribbean. The CVA evaluates each stock's overall climate vulnerability by combining two independently scored components — **sensitivity** (how biologically susceptible is the stock to climate change?) and **exposure** (how much is the stock's habitat projected to change?) — following the NOAA Fisheries Climate Vulnerability Assessment (FCVA) framework established by Morrison et al. (2015) and applied in recent regional CVAs. 

Analyses were conducted by Harris Analytics & Research LLC in support of [Isla Mar 501c3](https://www.islamar.org/).  

All code and materials are available under an open-access license, as per the Creative Commons CC0 1.0 license.  
We thank Dan Crear (ICATTC) and Tyler Loughran (NOAA) for their assistance in this work.

---

## Workflow Modules

The project is organized into nine numbered workflow modules. Each module has its own subdirectory with scripts and a `ReadMe` file.

| Module | Folder | Purpose | Key outputs |
|--------|--------|---------|-------------|
| 0 | `0-query-species-attributes-from-FishBase/` | Query biological traits and life-history attributes from FishBase via the `rfishbase` R package | `fishbase_species_attributes.csv` |
| 1 | `1-make-species-distribution-maps/` | Generate standardized PNG distribution maps for all 25 species from IUCN shapefiles | `outputs/disbribution-maps/*.png` |
| 2 | `2-exposure-anomalies/` | Calculate CMIP6-based standardized anomaly maps; produce 12-panel exposure-overlap figures reviewed by CVA experts; extract quantitative exposure scores | `outputs/exposure-overlap-12panel/`, `outputs/final-scores-compiled/quantitative-exposure-attribute-scores-all.csv` |
| 3 | `3-prelim-sensitivity-attribute-scoring/` | Extract and summarize preliminary sensitivity-attribute tallies from pre-workshop reviewer workbooks; generate per-stock LMHV summary plots for workshop preparation | `data/preliminary-scores/score_table_all.csv`, `outputs/prework/` |
| 4 | `4-final-attribute-exposure-scoring/` | Extract final reviewer scores, calculate attribute means, apply NOAA FCVA logic model to produce stock-level sensitivity, exposure, and overall vulnerability scores | `outputs/final-scores-compiled/overall-vulnerability-rankings/`, `figures/fig_overall_vulnerability.png` |
| 5 | `5-final-directional-effect-scoring/` | Extract directional-effect tallies (Positive / Neutral / Negative) from final workbooks; calculate stock-level directional-effect index | `outputs/final-scores-compiled/directional-effect/directional_effect_summary_by-stock.csv` |
| 6 | `6-final-data-quality-scoring/` | Extract reviewer data-quality scores (0–3) for each attribute; summarize and rank overall data quality per stock | `outputs/final-scores-compiled/data-quality/overall_data_quality_summary_by_stock.csv` |
| 7 | `7-scoring-distributions/` | Extract LMHV tally distributions from final workbooks; produce attribute score boxplots and tally-distribution figures | `outputs/final-tallies-long/`, `figures/fig_*_tally_distributions_by_stock.png`, `figures/fig_*_attribute_score_boxplot.png` |
| 8 | `8-uncertainty-analysis/` | Bootstrap resampling and leave-one-out influence analyses; quantify statistical robustness of final vulnerability rankings | `outputs/analyses/uncertainty-loo/`, `figures/fig_bootstrap_uncertainty.png`, `figures/fig_loo_bar_plots.png` |

---

## Directory Structure

```
Caribbean-CVA/
│
├── 0-query-species-attributes-from-FishBase/
│   ├── Query-species-attributes-from-FishBase.R
│   └── ReadMe.md
│
├── 1-make-species-distribution-maps/
│   ├── Make-species-distribution-maps.R
│   └── ReadMe.md
│
├── 2-exposure-anomalies/
│   ├── Exposure-anomalies.R
│   └── ReadMe.md
│
├── 3-prelim-sensitivity-attribute-scoring/
│   ├── 1-extract-scores.R
│   ├── 2-prework-scoring-summaries.R
│   └── ReadMe.md
│
├── 4-final-attribute-exposure-scoring/
│   ├── 1-extract-final-scores-from-all-reviewers.R
│   ├── 2-summarize-attribute-scores.R
│   ├── 3-calculate-overall-vulnerability-scores.R
│   ├── 4-plot-overall-vulnerability.R
│   └── ReadMe.MD
│
├── 5-final-directional-effect-scoring/
│   ├── 1-extract-directional-effect.R
│   ├── 2-summarize-directional-effect.R
│   └── ReadMe.MD
│
├── 6-final-data-quality-scoring/
│   ├── 1-extract-data-quality-scores.R
│   ├── 2-summarize-data-quality-scores.R
│   └── ReadMe.md
│
├── 7-scoring-distributions/
│   ├── 1-extract-tally-scores.R
│   ├── 2-finalize-tally-tables.R
│   ├── 3-plot-figures.R
│   └── ReadMe.md
│
├── 8-uncertainty-analysis/
│   ├── uncertainty-analyses.R
│   ├── 2-make-uncertainty-figures.R
│   └── ReadMe.MD
│
├── data/
│   ├── cmip6/                                       # CMIP6 NetCDF exposure files (*.nc)
│   ├── species-distribution-shapefiles/             # IUCN species range polygons (*.shp + sidecars)
│   ├── master-lists/                                # Reference lists
│   ├── preliminary-scores/                          # Pre-workshop reviewer workbooks
│   ├── final-scores/                                # Final reviewer workbooks (*.xlsx), one per reviewer
│   ├── attribute-list-rubric-completed.csv          # Expert rubric: which exposure factors apply to each stock
│   └── exposure-factor-filter-long.csv              # Long-form rubric (generated by Module 4 Script 3)
│
├── outputs/
│   ├── disbribution-maps/                           # Species distribution PNGs (Module 1)
│   ├── exposure-overlap/                            # Exposure-overlap figures, full layout
│   ├── exposure-overlap-12panel/                    # 12-panel exposure-overlap figures (Module 2)
│   │   └── <Species-Slug>/
│   │       ├── Distribution-Anomalies/              # Reference distribution + anomaly PNGs
│   │       └── Exposure-Overlap-12panel/            # 12-panel overlap PNGs
│   │
│   ├── final-scores-compiled/
│   │   ├── final-attribute-scores/                  # Compiled reviewer scores (Module 4 Script 1)
│   │   │   ├── table_final_attribute_scores_all.csv
│   │   │   └── quantitative-exposure-attribute-scores-all.csv
│   │   ├── overall-vulnerability-rankings/          # Final vulnerability scores (Module 4 Scripts 2–3)
│   │   │   ├── final_scores_uscar.csv               # Combined qual + quant scores, U.S. Caribbean
│   │   │   ├── attribute_means_uscar.csv            # Mean score per stock × attribute
│   │   │   ├── component_scores_uscar.csv           # Sensitivity + exposure component scores
│   │   │   └── overall_vulnerability_scores_uscar.csv
│   │   ├── directional-effect/                      # Directional effect summaries (Module 5)
│   │   │   ├── table_directional_effect_scores.csv
│   │   │   ├── directional_effect_wide_all.csv
│   │   │   └── directional_effect_summary_by-stock.csv
│   │   └── data-quality/                            # Data quality summaries (Module 6)
│   │       ├── table_data_quality_scores_extracted.csv
│   │       └── overall_data_quality_summary_by_stock.csv
│   │
│   ├── final-tallies-long/                          # LMHV tally tables (Module 7 Script 2)
│   │   ├── sensitivity_tallies_long.csv             # Per-reviewer, per-stock sensitivity tallies
│   │   ├── sensitivity_tallies_by_stock.csv
│   │   ├── directional_effect_tallies_long.csv
│   │   ├── directional_effect_tallies_by_stock.csv
│   │   ├── exposure_tallies_long.csv
│   │   └── exposure_tallies_by_stock.csv
│   │
│   ├── analyses/
│   │   ├── 1-inputs/                                # Raw tally extracts (Module 7 Script 1)
│   │   │   └── qa/                                  # QA diagnostic tables
│   │   └── uncertainty-loo/                         # Bootstrap and LOO outputs (Module 8)
│   │       ├── intermediate/                        # Validation and baseline-reproduction checks
│   │       └── final-tables/                        # Analysis-ready tables for figures
│   │
│   └── prework/                                     # Pre-workshop summaries and figures (Module 3)
│
├── figures/                                         # All final publication figures
│   ├── fig_overall_vulnerability.png                # Module 4
│   ├── fig_attribute_score_boxplot_combined.png     # Module 7
│   ├── fig_sensitivity_attribute_score_boxplot.png  # Module 7
│   ├── fig_exposure_attribute_score_boxplot.png     # Module 7
│   ├── fig_directional_effect_summary.png           # Module 7
│   ├── fig_sensitivity_tally_distributions_by_stock.png  # Module 7
│   ├── fig_exposure_tally_distributions_by_stock.png     # Module 7
│   ├── fig_loo_bar_plots.png                        # Module 8
│   └── fig_bootstrap_uncertainty.png               # Module 8
│
└── resources/
    └── HMS/                                         # Reference code from Loughran et al. 2025
```

---

## Data Pipeline

The diagram below shows the key file dependencies across modules. Modules in **bold** are the primary producers of cross-workflow data files.

```
 CMIP6 NetCDF files ─────────────────────────────────────────────┐
 Species shapefiles ─────────────────────────────────────────────┤
                                                                  │
                                                          Module 2 (Exposure-anomalies)
                                                                  │
                        ┌─────────────────────────────────────────┘
                        │
                        ▼
          quantitative-exposure-attribute-scores-all.csv
          exposure-overlap-12panel/ (expert review figures)
                        │
                        └──────────────────────────────┐
                                                        │
 data/final-scores/*.xlsx ──────────────────────────────┤
   (Final reviewer workbooks)                           │
                                                        ▼
                                           Module 4 (Final attribute scoring)
                                                        │
               ┌────────────────────────────────────────┤
               │                                        │
               ▼                                        ▼
  attribute_means_uscar.csv          overall_vulnerability_scores_uscar.csv
  exposure-factor-filter-long.csv    component_scores_uscar.csv
               │                     final_scores_uscar.csv
               │
               │         data/final-scores/*.xlsx ──────┐
               │                                        │
               │                                        ▼
               │                           Module 5 (Directional effect)
               │                                        │
               │                                        ▼
               │                    directional_effect_summary_by-stock.csv
               │
               │         data/final-scores/*.xlsx ──────┐
               │                                        │
               │                                        ▼
               │                           Module 6 (Data quality)
               │                                        │
               │                                        ▼
               │                    overall_data_quality_summary_by_stock.csv
               │
               │         data/final-scores/*.xlsx ──────┐
               │                                        │
               │                                        ▼
               │                      Module 7 (Scoring distributions)
               │                                        │
               │                ┌───────────────────────┘
               │                ▼
               │        outputs/final-tallies-long/
               │          sensitivity_tallies_long.csv
               │          directional_effect_tallies_long.csv
               │          exposure_tallies_long.csv
               │                │
               └────────────────┤
                                ▼
                      Module 8 (Uncertainty analysis)
                                │
                                ▼
                  outputs/analyses/uncertainty-loo/final-tables/
                  figures/fig_bootstrap_uncertainty.png
                  figures/fig_loo_bar_plots.png
```

---

## Workflow Module Details

### Module 0 — Query Species Attributes from FishBase

**Script:** `0-query-species-attributes-from-FishBase/Query-species-attributes-from-FishBase.R`

Queries biological traits and life-history parameters from [FishBase](https://www.fishbase.org) using the [`rfishbase`](https://github.com/ropensci/rfishbase) R package. Reads a `species-list.csv` with scientific names and retrieves species summaries, growth parameters (Von Bertalanffy K and L∞), reproductive mode, trophic level, and depth range. Outputs a single compiled `fishbase_species_attributes.csv`. Results are used as reference material during expert scoring.

---

### Module 1 — Make Species Distribution Maps

**Script:** `1-make-species-distribution-maps/Make-species-distribution-maps.R`

Loops through all species shapefiles in `data/species-distribution-shapefiles/` and produces standardized PNG distribution maps within a Caribbean bounding box (6°N–27.8°N, 92°W–57°W). All 25 maps use consistent symbology. Outputs are saved to `outputs/disbribution-maps/` (one PNG per species).

---

### Module 2 — Exposure Anomalies

**Script:** `2-exposure-anomalies/Exposure-anomalies.R`

Synthesizes CMIP6 multi-model ensemble projections with IUCN species range polygons to produce per-species, per-exposure-factor overlap analyses at three geographic scales (Western Atlantic, Caribbean Sea, U.S. Caribbean). Produces a 12-panel exposure-overlap figure for every stock × exposure factor combination (25 species × 13 factors = 325 figures) for expert review.

Also calculates and exports quantitative exposure scores as a weighted average (see Methods below), which feed into Module 4 as calculated exposure factor scores.

**Key output:** `outputs/final-scores-compiled/quantitative-exposure-attribute-scores-all.csv`

---

### Module 3 — Preliminary Sensitivity Attribute Scoring

**Scripts:** `3-prelim-sensitivity-attribute-scoring/1-extract-scores.R`, `2-prework-scoring-summaries.R`

Reads preliminary reviewer workbooks from `data/preliminary-scores/` and compiles LMHV tally scores into a master table. Calculates HMS-style weighted means and standard deviations per stock × attribute, and generates per-species stacked-bar panel plots (multi-page PDF) for use in reviewer orientation and pre-workshop preparation.

| File | Description |
|------|-------------|
| `data/preliminary-scores/score_table_all.csv` | Compiled long-format tally table from all preliminary reviewers |
| `outputs/prework/sp-x-att_score_summaries.csv` | Mean and SD per stock × attribute |
| `outputs/prework/prework_all_species.pdf` | Multi-page summary PDF for workshop preparation |

---

### Module 4 — Final Attribute and Exposure Scoring

**Scripts:** `4-final-attribute-exposure-scoring/1-extract-final-scores-from-all-reviewers.R`, `2-summarize-attribute-scores.R`, `3-calculate-overall-vulnerability-scores.R`, `4-plot-overall-vulnerability.R`

The core vulnerability scoring workflow. Extracts final reviewer-entered scores from completed workbooks, combines them with quantitative exposure scores from Module 2, and applies the NOAA FCVA logic model to assign overall sensitivity, exposure, and vulnerability ranks.

#### Script 1 — Extract final scores

Loops over all reviewer workbooks in `data/final-scores/`. Extracts final attribute scores (column K) from three workbook sections per stock sheet: Qualitative Exposure (rows 17–18), Sensitivity (rows 21–28), and Rigidity (rows 30–35). Rows 30–35 are labeled "Rigidity" in the workbooks but reclassified as Sensitivity in all downstream scripts.

**Output:** `outputs/final-scores-compiled/final-attribute-scores/table_final_attribute_scores_all.csv`

#### Script 2 — Summarize attribute scores

Standardizes and combines qualitative reviewer scores with the quantitative exposure scores from Module 2. Filters the combined table to the U.S. Caribbean region. Output columns follow the convention: `stock_name`, `region`, `attribute_type` (Sensitivity or Exposure), `score_type` (Qualitative or Calculated), `attribute_name`, `scorer`, `score`.

**Outputs:** `outputs/final-scores-compiled/overall-vulnerability-rankings/final_scores_uscar.csv`

#### Script 3 — Calculate overall vulnerability scores

Calculates attribute-level mean scores, applies the FCVA logic model to assign sensitivity and exposure component scores, and multiplies the two component scores to produce a final vulnerability rank. The expert exposure-factor rubric (`data/attribute-list-rubric-completed.csv`) controls which of the 15 exposure factors are included for each stock; all 14 sensitivity attributes are applied universally. The long-form rubric filter is written to `data/exposure-factor-filter-long.csv` for use by Modules 7 and 8.

**Outputs:**

| File | Description |
|------|-------------|
| `attribute_means_uscar.csv` | Mean score per stock × attribute (Sensitivity and Exposure) |
| `component_scores_uscar.csv` | Sensitivity and Exposure component scores and ranks per stock |
| `overall_vulnerability_scores_uscar.csv` | Final vulnerability score and rank per stock |
| `data/exposure-factor-filter-long.csv` | Long-form expert rubric: 375 rows (15 factors × 25 stocks), `include` column (TRUE/FALSE) |

#### Script 4 — Plot overall vulnerability

Produces the combined vulnerability summary figure (`figures/fig_overall_vulnerability.png`) showing exposure vs. sensitivity component scores, directional effect, and data quality rank.

---

### Module 5 — Final Directional Effect Scoring

**Scripts:** `5-final-directional-effect-scoring/1-extract-directional-effect.R`, `2-summarize-directional-effect.R`

Extracts reviewer directional-effect tallies (Positive / Neutral / Negative, rows 38–40, column M) from final workbooks. Four reviewers × 4 tallies per stock = 16 expected total tallies per stock. Calculates a stock-level directional-effect index as:

$$\text{Directional effect} = \frac{(-1 \times \text{Negative}) + (0 \times \text{Neutral}) + (1 \times \text{Positive})}{\text{Negative} + \text{Neutral} + \text{Positive}}$$

Classification thresholds: ≤ −0.333 = Negative, −0.333 to +0.333 = Neutral, ≥ +0.333 = Positive.

**Key output:** `outputs/final-scores-compiled/directional-effect/directional_effect_summary_by-stock.csv`

---

### Module 6 — Final Data Quality Scoring

**Scripts:** `6-final-data-quality-scoring/1-extract-data-quality-scores.R`, `2-summarize-data-quality-scores.R`

Extracts reviewer-assigned data-quality scores (0–3) from the final workbooks for each attribute and stock. Scores reflect the quality of the evidence underlying each attribute rating:

| Score | Label | Description |
|-------|-------|-------------|
| 3 | Adequate Data | Observed, modeled, or empirically measured for the species from a reputable source |
| 2 | Limited Data | Higher uncertainty; may be based on related species, data from outside the study area, or a less reliable source |
| 1 | Expert Judgment | Based on general knowledge of the species or ecosystem |
| 0 | No Data | No information available to support a score |

Overall data quality per stock is ranked by the proportion of scores ≥ 2: High (≥ 80%), Moderate (50–79%), Poor (< 50%).

**Key output:** `outputs/final-scores-compiled/data-quality/overall_data_quality_summary_by_stock.csv`

---

### Module 7 — Scoring Distributions

**Scripts:** `7-scoring-distributions/1-extract-tally-scores.R`, `2-finalize-tally-tables.R`, `3-plot-figures.R`

Extracts and visualizes the full LMHV tally distributions from the final reviewer workbooks. Unlike Module 4 (which uses final reviewer-entered scores), this module reads the tally columns (columns M–P) to capture the full distribution of reviewer votes across the four ordinal bins.

#### Script 1 — Extract tally scores

Reads tally columns (M–P) from three workbook sections per stock sheet: Qualitative Exposure (rows 17–19), Sensitivity (rows 21–28), and Rigidity (rows 30–35, relabelled as Sensitivity). Also reads Directional Effect tallies (rows 38–40, column M). Applies `standardize_attribute_names()` to normalize attribute labels. Writes long-format tally tables and QA diagnostics to `outputs/analyses/1-inputs/`.

#### Script 2 — Finalize tally tables

Recodes stock names to canonical form (see Stock Name Normalization below), binds qualitative and quantitative exposure tables, and writes the finalized long-format and grouped tally tables to `outputs/final-tallies-long/`.

#### Script 3 — Plot figures

Produces attribute score boxplots (using `attribute_means_uscar.csv`) and per-stock LMHV tally distribution figures. Exposure figures are filtered to expert-approved factor × stock pairs using `data/exposure-factor-filter-long.csv`.

**Key outputs (in `figures/`):** `fig_attribute_score_boxplot_combined.png`, `fig_directional_effect_summary.png`, `fig_sensitivity_tally_distributions_by_stock.png`, `fig_exposure_tally_distributions_by_stock.png`

---

### Module 8 — Uncertainty Analysis

**Scripts:** `8-uncertainty-analysis/uncertainty-analyses.R`, `2-make-uncertainty-figures.R`

Quantifies the statistical robustness of the final CVA vulnerability rankings using two complementary analyses following FCVA methods.

#### Bootstrap resampling

For each stock × sensitivity attribute, all reviewer tally votes are pooled (4 reviewers × 5 tallies = 20 votes per attribute). The pool is resampled with replacement 20 times and the mean is calculated; this is repeated for 10,000 iterations. Bootstrapped attribute means are passed through the FCVA logic model to produce a bootstrapped Sensitivity component rank. The Exposure component score is held fixed (derived from quantitative CMIP6 data, not reviewer tallies). A stock is flagged as **borderline** if its dominant (most frequent) rank accounts for fewer than 75% of iterations.

Directional effects are bootstrapped analogously from the directional-effect tally pool (4 reviewers × 4 tallies = 16 votes per stock; coded +1 / 0 / −1).

#### Leave-one-out (LOO) influence analysis

A deterministic analysis. For each stock × attribute (or factor), the attribute is removed, the FCVA logic model is re-run on the remaining attributes with the other component held fixed, and the new vulnerability rank is compared to baseline. The count of stocks that change rank when a given attribute is omitted is the influence measure for that attribute.

#### Baseline validation

As a QA step, `uncertainty-analyses.R` reproduces all 25 baseline vulnerability ranks from the tally inputs and halts if any disagree with `overall_vulnerability_scores_uscar.csv`. The `rank_threshold` setting must match the value used in Module 4 Script 3 (default: `attr_means_current` = 2, corresponding to the prior NOAA FCVA rule).

**Reads from:**
- `outputs/final-tallies-long/sensitivity_tallies_long.csv`
- `outputs/final-tallies-long/directional_effect_tallies_long.csv`
- `outputs/final-scores-compiled/overall-vulnerability-rankings/attribute_means_uscar.csv`
- `outputs/final-scores-compiled/overall-vulnerability-rankings/overall_vulnerability_scores_uscar.csv`
- `data/exposure-factor-filter-long.csv`

**Key outputs:** `outputs/analyses/uncertainty-loo/final-tables/`, `figures/fig_bootstrap_uncertainty.png`, `figures/fig_loo_bar_plots.png`

---

## Key Cross-Workflow Data Files

The table below lists the files consumed by more than one module.

| File | Produced by | Consumed by | Key columns |
|------|-------------|-------------|-------------|
| `data/final-scores/*.xlsx` | Reviewers | Modules 4, 5, 6, 7 | One tab per stock; scores in column K, tallies in columns M–P |
| `quantitative-exposure-attribute-scores-all.csv` | Module 2 | Module 4 Script 2, Module 7 Script 2 | `stock_name`, `attribute_name`, `spatial_extent`, `score` |
| `table_final_attribute_scores_all.csv` | Module 4 Script 1 | Module 4 Script 2 | `stock_name`, `Attribute_name`, `Attribute_type`, `Final_score`, `Scorer` |
| `final_scores_uscar.csv` | Module 4 Script 2 | Module 4 Script 3 | `stock_name`, `attribute_type`, `score_type`, `attribute_name`, `scorer`, `score` |
| `attribute_means_uscar.csv` | Module 4 Script 3 | Module 7 Script 3, Module 8 | `stock_name`, `attribute_type`, `score_type`, `attribute_name`, `attribute_mean` |
| `overall_vulnerability_scores_uscar.csv` | Module 4 Script 3 | Module 4 Script 4, Module 7 Script 3, Module 8 | `stock_name`, `Exp_score`, `Exp_rank`, `Sens_score`, `Sens_rank`, `Vuln_score`, `Vuln_rank` |
| `data/exposure-factor-filter-long.csv` | Module 4 Script 3 | Module 7 Script 3, Module 8 | `stock_name`, `attribute_name`, `include` (TRUE/FALSE) |
| `directional_effect_summary_by-stock.csv` | Module 5 | Module 4 Script 4, Module 7 Script 3 | `stock_name`, `Positive`, `Neutral`, `Negative`, `wt_avg`, `directional_effect` |
| `overall_data_quality_summary_by_stock.csv` | Module 6 | Module 4 Script 4 | `stock_name`, `prop_ge_2`, `data_quality_rank` |
| `sensitivity_tallies_long.csv` | Module 7 Script 2 | Module 8 | `reviewer_id`, `stock_name`, `attribute_name`, `tally_low`, `tally_moderate`, `tally_high`, `tally_very_high` |
| `directional_effect_tallies_long.csv` | Module 7 Script 2 | Module 8 | `reviewer_id`, `stock_name`, `effect_category`, `tally` |

---

## Species and Attribute Reference

### 25 Assessed Stocks

The following 25 stocks were assessed. **Canonical names** are sentence case and are the form used in all compiled output files. See Stock Name Normalization below for how these relate to names in reviewer workbooks and output figure directories.

| # | Canonical name (output files) |
|---|-------------------------------|
| 1 | Atlantic thread herring |
| 2 | Ballyhoo |
| 3 | Blue runner |
| 4 | Dolphinfish |
| 5 | Gray angelfish |
| 6 | Hogfish |
| 7 | King mackerel |
| 8 | Lane snapper |
| 9 | Long-spined sea urchin |
| 10 | Misty grouper |
| 11 | Mutton snapper |
| 12 | Nassau grouper |
| 13 | Queen conch |
| 14 | Queen snapper |
| 15 | Queen triggerfish |
| 16 | Rainbow parrotfish |
| 17 | Red grouper |
| 18 | Red hind |
| 19 | Sea cucumbers |
| 20 | Silk snapper |
| 21 | Spiny lobster |
| 22 | Stoplight parrotfish |
| 23 | White mullet |
| 24 | Yellowfin grouper |
| 25 | Yellowtail snapper |

### Stock Name Normalization

Reviewer workbooks use **Title Case** stock names (e.g., "Blue Runner"). Compiled output files use **sentence case** (e.g., "Blue runner"). Four stocks additionally have substantively different names between workbooks and outputs, requiring explicit manual overrides before the sentence-case conversion is applied:

| Workbook name | Canonical output name | Reason |
|---|---|---|
| Atlantic Herring | Atlantic thread herring | Different common name |
| Diadema | Long-spined sea urchin | Different common name |
| Redhind | Red hind | Missing space |
| Sea Cucumber | Sea cucumbers | Singular vs. plural |

The normalization is a two-step process applied in Modules 4, 7, and 8:
1. Apply the four manual overrides above.
2. Apply `stringr::str_to_sentence()` to convert remaining Title Case names to sentence case.

> **Note — output figure directories:** The `outputs/exposure-overlap-12panel/` subdirectories use Title Case hyphenated slugs (e.g., `King-Mackerel/`, `Red-Hind/`) derived from the original workbook names, not the canonical sentence-case names used in CSV outputs.

### Sensitivity Attributes

Fourteen sensitivity attributes are scored for all 25 stocks. Attributes from workbook rows 30–35 are labeled "Rigidity" in the reviewer workbooks but are reclassified as **Sensitivity** in all downstream scripts.

| # | Attribute name |
|---|---------------|
| 1 | Adult mobility |
| 2 | Complexity in reproductive strategy |
| 3 | Genetic diversity |
| 4 | Habitat specificity |
| 5 | Mobility and dispersal or early life stages |
| 6 | Other stressors |
| 7 | Population growth rate |
| 8 | Predation and competition dynamics |
| 9 | Prey specificity |
| 10 | Spawning characteristics |
| 11 | Species range |
| 12 | Specificity in early life history requirements |
| 13 | Stock Size Status |
| 14 | Tolerance to ocean acidification |

> **Naming note — "Stock Size Status":** This attribute appears in reviewer workbooks as "Stock Size/Status" (with slash). The `standardize_attribute_names()` function in Module 7 Script 1 removes the slash, producing "Stock Size Status" in the tally tables. However, `attribute_means_uscar.csv` (produced by Module 4 Script 3) stores the attribute as "Stock size/status" (sentence case, with slash). The figure scripts in Modules 7 and 8 handle this mismatch with an explicit `dplyr::recode()` call. The canonical form for lookups and figure labels is **"Stock Size Status"** (no slash).

### Exposure Factors

Up to 15 exposure factors are evaluated per stock — 13 quantitative factors derived from CMIP6 projections (Module 2) and 2 qualitative factors scored by expert reviewers. The subset of factors applied to each stock is determined by the expert rubric in `data/attribute-list-rubric-completed.csv`. All 15 factors are listed below.

**Quantitative (13) — CMIP6-derived:**

| Abbreviation | Full name |
|---|---|
| `bs` | Bottom salinity |
| `bt` | Bottom temperature |
| `chl` | Chlorophyll-a concentration |
| `mld` | Mixed layer depth |
| `msstg` | Mean sea surface temperature gradient |
| `o200` | Oxygen at 200m |
| `ph` | Surface pH |
| `pp` | Primary production |
| `precip` | Precipitation |
| `sso` | Sea surface oxygen |
| `sss` | Sea surface salinity |
| `sst` | Sea surface temperature |
| `swsm` | Surface wind speed magnitude |

**Qualitative (2) — expert-scored:**

| Name | Description |
|---|---|
| Sargassum influx | Projected changes in pelagic Sargassum influx to the U.S. Caribbean |
| Thermocline depth | Projected changes in thermocline depth |

---

## Acknowledgements

**References:**

- Morrison, W. E., Nelson, M. W., Howard, J. F., Stecher, H. A., Sheridan, P. F., & Resnick, M. L. (2015). Methodology for assessing the vulnerability of marine fish and shellfish species to a changing climate. *NOAA Technical Memorandum NMFS-OSF-3*. https://doi.org/10.7289/V5TM782J
- Hare, J. A., Morrison, W. E., Nelson, M. W., Stachura, M. M., Teeters, E. J., Griffis, R. B., Alexander, M. A., Scott, J. D., Alade, L., Bell, R. J., et al. (2016). A vulnerability assessment of fish and invertebrates to climate change on the Northeast U.S. continental shelf. *PLOS ONE*, 11(2), e0146756. https://doi.org/10.1371/journal.pone.0146756
- Loughran, C. E., Hazen, E. L., Brodie, S., Jacox, M. G., Whitney, F. A., Payne, M. R., et al. (2025). A climate vulnerability assessment of highly migratory species in the Northwest Atlantic Ocean. *PLOS Climate*, 4(8), e0000530. https://doi.org/10.1371/journal.pclm.0000530
- Craig, J. K., et al. (2025). A climate vulnerability assessment for managed species in the South Atlantic Bight. *PLOS Climate*. https://doi.org/10.1371/journal.pclm.0000543
