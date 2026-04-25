# Caribbean-CVA

## Project Overview

This repository contains the full analytical pipeline for a **Climate Vulnerability Assessment (CVA)** of 25 fish and invertebrate stocks managed in the U.S. Caribbean. The CVA evaluates each stock's overall climate vulnerability by combining two independently scored components — **sensitivity** (how biologically susceptible is the stock to climate change?) and **exposure** (how much is the stock's habitat projected to change?) — following the NOAA Fisheries Climate Vulnerability Assessment (FCVA) framework established by Morrison et al. (2015) and applied in recent regional CVAs including Loughran et al. (2025).

Analyses were conducted by Harris Analytics & Research LLC in support of [Isla Mar 501c3](https://www.islamar.org/).  
All code and materials are available under an open-access license, as per Creative Commons CC0 1.0.  
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

Quantifies the statistical robustness of the final CVA vulnerability rankings using two complementary analyses following Loughran et al. (2025):

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

## Methods

### Goals and Scope

The overall goal of the exposure factor analyses is to compare projected future ocean conditions against their past. To do so, we utilize spatially explicit ocean model projections (both historical and future ocean projections) to create a standardized anomaly map (i.e., a static comparison) and provide accompanying analyses. The values of the standardized anomalies are expressed in standard deviation units, which allows the result to be comparable across variables. For each exposure factor anomaly, we mapped the gridded cells of the standardized exposure anomaly and overlapped a given species distribution. Values from the overlapped cells were then compiled into frequency distributions (histograms) and categorized (bar plots).

The exposure overlap analyses below are presented within three geographic scopes:
1. A given species distribution within the Western Atlantic Ocean (5°S–72°N, 99°W–40°W),
2. The wider Caribbean region (6°N–28°N, 92°W–57°W), and
3. The U.S. Caribbean (16°N–20°N, 69°W–63°W).

The spatial scale of the analysis will have differing implications for a species' ecology and its management. All three spatial scales are presented for use at the discretion of the expert reviewers.

### Calculating Gridded Standardized Exposure Factor Anomalies

Similar to the HMS CVA, we utilized projections from the Coupled Model Intercomparison Project Phase 6 (CMIP6) multi-model ensemble under the SSP5-8.5 scenario, which is consistent with past NOAA CVAs. This scenario represents the highest greenhouse gas emissions pathway and likely represents the greatest changes that we can reasonably expect for the exposure factor; it should be considered an upper bound for CMIP6 projections.

Monthly outputs for each exposure factor (e.g., sea surface temperature) were obtained at a grid cell resolution of 1° latitude × 1° longitude. To calculate detrended standardized anomalies, a unitless value (z) was computed for each 1° × 1° grid cell as the difference between the future and historical values divided by the historical standard deviation:

z = (μ future − μ historical) / σ historical,

where μ future is the mean of a given grid cell for months during 2020–2049, μ historical is the mean of a given grid cell for months during 1985–2014, and σ is the interannual standard deviation during the same historical baseline period (1985–2014). Dividing by the historical baseline standard deviation converts z into standard deviation units (σ).

A positive standardized anomaly value (+σ) means the exposure factor is projected to increase compared to the historical baseline, while a negative value (–σ) indicates a decrease.

### Calculating Quantitative Exposure Scores

For a given species or stock and exposure factor, standardized anomalies from the gridded overlapped area were grouped into 0.25 standard deviation bins. These bins were then coded by their absolute value in the following categories (LMHV): Low (|σ| < 0.5), Moderate (0.5 ≤ |σ| < 1.5), High (1.5 ≤ |σ| < 2.0), Very High (|σ| ≥ 2.0). A weighted average score was calculated as:

Weighted Average = (1L + 2M + 3H + 4V) / (L + M + H + V)

where L, M, H, and V are the total counts of cells in each category. This score ranges from 1 (all Low) to 4 (all Very High) and is used as the calculated exposure score for each stock × factor in Modules 4 and 7.

### NOAA FCVA Logic Model

Overall climate vulnerability was calculated using the NOAA FCVA framework. For each stock, attribute-level means are calculated (as the arithmetic mean of non-missing reviewer scores for qualitative attributes, or the calculated weighted-average score for quantitative exposure factors). The logic model then assigns component scores based on how many attribute or factor means exceed specified thresholds:

| Overall component rank | Numeric score | Logic rule |
|---|---:|---|
| Very High | 4 | More than 3 attribute or factor means ≥ 3.5 |
| High | 3 | More than 2 attribute or factor means ≥ 3.0 |
| Moderate | 2 | More than 2 attribute or factor means ≥ 2.5 |
| Low | 1 | All other cases |

The final vulnerability score is the product of Sensitivity score × Exposure score, assigned to a rank as follows:

| Exposure × Sensitivity score | Overall vulnerability rank |
|---:|---|
| 1–3 | Low |
| 4–6 | Moderate |
| 8–9 | High |
| 12–16 | Very High |

---

## Exposure Overlap Analysis

<p align="center">
  <img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/outputs/exposure-overlap-12panel/King-Mackerel/Exposure-Overlap-12panel/King-Mackerel_Exposure-Overlap_o200.png?raw=true" 
       alt="King Mackerel – Exposure Overlap (o200)" 
       width="700"/>
</p>

#### Figure 1: Example of exposure overlap analysis for King Mackerel (*Scomberomorus cavalla*) and dissolved oxygen at 200 m (o200).

The exposure overlap analysis figure is a 3 × 4 panel grid organized in two dimensions. The three columns represent the three nested geographic extents.
- **Left column (A, D, G)** represents the species' range within the Western Atlantic for stock-wide context.
- **Middle column (B, E, H)** crops the map for the wider Caribbean and panels in this column only show values for grid cells within this region.
- **Right column (C, F, I)** further crops the map and only includes grid cells in the U.S. federally-managed waters offshore Puerto Rico and the USVI.

The four rows represent overlap analyses for a given standardized exposure factor anomaly (σ).
- **Row 1 (A–C): Spatial overlap maps.** Gridded maps of the standardized CMIP6 exposure factor anomaly grid cells (1°×1° cells) that overlap with a given species' distribution. Standardized anomalies are expressed in σ-units (change relative to the detrended interannual variability of the 1985–2014 baseline).
- **Row 2 (D–F): Categorical summations.** Total count of standardized exposure factor anomaly values within five signed categories (< −1.5σ, −1.5 to −0.5σ, −0.5 to +0.5σ, +0.5 to +1.5σ, > +1.5σ), with proportions indicated.
- **Row 3 (G–I): Frequency distributions.** Frequency histograms (y-axis = percentage) in 0.25σ bins, color-coded at four LMHV exposure levels based on absolute value.
- **Row 4 (J–L): LMHV summary and weighted averages.** Proportions of Low, Moderate, High, and Very High anomalies (absolute values) with weighted average score in the top-left corner.

Exposure overlap figures for all species are available here: https://github.com/holden-harris/Caribbean-CVA/tree/main/outputs/exposure-overlap-12panel

---

## Oceanographic Exposure Variables

For each species under CVA review, exposure analyses were conducted for the following 13 oceanographic exposure variables.

| Abbreviation | Full Variable Name | Description |
|---|---|---|
| **bs** | Bottom Salinity | Mean salinity near the seafloor. Important for benthic organisms sensitive to freshwater inputs or stratification. |
| **bt** | Bottom Temperature | Mean temperature near the seafloor, influencing demersal fish and benthic invertebrates. |
| **chl** | Chlorophyll-a Concentration | A proxy for phytoplankton biomass, indicating primary productivity and food availability at the base of the food web. |
| **mld** | Mixed Layer Depth | Depth of the upper, well-mixed surface ocean layer, affecting nutrient availability, light, and stratification. |
| **msstg** | Mean Sea Surface Temperature Gradient | Spatial temperature gradient at the surface; an indicator of thermal fronts, ocean circulation, and habitat boundaries. |
| **o200** | Oxygen at 200 m | Dissolved oxygen concentration at ~200 meters depth; reflects mid-water oxygen availability and deoxygenation trends. |
| **ph** | Surface pH | Measure of ocean acidity (linked to CO₂ uptake and ocean acidification). Lower values = more acidic. |
| **pp** | Primary Production | Gross primary productivity of phytoplankton; determines energy input to marine food webs. |
| **precip** | Precipitation | Rainfall over the ocean, relevant for freshwater input, stratification, and coastal salinity changes. |
| **sso** | Sea Surface Oxygen | Dissolved oxygen concentration at the surface, important for respiration of pelagic species. |
| **sss** | Sea Surface Salinity | Surface salt concentration, reflecting freshwater input, evaporation, and circulation. |
| **sst** | Sea Surface Temperature | Temperature of the upper ocean, widely used as a climate indicator and driver of species distributions. |
| **swsm** | Surface Wind Speed Magnitude | Intensity of winds at the ocean surface, a driver of mixing, upwelling, and surface currents. |

---

## Acknowledgements

Analyses were conducted by Harris Analytics & Research LLC in support of [Isla Mar 501c3](https://www.islamar.org/). These were built on the efforts from past CVAs.  
All code and materials are available under an open-access license, as per Creative Commons CC0 1.0.  
We thank Dan Crear (ICATTC) and Tyler Loughran (NOAA) for their assistance in this work.

**References:**

- Morrison, W. E., Nelson, M. W., Howard, J. F., Stecher, H. A., Sheridan, P. F., & Resnick, M. L. (2015). Methodology for assessing the vulnerability of marine fish and shellfish species to a changing climate. *NOAA Technical Memorandum NMFS-OSF-3*. https://doi.org/10.7289/V5TM782J
- Hare, J. A., Morrison, W. E., Nelson, M. W., Stachura, M. M., Teeters, E. J., Griffis, R. B., Alexander, M. A., Scott, J. D., Alade, L., Bell, R. J., et al. (2016). A vulnerability assessment of fish and invertebrates to climate change on the Northeast U.S. continental shelf. *PLOS ONE*, 11(2), e0146756. https://doi.org/10.1371/journal.pone.0146756
- Loughran, C. E., Hazen, E. L., Brodie, S., Jacox, M. G., Whitney, F. A., Payne, M. R., et al. (2025). A climate vulnerability assessment of highly migratory species in the Northwest Atlantic Ocean. *PLOS Climate*, 4(8), e0000530. https://doi.org/10.1371/journal.pclm.0000530
- Craig, J. K., et al. (2025). A climate vulnerability assessment for managed species in the South Atlantic Bight. *PLOS Climate*. https://doi.org/10.1371/journal.pclm.0000543
