# Caribbean-CVA

This repository contains the complete analytical pipeline for a NOAA Fisheries Climate Vulnerability Assessment (FCVA) of 25 fish and invertebrate stocks managed in the U.S. Caribbean. This workflow includes processing oceanographic projections and overlapped exposures (Modules 1–2); extracting data from workbooks filled by the expert reviewers (Modules 3 and 4.1); calculating results for overall vulnerability (Module 4.3), directional effect (Module 5), data quality (Module 6), and potential for distribution change (Module 9); evaluating uncertainty with bootstrap resampling (Modules 7-8) and leave-one-out analysis (Module 8); and, creating final figures and tables (Module 10).

> All data syntheses and analyses are in R. Each run writes figures to `figures/{run_label}/` and analysis outputs to `outputs/{run_label}/`; the active run is selected by `active_run` in `config.R` (see **Shared Configuration** below). To run the full pipeline, open `Caribbean-CVA.Rproj` in RStudio, set `active_run` in `config.R`, and source `run-all.R` from the project root.

All code and materials were developed by **Harris Analytics & Research LLC** in support of **[Isla Mar 501c3](https://www.islamar.org/)**. All code and data are available under an open-access Creative Commons CC0 1.0 license. Please cite if you use it. For more information, please contact holden@harris-analytics.com. 

## Workflow Modules

The project is organized into eleven numbered workflow modules. Each module has its own subdirectory with scripts and a `ReadMe` file. The pipeline is split into multiple phases: **Exposure Overlap (Modules 1-2)**; **Workbook Data Extraction (Modules 3-4)**; **Analyses (Modules 4–9)** produce all data outputs; **Figures (Module 10)** reads those outputs and produces all publication figures. Use `run-all.R` at the project root to execute final workbook extraction, final analyses, uncertainty analyses, and figure generation in order. 

| Module | Folder | Purpose | Key outputs |
|--------|--------|---------|-------------|
| 0 | `00-query-species-attributes-from-FishBase/` | Query biological traits and life-history attributes from FishBase via the `rfishbase` R package | `fishbase_species_attributes.csv` |
| 1 | `01-make-species-distribution-maps/` | Generate standardized PNG distribution maps for all 25 species from IUCN shapefiles | `outputs/disbribution-maps/*.png` |
| 2 | `02-exposure-anomalies/` | Calculate CMIP6-based standardized anomaly maps; produce 12-panel exposure-overlap figures reviewed by CVA experts; extract quantitative exposure scores | `outputs/exposure-overlap-12panel/`, `outputs/final-scores-compiled/quantitative-exposure-attribute-scores-all.csv` |
| 3 | `03-prelim-sensitivity-attribute-scoring/` | Extract and summarize preliminary sensitivity-attribute tallies from pre-workshop reviewer workbooks; generate per-stock LMHV summary plots for workshop preparation | `data/preliminary-scores/score_table_all.csv`, `outputs/prework/` |
| 4 | `04-final-attribute-exposure-scoring/` | Extract final reviewer scores, calculate attribute means, apply NOAA FCVA logic model to produce stock-level sensitivity, exposure, and overall vulnerability scores | `outputs/final-scores-compiled/overall-vulnerability-rankings/` |
| 5 | `05-final-directional-effect-scoring/` | Extract directional-effect tallies (Positive / Neutral / Negative) from final workbooks; calculate stock-level directional-effect index | `outputs/final-scores-compiled/directional-effect/directional_effect_summary_by-stock.csv` |
| 6 | `06-final-data-quality-scoring/` | Extract reviewer data-quality scores (0–3) for each attribute; summarize and rank overall data quality per stock | `outputs/final-scores-compiled/data-quality/overall_data_quality_summary_by_stock.csv` |
| 7 | `07-scoring-distributions/` | Extract LMHV tally distributions from final workbooks; finalize long-format tally tables for uncertainty analysis and figures | `outputs/final-tallies-long/` |
| 8 | `08-uncertainty-analysis/` | Bootstrap resampling and leave-one-out influence analyses; quantify statistical robustness of final vulnerability rankings | `outputs/{run_label}/analyses/uncertainty-loo/` |
| 9 | `09-distributional-change-potential/` | Calculate each stock's potential for distributional shift using four sensitivity attributes; bootstrap uncertainty via draw-pile resampling | `outputs/{run_label}/distribution-change-potential/` |
| 10 | `10-figures/` | Produce all publication figures from analysis outputs (Modules 4–9); no new data are generated here | `figures/{run_label}/fig_*.png` (13 publication figures per run) |

---

## Shared Configuration (`config.R`)

`config.R` at the project root is sourced by every analysis and figure script immediately after `rm(list = ls())`. It is the single authoritative location for all constants shared across the pipeline:

| Constant | Value / Description | Used in |
|---|---|---|
| `active_run` | `"broadened_distribution"` or `"cross_region_comparable"` — selects the active run configuration | All modules (sourced at top of every script) |
| `rank_threshold` | Derived from `active_run`; controls FCVA logic model cutoffs | FCVA logic model (Modules 4, 8, 9) |
| `sens_attrs_drop` | Derived from `active_run`; sensitivity attribute names excluded from this run | Module 4 Script 3, Module 8 |
| `run_label` | Derived from `active_run`; used as output subdirectory name (`"broadened-distribution"` or `"cross-region-comparable"`) | All modules that write outputs |
| `dir_eff_threshold` | `1/3` | Directional effect classification (Module 5, 8) |
| `borderline_prop` | `0.75` | Bootstrap borderline flag (Modules 8, 9) |
| `cert_very_high / cert_high / cert_moderate` | `0.95 / 0.90 / 0.67` | Bootstrap certainty bins (Module 10) |
| `rank_levels` | `c("Low","Moderate","High","Very High")` | Ordered factor levels across all scripts |
| `rank_colors` | green3 / yellow2 / orange2 / red3 | Standard rank palette (Module 10 Scripts 2 & 3) |
| `dir_levels` / `dir_colors` | Negative / Neutral / Positive | Directional effect palette (Module 10) |
| `stock_name_recode` | 27-entry lookup | Display-name normalization in all figure scripts |
| `attr_short_names` | 14 sensitivity attribute short labels | Axis labels in Module 10 Scripts 1 & 2 |
| `exp_attr_short_names` | 15 exposure factor short labels | Axis labels in Module 10 Scripts 1 & 2 |

Scripts with intentionally different palettes define local overrides that shadow the config values (Module 10 Script 1 uses lighter tally-bar colors; Module 10 Script 4 uses hex `vuln_colors` for the tile grid).

To change the logic model or threshold, edit the relevant constant in `config.R` and re-run `run-all.R`. No individual script needs to be touched.

### Adjustable logic model for overall vulnerability

The `rank_threshold` constant controls how many attributes must exceed each score cutoff before a Sensitivity or Exposure component is assigned a given rank. It is applied identically in Modules 04, 08, and 09. The logic model evaluates each condition in order; the first condition met wins:

| Condition | Score cutoff | Attributes required (at `rank_threshold = 1`) | Component rank |
|---|---|---|---|
| `n_attrs with mean ≥ 3.5  >  rank_threshold + 1` | ≥ 3.5 | ≥ 3 | Very High |
| `n_attrs with mean ≥ 3.0  >  rank_threshold`     | ≥ 3.0 | ≥ 2 | High |
| `n_attrs with mean ≥ 2.5  >  rank_threshold`     | ≥ 2.5 | ≥ 2 | Moderate |
| Otherwise | — | 0 or 1 attribute meets any threshold | Low |

The default value of `1` matches the standard NOAA FCVA methodology. Changing it shifts every rank boundary uniformly:

| `rank_threshold` | Named constant | Model character | Very High requires | High / Moderate require |
|---|---|---|---|---|
| `0L` | — | More permissive | ≥ 2 attributes with mean ≥ 3.5 | ≥ 1 attribute at respective cutoff |
| `1L` | `attr_means_current` | Standard NOAA FCVA (`cross_region_comparable` run) | ≥ 3 attributes with mean ≥ 3.5 | ≥ 2 attributes at respective cutoff |
| `2L` | `attr_means_plus1` | Revised / broader distribution (`broadened_distribution` run) | ≥ 4 attributes with mean ≥ 3.5 | ≥ 3 attributes at respective cutoff |

To switch runs, set `active_run` in `config.R` to one of the two named configurations and re-run `run-all.R`:

| `active_run` | FCVA logic model | Sensitivity attributes | Purpose |
|---|---|---|---|
| `broadened_distribution` | `attr_means_plus1` (rank_threshold = 2L) | All 14 | Broader relative distribution; conservation priority-setting for CFMC |
| `cross_region_comparable` | `attr_means_current` (rank_threshold = 1L) | 12 (drops Genetic diversity and Predation and competition dynamics) | Directly comparable to other NOAA FCVAs |

Outputs are written to `outputs/{run_label}/` and figures to `figures/{run_label}/` so both runs coexist on disk.

---

## Directory Structure

```
Caribbean-CVA/
│
├── config.R                                             # Shared constants; defines two named run configs
│                                                        #   (active_run switch, rank_threshold, sens_attrs_drop,
│                                                        #   run_label, colors, stock name recode, display labels)
├── run-all.R                                            # Master orchestration: Phase 1 analyses → Phase 2 figures
│
├── 00-query-species-attributes-from-FishBase/
│   ├── Query-species-attributes-from-FishBase.R
│   └── ReadMe.md
│
├── 01-make-species-distribution-maps/
│   ├── Make-species-distribution-maps.R
│   └── ReadMe.md
│
├── 02-exposure-anomalies/
│   ├── Exposure-anomalies.R
│   └── ReadMe.md
│
├── 03-prelim-sensitivity-attribute-scoring/
│   ├── 1-extract-scores.R
│   ├── 2-prework-scoring-summaries.R
│   └── ReadMe.md
│
├── 04-final-attribute-exposure-scoring/
│   ├── 1-extract-final-scores-from-all-reviewers.R
│   ├── 2-summarize-attribute-scores.R
│   ├── 3-calculate-overall-vulnerability-scores.R
│   └── ReadMe.MD
│
├── 05-final-directional-effect-scoring/
│   ├── 1-extract-directional-effect.R
│   ├── 2-summarize-directional-effect.R
│   └── ReadMe.MD
│
├── 06-final-data-quality-scoring/
│   ├── 1-extract-data-quality-scores.R
│   ├── 2-summarize-data-quality-scores.R
│   └── ReadMe.md
│
├── 07-scoring-distributions/
│   ├── 1-extract-tally-scores.R
│   ├── 2-finalize-tally-tables.R
│   └── ReadMe.md
│
├── 08-uncertainty-analysis/
│   ├── uncertainty-analyses.R
│   └── ReadMe.MD
│
├── 09-distributional-change-potential/
│   ├── 1-calculate-distributional-change-potential.R
│   ├── 2-bootstrap-distributional-change.R
│   └── ReadMe.md
│
├── 10-figures/
│   ├── 1-plot-scoring-distributions.R                  # Score boxplots, tally distributions, directional effect
│   ├── 2-plot-uncertainty-figures.R                    # LOO bar charts, bootstrap uncertainty panels
│   ├── 3-plot-distributional-change.R                  # DCP rank column chart, DCP vs. vulnerability cross-plot
│   └── 4-plot-overall-vulnerability.R                  # Overall vulnerability grid + directional effect panel
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
│   │   └── 1-inputs/                                # Raw tally extracts (Module 7 Script 1)
│   │       └── qa/                                  # QA diagnostic tables
│   │
│   ├── broadened-distribution/                      # Run outputs — attr_means_plus1, all 14 sensitivity attributes
│   │   ├── final-scores-compiled/
│   │   │   └── overall-vulnerability-rankings/      # attribute_means_uscar.csv, component_scores_uscar.csv,
│   │   │                                            #   overall_vulnerability_scores_uscar.csv, exposure_factor_qa.csv
│   │   ├── analyses/
│   │   │   └── uncertainty-loo/                     # Bootstrap and LOO outputs (Module 8)
│   │   │       ├── intermediate/                    # Validation and baseline-reproduction checks
│   │   │       └── final-tables/                    # Analysis-ready tables for figures
│   │   ├── distribution-change-potential/           # DCP scores and bootstrap (Module 9)
│   │   └── tables/                                  # Publication results tables (Module 10 Script 5)
│   │
│   ├── cross-region-comparable/                     # Run outputs — attr_means_current, 12 sensitivity attributes
│   │   └── [same structure as broadened-distribution/]
│   │
│   └── prework/                                     # Pre-workshop summaries and figures (Module 3)
│
├── figures/                                         # Publication figures — one subfolder per run (Module 10)
│   ├── broadened-distribution/                      # Figures for broadened_distribution run
│   │   ├── fig_overall_vulnerability.png            # Module 10 Script 4
│   │   ├── fig_attribute_score_boxplot_combined_*.png  # Module 10 Script 1 (horizontal + vertical)
│   │   ├── fig_sensitivity_attribute_score_boxplot.png # Module 10 Script 1
│   │   ├── fig_exposure_attribute_score_boxplot.png    # Module 10 Script 1
│   │   ├── fig_directional_effect_summary.png          # Module 10 Script 1
│   │   ├── fig_sensitivity_tally_distributions_by_stock.png  # Module 10 Script 1
│   │   ├── fig_exposure_tally_distributions_by_stock.png     # Module 10 Script 1
│   │   ├── fig_reviewer_stock_coverage.png              # Module 10 Script 1 (QA)
│   │   ├── fig_loo_bar_plots.png                        # Module 10 Script 2
│   │   ├── fig_bootstrap_uncertainty.png                # Module 10 Script 2
│   │   ├── fig_distributional_change_ranks.png          # Module 10 Script 3
│   │   └── fig_distributional_change_vs_vulnerability.png  # Module 10 Script 3
│   └── cross-region-comparable/                     # Figures for cross_region_comparable run
│       └── [same 13 figures as broadened-distribution/]
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
                  outputs/{run_label}/analyses/uncertainty-loo/final-tables/

 outputs/{run_label}/attribute_means_uscar.csv ─────┐
 sensitivity_tallies_long.csv ──────────────────────┤
 outputs/{run_label}/overall_vulnerability_scores_uscar.csv ─┤
                                                     ▼
                                         Module 9 (Distributional change potential)
                                                     │
                                                     ▼
                         outputs/{run_label}/distribution-change-potential/
                           distributional_change_potential_uscar.csv
                           distributional_change_bootstrap_uscar.csv

 ════════════════════════════════════════════════════════════════════
  PHASE 2 — Figures (Module 10)
  Reads run-specific analysis outputs; writes figures/{run_label}/
 ════════════════════════════════════════════════════════════════════

 outputs/final-scores-compiled/ ────────────────────┐  (shared, run-independent)
 outputs/final-tallies-long/ ───────────────────────┤
 outputs/{run_label}/analyses/uncertainty-loo/ ─────┤
 outputs/{run_label}/distribution-change-potential/ ┤
                                                     ▼
                                         Module 10 (Figures)
                                                     │
                    ┌────────────────────────────────┤
                    │                                │
                    ▼                                ▼
            figures/{run_label}/fig_overall_vulnerability.png    figures/{run_label}/fig_distributional_change_*.png
            figures/{run_label}/fig_*_score_boxplot_*.png        figures/{run_label}/fig_bootstrap_uncertainty.png
            figures/{run_label}/fig_directional_effect_*.png     figures/{run_label}/fig_loo_bar_plots.png
            figures/{run_label}/fig_*_tally_distributions_*.png
```

---

## Workflow Module Details

### Module 0 — Query Species Attributes from FishBase

**Script:** `00-query-species-attributes-from-FishBase/Query-species-attributes-from-FishBase.R`

Queries biological traits and life-history parameters from [FishBase](https://www.fishbase.org) using the [`rfishbase`](https://github.com/ropensci/rfishbase) R package. Reads a `species-list.csv` with scientific names and retrieves species summaries, growth parameters (Von Bertalanffy K and L∞), reproductive mode, trophic level, and depth range. Outputs a single compiled `fishbase_species_attributes.csv`. Results are used as reference material during expert scoring.

---

### Module 1 — Make Species Distribution Maps

**Script:** `01-make-species-distribution-maps/Make-species-distribution-maps.R`

Loops through all species shapefiles in `data/species-distribution-shapefiles/` and produces standardized PNG distribution maps within a Caribbean bounding box (6°N–27.8°N, 92°W–57°W). All 25 maps use consistent symbology. Outputs are saved to `outputs/disbribution-maps/` (one PNG per species).

---

### Module 2 — Exposure Anomalies

**Script:** `02-exposure-anomalies/Exposure-anomalies.R`

Synthesizes CMIP6 multi-model ensemble projections with IUCN species range polygons to produce per-species, per-exposure-factor overlap analyses at three geographic scales (Western Atlantic, Caribbean Sea, U.S. Caribbean). Produces a 12-panel exposure-overlap figure for every stock × exposure factor combination (25 species × 13 factors = 325 figures) for expert review.

Also calculates and exports quantitative exposure scores as a weighted average (see Methods below), which feed into Module 4 as calculated exposure factor scores.

**Key output:** `outputs/final-scores-compiled/quantitative-exposure-attribute-scores-all.csv`

---

### Module 3 — Preliminary Sensitivity Attribute Scoring

**Scripts:** `03-prelim-sensitivity-attribute-scoring/1-extract-scores.R`, `2-prework-scoring-summaries.R`

Reads preliminary reviewer workbooks from `data/preliminary-scores/` and compiles LMHV tally scores into a master table. Calculates HMS-style weighted means and standard deviations per stock × attribute, and generates per-species stacked-bar panel plots (multi-page PDF) for use in reviewer orientation and pre-workshop preparation.

| File | Description |
|------|-------------|
| `data/preliminary-scores/score_table_all.csv` | Compiled long-format tally table from all preliminary reviewers |
| `outputs/prework/sp-x-att_score_summaries.csv` | Mean and SD per stock × attribute |
| `outputs/prework/prework_all_species.pdf` | Multi-page summary PDF for workshop preparation |

---

### Module 4 — Final Attribute and Exposure Scoring

**Scripts:** `04-final-attribute-exposure-scoring/1-extract-final-scores-from-all-reviewers.R`, `2-summarize-attribute-scores.R`, `3-calculate-overall-vulnerability-scores.R`

The core vulnerability scoring workflow. Extracts final reviewer-entered scores from completed workbooks, combines them with quantitative exposure scores from Module 2, and applies the NOAA FCVA logic model to assign overall sensitivity, exposure, and vulnerability ranks.

#### Script 1 — Extract final scores

Loops over all reviewer workbooks in `data/final-scores/`. Extracts final attribute scores (column K) from three workbook sections per stock sheet: Qualitative Exposure (rows 17–18), Sensitivity (rows 21–28), and Rigidity (rows 30–35). Rows 30–35 are labeled "Rigidity" in the workbooks but reclassified as Sensitivity in all downstream scripts.

**Output:** `outputs/final-scores-compiled/final-attribute-scores/table_final_attribute_scores_all.csv`

#### Script 2 — Summarize attribute scores

Standardizes and combines qualitative reviewer scores with the quantitative exposure scores from Module 2. Filters the combined table to the U.S. Caribbean region. Output columns follow the convention: `stock_name`, `region`, `attribute_type` (Sensitivity or Exposure), `score_type` (Qualitative or Calculated), `attribute_name`, `scorer`, `score`.

**Outputs:** `outputs/final-scores-compiled/overall-vulnerability-rankings/final_scores_uscar.csv`

#### Script 3 — Calculate overall vulnerability scores

Calculates attribute-level mean scores, applies the FCVA logic model to assign sensitivity and exposure component scores, and multiplies the two component scores to produce a final vulnerability rank. The expert exposure-factor rubric (`data/attribute-list-rubric-completed.csv`) controls which of the 15 exposure factors are included for each stock. The sensitivity attributes applied are controlled by `sens_attrs_drop` in `config.R`: all 14 for `broadened_distribution`; 12 for `cross_region_comparable` (drops Genetic diversity and Predation and competition dynamics). The long-form rubric filter is written to `data/exposure-factor-filter-long.csv` for use by Module 8.

**Outputs** (written to `outputs/{run_label}/final-scores-compiled/overall-vulnerability-rankings/`):

| File | Description |
|------|-------------|
| `attribute_means_uscar.csv` | Mean score per stock × attribute (Sensitivity and Exposure) |
| `component_scores_uscar.csv` | Sensitivity and Exposure component scores and ranks per stock |
| `overall_vulnerability_scores_uscar.csv` | Final vulnerability score and rank per stock |
| `data/exposure-factor-filter-long.csv` | Long-form expert rubric: 375 rows (15 factors × 25 stocks), `include` column (TRUE/FALSE) — written to `data/`, run-independent |

---

### Module 5 — Final Directional Effect Scoring

**Scripts:** `05-final-directional-effect-scoring/1-extract-directional-effect.R`, `2-summarize-directional-effect.R`

Extracts reviewer directional-effect tallies (Positive / Neutral / Negative, rows 38–40, column M) from final workbooks. Four reviewers × 4 tallies per stock = 16 expected total tallies per stock. Calculates a stock-level directional-effect index as:

$$\text{Directional effect} = \frac{(-1 \times \text{Negative}) + (0 \times \text{Neutral}) + (1 \times \text{Positive})}{\text{Negative} + \text{Neutral} + \text{Positive}}$$

Classification thresholds: ≤ −0.333 = Negative, −0.333 to +0.333 = Neutral, ≥ +0.333 = Positive.

**Key output:** `outputs/final-scores-compiled/directional-effect/directional_effect_summary_by-stock.csv`

---

### Module 6 — Final Data Quality Scoring

**Scripts:** `06-final-data-quality-scoring/1-extract-data-quality-scores.R`, `2-summarize-data-quality-scores.R`

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

**Scripts:** `07-scoring-distributions/1-extract-tally-scores.R`, `2-finalize-tally-tables.R`

Extracts the full LMHV tally distributions from the final reviewer workbooks. Unlike Module 4 (which uses final reviewer-entered scores), this module reads the tally columns (columns M–P) to capture the full distribution of reviewer votes across the four ordinal bins. Figures from these tables are produced by Module 10 Script 1.

#### Script 1 — Extract tally scores

Reads tally columns (M–P) from three workbook sections per stock sheet: Qualitative Exposure (rows 17–19), Sensitivity (rows 21–28), and Rigidity (rows 30–35, relabelled as Sensitivity). Also reads Directional Effect tallies (rows 38–40, column M). Applies `standardize_attribute_names()` to normalize attribute labels. Writes long-format tally tables and QA diagnostics to `outputs/analyses/1-inputs/`.

#### Script 2 — Finalize tally tables

Recodes stock names to canonical form (see Stock Name Normalization below), binds qualitative and quantitative exposure tables, and writes the finalized long-format and grouped tally tables to `outputs/final-tallies-long/`.

**Key outputs (in `outputs/final-tallies-long/`):** `sensitivity_tallies_long.csv`, `directional_effect_tallies_long.csv`, `exposure_tallies_long.csv`, and corresponding grouped-by-stock versions.

---

### Module 8 — Uncertainty Analysis

**Script:** `08-uncertainty-analysis/uncertainty-analyses.R`

Quantifies the statistical robustness of the final CVA vulnerability rankings using two complementary analyses following FCVA methods. Figures from these outputs are produced by Module 10 Script 2.

#### Bootstrap resampling

For each stock × sensitivity attribute, all reviewer tally votes are pooled (4 reviewers × 5 tallies = 20 votes per attribute). The pool is resampled with replacement 20 times and the mean is calculated; this is repeated for 10,000 iterations. Bootstrapped attribute means are passed through the FCVA logic model to produce a bootstrapped Sensitivity component rank. The Exposure component score is held fixed (derived from quantitative CMIP6 data, not reviewer tallies). A stock is flagged as **borderline** if its dominant (most frequent) rank accounts for fewer than 75% of iterations.

Directional effects are bootstrapped analogously from the directional-effect tally pool (4 reviewers × 4 tallies = 16 votes per stock; coded +1 / 0 / −1).

#### Leave-one-out (LOO) influence analysis

A deterministic analysis. For each stock × attribute (or factor), the attribute is removed, the FCVA logic model is re-run on the remaining attributes with the other component held fixed, and the new vulnerability rank is compared to baseline. The count of stocks that change rank when a given attribute is omitted is the influence measure for that attribute.

#### Baseline validation

As a QA step, `uncertainty-analyses.R` reproduces all 25 baseline vulnerability ranks from the tally inputs and halts if any disagree with `overall_vulnerability_scores_uscar.csv`. The `rank_threshold` and `sens_attrs_drop` are sourced from `config.R` via `active_run` and must match the settings used in Module 4 Script 3.

**Reads from:**
- `outputs/final-tallies-long/sensitivity_tallies_long.csv`
- `outputs/final-tallies-long/directional_effect_tallies_long.csv`
- `outputs/{run_label}/final-scores-compiled/overall-vulnerability-rankings/attribute_means_uscar.csv`
- `outputs/{run_label}/final-scores-compiled/overall-vulnerability-rankings/overall_vulnerability_scores_uscar.csv`
- `data/exposure-factor-filter-long.csv`

**Key outputs:** `outputs/{run_label}/analyses/uncertainty-loo/final-tables/`

---

### Module 9 — Potential for Distributional Change

**Scripts:** `09-distributional-change-potential/1-calculate-distributional-change-potential.R`, `2-bootstrap-distributional-change.R`

Calculates each stock's potential for distributional shift under changing environmental conditions, following the methodology of prior NOAA CVAs (HMS: Loughran et al. 2025; South Atlantic: Craig et al. 2025; GoM: Quinlan et al. 2023). Four sensitivity attributes are used; three movement-related attributes are inverted (`5 − mean`) before the FCVA logic model is applied. Figures from these outputs are produced by Module 10 Script 3.

**Attribute set and inversion logic:**

| Attribute | Direction |
|-----------|-----------|
| Adult mobility | Inverted |
| Habitat specificity | Inverted |
| Mobility and dispersal or early life stages | Inverted |
| Species range | Not inverted (analog for Sensitivity to Temperature) |

> **Open decision:** `Species range` substitutes for "Sensitivity to Temperature" used in prior CVAs. See `09-distributional-change-potential/ReadMe.md` for rationale and alternatives.

#### Script 1 — Baseline DCP scores

Reads `attribute_means_uscar.csv`, applies inversion, and applies the FCVA logic model (same `rank_threshold` as Modules 4 and 8, sourced from `config.R`) to produce a DCP rank for each stock.

#### Script 2 — Bootstrap uncertainty

Mirrors Module 8: builds 20-vote draw piles from `sensitivity_tallies_long.csv` (swapping tally counts for inverted attributes), runs a baseline reproduction gate, then executes 10,000 bootstrap iterations with `bootstrap_seed = 99`. Stocks are flagged borderline if the dominant rank accounts for fewer than 75% of iterations. 

**Reads from:**
- `outputs/{run_label}/final-scores-compiled/overall-vulnerability-rankings/attribute_means_uscar.csv`
- `outputs/final-tallies-long/sensitivity_tallies_long.csv`

**Key outputs:** `outputs/{run_label}/distribution-change-potential/distributional_change_potential_uscar.csv`, `outputs/{run_label}/distribution-change-potential/distributional_change_bootstrap_uscar.csv`

---

### Module 10 — Figures

**Scripts:** `10-figures/1-plot-scoring-distributions.R`, `2-plot-uncertainty-figures.R`, `3-plot-distributional-change.R`, `4-plot-overall-vulnerability.R`

Produces all 13 publication figures. This module has no analysis logic — it reads finalized outputs from Modules 4–9 and writes PNG files to `figures/{run_label}/`. All five scripts source `config.R` for the active run label, shared color palettes, display labels, and stock name recoding.

#### Script 1 — Score distributions

Produces attribute score boxplots (using `attribute_means_uscar.csv`) and per-stock LMHV tally distribution figures. Exposure figures are filtered to expert-approved factor × stock pairs using `data/exposure-factor-filter-long.csv`. Also produces the directional effect summary bar chart.

**Key outputs:** `fig_attribute_score_boxplot_combined.png`, `fig_sensitivity_attribute_score_boxplot.png`, `fig_exposure_attribute_score_boxplot.png`, `fig_directional_effect_summary.png`, `fig_sensitivity_tally_distributions_by_stock.png`, `fig_exposure_tally_distributions_by_stock.png`

#### Script 2 — Uncertainty figures

Produces LOO influence bar charts and bootstrap uncertainty stacked-bar panels for both sensitivity ranks and directional effects.

**Key outputs:** `fig_loo_bar_plots.png`, `fig_bootstrap_uncertainty.png`

#### Script 3 — Distributional change figures

Produces two publication figures from Module 9 outputs:
- **Figure A** (`fig_distributional_change_ranks.png`) — stacked column showing stocks by DCP rank category with certainty-encoded font
- **Figure B** (`fig_distributional_change_vs_vulnerability.png`) — 4×4 cross-plot of DCP rank vs. overall climate vulnerability; the High/VH vulnerability + Low/Moderate DCP quadrant identifies stocks of highest management concern

#### Script 4 — Overall vulnerability

Produces the combined vulnerability summary figure (`fig_overall_vulnerability.png`) as a two-panel layout:
- **Panel A** — 4×4 tile grid: Climate Exposure rank (x) × Biological Sensitivity rank (y); stock labels colored and styled by bootstrap certainty
- **Panel B** — Directional effect column chart; stock labels styled by bootstrap certainty

**Reads from:**
- `outputs/{run_label}/final-scores-compiled/overall-vulnerability-rankings/overall_vulnerability_scores_uscar.csv`
- `outputs/final-tallies-long/directional_effect_tallies_by_stock.csv`
- `outputs/{run_label}/analyses/uncertainty-loo/final-tables/table_bootstrap_uncertainty_stock.csv`
- `outputs/{run_label}/analyses/uncertainty-loo/final-tables/table_directional_effect_bootstrap.csv`

---

## Key Cross-Workflow Data Files

The table below lists the files consumed by more than one module.

| File | Produced by | Consumed by | Key columns |
|------|-------------|-------------|-------------|
| `data/final-scores/*.xlsx` | Reviewers | Modules 4, 5, 6, 7 | One tab per stock; scores in column K, tallies in columns M–P |
| `quantitative-exposure-attribute-scores-all.csv` | Module 2 | Module 4 Script 2, Module 7 Script 2 | `stock_name`, `attribute_name`, `spatial_extent`, `score` |
| `table_final_attribute_scores_all.csv` | Module 4 Script 1 | Module 4 Script 2 | `stock_name`, `Attribute_name`, `Attribute_type`, `Final_score`, `Scorer` |
| `final_scores_uscar.csv` | Module 4 Script 2 | Module 4 Script 3 | `stock_name`, `attribute_type`, `score_type`, `attribute_name`, `scorer`, `score` |
| `outputs/{run_label}/…/attribute_means_uscar.csv` | Module 4 Script 3 | Module 8, Module 9 Script 1, Module 10 Scripts 1 & 4 | `stock_name`, `attribute_type`, `score_type`, `attribute_name`, `attribute_mean` |
| `outputs/{run_label}/…/overall_vulnerability_scores_uscar.csv` | Module 4 Script 3 | Module 8, Module 10 Scripts 3 & 4 | `stock_name`, `Exp_score`, `Exp_rank`, `Sens_score`, `Sens_rank`, `Vuln_score`, `Vuln_rank` |
| `data/exposure-factor-filter-long.csv` | Module 4 Script 3 | Module 8, Module 10 Script 1 | `stock_name`, `attribute_name`, `include` (TRUE/FALSE) |
| `directional_effect_summary_by-stock.csv` | Module 5 | Not consumed downstream (standalone supplementary output) | `stock_name`, `Positive`, `Neutral`, `Negative`, `wt_avg`, `directional_effect` |
| `overall_data_quality_summary_by_stock.csv` | Module 6 | Not consumed downstream (standalone supplementary output) | `stock_name`, `prop_ge_2`, `data_quality_rank` |
| `sensitivity_tallies_long.csv` | Module 7 Script 2 | Module 8, Module 9 Script 2, Module 10 Script 1 | `reviewer_id`, `stock_name`, `attribute_name`, `tally_L`, `tally_M`, `tally_H`, `tally_VH` |
| `directional_effect_tallies_long.csv` | Module 7 Script 2 | Module 8 | `reviewer_id`, `stock_name`, `effect_category`, `tally` |
| `directional_effect_tallies_by_stock.csv` | Module 7 Script 2 | Module 10 Script 4 | `stock_name`, `tally_neg`, `tally_neut`, `tally_pos`, `n_tallies` |
| `outputs/{run_label}/…/table_bootstrap_uncertainty_stock.csv` | Module 8 | Module 10 Scripts 2 & 4 | `stock_name`, `vuln_rank`, `prop` |
| `outputs/{run_label}/…/table_directional_effect_bootstrap.csv` | Module 8 | Module 10 Scripts 2 & 4 | `stock_name`, `dir_rank`, `prop` |
| `outputs/{run_label}/…/distributional_change_potential_uscar.csv` | Module 9 Script 1 | Module 9 Script 2, Module 10 Script 3 | `stock_name`, `dcp_rank`, `dcp_numeric` |
| `outputs/{run_label}/…/distributional_change_bootstrap_uscar.csv` | Module 9 Script 2 | Module 10 Script 3 | `stock_name`, `dominant_rank`, `dominant_prop`, `borderline` |

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

Stock names appear in multiple forms across the pipeline:
- **Compiled analysis outputs** (e.g., `overall_vulnerability_scores_uscar.csv`) use **sentence case** (e.g., `"Red hind"`, `"Atlantic thread herring"`).
- **Reviewer workbooks and tally files** use **Title Case** (e.g., `"Red Hind"`, `"Redhind"`, `"Atlantic Herring"`).

All figure scripts (Module 10) normalize stock names to a consistent **Title Case display form** using the canonical `stock_name_recode` lookup defined in `config.R`. This lookup handles both simple capitalization differences and the four stocks whose workbook names differ substantively from their assessment names:

| Input name (analysis outputs or tallies) | Display name (figures) | Reason |
|---|---|---|
| `Atlantic thread herring` | `Atlantic Herring` | Different common name |
| `Long-spined sea urchin` | `Diadema` | Different common name |
| `Red hind` | `Red Hind` | Capitalization + two-word form |
| `Redhind` | `Red Hind` | Workbook single-word form |
| `Sea cucumbers` | `Sea Cucumber` | Plural vs. singular |

The analysis modules (4, 7, 8, 9) use a complementary normalization: four manual overrides followed by `stringr::str_to_sentence()`, which converts workbook Title Case names into the sentence-case form used in compiled outputs.

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

> **Run configuration note:** The `cross_region_comparable` run drops attributes #3 (Genetic diversity) and #8 (Predation and competition dynamics), using the same 12-attribute set as all other published NOAA FCVAs. The `broadened_distribution` run retains all 14 attributes. The dropped attributes are specified by `sens_attrs_drop` in `config.R` and are excluded in Module 4 Script 3 and Module 8 before any calculations.

> **Naming note — "Stock Size Status":** This attribute appears in reviewer workbooks as "Stock Size/Status" (with slash). The `standardize_attribute_names()` function in Module 7 Script 1 removes the slash, producing "Stock Size Status" in the tally tables. However, `attribute_means_uscar.csv` (produced by Module 4 Script 3) stores the attribute as "Stock size/status" (sentence case, with slash). The figure scripts in Module 10 handle this mismatch with an explicit `dplyr::recode()` call. The canonical form for lookups and figure labels is **"Stock Size Status"** (no slash).

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

We thank Tyler Loughran (NOAA) and Dan Crear (ICATTC) for their assistance in this work.

**References:**

- Morrison, W. E., Nelson, M. W., Howard, J. F., Stecher, H. A., Sheridan, P. F., & Resnick, M. L. (2015). Methodology for assessing the vulnerability of marine fish and shellfish species to a changing climate. *NOAA Technical Memorandum NMFS-OSF-3*. https://doi.org/10.7289/V5TM782J
- Hare, J. A., Morrison, W. E., Nelson, M. W., Stachura, M. M., Teeters, E. J., Griffis, R. B., Alexander, M. A., Scott, J. D., Alade, L., Bell, R. J., et al. (2016). A vulnerability assessment of fish and invertebrates to climate change on the Northeast U.S. continental shelf. *PLOS ONE*, 11(2), e0146756. https://doi.org/10.1371/journal.pone.0146756
- Loughran, C. E., Hazen, E. L., Brodie, S., Jacox, M. G., Whitney, F. A., Payne, M. R., et al. (2025). A climate vulnerability assessment of highly migratory species in the Northwest Atlantic Ocean. *PLOS Climate*, 4(8), e0000530. https://doi.org/10.1371/journal.pclm.0000530
- Craig, J. K., et al. (2025). A climate vulnerability assessment for managed species in the South Atlantic Bight. *PLOS Climate*. https://doi.org/10.1371/journal.pclm.0000543
