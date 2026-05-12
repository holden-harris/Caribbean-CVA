# Table Inventory — Caribbean CVA Pipeline

This document catalogs all 51 CSV tables produced across Modules 00–10. Use it as a data dictionary and pipeline map.

**51 tables total:** 40 main analysis tables + 11 QA/validation tables.

**Column definitions:**

| Column | Meaning |
|---|---|
| `table` | CSV filename (no path) |
| `description` | Purpose of the table |
| `output-dir` | Path relative to repo root |
| `module` | Numbered module folder |
| `script` | Script that writes the file |
| `rows` | What one row represents |
| `columns` | All column names |
| `consumed-by` | Downstream module/script(s) that `read_csv` this file; `—` = publication output only |
| `type` | `intermediate` / `final` / `figure-input` / `QA` |
| `format` | `long` / `wide` |

---

## Main analysis tables

### Module 00 — FishBase query

| table | description | output-dir | module | script | rows | columns | consumed-by | type | format |
|---|---|---|---|---|---|---|---|---|---|
| `fishbase_species_attributes.csv` | Life-history and ecology traits pulled from FishBase API for each assessed species | `outputs/` | `00-query-species-attributes-from-FishBase` | `Query-species-attributes-from-FishBase.R` | 1 per species | `scientific_name, Species, FBname, BodyShapeI, DepthRangeShallow, DepthRangeDeep, DepthRangeComShallow, DepthRangeComDeep, LongevityWild, LongevityCaptive, Vulnerability, VulnerabilityClimate, Length, CommonLength, Weight, n_growth_studies, K_avg, Loo_avg, ReproMode, Spawning, RepGuild1, RepGuild2, AddInfos, DietTroph, DietSeTroph, DietRemark, AddRems` | — | final | wide |

### Module 02 — Exposure anomalies

| table | description | output-dir | module | script | rows | columns | consumed-by | type | format |
|---|---|---|---|---|---|---|---|---|---|
| `quantitative-exposure-attribute-scores-all.csv` | CMIP6-derived quantitative exposure scores for each stock × factor × spatial extent combination | `outputs/final-scores-compiled/` | `02-exposure-anomalies` | `exposure-factor-scores.R` | 1 per stock × exposure_factor × spatial_extent (~975 rows) | `stock_name, quantitative_exposure_factor, full_names, spatial_extent, attribute_score, tally_L, tally_M, tally_H, tally_VH, n_tallies` | `04/1-extract-final-scores-from-all-reviewers.R`, `07/1-extract-tally-scores.R` | intermediate | long |

### Module 03 — Preliminary sensitivity scoring

| table | description | output-dir | module | script | rows | columns | consumed-by | type | format |
|---|---|---|---|---|---|---|---|---|---|
| `score_table_all_clean.csv` | Pre-workshop sensitivity attribute scores from all pre-work reviewers | `data/preliminary-scores/` | `03-prelim-sensitivity-attribute-scoring` | `extract-preliminary-scores-from-all-reviewers.R` | 1 per scorer × stock × attribute | `SourceFile, Scorer, stock_name, row_idx, Attribute_name, Data_quality, Scoring_rank_1, Scoring_rank_2, Scoring_rank_3, Scoring_rank_4, Attribute_type` | — | final | long |

### Module 04 — Final attribute-exposure scoring

**Script 1:** `1-extract-final-scores-from-all-reviewers.R`

| table | description | output-dir | module | script | rows | columns | consumed-by | type | format |
|---|---|---|---|---|---|---|---|---|---|
| `table_final_attribute_scores_all.csv` | Raw final sensitivity and qualitative exposure scores extracted from each reviewer's workbook | `outputs/final-scores-compiled/final-attribute-scores/` | `04-final-attribute-exposure-scoring` | `1-extract-final-scores-from-all-reviewers.R` | 1 per reviewer × stock × attribute | `Attribute_name, Final_score, SourceFile, Scorer, stock_name, row_idx, Attribute_type` | `04/2-summarize-attribute-scores.R` | intermediate | long |

**Script 2:** `2-summarize-attribute-scores.R`

| table | description | output-dir | module | script | rows | columns | consumed-by | type | format |
|---|---|---|---|---|---|---|---|---|---|
| `final_scores_compiled.csv` | Reviewer-level final scores for all stocks and regions combined (sensitivity + qualitative + quantitative exposure) | `outputs/final-scores-compiled/overall-vulnerability-rankings/` | `04-final-attribute-exposure-scoring` | `2-summarize-attribute-scores.R` | 1 per reviewer × stock × attribute (all regions) | `stock_name, region, attribute_type, score_type, attribute_name, scorer, score` | — | intermediate | long |
| `final_scores_uscar.csv` | Same as `final_scores_compiled.csv` filtered to U.S. Caribbean region only | `outputs/final-scores-compiled/overall-vulnerability-rankings/` | `04-final-attribute-exposure-scoring` | `2-summarize-attribute-scores.R` | 1 per reviewer × stock × attribute (U.S. Caribbean only) | `stock_name, region, attribute_type, score_type, attribute_name, scorer, score` | `04/3-calculate-overall-vulnerability-scores.R` | intermediate | long |

**Script 3:** `3-calculate-overall-vulnerability-scores.R`

| table | description | output-dir | module | script | rows | columns | consumed-by | type | format |
|---|---|---|---|---|---|---|---|---|---|
| `exposure-factor-filter-long.csv` | Per-stock exposure factor inclusion flags; marks which of the 13 CMIP6 factors apply to each stock | `data/` | `04-final-attribute-exposure-scoring` | `3-calculate-overall-vulnerability-scores.R` | 1 per stock × exposure_factor | `stock_name, attribute_name, include` | `08/uncertainty-analyses.R` | intermediate | long |
| `attribute_means_uscar.csv` | Mean (and SD) final score per stock × attribute across all reviewers, for both sensitivity and exposure attributes | `outputs/final-scores-compiled/overall-vulnerability-rankings/` | `04-final-attribute-exposure-scoring` | `3-calculate-overall-vulnerability-scores.R` | 1 per stock × attribute | `stock_name, attribute_type, score_type, attribute_name, attribute_mean, attribute_sd, n_scores, attribute_mean_lmhv` | `08/uncertainty-analyses.R`, `09/1-calculate-distributional-change-potential.R` | final | long |
| `component_scores_uscar.csv` | Sensitivity and exposure component scores per stock, including attribute-count thresholds used by the FCVA logic model | `outputs/final-scores-compiled/overall-vulnerability-rankings/` | `04-final-attribute-exposure-scoring` | `3-calculate-overall-vulnerability-scores.R` | 1 per stock × attribute_type | `stock_name, attribute_type, mean_attribute_score, sd_attribute_score, n_attributes, n_ge_2_5, n_ge_3_0, n_ge_3_5, component_score, component_rank` | `10/4-plot-overall-vulnerability.R` | final | long |
| `overall_vulnerability_scores_uscar.csv` | Baseline overall vulnerability scores and ranks (Sensitivity, Exposure, Overall) for each stock | `outputs/final-scores-compiled/overall-vulnerability-rankings/` | `04-final-attribute-exposure-scoring` | `3-calculate-overall-vulnerability-scores.R` | 1 per stock | `stock_name, Exp_score, Exp_rank, Sens_score, Sens_rank, n_Exp_fact, Vuln_score, Vuln_rank` | `08/uncertainty-analyses.R`, `09/1-calculate-distributional-change-potential.R`, `10/4-plot-overall-vulnerability.R` | final | wide |

### Module 05 — Directional effect scoring

**Script 1:** `1-extract-directional-effect.R`

| table | description | output-dir | module | script | rows | columns | consumed-by | type | format |
|---|---|---|---|---|---|---|---|---|---|
| `table_directional_effect.csv` | Raw directional effect (Positive/Neutral/Negative) tallies extracted from each reviewer's workbook | `outputs/final-scores-compiled/directional-effect/` | `05-final-directional-effect-scoring` | `1-extract-directional-effect.R` | 1 per reviewer × stock | `SourceFile, Scorer, stock_name, row_idx, Directional_effect, Directional_score` | `05/2-summarize-directional-effect.R` | intermediate | long |

**Script 2:** `2-summarize-directional-effect.R`

| table | description | output-dir | module | script | rows | columns | consumed-by | type | format |
|---|---|---|---|---|---|---|---|---|---|
| `directional_effect_wide_all.csv` | Reviewer-level directional effect tally counts and weighted average, pivoted wide | `outputs/final-scores-compiled/directional-effect/` | `05-final-directional-effect-scoring` | `2-summarize-directional-effect.R` | 1 per reviewer × stock | `Scorer, stock_name, Positive, Neutral, Negative, n_tallies, wt_avg, overall` | — | intermediate | wide |
| `directional_effect_summary_by-stock.csv` | Stock-level directional effect summary: pooled tally counts, weighted average, and final classification | `outputs/final-scores-compiled/directional-effect/` | `05-final-directional-effect-scoring` | `2-summarize-directional-effect.R` | 1 per stock | `stock_name, n_positive, n_neutral, n_negative, n_tallies, n_reviewers, wt_avg, overall` | `10/5-produce-results-tables.R` | final | wide |

### Module 06 — Data quality scoring

**Script 1:** `1-extract-data-quality-scores.R`

| table | description | output-dir | module | script | rows | columns | consumed-by | type | format |
|---|---|---|---|---|---|---|---|---|---|
| `table_data_quality_scores_extracted.csv` | Raw data quality scores (0–3 scale) extracted from each reviewer's workbook for every stock × attribute | `outputs/final-scores-compiled/data-quality/` | `06-final-data-quality-scoring` | `1-extract-data-quality-scores.R` | 1 per reviewer × stock × attribute | `SourceFile, Scorer, stock_name, row_idx, Attribute_type, Attribute_name, Data_quality_score, score_cell, label_cell, Data_quality_label, score_ge_2` | `06/2-summarize-data-quality-scores.R` | intermediate | long |

**Script 2:** `2-summarize-data-quality-scores.R`

| table | description | output-dir | module | script | rows | columns | consumed-by | type | format |
|---|---|---|---|---|---|---|---|---|---|
| `overall_data_quality_summary_by_stock.csv` | Stock-level data quality summary: count of each quality tier, proportion ≥ 2, mean score, and final rank | `outputs/final-scores-compiled/data-quality/` | `06-final-data-quality-scoring` | `2-summarize-data-quality-scores.R` | 1 per stock | `stock_name, n_3, n_2, n_1, n_0, n_ge_2, prop_ge_2, mean_score, sd_score, data_quality_rank` | `10/5-produce-results-tables.R` | final | wide |

### Module 07 — Scoring distributions

**Script 1:** `1-extract-tally-scores.R` — raw workbook extracts

| table | description | output-dir | module | script | rows | columns | consumed-by | type | format |
|---|---|---|---|---|---|---|---|---|---|
| `table_sensitivity_tallies_long.csv` | Reviewer-level tally counts (L/M/H/VH) for all sensitivity and rigidity attributes, extracted directly from workbooks | `outputs/analyses/1-inputs/` | `07-scoring-distributions` | `1-extract-tally-scores.R` | 1 per reviewer × stock × sensitivity attribute | `source_file, reviewer_id, stock_name, sheet_name, attribute_type, attribute_name, row_num, final_score, data_quality_index, tally_L, tally_M, tally_H, tally_VH, n_tallies` | `07/2-finalize-tally-tables.R` | intermediate | long |
| `table_qualitative_exposure_tallies_long.csv` | Reviewer-level tally counts for all qualitative exposure attributes, extracted directly from workbooks | `outputs/analyses/1-inputs/` | `07-scoring-distributions` | `1-extract-tally-scores.R` | 1 per reviewer × stock × exposure attribute | `source_file, reviewer_id, stock_name, sheet_name, attribute_type, attribute_name, row_num, final_score, data_quality_index, tally_L, tally_M, tally_H, tally_VH, n_tallies` | `07/2-finalize-tally-tables.R` | intermediate | long |
| `table_directional_effect_tallies_long.csv` | Reviewer-level directional effect tally counts, extracted directly from workbooks | `outputs/analyses/1-inputs/` | `07-scoring-distributions` | `1-extract-tally-scores.R` | 1 per reviewer × stock × effect_category | `source_file, reviewer_id, stock_name, sheet_name, effect_category, row_num, tally, n_tallies` | `07/2-finalize-tally-tables.R` | intermediate | long |
| `table_all_qualitative_tallies_long.csv` | Sensitivity and qualitative exposure tallies combined into a single long table | `outputs/analyses/1-inputs/` | `07-scoring-distributions` | `1-extract-tally-scores.R` | 1 per reviewer × stock × attribute (sensitivity + qual. exposure) | `source_file, reviewer_id, stock_name, sheet_name, attribute_type, attribute_name, row_num, final_score, data_quality_index, tally_L, tally_M, tally_H, tally_VH, n_tallies` | — | intermediate | long |

**Script 2:** `2-finalize-tally-tables.R` — analysis-ready final tables

| table | description | output-dir | module | script | rows | columns | consumed-by | type | format |
|---|---|---|---|---|---|---|---|---|---|
| `sensitivity_tallies_long.csv` | Finalized reviewer-level sensitivity tally table (`source_file` and `sheet_name` dropped; analysis-ready) | `outputs/final-tallies-long/` | `07-scoring-distributions` | `2-finalize-tally-tables.R` | 1 per reviewer × stock × sensitivity attribute | `reviewer_id, stock_name, attribute_type, attribute_name, row_num, final_score, data_quality_index, tally_L, tally_M, tally_H, tally_VH, n_tallies` | `08/uncertainty-analyses.R`, `09/2-bootstrap-distributional-change.R` | final | long |
| `sensitivity_tallies_by_stock.csv` | Sensitivity tallies pooled across reviewers per stock × attribute, with reviewer count | `outputs/final-tallies-long/` | `07-scoring-distributions` | `2-finalize-tally-tables.R` | 1 per stock × sensitivity attribute | `stock_name, attribute_name, tally_L, tally_M, tally_H, tally_VH, n_scorers, n_tallies` | `10/1-plot-scoring-distributions.R` | final | wide |
| `directional_effect_tallies_long.csv` | Finalized reviewer-level directional effect tally table | `outputs/final-tallies-long/` | `07-scoring-distributions` | `2-finalize-tally-tables.R` | 1 per reviewer × stock × effect_category | `reviewer_id, stock_name, effect_category, tally` | `08/uncertainty-analyses.R` | final | long |
| `directional_effect_tallies_by_stock.csv` | Directional effect tallies pooled across reviewers per stock | `outputs/final-tallies-long/` | `07-scoring-distributions` | `2-finalize-tally-tables.R` | 1 per stock | `stock_name, tally_neg, tally_neut, tally_pos, n_scorers, n_tallies` | `10/1-plot-scoring-distributions.R` | final | wide |
| `exposure_tallies_long.csv` | Finalized reviewer-level exposure tally table combining qualitative and quantitative exposure attributes | `outputs/final-tallies-long/` | `07-scoring-distributions` | `2-finalize-tally-tables.R` | 1 per reviewer × stock × exposure attribute | `reviewer_id, stock_name, attribute_type, attribute_name, tally_L, tally_M, tally_H, tally_VH, n_tallies` | `08/uncertainty-analyses.R`, `10/1-plot-scoring-distributions.R` | final | long |
| `exposure_tallies_by_stock.csv` | Exposure tallies pooled across reviewers per stock × attribute | `outputs/final-tallies-long/` | `07-scoring-distributions` | `2-finalize-tally-tables.R` | 1 per stock × exposure attribute | `stock_name, attribute_type, attribute_name, tally_L, tally_M, tally_H, tally_VH, n_tallies` | `10/1-plot-scoring-distributions.R` | final | wide |

### Module 08 — Uncertainty analysis

| table | description | output-dir | module | script | rows | columns | consumed-by | type | format |
|---|---|---|---|---|---|---|---|---|---|
| `table_bootstrap_uncertainty_stock.csv` | Bootstrap rank distribution for overall vulnerability; one row per stock × possible rank | `outputs/analyses/uncertainty-loo/final-tables/` | `08-uncertainty-analysis` | `uncertainty-analyses.R` | 1 per stock × vuln_rank | `stock_name, vuln_rank, n, prop, borderline` | `10/2-plot-uncertainty-figures.R` | final | long |
| `table_bootstrap_uncertainty_sensitivity_component.csv` | Bootstrap rank distribution for the sensitivity component only | `outputs/analyses/uncertainty-loo/final-tables/` | `08-uncertainty-analysis` | `uncertainty-analyses.R` | 1 per stock × sens_rank | `stock_name, sens_rank, n, prop` | `10/2-plot-uncertainty-figures.R` | final | long |
| `table_directional_effect_bootstrap.csv` | Bootstrap rank distribution for directional effect, including baseline rank and weighted mean | `outputs/analyses/uncertainty-loo/final-tables/` | `08-uncertainty-analysis` | `uncertainty-analyses.R` | 1 per stock × dir_rank | `stock_name, dir_rank, n, prop, baseline_dir_rank, baseline_w_mean` | `10/5-produce-results-tables.R` | final | long |
| `table_leave_one_out_sensitivity_long.csv` | LOO influence results for sensitivity attributes: baseline vs. new ranks when each attribute is omitted, per stock | `outputs/analyses/uncertainty-loo/final-tables/` | `08-uncertainty-analysis` | `uncertainty-analyses.R` | 1 per stock × attribute_omitted | `stock_name, attribute_omitted, baseline_sens_rank, new_sens_rank, baseline_vuln_rank, new_vuln_rank, baseline_vuln_score, new_vuln_score, rank_changed, rank_change_direction` | `10/2-plot-uncertainty-figures.R` | final | long |
| `table_leave_one_out_sensitivity_summary.csv` | LOO sensitivity influence summary: how many stocks change rank when each attribute is omitted | `outputs/analyses/uncertainty-loo/final-tables/` | `08-uncertainty-analysis` | `uncertainty-analyses.R` | 1 per attribute_omitted | `attribute_omitted, n_stocks_tested, n_rank_changed, n_rank_lower, n_rank_higher, prop_rank_changed` | `10/2-plot-uncertainty-figures.R` | final | wide |
| `table_leave_one_out_exposure_long.csv` | LOO influence results for exposure factors: baseline vs. new ranks when each factor is omitted, per stock | `outputs/analyses/uncertainty-loo/final-tables/` | `08-uncertainty-analysis` | `uncertainty-analyses.R` | 1 per stock × factor_omitted | `stock_name, factor_omitted, baseline_exp_rank, new_exp_rank, baseline_vuln_rank, new_vuln_rank, baseline_vuln_score, new_vuln_score, rank_changed, rank_change_direction` | `10/2-plot-uncertainty-figures.R` | final | long |
| `table_leave_one_out_exposure_summary.csv` | LOO exposure influence summary: how many stocks change rank when each exposure factor is omitted | `outputs/analyses/uncertainty-loo/final-tables/` | `08-uncertainty-analysis` | `uncertainty-analyses.R` | 1 per factor_omitted | `factor_omitted, n_stocks_tested, n_rank_changed, n_rank_lower, n_rank_higher, prop_rank_changed` | `10/2-plot-uncertainty-figures.R` | final | wide |
| `table_exposure_factor_scores_for_plot.csv` | Mean exposure factor scores with component rank, formatted for figure generation | `outputs/analyses/uncertainty-loo/final-tables/` | `08-uncertainty-analysis` | `uncertainty-analyses.R` | 1 per stock × exposure_factor | `stock_name, exposure_factor, mean_score, component_rank_baseline` | `10/2-plot-uncertainty-figures.R` | figure-input | long |
| `table_sensitivity_attribute_scores_for_plot.csv` | Mean sensitivity attribute scores with component rank, formatted for figure generation | `outputs/analyses/uncertainty-loo/final-tables/` | `08-uncertainty-analysis` | `uncertainty-analyses.R` | 1 per stock × attribute | `stock_name, attribute_name, mean_score, component_rank_baseline` | `10/2-plot-uncertainty-figures.R` | figure-input | long |
| `table_bootstrap_final_summary.csv` | Publication-ready bootstrap summary: dominant rank and rank proportion columns for all vulnerability components | `outputs/analyses/uncertainty-loo/final-tables/` | `08-uncertainty-analysis` | `uncertainty-analyses.R` | 1 per stock | `Stock, N Exp factors, Exp, Sens, Vul, L, M, H, VH` | `10/5-produce-results-tables.R` | final | wide |

### Module 09 — Distributional change potential

| table | description | output-dir | module | script | rows | columns | consumed-by | type | format |
|---|---|---|---|---|---|---|---|---|---|
| `distributional_change_potential_uscar.csv` | Baseline distributional change potential (DCP) scores and ranks derived from 4 BSAs with mobility attributes inverted | `outputs/distribution-change-potential/` | `09-distributional-change-potential` | `1-calculate-distributional-change-potential.R` | 1 per stock | `stock_name, adult_mobility_inverted, habitat_specificity_inverted, early_life_dispersal_inverted, species_range, dcp_rank, dcp_numeric` | `09/2-bootstrap-distributional-change.R`, `10/3-plot-distributional-change.R` | final | wide |
| `distributional_change_bootstrap_uscar.csv` | Bootstrap rank distribution for DCP scores | `outputs/distribution-change-potential/` | `09-distributional-change-potential` | `2-bootstrap-distributional-change.R` | 1 per stock | `stock_name, prop_L, prop_M, prop_H, prop_VH, dominant_rank, dominant_prop, borderline` | `10/3-plot-distributional-change.R` | final | wide |
| `distributional_change_full_uscar.csv` | Combined baseline DCP scores and bootstrap uncertainty results in a single table | `outputs/distribution-change-potential/` | `09-distributional-change-potential` | `2-bootstrap-distributional-change.R` | 1 per stock | `stock_name, adult_mobility_inverted, habitat_specificity_inverted, early_life_dispersal_inverted, species_range, dcp_rank, dcp_numeric, prop_L, prop_M, prop_H, prop_VH, dominant_rank, dominant_prop, borderline` | `10/3-plot-distributional-change.R` | final | wide |

### Module 10 — Figures

| table | description | output-dir | module | script | rows | columns | consumed-by | type | format |
|---|---|---|---|---|---|---|---|---|---|
| `table_directional_effect_results.csv` | Publication results table: directional effect category, weighted mean, tally counts, and bootstrap proportions by stock | `outputs/tables/` | `10-figures` | `5-produce-results-tables.R` | 1 per stock | `Stock, Vuln_rank, Dir_effect, Wt_mean, N_negative, N_neutral, N_positive, N_tallies, Boot_Negative, Boot_Neutral, Boot_Positive, Dominant_prop, Borderline` | — | final | wide |
| `table_data_quality_results.csv` | Publication results table: data quality rank, proportion ≥ 2, mean score, and tier counts by stock | `outputs/tables/` | `10-figures` | `5-produce-results-tables.R` | 1 per stock | `Stock, Vuln_rank, Data_quality_rank, Prop_ge_2, Mean_score, N_adequate, N_limited, N_expert, N_nodata` | — | final | wide |

---

## QA / validation tables

These files are written by the pipeline for data-integrity checks and are not consumed by downstream analysis scripts.

| table | description | output-dir | module | script | rows | columns |
|---|---|---|---|---|---|---|
| `n_reviews_final_scores.csv` | Number of reviewer workbooks successfully read per stock | `outputs/final-scores-compiled/final-attribute-scores/` | `04-final-attribute-exposure-scoring` | `1-extract-final-scores-from-all-reviewers.R` | 1 per stock | `stock_name, n_reviews, reviewers` |
| `exposure_factor_qa.csv` | Number and list of exposure factors included per stock after expert filtering | `outputs/final-scores-compiled/overall-vulnerability-rankings/` | `04-final-attribute-exposure-scoring` | `3-calculate-overall-vulnerability-scores.R` | 1 per stock | `N-Exp-Fact, Exposure-Factors, stock_name` |
| `n_reviews_directional_effect.csv` | Number of reviewer workbooks with directional effect data per stock | `outputs/final-scores-compiled/directional-effect/` | `05-final-directional-effect-scoring` | `1-extract-directional-effect.R` | 1 per stock | `stock_name, n_reviews, reviewers` |
| `qa_data_quality_by_scorers.csv` | Number of reviewer workbooks with data quality data per stock | `outputs/final-scores-compiled/data-quality/` | `06-final-data-quality-scoring` | `1-extract-data-quality-scores.R` | 1 per stock | `stock_name, n_reviews, reviewers` |
| `qa_data_quality_by_stock.csv` | Row-count check for extracted data quality scores per stock | `outputs/final-scores-compiled/data-quality/` | `06-final-data-quality-scoring` | `1-extract-data-quality-scores.R` | 1 per stock | `stock_name, n_scorers, n_rows` |
| `qa_tally_row_checks.csv` | Row-level tally integrity checks: flags rows where tallies don't sum to 5 or contain NAs | `outputs/analyses/1-inputs/qa/` | `07-scoring-distributions` | `1-extract-tally-scores.R` | 1 per reviewer × stock × attribute | `stock_name, reviewer_id, attribute_name, row_num, attribute_type, final_score, data_quality_index, tally_L, tally_M, tally_H, tally_VH, n_tallies, n_tallies_ok, has_na_tally, any_row_issue` |
| `qa_reviewer_stock_coverage.csv` | Reviewer × stock coverage matrix showing attribute count per reviewer per stock | `outputs/analyses/1-inputs/qa/` | `07-scoring-distributions` | `1-extract-tally-scores.R` | 1 per stock (reviewers as columns) | `stock_name, [reviewer_id columns]` |
| `qa_attribute_tally_summary.csv` | Pooled tally sums per stock × attribute compared against expected values; flags mismatches | `outputs/analyses/1-inputs/qa/` | `07-scoring-distributions` | `1-extract-tally-scores.R` | 1 per stock × attribute | `stock_name, attribute_type, attribute_name, n_reviewers, tally_L, tally_M, tally_H, tally_VH, pooled_tally_sum, expected_tally_sum, tally_sum_ok, pooled_mean_score` |
| `qa_directional_effect_summary.csv` | Directional effect tally count check per stock | `outputs/analyses/1-inputs/qa/` | `07-scoring-distributions` | `1-extract-tally-scores.R` | 1 per stock | `stock_name, n_reviewers, total_tally` |
| `01_preanalysis_check_log.csv` | Pre-analysis validation log for the uncertainty module: named checks with pass/fail status and detail | `outputs/analyses/uncertainty-loo/intermediate/` | `08-uncertainty-analysis` | `uncertainty-analyses.R` | 1 per check | `check_name, status, detail` |
| `02_baseline_reproduced_scores.csv` | Reproduced baseline vulnerability scores computed from tally tables; verifies pipeline reproducibility before bootstrap runs | `outputs/analyses/uncertainty-loo/intermediate/` | `08-uncertainty-analysis` | `uncertainty-analyses.R` | 1 per stock | `stock_name, sensitivity_rank_repro, sensitivity_score_numeric_repro, exposure_rank_repro, exposure_score_numeric_repro, vulnerability_score_numeric_repro, vulnerability_rank_repro` |
