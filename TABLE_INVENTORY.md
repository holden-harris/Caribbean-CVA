# Table Inventory — Caribbean CVA Pipeline

51 CSV outputs across Modules 00–10, in pipeline order. **Type key:** `intermediate` = consumed within the pipeline; `final` = pipeline endpoint or publication supplement; `figure-input` = pre-formatted for a specific figure script; `QA` = data integrity check, not consumed downstream.

<table>
<colgroup>
<col width="11%">
<col width="19%">
<col width="10%">
<col width="4%">
<col width="9%">
<col width="8%">
<col width="25%">
<col width="8%">
<col width="4%">
<col width="2%">
</colgroup>
<thead>
<tr>
<th>table</th>
<th>description</th>
<th>output-dir</th>
<th>module</th>
<th>script</th>
<th>rows</th>
<th>columns</th>
<th>consumed-by</th>
<th>type</th>
<th>format</th>
</tr>
</thead>
<tbody>
<tr>
<td><code>fishbase_species_attributes.csv</code></td>
<td>Life-history and ecology traits pulled from FishBase API for each assessed species</td>
<td><code>outputs/</code></td>
<td>00</td>
<td><code>Query-species-attributes-from-FishBase.R</code></td>
<td>1 per species</td>
<td><code>scientific_name, Species, FBname, BodyShapeI, DepthRangeShallow, DepthRangeDeep, DepthRangeComShallow, DepthRangeComDeep, LongevityWild, LongevityCaptive, Vulnerability, VulnerabilityClimate, Length, CommonLength, Weight, n_growth_studies, K_avg, Loo_avg, ReproMode, Spawning, RepGuild1, RepGuild2, AddInfos, DietTroph, DietSeTroph, DietRemark, AddRems</code></td>
<td>—</td>
<td>final</td>
<td>wide</td>
</tr>
<tr>
<td><code>quantitative-exposure-attribute-scores-all.csv</code></td>
<td>CMIP6-derived quantitative exposure scores for each stock × factor × spatial extent combination</td>
<td><code>outputs/final-scores-compiled/</code></td>
<td>02</td>
<td><code>exposure-factor-scores.R</code></td>
<td>1 per stock × exposure_factor × spatial_extent (~975)</td>
<td><code>stock_name, quantitative_exposure_factor, full_names, spatial_extent, attribute_score, tally_L, tally_M, tally_H, tally_VH, n_tallies</code></td>
<td><code>04/1-extract-final-scores-from-all-reviewers.R</code>, <code>07/1-extract-tally-scores.R</code></td>
<td>intermediate</td>
<td>long</td>
</tr>
<tr>
<td><code>score_table_all_clean.csv</code></td>
<td>Pre-workshop sensitivity attribute scores from all pre-work reviewers</td>
<td><code>data/preliminary-scores/</code></td>
<td>03</td>
<td><code>extract-preliminary-scores-from-all-reviewers.R</code></td>
<td>1 per scorer × stock × attribute</td>
<td><code>SourceFile, Scorer, stock_name, row_idx, Attribute_name, Data_quality, Scoring_rank_1, Scoring_rank_2, Scoring_rank_3, Scoring_rank_4, Attribute_type</code></td>
<td>—</td>
<td>final</td>
<td>long</td>
</tr>
<tr>
<td><code>n_reviews_final_scores.csv</code></td>
<td>Number of reviewer workbooks successfully read per stock</td>
<td><code>outputs/final-scores-compiled/final-attribute-scores/</code></td>
<td>04</td>
<td><code>1-extract-final-scores-from-all-reviewers.R</code></td>
<td>1 per stock</td>
<td><code>stock_name, n_reviews, reviewers</code></td>
<td>—</td>
<td>QA</td>
<td>wide</td>
</tr>
<tr>
<td><code>table_final_attribute_scores_all.csv</code></td>
<td>Raw final sensitivity and qualitative exposure scores extracted from each reviewer's workbook</td>
<td><code>outputs/final-scores-compiled/final-attribute-scores/</code></td>
<td>04</td>
<td><code>1-extract-final-scores-from-all-reviewers.R</code></td>
<td>1 per reviewer × stock × attribute</td>
<td><code>Attribute_name, Final_score, SourceFile, Scorer, stock_name, row_idx, Attribute_type</code></td>
<td><code>04/2-summarize-attribute-scores.R</code></td>
<td>intermediate</td>
<td>long</td>
</tr>
<tr>
<td><code>final_scores_compiled.csv</code></td>
<td>Reviewer-level final scores for all stocks and regions combined (sensitivity + qualitative + quantitative exposure)</td>
<td><code>outputs/final-scores-compiled/overall-vulnerability-rankings/</code></td>
<td>04</td>
<td><code>2-summarize-attribute-scores.R</code></td>
<td>1 per reviewer × stock × attribute (all regions)</td>
<td><code>stock_name, region, attribute_type, score_type, attribute_name, scorer, score</code></td>
<td>—</td>
<td>intermediate</td>
<td>long</td>
</tr>
<tr>
<td><code>final_scores_uscar.csv</code></td>
<td>Same as <code>final_scores_compiled.csv</code> filtered to U.S. Caribbean region only</td>
<td><code>outputs/final-scores-compiled/overall-vulnerability-rankings/</code></td>
<td>04</td>
<td><code>2-summarize-attribute-scores.R</code></td>
<td>1 per reviewer × stock × attribute (U.S. Caribbean only)</td>
<td><code>stock_name, region, attribute_type, score_type, attribute_name, scorer, score</code></td>
<td><code>04/3-calculate-overall-vulnerability-scores.R</code></td>
<td>intermediate</td>
<td>long</td>
</tr>
<tr>
<td><code>exposure-factor-filter-long.csv</code></td>
<td>Per-stock exposure factor inclusion flags; marks which of the 13 CMIP6 factors apply to each stock</td>
<td><code>data/</code></td>
<td>04</td>
<td><code>3-calculate-overall-vulnerability-scores.R</code></td>
<td>1 per stock × exposure_factor</td>
<td><code>stock_name, attribute_name, include</code></td>
<td><code>08/uncertainty-analyses.R</code></td>
<td>intermediate</td>
<td>long</td>
</tr>
<tr>
<td><code>attribute_means_uscar.csv</code></td>
<td>Mean (and SD) final score per stock × attribute across all reviewers, for both sensitivity and exposure attributes</td>
<td><code>outputs/final-scores-compiled/overall-vulnerability-rankings/</code></td>
<td>04</td>
<td><code>3-calculate-overall-vulnerability-scores.R</code></td>
<td>1 per stock × attribute</td>
<td><code>stock_name, attribute_type, score_type, attribute_name, attribute_mean, attribute_sd, n_scores, attribute_mean_lmhv</code></td>
<td><code>08/uncertainty-analyses.R</code>, <code>09/1-calculate-distributional-change-potential.R</code></td>
<td>final</td>
<td>long</td>
</tr>
<tr>
<td><code>component_scores_uscar.csv</code></td>
<td>Sensitivity and exposure component scores per stock, including attribute-count thresholds used by the FCVA logic model</td>
<td><code>outputs/final-scores-compiled/overall-vulnerability-rankings/</code></td>
<td>04</td>
<td><code>3-calculate-overall-vulnerability-scores.R</code></td>
<td>1 per stock × attribute_type</td>
<td><code>stock_name, attribute_type, mean_attribute_score, sd_attribute_score, n_attributes, n_ge_2_5, n_ge_3_0, n_ge_3_5, component_score, component_rank</code></td>
<td><code>10/4-plot-overall-vulnerability.R</code></td>
<td>final</td>
<td>long</td>
</tr>
<tr>
<td><code>overall_vulnerability_scores_uscar.csv</code></td>
<td>Baseline overall vulnerability scores and ranks (Sensitivity, Exposure, Overall) for each stock</td>
<td><code>outputs/final-scores-compiled/overall-vulnerability-rankings/</code></td>
<td>04</td>
<td><code>3-calculate-overall-vulnerability-scores.R</code></td>
<td>1 per stock</td>
<td><code>stock_name, Exp_score, Exp_rank, Sens_score, Sens_rank, n_Exp_fact, Vuln_score, Vuln_rank</code></td>
<td><code>08/uncertainty-analyses.R</code>, <code>09/1-calculate-distributional-change-potential.R</code>, <code>10/4-plot-overall-vulnerability.R</code></td>
<td>final</td>
<td>wide</td>
</tr>
<tr>
<td><code>exposure_factor_qa.csv</code></td>
<td>Number and list of exposure factors included per stock after expert filtering</td>
<td><code>outputs/final-scores-compiled/overall-vulnerability-rankings/</code></td>
<td>04</td>
<td><code>3-calculate-overall-vulnerability-scores.R</code></td>
<td>1 per stock</td>
<td><code>N-Exp-Fact, Exposure-Factors, stock_name</code></td>
<td>—</td>
<td>QA</td>
<td>wide</td>
</tr>
<tr>
<td><code>n_reviews_directional_effect.csv</code></td>
<td>Number of reviewer workbooks with directional effect data per stock</td>
<td><code>outputs/final-scores-compiled/directional-effect/</code></td>
<td>05</td>
<td><code>1-extract-directional-effect.R</code></td>
<td>1 per stock</td>
<td><code>stock_name, n_reviews, reviewers</code></td>
<td>—</td>
<td>QA</td>
<td>wide</td>
</tr>
<tr>
<td><code>table_directional_effect.csv</code></td>
<td>Raw directional effect (Positive/Neutral/Negative) tallies extracted from each reviewer's workbook</td>
<td><code>outputs/final-scores-compiled/directional-effect/</code></td>
<td>05</td>
<td><code>1-extract-directional-effect.R</code></td>
<td>1 per reviewer × stock</td>
<td><code>SourceFile, Scorer, stock_name, row_idx, Directional_effect, Directional_score</code></td>
<td><code>05/2-summarize-directional-effect.R</code></td>
<td>intermediate</td>
<td>long</td>
</tr>
<tr>
<td><code>directional_effect_wide_all.csv</code></td>
<td>Reviewer-level directional effect tally counts and weighted average, pivoted wide</td>
<td><code>outputs/final-scores-compiled/directional-effect/</code></td>
<td>05</td>
<td><code>2-summarize-directional-effect.R</code></td>
<td>1 per reviewer × stock</td>
<td><code>Scorer, stock_name, Positive, Neutral, Negative, n_tallies, wt_avg, overall</code></td>
<td>—</td>
<td>intermediate</td>
<td>wide</td>
</tr>
<tr>
<td><code>directional_effect_summary_by-stock.csv</code></td>
<td>Stock-level directional effect summary: pooled tally counts, weighted average, and final classification</td>
<td><code>outputs/final-scores-compiled/directional-effect/</code></td>
<td>05</td>
<td><code>2-summarize-directional-effect.R</code></td>
<td>1 per stock</td>
<td><code>stock_name, n_positive, n_neutral, n_negative, n_tallies, n_reviewers, wt_avg, overall</code></td>
<td><code>10/5-produce-results-tables.R</code></td>
<td>final</td>
<td>wide</td>
</tr>
<tr>
<td><code>qa_data_quality_by_scorers.csv</code></td>
<td>Number of reviewer workbooks with data quality data per stock</td>
<td><code>outputs/final-scores-compiled/data-quality/</code></td>
<td>06</td>
<td><code>1-extract-data-quality-scores.R</code></td>
<td>1 per stock</td>
<td><code>stock_name, n_reviews, reviewers</code></td>
<td>—</td>
<td>QA</td>
<td>wide</td>
</tr>
<tr>
<td><code>qa_data_quality_by_stock.csv</code></td>
<td>Row-count check for extracted data quality scores per stock</td>
<td><code>outputs/final-scores-compiled/data-quality/</code></td>
<td>06</td>
<td><code>1-extract-data-quality-scores.R</code></td>
<td>1 per stock</td>
<td><code>stock_name, n_scorers, n_rows</code></td>
<td>—</td>
<td>QA</td>
<td>wide</td>
</tr>
<tr>
<td><code>table_data_quality_scores_extracted.csv</code></td>
<td>Raw data quality scores (0–3 scale) extracted from each reviewer's workbook for every stock × attribute</td>
<td><code>outputs/final-scores-compiled/data-quality/</code></td>
<td>06</td>
<td><code>1-extract-data-quality-scores.R</code></td>
<td>1 per reviewer × stock × attribute</td>
<td><code>SourceFile, Scorer, stock_name, row_idx, Attribute_type, Attribute_name, Data_quality_score, score_cell, label_cell, Data_quality_label, score_ge_2</code></td>
<td><code>06/2-summarize-data-quality-scores.R</code></td>
<td>intermediate</td>
<td>long</td>
</tr>
<tr>
<td><code>overall_data_quality_summary_by_stock.csv</code></td>
<td>Stock-level data quality summary: count of each quality tier, proportion ≥ 2, mean score, and final rank</td>
<td><code>outputs/final-scores-compiled/data-quality/</code></td>
<td>06</td>
<td><code>2-summarize-data-quality-scores.R</code></td>
<td>1 per stock</td>
<td><code>stock_name, n_3, n_2, n_1, n_0, n_ge_2, prop_ge_2, mean_score, sd_score, data_quality_rank</code></td>
<td><code>10/5-produce-results-tables.R</code></td>
<td>final</td>
<td>wide</td>
</tr>
<tr>
<td><code>table_sensitivity_tallies_long.csv</code></td>
<td>Reviewer-level tally counts (L/M/H/VH) for all sensitivity and rigidity attributes, extracted directly from workbooks</td>
<td><code>outputs/analyses/1-inputs/</code></td>
<td>07</td>
<td><code>1-extract-tally-scores.R</code></td>
<td>1 per reviewer × stock × sensitivity attribute</td>
<td><code>source_file, reviewer_id, stock_name, sheet_name, attribute_type, attribute_name, row_num, final_score, data_quality_index, tally_L, tally_M, tally_H, tally_VH, n_tallies</code></td>
<td><code>07/2-finalize-tally-tables.R</code></td>
<td>intermediate</td>
<td>long</td>
</tr>
<tr>
<td><code>table_qualitative_exposure_tallies_long.csv</code></td>
<td>Reviewer-level tally counts for all qualitative exposure attributes, extracted directly from workbooks</td>
<td><code>outputs/analyses/1-inputs/</code></td>
<td>07</td>
<td><code>1-extract-tally-scores.R</code></td>
<td>1 per reviewer × stock × exposure attribute</td>
<td><code>source_file, reviewer_id, stock_name, sheet_name, attribute_type, attribute_name, row_num, final_score, data_quality_index, tally_L, tally_M, tally_H, tally_VH, n_tallies</code></td>
<td><code>07/2-finalize-tally-tables.R</code></td>
<td>intermediate</td>
<td>long</td>
</tr>
<tr>
<td><code>table_directional_effect_tallies_long.csv</code></td>
<td>Reviewer-level directional effect tally counts, extracted directly from workbooks</td>
<td><code>outputs/analyses/1-inputs/</code></td>
<td>07</td>
<td><code>1-extract-tally-scores.R</code></td>
<td>1 per reviewer × stock × effect_category</td>
<td><code>source_file, reviewer_id, stock_name, sheet_name, effect_category, row_num, tally, n_tallies</code></td>
<td><code>07/2-finalize-tally-tables.R</code></td>
<td>intermediate</td>
<td>long</td>
</tr>
<tr>
<td><code>table_all_qualitative_tallies_long.csv</code></td>
<td>Sensitivity and qualitative exposure tallies combined into a single long table</td>
<td><code>outputs/analyses/1-inputs/</code></td>
<td>07</td>
<td><code>1-extract-tally-scores.R</code></td>
<td>1 per reviewer × stock × attribute (sensitivity + qual. exposure)</td>
<td><code>source_file, reviewer_id, stock_name, sheet_name, attribute_type, attribute_name, row_num, final_score, data_quality_index, tally_L, tally_M, tally_H, tally_VH, n_tallies</code></td>
<td>—</td>
<td>intermediate</td>
<td>long</td>
</tr>
<tr>
<td><code>qa_tally_row_checks.csv</code></td>
<td>Row-level tally integrity checks: flags rows where tallies don't sum to 5 or contain NAs</td>
<td><code>outputs/analyses/1-inputs/qa/</code></td>
<td>07</td>
<td><code>1-extract-tally-scores.R</code></td>
<td>1 per reviewer × stock × attribute</td>
<td><code>stock_name, reviewer_id, attribute_name, row_num, attribute_type, final_score, data_quality_index, tally_L, tally_M, tally_H, tally_VH, n_tallies, n_tallies_ok, has_na_tally, any_row_issue</code></td>
<td>—</td>
<td>QA</td>
<td>long</td>
</tr>
<tr>
<td><code>qa_reviewer_stock_coverage.csv</code></td>
<td>Reviewer × stock coverage matrix showing attribute count scored per reviewer per stock</td>
<td><code>outputs/analyses/1-inputs/qa/</code></td>
<td>07</td>
<td><code>1-extract-tally-scores.R</code></td>
<td>1 per stock (reviewers as columns)</td>
<td><code>stock_name, [reviewer_id columns]</code></td>
<td>—</td>
<td>QA</td>
<td>wide</td>
</tr>
<tr>
<td><code>qa_attribute_tally_summary.csv</code></td>
<td>Pooled tally sums per stock × attribute compared against expected values; flags mismatches</td>
<td><code>outputs/analyses/1-inputs/qa/</code></td>
<td>07</td>
<td><code>1-extract-tally-scores.R</code></td>
<td>1 per stock × attribute</td>
<td><code>stock_name, attribute_type, attribute_name, n_reviewers, tally_L, tally_M, tally_H, tally_VH, pooled_tally_sum, expected_tally_sum, tally_sum_ok, pooled_mean_score</code></td>
<td>—</td>
<td>QA</td>
<td>wide</td>
</tr>
<tr>
<td><code>qa_directional_effect_summary.csv</code></td>
<td>Directional effect tally count check per stock</td>
<td><code>outputs/analyses/1-inputs/qa/</code></td>
<td>07</td>
<td><code>1-extract-tally-scores.R</code></td>
<td>1 per stock</td>
<td><code>stock_name, n_reviewers, total_tally</code></td>
<td>—</td>
<td>QA</td>
<td>wide</td>
</tr>
<tr>
<td><code>sensitivity_tallies_long.csv</code></td>
<td>Finalized reviewer-level sensitivity tally table (<code>source_file</code> and <code>sheet_name</code> dropped; analysis-ready)</td>
<td><code>outputs/final-tallies-long/</code></td>
<td>07</td>
<td><code>2-finalize-tally-tables.R</code></td>
<td>1 per reviewer × stock × sensitivity attribute</td>
<td><code>reviewer_id, stock_name, attribute_type, attribute_name, row_num, final_score, data_quality_index, tally_L, tally_M, tally_H, tally_VH, n_tallies</code></td>
<td><code>08/uncertainty-analyses.R</code>, <code>09/2-bootstrap-distributional-change.R</code></td>
<td>final</td>
<td>long</td>
</tr>
<tr>
<td><code>sensitivity_tallies_by_stock.csv</code></td>
<td>Sensitivity tallies pooled across reviewers per stock × attribute, with reviewer count</td>
<td><code>outputs/final-tallies-long/</code></td>
<td>07</td>
<td><code>2-finalize-tally-tables.R</code></td>
<td>1 per stock × sensitivity attribute</td>
<td><code>stock_name, attribute_name, tally_L, tally_M, tally_H, tally_VH, n_scorers, n_tallies</code></td>
<td><code>10/1-plot-scoring-distributions.R</code></td>
<td>final</td>
<td>wide</td>
</tr>
<tr>
<td><code>directional_effect_tallies_long.csv</code></td>
<td>Finalized reviewer-level directional effect tally table</td>
<td><code>outputs/final-tallies-long/</code></td>
<td>07</td>
<td><code>2-finalize-tally-tables.R</code></td>
<td>1 per reviewer × stock × effect_category</td>
<td><code>reviewer_id, stock_name, effect_category, tally</code></td>
<td><code>08/uncertainty-analyses.R</code></td>
<td>final</td>
<td>long</td>
</tr>
<tr>
<td><code>directional_effect_tallies_by_stock.csv</code></td>
<td>Directional effect tallies pooled across reviewers per stock</td>
<td><code>outputs/final-tallies-long/</code></td>
<td>07</td>
<td><code>2-finalize-tally-tables.R</code></td>
<td>1 per stock</td>
<td><code>stock_name, tally_neg, tally_neut, tally_pos, n_scorers, n_tallies</code></td>
<td><code>10/1-plot-scoring-distributions.R</code></td>
<td>final</td>
<td>wide</td>
</tr>
<tr>
<td><code>exposure_tallies_long.csv</code></td>
<td>Finalized reviewer-level exposure tally table combining qualitative and quantitative exposure attributes</td>
<td><code>outputs/final-tallies-long/</code></td>
<td>07</td>
<td><code>2-finalize-tally-tables.R</code></td>
<td>1 per reviewer × stock × exposure attribute</td>
<td><code>reviewer_id, stock_name, attribute_type, attribute_name, tally_L, tally_M, tally_H, tally_VH, n_tallies</code></td>
<td><code>08/uncertainty-analyses.R</code>, <code>10/1-plot-scoring-distributions.R</code></td>
<td>final</td>
<td>long</td>
</tr>
<tr>
<td><code>exposure_tallies_by_stock.csv</code></td>
<td>Exposure tallies pooled across reviewers per stock × attribute</td>
<td><code>outputs/final-tallies-long/</code></td>
<td>07</td>
<td><code>2-finalize-tally-tables.R</code></td>
<td>1 per stock × exposure attribute</td>
<td><code>stock_name, attribute_type, attribute_name, tally_L, tally_M, tally_H, tally_VH, n_tallies</code></td>
<td><code>10/1-plot-scoring-distributions.R</code></td>
<td>final</td>
<td>wide</td>
</tr>
<tr>
<td><code>01_preanalysis_check_log.csv</code></td>
<td>Pre-analysis validation log for the uncertainty module: named checks with pass/fail status and detail</td>
<td><code>outputs/analyses/uncertainty-loo/intermediate/</code></td>
<td>08</td>
<td><code>uncertainty-analyses.R</code></td>
<td>1 per check</td>
<td><code>check_name, status, detail</code></td>
<td>—</td>
<td>QA</td>
<td>wide</td>
</tr>
<tr>
<td><code>02_baseline_reproduced_scores.csv</code></td>
<td>Reproduced baseline vulnerability scores computed from tally tables; verifies pipeline reproducibility before bootstrap runs</td>
<td><code>outputs/analyses/uncertainty-loo/intermediate/</code></td>
<td>08</td>
<td><code>uncertainty-analyses.R</code></td>
<td>1 per stock</td>
<td><code>stock_name, sensitivity_rank_repro, sensitivity_score_numeric_repro, exposure_rank_repro, exposure_score_numeric_repro, vulnerability_score_numeric_repro, vulnerability_rank_repro</code></td>
<td>—</td>
<td>QA</td>
<td>wide</td>
</tr>
<tr>
<td><code>table_bootstrap_uncertainty_stock.csv</code></td>
<td>Bootstrap rank distribution for overall vulnerability; one row per stock × possible rank</td>
<td><code>outputs/analyses/uncertainty-loo/final-tables/</code></td>
<td>08</td>
<td><code>uncertainty-analyses.R</code></td>
<td>1 per stock × vuln_rank</td>
<td><code>stock_name, vuln_rank, n, prop, borderline</code></td>
<td><code>10/2-plot-uncertainty-figures.R</code></td>
<td>final</td>
<td>long</td>
</tr>
<tr>
<td><code>table_bootstrap_uncertainty_sensitivity_component.csv</code></td>
<td>Bootstrap rank distribution for the sensitivity component only</td>
<td><code>outputs/analyses/uncertainty-loo/final-tables/</code></td>
<td>08</td>
<td><code>uncertainty-analyses.R</code></td>
<td>1 per stock × sens_rank</td>
<td><code>stock_name, sens_rank, n, prop</code></td>
<td><code>10/2-plot-uncertainty-figures.R</code></td>
<td>final</td>
<td>long</td>
</tr>
<tr>
<td><code>table_directional_effect_bootstrap.csv</code></td>
<td>Bootstrap rank distribution for directional effect, including baseline rank and weighted mean</td>
<td><code>outputs/analyses/uncertainty-loo/final-tables/</code></td>
<td>08</td>
<td><code>uncertainty-analyses.R</code></td>
<td>1 per stock × dir_rank</td>
<td><code>stock_name, dir_rank, n, prop, baseline_dir_rank, baseline_w_mean</code></td>
<td><code>10/5-produce-results-tables.R</code></td>
<td>final</td>
<td>long</td>
</tr>
<tr>
<td><code>table_leave_one_out_sensitivity_long.csv</code></td>
<td>LOO influence results for sensitivity attributes: baseline vs. new ranks when each attribute is omitted, per stock</td>
<td><code>outputs/analyses/uncertainty-loo/final-tables/</code></td>
<td>08</td>
<td><code>uncertainty-analyses.R</code></td>
<td>1 per stock × attribute_omitted</td>
<td><code>stock_name, attribute_omitted, baseline_sens_rank, new_sens_rank, baseline_vuln_rank, new_vuln_rank, baseline_vuln_score, new_vuln_score, rank_changed, rank_change_direction</code></td>
<td><code>10/2-plot-uncertainty-figures.R</code></td>
<td>final</td>
<td>long</td>
</tr>
<tr>
<td><code>table_leave_one_out_sensitivity_summary.csv</code></td>
<td>LOO sensitivity influence summary: how many stocks change rank when each attribute is omitted</td>
<td><code>outputs/analyses/uncertainty-loo/final-tables/</code></td>
<td>08</td>
<td><code>uncertainty-analyses.R</code></td>
<td>1 per attribute_omitted</td>
<td><code>attribute_omitted, n_stocks_tested, n_rank_changed, n_rank_lower, n_rank_higher, prop_rank_changed</code></td>
<td><code>10/2-plot-uncertainty-figures.R</code></td>
<td>final</td>
<td>wide</td>
</tr>
<tr>
<td><code>table_leave_one_out_exposure_long.csv</code></td>
<td>LOO influence results for exposure factors: baseline vs. new ranks when each factor is omitted, per stock</td>
<td><code>outputs/analyses/uncertainty-loo/final-tables/</code></td>
<td>08</td>
<td><code>uncertainty-analyses.R</code></td>
<td>1 per stock × factor_omitted</td>
<td><code>stock_name, factor_omitted, baseline_exp_rank, new_exp_rank, baseline_vuln_rank, new_vuln_rank, baseline_vuln_score, new_vuln_score, rank_changed, rank_change_direction</code></td>
<td><code>10/2-plot-uncertainty-figures.R</code></td>
<td>final</td>
<td>long</td>
</tr>
<tr>
<td><code>table_leave_one_out_exposure_summary.csv</code></td>
<td>LOO exposure influence summary: how many stocks change rank when each exposure factor is omitted</td>
<td><code>outputs/analyses/uncertainty-loo/final-tables/</code></td>
<td>08</td>
<td><code>uncertainty-analyses.R</code></td>
<td>1 per factor_omitted</td>
<td><code>factor_omitted, n_stocks_tested, n_rank_changed, n_rank_lower, n_rank_higher, prop_rank_changed</code></td>
<td><code>10/2-plot-uncertainty-figures.R</code></td>
<td>final</td>
<td>wide</td>
</tr>
<tr>
<td><code>table_exposure_factor_scores_for_plot.csv</code></td>
<td>Mean exposure factor scores with component rank, pre-formatted for figure generation</td>
<td><code>outputs/analyses/uncertainty-loo/final-tables/</code></td>
<td>08</td>
<td><code>uncertainty-analyses.R</code></td>
<td>1 per stock × exposure_factor</td>
<td><code>stock_name, exposure_factor, mean_score, component_rank_baseline</code></td>
<td><code>10/2-plot-uncertainty-figures.R</code></td>
<td>figure-input</td>
<td>long</td>
</tr>
<tr>
<td><code>table_sensitivity_attribute_scores_for_plot.csv</code></td>
<td>Mean sensitivity attribute scores with component rank, pre-formatted for figure generation</td>
<td><code>outputs/analyses/uncertainty-loo/final-tables/</code></td>
<td>08</td>
<td><code>uncertainty-analyses.R</code></td>
<td>1 per stock × attribute</td>
<td><code>stock_name, attribute_name, mean_score, component_rank_baseline</code></td>
<td><code>10/2-plot-uncertainty-figures.R</code></td>
<td>figure-input</td>
<td>long</td>
</tr>
<tr>
<td><code>table_bootstrap_final_summary.csv</code></td>
<td>Publication-ready bootstrap summary: dominant rank and rank proportion columns for all vulnerability components</td>
<td><code>outputs/analyses/uncertainty-loo/final-tables/</code></td>
<td>08</td>
<td><code>uncertainty-analyses.R</code></td>
<td>1 per stock</td>
<td><code>Stock, N Exp factors, Exp, Sens, Vul, L, M, H, VH</code></td>
<td><code>10/5-produce-results-tables.R</code></td>
<td>final</td>
<td>wide</td>
</tr>
<tr>
<td><code>distributional_change_potential_uscar.csv</code></td>
<td>Baseline distributional change potential (DCP) scores and ranks derived from 4 BSAs with mobility attributes inverted</td>
<td><code>outputs/distribution-change-potential/</code></td>
<td>09</td>
<td><code>1-calculate-distributional-change-potential.R</code></td>
<td>1 per stock</td>
<td><code>stock_name, adult_mobility_inverted, habitat_specificity_inverted, early_life_dispersal_inverted, species_range, dcp_rank, dcp_numeric</code></td>
<td><code>09/2-bootstrap-distributional-change.R</code>, <code>10/3-plot-distributional-change.R</code></td>
<td>final</td>
<td>wide</td>
</tr>
<tr>
<td><code>distributional_change_bootstrap_uscar.csv</code></td>
<td>Bootstrap rank distribution for DCP scores</td>
<td><code>outputs/distribution-change-potential/</code></td>
<td>09</td>
<td><code>2-bootstrap-distributional-change.R</code></td>
<td>1 per stock</td>
<td><code>stock_name, prop_L, prop_M, prop_H, prop_VH, dominant_rank, dominant_prop, borderline</code></td>
<td><code>10/3-plot-distributional-change.R</code></td>
<td>final</td>
<td>wide</td>
</tr>
<tr>
<td><code>distributional_change_full_uscar.csv</code></td>
<td>Combined baseline DCP scores and bootstrap uncertainty results in a single table</td>
<td><code>outputs/distribution-change-potential/</code></td>
<td>09</td>
<td><code>2-bootstrap-distributional-change.R</code></td>
<td>1 per stock</td>
<td><code>stock_name, adult_mobility_inverted, habitat_specificity_inverted, early_life_dispersal_inverted, species_range, dcp_rank, dcp_numeric, prop_L, prop_M, prop_H, prop_VH, dominant_rank, dominant_prop, borderline</code></td>
<td><code>10/3-plot-distributional-change.R</code></td>
<td>final</td>
<td>wide</td>
</tr>
<tr>
<td><code>table_directional_effect_results.csv</code></td>
<td>Publication results table: directional effect category, weighted mean, tally counts, and bootstrap proportions by stock</td>
<td><code>outputs/tables/</code></td>
<td>10</td>
<td><code>5-produce-results-tables.R</code></td>
<td>1 per stock</td>
<td><code>Stock, Vuln_rank, Dir_effect, Wt_mean, N_negative, N_neutral, N_positive, N_tallies, Boot_Negative, Boot_Neutral, Boot_Positive, Dominant_prop, Borderline</code></td>
<td>—</td>
<td>final</td>
<td>wide</td>
</tr>
<tr>
<td><code>table_data_quality_results.csv</code></td>
<td>Publication results table: data quality rank, proportion ≥ 2, mean score, and tier counts by stock</td>
<td><code>outputs/tables/</code></td>
<td>10</td>
<td><code>5-produce-results-tables.R</code></td>
<td>1 per stock</td>
<td><code>Stock, Vuln_rank, Data_quality_rank, Prop_ge_2, Mean_score, N_adequate, N_limited, N_expert, N_nodata</code></td>
<td>—</td>
<td>final</td>
<td>wide</td>
</tr>
</tbody>
</table>
