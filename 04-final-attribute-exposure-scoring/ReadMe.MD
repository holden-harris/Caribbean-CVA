# 4-final-attribute-exposure-scoring

## Overview

This folder compiles and summarizes final attribute scores from individually scored expert review workbooks for the Caribbean Climate Vulnerability Assessment (CVA). The workflow first extracts final reviewer scores into a single compiled table, then summarizes those scores across attribute types and joins them with calculated quantitative exposure scores for downstream vulnerability scoring.

# `1-extract-final-scores-from-all-reviewers.R`

This script compiles final reviewer-entered CVA attribute scores from all completed scoring workbooks into a single standardized, long-format table for downstream analysis.

#### Purpose

The script reads the final reviewer workbooks stored in:

`./data/final-scores/`

and extracts the final attribute scores recorded for each stock and attribute. These outputs are combined into a single analysis-ready table used by the downstream summary script.

### Input Directory Structure
```text
data/
  └── final-scores/              ## Reviewer workbooks (*.xlsx)
outputs/
  ├── final-attribute-scores/    ## Compiled and summary tables
  └── overall-vulnerability-rankings
4-final-attribute-exposure-scoring/
├── 1-extract-final-scores-from-all-reviewers.R
├── 2-summarize-attribute-scores.R
├── 3-calculate-overall-vulnerability-scores
└── README.md
```

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
   - scorer  
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
- Preserve raw reviewer workbooks in `./data/final-scores/`
- Avoid editing compiled CSV outputs by hand
- Rerun the extraction step whenever reviewer workbooks are added or revised

## Outputs

| File | Description |
|------|-------------|
| `table_final_attribute_scores_all.csv` | Long-format compiled table of final attribute scores from all reviewer workbooks |
| `quantitative-exposure-attribute-scores-all.csv` | Quantitative exposure attribute scores used downstream in the summary step |
| Additional summary outputs | Created by `2-summarize-attribute-scores.R` in the same output directory |

```
outputs/
└── final-scores-compiled/
    └── final-attribute-scores/ 
        ├── table_final_attribute_scores_all.csv
        └── quantitative-exposure-attribute-scores-all.csv
```

# `2-summarize-attribute-scores.R`

### Coding Workflow 

This script compiles reviewer-based qualitative CVA scores and calculates quantitative exposure scores into a standardized long-format table, filters the dataset to the U.S. Caribbean, performs quality-assurance checks on score completeness, and then applies the NOAA FCVA logic model to calculate stock-level exposure, sensitivity, and overall vulnerability rankings.

The workflow has the following steps:

1. **Read and standardize qualitative reviewer scores**
   - Load the compiled table of final reviewer-entered attribute scores.
   - Recode `Rigidity` as `Sensitivity`.
   - Recode `Qualitative Exposure Factors` as `Exposure`.
   - Standardize fields to a common structure:
     - `stock_name`
     - `region` (i.e., U.S. Caribbean, Wider Caribbean, and Western Atlantic)
     - `attribute_type` (i.e., exposure factor or biological sensitivity attribute)
     - `score_type` (i.e., qualitative or calculated)
     - `attribute_name`
     - `scorer`
     - `score` (i.e., specific value 1-4)

2. **Read and standardize calculated exposure scores**
   - Load the quantitative exposure attribute scores generated from the exposure-anomalies workflow.
   - Harmonize stock names to match the qualitative table.
   - Assign region labels from the exposure analysis extent:
     - `U.S. Caribbean`
     - `Wider Caribbean`
     - `Western Atlantic`
   - Standardize fields to the same long-format structure as the qualitative table.
   - Assign:
     - `attribute_type = "Exposure"`
     - `score_type = "Calculated"`
     - `scorer = "Calculated"`

3. **Combine qualitative and calculated scores**
   - Append the qualitative and calculated tables into one compiled dataset.
   - Sort rows by stock, attribute type, score type, region, attribute name, and scorer.

4. **Filter to the U.S. Caribbean**
   - Remove rows where `region` is `Wider Caribbean` or `Western Atlantic`.
   - Retain:
     - all qualitative scores (`region = NA`)
     - calculated exposure scores for `U.S. Caribbean`
   - Write the U.S. Caribbean subset to disk.

5. **Run QA checks on score completeness**
   - Count the number of rows per `stock_name × attribute_name × score_type`.
   - Confirm that:
     - calculated scores have exactly **1** score per stock × attribute
     - qualitative scores have exactly **4** scores per stock × attribute
   - Identify any problem rows and list the associated scorers.

6. **Write the full compiled tables to disk.**
   - `table_final_attribute_exposure_uscar_all.csv`
   - `table_final_attribute_exposure_uscar_call.csv`

| File | Description |
|------|-------------|
| `final_scores_all.csv` | Long-format compiled table combining qualitative reviewer scores and calculated quantitative exposure scores |
| `final_scores_uscar.csv` | U.S. Caribbean subset used in the final FCVA workflow, i.e., subset without exposure factors calculated for Wider Caribbean or Western Atlantic |

# 3-calculate-overall-vulnerability-scores.R

Overall climate vulnerability was calculated using the NOAA Fisheries Climate Vulnerability Assessment (FCVA) framework, in which sensitivity and exposure are first summarized separately and then combined into a final vulnerability score and rank. The original NOAA FCVA methodology uses three steps: (1) calculate mean scores for individual sensitivity attributes and exposure factors, (2) assign an overall sensitivity score and an overall exposure score using a logic model based on the number of attribute or factor means above threshold values, and (3) multiply the two component scores to assign a final vulnerability rank. This framework was established by Morrison et al. and has been retained in more recent regional assessments, including the South Atlantic, California marine fisheries, and highly migratory species CVAs.

For the Caribbean CVA, a compiled scoring table was first assembled for each stock. This table contained both reviewer-derived qualitative scores and calculated quantitative exposure scores. Each row in the compiled table included the stock name, attribute type, score type, attribute name, scorer, score, and region. Attribute type was standardized to either **Sensitivity** or **Exposure**, and score type was standardized to either **Qualitative** or **Calculated**.

## Step 1. Calculate mean scores for individual attributes and factors

For each stock, an average score was calculated for every individual sensitivity attribute and exposure factor. In the original NOAA FCVA framework, attribute and factor means are calculated as weighted averages of expert tallies distributed across four ordinal bins: Low = 1, Moderate = 2, High = 3, and Very High = 4. The weighted average is:

$$
\text{Attribute or Factor Mean} =
\frac{(L \times 1) + (M \times 2) + (H \times 3) + (VH \times 4)}
{L + M + H + VH}
$$

where *L*, *M*, *H*, and *VH* are the number of tallies assigned to the four scoring bins.

In the present Caribbean workflow, the compiled qualitative table already contained final reviewer-entered scores on the same 1 to 4 ordinal scale, rather than raw tallies. Therefore, for each stock × sensitivity attribute combination, the attribute mean was calculated as the arithmetic mean of all non-missing reviewer scores. For calculated exposure factors, each stock × factor combination in the U.S. Caribbean had a single derived value from the quantitative exposure analysis, so the factor mean equaled that score. In total, 14 sensitivity attributes were used in the sensitivity component scoring for all stocks. The number of exposure factors used in the exposure component scoring varies by stock and is determined by an expert rubric (`data/attribute-list-rubric-completed.csv`), with a maximum of 15 factors (13 quantitative + 2 qualitative); each stock uses only the subset identified as relevant by the CVA reviewers.

## Step 2. Assign overall sensitivity and exposure component scores using the logic model

After calculating individual attribute and factor means, the NOAA FCVA logic model was applied separately to sensitivity and exposure. Under this framework, the overall component score depends on how many attribute or factor means exceed specified thresholds, rather than on the average of all attribute means. Morrison et al. adopted this approach because simple averaging can downweight a small number of very high-risk attributes or factors; the logic model, in contrast, highlights cases in which a stock has multiple sources of elevated sensitivity or exposure.

For each stock, the number of sensitivity attribute means and exposure factor means that were greater than or equal to 2.5, 3.0, and 3.5 was counted. These counts were then translated into component scores and categories using the NOAA FCVA logic rules shown in Table 1.

### Table 1. Logic model used to assign overall sensitivity and exposure component scores

| Overall sensitivity or exposure category | Numeric score | Logic rule |
|---|---:|---|
| Very High | 4 | More than 3 attribute or factor means ≥ 3.5 |
| High | 3 | 2 or more attribute or factor means ≥ 3.0 |
| Moderate | 2 | 2 or more attribute or factor means ≥ 2.5 |
| Low | 1 | All other cases |

Recent regional CVAs have described these rules as “more than 3” attributes or factors with a mean ≥ 3.5 for a Very High component score, “more than 2” with mean ≥ 3.0 for High, and “more than 2” with mean ≥ 2.5 for Moderate; all other cases are Low. These thresholds correspond to `rank_threshold = attr_means_current = 2` in the R script, consistent with the original NOAA FCVA methodology.

## Step 3. Calculate the overall vulnerability score and rank

The final climate vulnerability score for each stock was calculated as the product of the exposure and sensitivity component scores. Component categories Low, Moderate, High, and Very High were assigned numeric values of 1, 2, 3, and 4, respectively. The product of these two component scores was then converted to the final overall vulnerability rank using the standard NOAA FCVA thresholds shown in Table 2. Specifically, products of 1–3 were classified as Low, 4–6 as Moderate, 8–9 as High, and 12–16 as Very High. This product-based classification follows the original NOAA FCVA methodology and has been used in subsequent U.S. regional CVAs.

### Table 2. Combined score thresholds used to assign the final overall vulnerability rank

| Exposure score × Sensitivity score | Overall vulnerability rank |
|---:|---|
| 1–3 | Low |
| 4–6 | Moderate |
| 8–9 | High |
| 12–16 | Very High |

## Interpretation

Under this framework, a stock must be both sensitive and exposed to receive a high overall vulnerability ranking. Thus, the final vulnerability rank generally reflects the lesser of the two component scores, although the multiplication step preserves the distinction between combinations such as Moderate × Very High and High × High. Morrison et al. note that this framework is intended to identify stocks with multiple, concurrent sources of risk rather than to estimate the absolute magnitude of climate impacts.

## Code Workflow

The workflow produces three main outputs:

- **Attribute means table**: one row per stock × attribute, containing the mean score, standard deviation, number of non-missing scores, and optional LMHV bin.
- **Component scores table**: one row per stock × component, containing the number of attributes exceeding each threshold and the resulting exposure or sensitivity score and rank.
- **Overall vulnerability table**: one row per stock, containing exposure score, sensitivity score, combined vulnerability score, and final vulnerability rank.

1. **Calculate attribute-level mean scores**
   - For each `stock_name × attribute_type × score_type × attribute_name`, calculate:
     - `attribute_mean`
     - `attribute_sd`
     - `n_scores`
   - Categorize each `attribute_mean` into an LMHV bin for QA:
     - `Low`
     - `Moderate`
     - `High`
     - `Very High`
   - Write the attribute means table to disk.

2. **Apply the NOAA FCVA logic model**
   - Restrict the logic-model step to:
     - qualitative `Sensitivity` attributes
     - calculated `Exposure` attributes
   - **Apply expert exposure factor filter**
     - Read `data/attribute-list-rubric-completed.csv` and pivot to long form
       (`stock_name`, `attribute_name`, `include`). Write the result to
       `data/exposure-factor-filter-long.csv` for use by downstream scripts.
     - Retain only exposure factor rows where `include == TRUE` in the rubric.
       Sensitivity attribute rows are unaffected — all 14 attributes apply to all stocks.
     - A name-alignment guard confirms all exposure factor names in the attribute means
       match the rubric before filtering proceeds; mismatches stop execution.
   - For each `stock_name × attribute_type`, calculate:
     - `mean_attribute_score`
     - `sd_attribute_score`
     - `n_attributes`
     - number of attribute means ≥ 2.5
     - number of attribute means ≥ 3.0
     - number of attribute means ≥ 3.5
   - Assign the component score and rank using the FCVA logic model:
     - **Very High (4)** = more than 3 attribute means ≥ 3.5
     - **High (3)** = more than 2 attribute means ≥ 3.0
     - **Moderate (2)** = more than 2 attribute means ≥ 2.5
     - **Low (1)** = all other cases
   - Write the component scores table to disk.

3. **Calculate overall vulnerability**
   - Pivot the component scores wider so each stock has one row.
   - Create stock-level columns for:
     - `Exp_score`
     - `Exp_rank`
     - `Sens_score`
     - `Sens_rank`
   - Calculate:
     - `Vuln_score = Exp_score × Sens_score`
   - Assign the final vulnerability rank:
     - `1–3 = Low`
     - `4–6 = Moderate`
     - `8–9 = High`
     - `12–16 = Very High`
   - Write the final stock-level vulnerability table to disk.

### Final outputs

| File | Description |
|------|-------------|
| `attribute_means_uscar.csv` | Attribute-level mean scores and LMHV bins for all attributes and factors (complete; unfiltered reference) |
| `component_scores_uscar.csv` | Exposure and Sensitivity component scores and ranks for each stock; `n_attributes` for the Exposure component now varies by stock reflecting the per-stock expert filter |
| `overall_vulnerability_scores_uscar.csv` | Final stock-level vulnerability scores and ranks |
| `data/exposure-factor-filter-long.csv` | Long-form expert rubric filter: one row per stock × exposure factor (375 rows = 15 factors × 25 stocks), with `include` (TRUE/FALSE) indicating expert-approved applicability. Also read by Script 8 (`uncertainty-analyses.R`). |

### Output directory structure
```text
outputs/
└── final-scores-compiled/
    └── final-attribute-scores/
        ├── table_final_attribute_scores_compiled.csv
        ├── final_scores_uscar.csv
        ├── attribute_means_uscar.csv
        ├── component_scores_uscar.csv
        └── overall_vulnerability_scores_uscar.csv
```

# `4-plot-overall-vulnerability.R`

Produces the two-panel publication figure summarizing overall climate vulnerability (Panel A) and dominant directional effect (Panel B) for all 25 Caribbean CVA stocks.

## Inputs

| File | Location | Description |
|------|----------|-------------|
| `overall_vulnerability_scores_uscar.csv` | `outputs/final-scores-compiled/overall-vulnerability-rankings/` | Stock-level exposure, sensitivity, and vulnerability scores and ranks |
| `table_bootstrap_uncertainty_stock.csv` | `outputs/analyses/uncertainty-loo/final-tables/` | Bootstrap certainty by stock and vulnerability rank (proportion of 10,000 iterations matching assigned rank) |
| `table_directional_effect_bootstrap.csv` | `outputs/analyses/uncertainty-loo/final-tables/` | Bootstrap certainty by stock and directional effect category |
| `directional_effect_tallies_by_stock.csv` | `outputs/final-tallies-long/` | Reviewer tally counts (negative, neutral, positive) per stock |

## Output

| File | Location | Description |
|------|----------|-------------|
| `fig_overall_vulnerability.png` | `figures/` | Two-panel figure (Panel A: vulnerability grid; Panel B: directional effect) |

## Panel A — Overall Vulnerability Grid

Panel A is a 4×4 colored tile matrix where the x-axis is **Climate Exposure rank** (Low → Very High) and the y-axis is **Biological Sensitivity rank** (Low → Very High). Each cell is filled with the overall vulnerability color derived from the FCVA scoring logic (Vuln_score = Exp_score × Sens_score; thresholds: ≤3 = Low, ≤6 = Moderate, ≤9 = High, >9 = Very High). Stock names are placed in the cell corresponding to their assigned Exposure and Sensitivity ranks. Within each occupied cell, stocks are stacked vertically in descending order of bootstrap certainty (most certain at top).

The Very High sensitivity row is rendered at reduced height because no stocks fall in that row.

## Panel B — Directional Effect

Panel B is a three-column colored tile table (Negative, Neutral, Positive) showing the dominant directional effect assigned to each stock by plurality of reviewer tallies. Stocks are stacked within each column in descending order of directional-effect bootstrap certainty.

## Font Encoding

Text face and color encode bootstrap certainty consistently across both panels:

| Certainty | Bootstrap proportion | Font face | Text color |
|-----------|---------------------|-----------|------------|
| Very High | > 95% | Bold | Black |
| High | 90–95% | Italic | Black |
| Moderate | 67–89% | Bold | White |
| Low | ≤ 66% | Italic | White |

Panel A certainty is drawn from `table_bootstrap_uncertainty_stock.csv` (vulnerability rank bootstrap). Panel B certainty is drawn from `table_directional_effect_bootstrap.csv` (directional effect bootstrap).

## Stock Name Recoding

Stock names in the vulnerability and bootstrap tables use sentence-case canonical names (e.g., "Atlantic thread herring"). The `stock_name_recode` lookup converts these to the Title Case display names used in the tally tables and figure labels (e.g., "Atlantic Herring"). Four names require non-trivial remapping: Atlantic Herring → Atlantic thread herring, Diadema → Long-spined sea urchin, Redhind → Red hind, Sea Cucumber → Sea cucumbers.
