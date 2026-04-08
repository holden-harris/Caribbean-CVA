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

# `1-extract-final-scores-from-all-reviewers.R`

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

from pathlib import Path

# `2-summarize-attribute-scores.R`

Overall climate vulnerability was calculated using the NOAA Fisheries Climate Vulnerability Assessment (FCVA) framework, in which sensitivity and exposure are first summarized separately and then combined into a final vulnerability score and rank. The original NOAA FCVA methodology uses three steps: (1) calculate mean scores for individual sensitivity attributes and exposure factors, (2) assign an overall sensitivity score and an overall exposure score using a logic model based on the number of attribute or factor means above threshold values, and (3) multiply the two component scores to assign a final vulnerability rank. This framework was established by Morrison et al. and has been retained in more recent regional assessments, including the South Atlantic, California marine fisheries, and highly migratory species CVAs.

For the Caribbean CVA, a compiled scoring table was first assembled for each stock. This table contained both reviewer-derived qualitative scores and calculated quantitative exposure scores. Each row in the compiled table included the stock name, attribute type, score type, attribute name, scorer, score, and region. Attribute type was standardized to either **Sensitivity** or **Exposure**, and score type was standardized to either **Qualitative** or **Calculated**.

## Workflow summary

1. Compile qualitative reviewer scores and calculated exposure scores into a single standardized table.
2. Calculate one mean score for each stock × attribute combination.
3. Apply the NOAA FCVA logic model separately to the sensitivity and exposure attribute means.
4. Convert the resulting component categories to numeric scores.
5. Multiply the exposure and sensitivity scores to calculate the final vulnerability score.
6. Assign the final vulnerability rank using the combined-score thresholds.

## Suggested outputs

The workflow produces three main outputs:

- **Attribute means table**: one row per stock × attribute, containing the mean score, standard deviation, number of non-missing scores, and optional LMHV bin.
- **Component scores table**: one row per stock × component, containing the number of attributes exceeding each threshold and the resulting exposure or sensitivity score and rank.
- **Overall vulnerability table**: one row per stock, containing exposure score, sensitivity score, combined vulnerability score, and final vulnerability rank.

## Step 1. Calculate mean scores for individual attributes and factors

For each stock, an average score was calculated for every individual sensitivity attribute and exposure factor. In the original NOAA FCVA framework, attribute and factor means are calculated as weighted averages of expert tallies distributed across four ordinal bins: Low = 1, Moderate = 2, High = 3, and Very High = 4. The weighted average is:

\\[
\\text{Attribute or Factor Mean} =
\\frac{(L \\times 1) + (M \\times 2) + (H \\times 3) + (VH \\times 4)}
{L + M + H + VH}
\\]

where *L*, *M*, *H*, and *VH* are the number of tallies assigned to the four scoring bins.

In the present Caribbean workflow, the compiled qualitative table already contained final reviewer-entered scores on the same 1 to 4 ordinal scale, rather than raw tallies. Therefore, for each stock × sensitivity attribute combination, the attribute mean was calculated as the arithmetic mean of all non-missing reviewer scores. For calculated exposure factors, each stock × factor combination in the U.S. Caribbean had a single derived value from the quantitative exposure analysis, so the factor mean equaled that score. In total, 14 sensitivity attributes and 13 quantitative exposure factors were used in the final component scoring.

## Step 2. Assign overall sensitivity and exposure component scores using the logic model

After calculating individual attribute and factor means, the NOAA FCVA logic model was applied separately to sensitivity and exposure. Under this framework, the overall component score depends on how many attribute or factor means exceed specified thresholds, rather than on the average of all attribute means. Morrison et al. adopted this approach because simple averaging can downweight a small number of very high-risk attributes or factors; the logic model, in contrast, highlights cases in which a stock has multiple sources of elevated sensitivity or exposure.

For each stock, the number of sensitivity attribute means and exposure factor means that were greater than or equal to 2.5, 3.0, and 3.5 was counted. These counts were then translated into component scores and categories using the NOAA FCVA logic rules shown in Table 1.

### Table 1. Logic model used to assign overall sensitivity and exposure component scores

| Overall sensitivity or exposure category | Numeric score | Logic rule |
|---|---:|---|
| Very High | 4 | 3 or more attribute or factor means ≥ 3.5 |
| High | 3 | 2 or more attribute or factor means ≥ 3.0 |
| Moderate | 2 | 2 or more attribute or factor means ≥ 2.5 |
| Low | 1 | All other cases |

Recent regional CVAs have described these rules as “3 or more” attributes or factors with a mean ≥ 3.5 for a Very High component score, “2 or more” with mean ≥ 3.0 for High, and “2 or more” with mean ≥ 2.5 for Moderate; all other cases are Low.

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
