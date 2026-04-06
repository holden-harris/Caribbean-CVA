# 4-final-directional-effect-scoring

## Overview
This folder contains scripts used to compile and summarize final **Directional Effect** scores from completed Caribbean CVA reviewer workbooks. The workflow extracts reviewer tallies for **Positive**, **Neutral**, and **Negative** directional effects from each species/stock worksheet, combines them into a single long-format table, and then summarizes those tallies at both the **reviewer × stock** level and the **stock-level across all reviewers**.

This folder supports transparency and reproducibility by documenting how final directional-effect tallies were pulled directly from the scoring workbooks and converted into analysis-ready summary tables.

## Final products
Scripts in this folder produce the following final products:

- `table_directional_effect_scores.csv`  
  Long-format compiled table of directional-effect tallies by reviewer and stock.

- `n_reviews_directional-effect.csv`  
  Quality-assurance table listing the number of reviewers contributing directional-effect scores for each stock.

- `directional_effect_wide_all.csv`  
  Reviewer-level wide summary table showing Positive, Neutral, and Negative tallies and the resulting weighted-average directional effect for each `stock_name × Scorer`.

- `directional_effect_summary_by-stock.csv`  
  Final stock-level summary table aggregating raw tallies across all reviewers and calculating one overall directional-effect category per stock.

## Input directory structure
The scripts in this folder expect inputs organized approximately as follows:

```text
project/
├─ data/
│  └─ final-scores/                    # Final Caribbean CVA reviewer Excel workbooks
│     ├─ [workbook_1].xlsx
│     ├─ [workbook_2].xlsx
│     └─ ...
├─ 4-final-directional-effect-scoring/
│  ├─ 1-extract-directional-effect.R
│  └─ 2-summarize-directional-effect.R
└─ outputs/
   └─ final-scores-compiled/           # Compiled tables written by this workflow
