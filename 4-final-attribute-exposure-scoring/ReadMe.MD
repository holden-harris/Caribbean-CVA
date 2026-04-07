# 4-final-attribute-exposure-scoring

## Overview

This folder compiles and summarizes final **attribute scores** from individually scored expert review workbooks for the Caribbean Climate Vulnerability Assessment (CVA). The workflow first extracts final reviewer scores into a single compiled table, then summarizes those scores across attribute types and joins them with calculated quantitative exposure scores for downstream vulnerability scoring.

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
