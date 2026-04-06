# 4-final-directional-effect-scoring

## Overview

This folder compiles **directional effect scores** (Positive, Neutral, Negative) from individually scored expert review workbooks as part of the Caribbean Climate Vulnerability Assessment (CVA). Each reviewer independently scored fish stocks using a standardized Excel template, and this code extracts and consolidates those scores into a single, analysis-ready table.

Directional effect scores characterize whether climate change is expected to have a positive, neutral, or negative impact on a given fish stock. Each score is read from a fixed cell range (M38:M40) within species-specific tabs of each reviewer's workbook.

## Final Products

| File | Description |
|------|-------------|
| `table_directional_effect_scores.csv` | Long-format table of all directional effect scores by stock, reviewer, and effect direction |
| `n_reviews_directional-effect.csv` | Summary of how many reviewers scored each stock, with reviewer names listed |

## Setup

### Dependencies

- R (≥ 4.0)
- R packages: `readxl`, `xml2`, `dplyr`, `stringr`, `tidyr`, `tibble`, `openxlsx`

### Input Directory Structure
