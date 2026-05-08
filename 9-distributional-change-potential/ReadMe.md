# Module 9 — Potential for Distributional Change

## Purpose

Calculates and visualizes each stock's potential to shift its geographic distribution in response to changing environmental conditions. Stocks with mobile adults, dispersive early life stages, generalist habitat use, and high temperature sensitivity are classified as having high distributional change potential (DCP); stocks with the opposite traits are classified as low.

This analysis follows the same approach used in prior NOAA CVAs (HMS: Loughran et al. 2025; South Atlantic: Craig et al. 2025; GoM: Quinlan et al. 2023), applying the same FCVA logic model to a subset of biological sensitivity attributes (BSAs) scored in Module 4. Bootstrap uncertainty is quantified using the same draw-pile resampling method as Module 8.

---

## Attribute selection and inversion logic

Four BSAs are used. Three movement-related attributes are **inverted** (`5 − mean`) so that a high transformed score consistently indicates greater propensity to shift:

| BSA (exact CSV name) | Equivalent in prior CVAs | Direction |
|---|---|---|
| `Adult mobility` | Adult Mobility | Inverted |
| `Habitat specificity` | Habitat Specificity | Inverted |
| `Mobility and dispersal or early life stages` | Mobility and Dispersal of Early Life Stages | Inverted |
| `Species range` | Sensitivity to Temperature (analog) | Not inverted |

**Open methodological decision — Species Range substitution:**

Prior NOAA CVAs use a "Sensitivity to Temperature" attribute derived from known temperature of occurrence and latitudinal range. The Caribbean CVA does not score this attribute directly. `Species range` is used as an analog here because both attributes use latitudinal extent as a proxy for thermal tolerance: a wider latitudinal range implies broader thermal tolerance (lower temperature sensitivity, lower propensity to shift).

> **This substitution should be confirmed by the project lead before the analysis is finalized.** To use only the three movement-related attributes, remove `"Species range"` from `target_attributes` at the top of Scripts 1 and 2. No other code changes are required.

---

## Inputs

| File | Produced by | Used in |
|------|-------------|---------|
| `outputs/final-scores-compiled/overall-vulnerability-rankings/attribute_means_uscar.csv` | Module 4 Script 3 | Script 1 |
| `outputs/final-tallies-long/sensitivity_tallies_long.csv` | Module 7 Script 2 | Script 2 |
| `outputs/final-scores-compiled/overall-vulnerability-rankings/overall_vulnerability_scores_uscar.csv` | Module 4 Script 3 | Script 3 |

---

## Outputs

| File | Description |
|------|-------------|
| `outputs/distribution-change-potential/distributional_change_potential_uscar.csv` | Baseline DCP scores and ranks for all 25 stocks |
| `outputs/distribution-change-potential/distributional_change_bootstrap_uscar.csv` | Bootstrap rank-proportion distribution and borderline flags |
| `outputs/distribution-change-potential/distributional_change_full_uscar.csv` | Scripts 1 + 2 joined; primary input for Script 3 |
| `figures/fig_distributional_change_ranks.png` | Figure A: stocks by DCP rank category |
| `figures/fig_distributional_change_vs_vulnerability.png` | Figure B: DCP vs. overall climate vulnerability cross-plot |

---

## Execution order

Run scripts in numbered order from the project root (`.Rproj` file open in RStudio):

1. `1-calculate-distributional-change-potential.R` — baseline DCP scores
2. `2-bootstrap-distributional-change.R` — uncertainty quantification
3. `3-plot-distributional-change.R` — publication figures

Script 2 requires the output of Script 1. Script 3 requires the output of Script 2.

---

## Methods

**Inversion rule.** For movement attributes, the transformed score is `5 − original_mean`, mapping Low (1) → Very High (4), Moderate (2) → High (3), High (3) → Moderate (2), Very High (4) → Low (1). After inversion, a high value indicates high propensity to shift across all four attributes.

**FCVA logic model.** Identical to Modules 4 and 8. With `rank_threshold = 2`:

| DCP Rank | Condition |
|----------|-----------|
| Very High | > 3 transformed means ≥ 3.5 |
| High | > 2 transformed means ≥ 3.0 |
| Moderate | > 2 transformed means ≥ 2.5 |
| Low | all other cases |

**Bootstrap.** For each stock × attribute, reviewer tally votes are pooled across 4 reviewers (4 × 5 = 20 votes). For inverted attributes, tally counts are swapped (L↔VH, M↔H) before the draw pile is built, producing the same result as inverting each individual vote. 10,000 bootstrap iterations resample each draw pile with replacement. `bootstrap_seed = 99` matches Module 8.

**Borderline flag.** A stock is flagged borderline if the dominant rank accounts for fewer than 75% of iterations (`borderline_threshold = 0.25`). Matches Module 8.

**Baseline reproduction gate.** Before bootstrapping, Script 2 verifies that draw-pile means reproduce the Script 1 DCP ranks for all 25 stocks. The script halts with an informative error if any stock disagrees.

**Figure font encoding.** Bootstrap certainty encoded as font face and color (matches Module 4 Fig. 1):

| Dominant prop | Certainty | Face | Color |
|---|---|---|---|
| > 0.95 | Very High | bold | black |
| 0.90–0.95 | High | italic | black |
| 0.67–0.89 | Moderate | bold | white |
| ≤ 0.66 | Low | italic | white |
