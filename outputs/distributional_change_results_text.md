# Distributional Change Potential — Methods and Results Text

> **Status:** Methods paragraph is finalized from script review. Results paragraph contains bracketed placeholders — fill in after running Scripts 1–3 and examining console output and figures.

---

## Methods

### Conceptual framework

The potential for distributional change (DCP) quantifies each stock's intrinsic capacity to shift its geographic range in response to changing ocean conditions, independent of whether climate change is projected to force such a shift. Stocks with high adult mobility, dispersive early life stages, generalist habitat use, and broad latitudinal ranges are assumed to have greater ability to track favorable conditions as they move; stocks with the opposite traits are more likely to remain within their current range even as conditions deteriorate. This component was assessed for all 25 U.S. Caribbean stocks using the framework established by prior NOAA Fisheries Climate Vulnerability Assessments (CVAs) for the Gulf of Mexico (Quinlan et al. 2023), South Atlantic (Craig et al. 2025), and Highly Migratory Species (Loughran et al. 2025).

### Attribute selection

Four Biological Sensitivity Attributes (BSAs) scored in the Module 4 reviewer process were used:

| Attribute | Role | Score direction |
|---|---|---|
| Adult mobility | Movement | Inverted |
| Habitat specificity | Movement | Inverted |
| Mobility and dispersal or early life stages | Movement | Inverted |
| Species range | Temperature sensitivity analog | Not inverted |

The first three attributes relate to the physical capacity for movement or dispersal and were inverted before analysis (see below). The fourth attribute — Species range — serves as a Caribbean-specific analog for the "Sensitivity to Temperature" attribute used in prior NOAA CVAs. Both attributes use latitudinal extent as a proxy for thermal tolerance: stocks with wide latitudinal ranges occupy broad climate envelopes and are assumed to have lower temperature sensitivity and therefore lower innate propensity to shift in response to warming. Species range was therefore used without inversion, consistent with prior CVA convention.

> **Open methodological decision:** The Caribbean CVA does not directly score "Sensitivity to Temperature." The substitution of Species range as an analog is a methodological choice that should be confirmed by the project lead before this analysis is finalized. The analysis scripts are parameterized so that Species range can be excluded by removing it from the target attribute list without any other code changes.

### Score transformation (inversion)

Because a low score on movement-related attributes indicates *high* capacity to shift distribution (e.g., a score of Low (1) on "Adult mobility" means the stock is highly mobile, and a score of Very High (4) means it is essentially sessile), the three movement attributes were inverted prior to analysis using the transformation:

**Transformed score = 5 − original mean score**

This maps the original ordinal scale (Low = 1, Moderate = 2, High = 3, Very High = 4) to a reversed scale where higher values consistently indicate greater propensity to shift:

| Original score | Original category | Transformed score | Transformed category |
|---|---|---|---|
| 1 | Low | 4 | Very High |
| 2 | Moderate | 3 | High |
| 3 | High | 2 | Moderate |
| 4 | Very High | 1 | Low |

Species range was not transformed; a higher original score (wider range) was treated directly as greater distributional change potential.

The transformation was applied to attribute means derived from Module 4. Attribute means reflect the consensus of four expert reviewers, each providing tally votes across a standardized 1–5 vote pool (five total votes per reviewer per attribute). The mean of all reviewer votes is stored in `attribute_means_uscar.csv` and used as the transformed-mean input to the DCP logic model.

### DCP classification (FCVA logic model)

Each stock was assigned a DCP rank using the same cascading threshold logic model applied to overall climate vulnerability in Modules 4 and 8. With a rank threshold of 2, the rules are:

| DCP Rank | Condition |
|---|---|
| Very High | More than 3 of the 4 transformed attribute means ≥ 3.5 |
| High | More than 2 of the 4 transformed attribute means ≥ 3.0 |
| Moderate | More than 2 of the 4 transformed attribute means ≥ 2.5 |
| Low | All other cases |

The thresholds (3.5, 3.0, 2.5) and rank threshold (2) are identical to those used in Modules 4 and 8, ensuring methodological consistency across all CVA components. The resulting DCP rank is ordinal (Low < Moderate < High < Very High) and represents each stock's potential to shift distribution in response to climate change.

### Bootstrap uncertainty quantification

Bootstrap uncertainty was quantified using the same draw-pile resampling method applied in Module 8 (overall vulnerability uncertainty analysis), adapted for the DCP attribute set and inversion structure.

**Draw pile construction.** For each stock × attribute combination, reviewer tally votes were pooled across all four reviewers to form a draw pile of 20 individual votes (4 reviewers × 5 votes each). Tally counts were drawn from `sensitivity_tallies_long.csv`, which records the number of votes each reviewer assigned to each ordinal category (Low, Moderate, High, Very High) for each attribute.

**Tally-level inversion.** For the three movement attributes, inversion was applied at the tally count level before the draw piles were constructed by swapping the counts of opposite categories:

- New tally_L = original tally_VH
- New tally_M = original tally_H  
- New tally_H = original tally_M
- New tally_VH = original tally_L

This swap is mathematically equivalent to applying the transformation `5 − vote` to each individual tally vote before pooling, and produces draw piles with means equal to `5 − (original draw-pile mean)`. Performing the swap at the count level — rather than inverting each draw at sample time — is computationally equivalent and avoids repeated per-iteration conditional logic.

**Baseline reproduction gate.** Before bootstrapping, the mean of each inverted draw pile was computed and the FCVA logic model was applied to verify that the resulting DCP ranks reproduce the baseline ranks computed from attribute means in Script 1. Any stock for which the draw-pile-derived rank differed from the Script 1 baseline rank caused the script to halt with an informative error message, preventing bootstrap results from being produced from an inconsistent parameterization. This gate passed for all 25 stocks.

**Bootstrap procedure.** Ten thousand (10,000) bootstrap iterations were run per stock (`set.seed(99)`, matching Module 8). In each iteration, all four draw piles for that stock were independently resampled with replacement (20 draws per pile), the resampled mean was computed for each attribute, and the FCVA logic model was applied to produce a bootstrap DCP rank. The proportion of iterations yielding each rank (prop_L, prop_M, prop_H, prop_VH) was computed from the 10,000-iteration distribution.

**Dominant rank and borderline flagging.** The dominant rank is the rank that occurred most frequently across bootstrap iterations; the dominant proportion (dominant_prop) is the fraction of iterations in which this rank was assigned. Stocks were flagged as borderline if their dominant proportion was less than 0.75 (i.e., the dominant rank was not returned in at least 75% of iterations), consistent with the Module 8 borderline threshold.

**Certainty encoding.** Bootstrap certainty was encoded visually in the figures using font face and color, following the scheme established in Module 4 (Figure 1) and applied consistently throughout the CVA:

| Dominant proportion | Certainty level | Font face | Font color |
|---|---|---|---|
| > 0.95 | Very High | Bold | Black |
| 0.90–0.95 | High | Italic | Black |
| 0.67–0.89 | Moderate | Bold | White |
| ≤ 0.66 | Low | Italic | White |

---

## Results paragraph

> **Status: FILLED IN** — values derived from Scripts 1–3 run on 2026-05-07.

None of the 25 assessed stocks were classified as having **Very High** potential for distributional change; 3 were classified as **High**, 11 as **Moderate**, and 11 as **Low** (Figure A). Bootstrap uncertainty was substantial across the Moderate and Low ranks: 10 of the 25 stocks were borderline (dominant DCP rank assigned in fewer than 75% of bootstrap iterations), including seven stocks classified as Moderate (Atlantic thread herring, 57.1%; Ballyhoo, 55.6%; Red grouper, 56.0%; Yellowtail snapper, 59.0%; King mackerel, 68.8%; Queen conch, 70.3%; Silk snapper, 70.6%) and three classified as Low (Long-spined sea urchin, 63.8%; Queen snapper, 63.9%; White mullet, 71.9%). The three stocks ranked High — Blue runner, Dolphinfish, and Spiny lobster — were all assigned with moderate to very high certainty (dominant proportions of 85.1%, 79.6%, and 100%, respectively) and share traits of high adult mobility, broad latitudinal ranges, and relatively dispersive early life stages. Stocks classified as Low DCP were predominantly reef-associated or benthic invertebrates with sessile or low-mobility life histories, limited habitat generalism, and narrow depth or habitat associations (e.g., Gray angelfish, Hogfish, Misty grouper, Queen triggerfish, Sea cucumbers, Stoplight parrotfish, Rainbow parrotfish).

The cross-plot of DCP against overall climate vulnerability (Figure B) revealed a striking pattern: no stock with High or Very High overall vulnerability was assigned a DCP rank higher than Moderate, and 11 of the 25 stocks fell into the highest-concern quadrant — classified as both **High or Very High** in overall climate vulnerability and **Low or Moderate** in distributional change potential. Three of these stocks carry the highest climate vulnerability rank (Very High): Long-spined sea urchin (Diadema), Mutton snapper, and Rainbow parrotfish, all of which were classified as Low DCP, indicating limited capacity to seek thermal refuge through range shift. The remaining eight stocks in this quadrant — Gray angelfish, Hogfish, Misty grouper, Queen triggerfish, Sea cucumbers, and Stoplight parrotfish (all High vulnerability, Low DCP), and Lane snapper and Nassau grouper (High vulnerability, Moderate DCP) — face high climate exposure and sensitivity while also lacking the mobility traits that would allow distributional tracking of preferred environmental conditions. None of the three highest-DCP stocks (Blue runner, Dolphinfish, Spiny lobster) were ranked High or Very High in overall vulnerability, suggesting that the stocks most capable of distributional response are not those facing the greatest climate risk in the Caribbean context.

---

## Spot-check validation

Before finalizing, spot-check at least three stocks by computing transformed means and applying the logic model by hand:

1. Choose one Very High, one High, and one Low or Moderate stock.
2. Retrieve the four attribute means from `attribute_means_uscar.csv` for each stock.
3. Apply inversion to the three movement attributes: transformed mean = 5 − original mean.
4. Count how many transformed means meet each threshold (≥ 3.5, ≥ 3.0, ≥ 2.5).
5. Apply the logic model: if > 3 meet ≥ 3.5 → VH; else if > 2 meet ≥ 3.0 → H; else if > 2 meet ≥ 2.5 → M; else L.
6. Confirm the result matches `dcp_rank` in `distributional_change_potential_uscar.csv`.

| Stock | Adult mobility (inv.) | Habitat spec. (inv.) | Early life disp. (inv.) | Species range | Counts (≥ 3.5 / ≥ 3.0 / ≥ 2.5) | Hand rank | Script rank | Match? |
|---|---|---|---|---|---|---|---|---|
| Blue runner (High) | 5 − 1.35 = **3.65** | 5 − 1.25 = **3.75** | 5 − 1.45 = **3.55** | 1.70 | 3 / 3 / 3 | n_ge35=3, not >3 → not VH; n_ge30=3 >2 → **High** | High | ✓ |
| Nassau grouper (Moderate) | 5 − 2.10 = **2.90** | 5 − 2.35 = **2.65** | 5 − 2.20 = **2.80** | 2.70 | 0 / 0 / 4 | n_ge35=0, n_ge30=0; n_ge25=4 >2 → **Moderate** | Moderate | ✓ |
| Gray angelfish (Low) | 5 − 3.20 = **1.80** | 5 − 2.95 = **2.05** | 5 − 2.10 = **2.90** | 1.95 | 0 / 0 / 1 | n_ge35=0, n_ge30=0; n_ge25=1, not >2 → **Low** | Low | ✓ |

---

## References cited in this section

- Craig, J.K., et al. (2025). South Atlantic Climate Vulnerability Assessment. NOAA Technical Memorandum.
- Loughran, T., et al. (2025). Highly Migratory Species Climate Vulnerability Assessment. NOAA Technical Memorandum.
- Morrison, W.E., et al. (2015). A Methodology for Assessing Vulnerability of Marine Fish and Invertebrates to a Changing Climate. NOAA Technical Memorandum NMFS-OSF-2.
- Quinlan, J.A., et al. (2023). Gulf of Mexico Climate Vulnerability Assessment. NOAA Technical Memorandum.
