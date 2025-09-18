# Caribbean-CVA

## Overview
The general methods for the Caribbean CVA exposure factor analyses were adapted from [Loughran et al. (2025)](https://doi.org/10.1371/journal.pclm.0000530) for their CVA of U.S. highly migratory species, which adapted and used methods from over a decade of conducting CVAs by NOAA (e.g., [Morrison et al. 2015](https://doi.org/10.7289/V5TM782J), [Hare et al. 2016](https://doi.org/10.1371/journal.pone.0146756), [Craig et al. 2025](https://doi.org/10.1371/journal.pclm.0000543)). 

The overall goal of the exposure factor analyses is to compare projected future ocean conditions against their past. To do so, we utilize spatially explicit ocean model projections (both historical and future ocean projections) to create a standardized anomaly map (i.e., a static comparison) and provide accompanying analyses. The values of the standardized anomalies are expressed in standard deviation units, which allows the result to be comparable across variables. For each exposure factor anomaly, we mapped the gridded cells of the standardized exposure anomaly and overlapped a given species distribution. Values from the overlapped cells were then compiled into frequency distributions (histograms) and categorized (bar plots). 

The exposure overlap analyses below are presented within three geographic scopes: 
1. A given species distribution within the Western Atlantic Ocean (5°S–72°N, 99°W–40°W),
2. The wider Caribbean region (6°N–28°N, 92°W–57°W), and
3. The U.S Caribbean (16°N–20°N, 69°W–63°W)
The spatial scale of the analysis will have differing implications for a species' ecology and its management. All three spatial scales are presented for use at the discretion of the expert reviewers. 

---

## Methods for Calculating Gridded Standardized Exposure Factor Anomalies
Similar to the HMS CVA, we utilized projections from the Coupled Model Intercomparison Project Phase 6 (CMIP6) multi-model ensemble under the SSP5-8.5 scenario, which is consistent with past NOAA CVAs. This scenario represents the highest greenhouse gas emissions pathway and likely represents the greatest changes that we can reasonably expect for the exposure factor; it should be considered an upper bound for CMIP6 projections. 

Monthly outputs for each exposure factor (e.g., sea surface temperature) were obtained at a grid cell resolution of 1° latitude × 1° longitude. To calculate detrended standardized anomalies, a unitless value (z) was computed for each 1° × 1 ° grid cell as the difference between the future and historical values divided by the historical standard deviation:

z = (μ future - μ historical )  / 𝝈 historical,

where μ future is the mean of a given grid cell for months during 2020–2049, μ historical is the mean of a given grid cell for months during 1985–2014, and 𝝈 is the interannual standard deviation during the same historical baseline period (1985–2014). Dividing by the historical baseline standard deviation converts z into standard deviation units (σ). 

These standardized anomaly values represent the size of projected change relative to typical historical year-to-year variability. In other words, the size of the calculated standardized anomaly is determined by both the difference between the future and historical monthly average (a larger difference proportional to a larger standardized anomaly) and the size of historical variability (a larger historical standard deviation is proportional to a smaller standardized anomaly). 

The calculated standardized anomaly values can be positive or negative. A positive standardized anomaly value (+σ) means the exposure factor is projected to increase compared to the historical baseline, while a negative value (–σ) indicates a decrease. For example, a grid cell with a +2.0σ anomaly in sea surface temperature (SST) means that the average future SST is projected to be about two standard deviations higher than the historical mean, which would represent a shift twice as large as the typical year-to-year variability observed during 1985–2014. Conversely, a grid cell with a –2.0σ anomaly in SST would indicate that the average future SST is projected to be two standard deviations lower than the historical mean, likewise representing a shift twice the size of normal historical variability, but in the opposite direction.

The exposure factor analyses (example below) were then calculated as the overlap between the standardized exposure anomaly and a species distribution. Species ranges were determined from shape files for the International Union for Conservation of Nature (IUCN) distributions as the baseline geographical limits for each species. 

A weighted average was calculated consistent with past NOAA CVAs. For a given species or stock and exposure factor, standardized anomalies from the gridded overlapped area were grouped into 0.25 standard deviation bins. These bins were then coded by their absolute value in the following categories (LMHV): Low (|σ| < 0.5), Moderate (0.5 ≤ |σ| < 1.5), High (1.5 ≤ |σ| < 2.0), Very High (|σ| ≥ 2.0). For each exposure factor and geographic extent, a weighted average score of absolute values was calculated with the following:

Weighted Average = (1L + 2M +3H + 4V) / (L + M + H + V) 

where L, M, H, and V are the sum number values in each category (low, moderate, high, or very high, respectively) with corresponding multipliers of Lx1, Mx2, Hx3, or Vx4 used to calculate the weighted value. This weighted average score ranges from a minimum of 1 to a maximum of 4. For example, if all of the σ values were within -0.5 to +0.5, the calculated weighted average score would be 1.0; if 100% of the σ values were less than -2 or greater than +2, its weighted average score would be 4.0. 

---

## Exposure Overlap Analysis
<p align="center">
  <img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/outputs/exposure-overlap-12panel/King-Mackerel/Exposure-Overlap-12panel/King-Mackerel_Exposure-Overlap_o200.png?raw=true" 
       alt="King Mackerel – Exposure Overlap (o200)" 
       width="700"/>
</p>

#### Figure 1: Example of exposure overlap analysis for King Mackerel (Scomberomorus cavalla) and dissolved oxygen at 200 m (o200).

The exposure overlap analysis figure is a 3 × 4 panel grid organized in two dimensions. The three columns represent the three nested geographic extents.  
- Left column (A, D, G) represents the species’ range within the Western Atlantic for stock-wide context. 
- Middle column (B, E, H) crops the map for the wider Caribbean (center) and panels in this column only show values for grid cells within this region. 
- Right column (C, F, I) further crops the map and only includes grid cells in the U.S. federally-managed waters offshore Puerto Rico and the USVI. 

The four rows represent overlap analyses for a given standardized exposure factor anomaly (σ).
- **Row 1 (A–C): Spatial overlap maps.** 
  - Gridded maps of the standardized CMIP6 exposure factor anomaly grid cells (1°×1° cells)  that overlap with a given species’ distribution.
  - Standardized anomalies are expressed in σ-units (change relative to the detrended interannual variability of the 1985–2014 baseline).
  - The legend coloring is consistent across the three plots and determined by the range in values for the W. Atlantic species range (panel A).
- **Row 2 (D–F): Categorical summations.**
  - These panels show the total count (y-axis) of standardized exposure factor anomaly values within the following categories: less than -1.5σ, -1.5σ to -0.5σ, –0.5σ to +0.5σ, +0.5 to +1.5, greater than +1.5σ. The top of each bar indicates the proportion of each category, which sums to 100%.
  - Unlike the categorical summary plots in row 4, the sign of the anomaly is unchanged and so represents both magnitude and direction. For example, a positive sign for SST would indicate warming, and a negative sign would indicate cooling. 
- **Row 3 (G–I): Frequency distributions.**
  - Frequency histograms (y-axis = percentage) of standardized anomalies from the gridded overlapped 0.25σ bins ranging from the minimum and maximum standardized anomalies. 
  - These are color-coded at four exposure levels based on their absolute value and correspond to the ‘LMHV’ categorical plots in Row 4:  
     - L: Low (|σ| < 0.5)
     - M: Moderate (0.5 ≤ |σ| < 1.5)
     - H: High (1.5 ≤ |σ| < 2.0)
     - V: Very High (|σ| ≥ 2.0)
- **Row 4 (J–L): LMHV summary and weighted averages.**
  - The categorical bins are collapsed into proportions of Low, Moderate, High, and Very High anomalies (LMHV), which are congruent with past CVAs.
  - The LMHV categories represent the absolute magnitude of the standardized anomalies |σ| and removes directionality. For example, values less than –2σ and greater than +2σ anomalies both contribute to “Very High.”
  - The weighted average score is indicated in the top-left corner. 

The example plot shows the overlap for King Mackerel and the standardized anomaly values for dissolved oxygen at 200 m (o200), which indicate how much o200 is projected to change relative to historical year-to-year variability. A grid cell with a +1σ anomaly indicates that oxygen at 200 m is projected to be one standard deviation higher than the historical mean, suggesting an increase in oxygen availability beyond what would normally occur due to natural variability. Conversely, a –1σ anomaly indicates oxygen is projected to be one standard deviation lower than the historical mean, reflecting a comparable decrease in oxygen relative to typical variability. 

The contribution of these values to a species’ vulnerability assessment will be determined by expert reviewers. For example, an increase in oxygen (+σ) may alleviate physiological stress and expand suitable habitat, while a decrease (–σ), indicating deoxygenation, could compress available habitat and increase stress for pelagic predators such as King Mackerel. The absolute values of the standardized anomalies are shown in the LMHV plots on the final row. 

Exposure overlap for all species is available here: https://github.com/holden-harris/Caribbean-CVA/tree/main/outputs/exposure-overlap-12panel 

---

## Reference Figures: Species Distribution and Exposure Factor Anomaly
Plots for the species distribution and exposure factor, which were used to create the overlap exposure analyses, are also provided for reference. 

<p align="center">
  <img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/outputs/exposure-overlap-12panel/King-Mackerel/Distribution-Anomalies/King-Mackerel_Distribution-Anomalies_o200.png?raw=true" 
       alt="King Mackerel – Distribution Anomalies (o200)" 
       width="700"/>
</p>


#### Figure 2: Example of species distribution and standardized exposure anomaly for King Mackerel (Scomberomorus cavalla) and dissolved oxygen at 200 m (o200).

The **first row, panels A–B**, shows the species’ distribution: including (A) their entire distribution within the Western Atlantic and (B) Caribbean extents. 

**Rows 2 and 3 (panels C–F)** show the CMIP6 ensemble standardized anomaly for o200 (future minus baseline, divided by baseline interannual SD) at four nested domains: (C) global, (D) Western Atlantic, (E) Caribbean Sea, and (F) U.S. Caribbean. Warm colors denote higher o200 relative to the baseline; cool colors denote lower o200. 

---

## Oceanographic Exposure Variables
For each species under CVA review, exposure analyses were conducted for the following 13 oceanographic exposure variables. Note that all exposure overlap analyses were conducted, even if they may not be relevant for the species. For example, changes in bottom temperature may not be useful to consider for our case of King Mackerel, given their epipelagic distribution. All analyses are provided as reference, and they can be used or not at the discretion of the expert reviewers. 

| Abbreviation | Full Variable Name                    | Description                                                                                                                             |
| ------------ | ------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------- |
| **bs**       | Bottom Salinity                       | Mean salinity (salt concentration) near the seafloor. Important for benthic organisms sensitive to freshwater inputs or stratification. |
| **bt**       | Bottom Temperature                    | Mean temperature near the seafloor, influencing demersal fish and benthic invertebrates.                                                |
| **chl**      | Chlorophyll-a Concentration           | A proxy for phytoplankton biomass, indicating primary productivity and food availability at the base of the food web.                   |
| **mld**      | Mixed Layer Depth                     | Depth of the upper, well-mixed surface ocean layer, affecting nutrient availability, light, and stratification.                         |
| **msstg**    | Mean Sea Surface Temperature Gradient | Spatial temperature gradient at the surface; an indicator of thermal fronts, ocean circulation, and habitat boundaries.                 |
| **o200**     | Oxygen at 200 m                       | Dissolved oxygen concentration at \~200 meters depth; reflects mid-water oxygen availability and deoxygenation trends.                  |
| **ph**       | Surface pH                            | Measure of ocean acidity (linked to CO₂ uptake and ocean acidification). Lower values = more acidic.                                    |
| **pp**       | Primary Production                    | Gross primary productivity of phytoplankton; determines energy input to marine food webs.                                               |
| **precip**   | Precipitation                         | Rainfall over the ocean, relevant for freshwater input, stratification, and coastal salinity changes.                                   |
| **sso**      | Sea Surface Oxygen                    | Dissolved oxygen concentration at the surface, important for respiration of pelagic species.                                            |
| **sss**      | Sea Surface Salinity                  | Surface salt concentration, reflecting freshwater input, evaporation, and circulation.                                                  |
| **sst**      | Sea Surface Temperature               | Temperature of the upper ocean, widely used as a climate indicator and driver of species distributions.                                 |
| **swsm**     | Surface Wind Speed Magnitude          | Intensity of winds at the ocean surface, a driver of mixing, upwelling, and surface currents.                                           |
