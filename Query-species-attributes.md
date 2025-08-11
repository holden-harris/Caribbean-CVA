# Query FishBase Environmental and Trait Data
This script automates the extraction of species-level traits and environmental attributes from [FishBase](https://www.fishbase.org) using the [`rfishbase`](https://github.com/ropensci/rfishbase) R package as part of the **Caribbean Climate Vulnerability Assessment (CVA)** 

## Final Column Descriptions
| Column Name          | Description                                    |
| -------------------- | ---------------------------------------------- |
| scientific\_name     | Scientific name (from input list)              |
| FBname               | FishBase common name                           |
| BodyShapeI           | Body shape classification                      |
| DepthRangeShallow    | Minimum recorded depth (m)                     |
| DepthRangeDeep       | Maximum recorded depth (m)                     |
| DepthRangeComShallow | Common shallow depth (if reported)             |
| DepthRangeComDeep    | Common deep depth (if reported)                |
| LongevityWild        | Lifespan in the wild (years)                   |
| LongevityCaptive     | Lifespan in captivity (years)                  |
| Vulnerability        | Intrinsic vulnerability index from FishBase    |
| VulnerabilityClimate | Climate-specific vulnerability score           |
| Length               | Maximum length (cm)                            |
| CommonLength         | Typical/common length (cm)                     |
| Weight               | Maximum recorded weight (g)                    |
| n\_growth\_studies   | Number of growth studies retrieved             |
| K\_avg               | Average growth coefficient (Von Bertalanffy K) |
| Loo\_avg             | Average asymptotic length (L∞)                 |
| ReproMode            | Reproductive mode (e.g., oviparous)            |
| Spawning             | Spawning behavior                              |
| RepGuild1, RepGuild2 | Reproductive guild classifications             |
| AddInfos             | Additional reproductive notes                  |
| DietTroph            | Estimated trophic level                        |
| DietSeTroph          | Standard error of trophic level estimate       |
| DietRemark           | Qualitative diet remarks                       |
| AddRems              | Additional ecological notes                    |

## Code Explanation
The R script follows a modular structure that performs the following steps:

### 1. Load species list
The script begins by reading a CSV file named `species-list.csv`, which should include a column `scientific_name` containing the Latin names of fish species of interest. This list serves as the input to all subsequent queries.

```r
sp_list <- read.csv("./species-list.csv")
```
---

### 2. Query FishBase data tables
The script then uses the `rfishbase` package to access multiple FishBase tables via API calls. Each function returns a dataframe containing biological and ecological data for the species in the list:

- `species()` retrieves general species information (e.g., depth range, lifespan, body shape)
- `popgrowth()` returns growth parameters like Von Bertalanffy K and L∞
- `ecology()` provides trophic level and diet information
- `reproduction()` includes reproductive mode and guild classification

```r
spp_raw     <- species(sp_list$scientific_name)
growth_raw  <- popgrowth(sp_list$scientific_name)
ecology_raw <- ecology(sp_list$scientific_name)
repro_raw   <- reproduction(sp_list$scientific_name)
```
> ⚠️ These API calls may take ~1 minute to complete, depending on the number of species.
---

### 3. Join species-level attributes
For each data table, only a subset of relevant columns is retained using `select()`, and the data are joined to the original species list using `left_join()`. The joining is performed on `scientific_name` from the species list and `Species` from FishBase.

#### 3a. Join species summary data
This includes depth ranges, body shape, length, weight, vulnerability, and longevity fields.

```r
sp_attributes <- sp_list %>%
  left_join(select(spp_raw, ...), by = c("scientific_name" = "Species"))
```

#### 3b. Summarize and join growth data
Because a species may have multiple entries in FishBase's growth table, the script calculates the **mean K** and **mean L∞ (Loo)** per species, and joins these summary values.
```r
growth_summary <- growth_raw %>%
  group_by(Species) %>%
  summarise(K_avg = mean(K, na.rm = TRUE), Loo_avg = mean(Loo, na.rm = TRUE), ...)
```

```r
sp_attributes <- sp_attributes %>%
  left_join(growth_summary, by = c("scientific_name" = "Species"))
```

#### 3c. Join reproduction data
Adds reproductive mode, spawning behavior, and guild categories.
```r
left_join(select(repro_raw, ...), by = c("scientific_name" = "Species"))
```

#### 3d. Join ecological traits
Adds trophic level, diet standard error, and remarks from the ecology table.

```r
left_join(select(ecology_raw, ...), by = c("scientific_name" = "Species"))
```

---

### 4. Save final output
The final merged dataframe, `sp_attributes`, contains a comprehensive trait profile for each species.
```r
write.csv(sp_attributes, "fishbase_species_attributes.csv", row.names = FALSE)
```
