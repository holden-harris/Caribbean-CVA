# Query FishBase Environmental and Trait Data
This script automates the extraction of species-level traits and environmental attributes from [FishBase](https://www.fishbase.org) using the [`rfishbase`](https://github.com/ropensci/rfishbase) R package as part of the **Caribbean Climate Vulnerability Assessment (CVA)** led by Isla Mar and Harris Analytics & Research LLC.

## 📥 Input
The script expects a CSV file named `species-list.csv` in the working directory. It must contain scientific names. 

## Column descriptions
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
