# Caribbean CVA: Exposure Analyses
This script creates figures to evaluate climate exposures for federally managed and ecologically important species in the U.S. Caribbean.

---
## Setup
### Packages needed
```r
library(dplyr)        # Data wrangling
library(raster)       # Raster data handling
library(terra)        # Modern raster/spatial package
library(sp)           # Spatial data classes
library(sf)           # Modern simple features (vector data)
library(viridisLite)  # Color palettes
library(ggplot2)      # Data visualization
library(rnaturalearth) # Base maps and country outlines
library(scales)       # Axis formatting
library(ggplotify)    # Convert plots to grobs
library(patchwork)    # Combine plots
```

### Set directory paths
The directory structure for this project is organized as follows:
```r
project/
├─ data/
│  ├─ cmip6/                              # NetCDF exposure files (*.nc)
│  └─ species-distribution-shapefiles/    # Species polygons (*.shp + sidecars)
└─ outputs/
   └─ exposure-overlap/                   # PDFs & PNGs will be written here
```

### Set geographic extent 
Define spatial domains used in plotting and analysis. 
Bounding box: W, E, S, N

```r
## Caribbean Sea
xlim_carib <- c(-92, -57)
ylim_carib <- c(  6,  28)

## U.S. Caribbean
xlim_uscar <- c(-69, -63.0)
ylim_uscar <- c( 16,  20.0)

## W. Atlantic Ocean
xlim_nwa <- c(-99, -40)  # W Atlantic: 99°W to 40°W
ylim_nwa <- c( -5,  72)
```
