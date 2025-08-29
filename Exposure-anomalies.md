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
## Utility Functions
Helper functions used throughout the exposure–overlap analysis. They make it easier to extract clean names, handle bounding boxes, and safely open PDF devices.

### name_from_shp
Extracts a **clean species name** from a shapefile path.  
- Removes file extensions (`.shp`)  
- Inserts spaces between CamelCase (e.g., `AtlanticHerring` → `Atlantic Herring`)  
- Replaces underscores, dashes, and dots with spaces. Trims and normalizes spacing  

### exp_name_from_nc
- Extracts the exposure factor code from a NetCDF filename. Everything before the first underscore (_) is returned.

### bbox_with_pad
- Calculates a bounding box (xlim, ylim) for an sf object (species polygon).
- Adds a proportional padding (default = 5%) so plots have margins around the species range

### safe_open_pdf
Opens a PDF graphics device safely.
- Uses cairo_pdf if available (better text rendering, avoids font issues).
- Falls back to base pdf() if Cairo is unavailable.
- Needed for writing multi-page PDFs of figures.

### species_name_clean
- Formats species names for filenames by replacing spaces with hyphens.

---

## Plotting Functions
Documents the plotting utilities used to visualize species ranges, climate anomalies, and their overlap. Each function is designed to be composable so figures can be combined with **patchwork** or embedded into multi-panel PDFs.

### `plot_distribution()`
**Description:** Draws a clean species distribution map for a given spatial domain.
#### Arguments:
- **sp** *(sf)*: Species polygon(s) in EPSG:4326 (longitude/latitude).  
- **xlim, ylim** *(numeric length-2)*: Map window bounds given as West–East and South–North.  
- **title** *(character, optional)*: Plot title. Typically, a species name.  
#### Notes:
- Uses **`rnaturalearth`** (medium scale) to provide land polygons for geographic context.
- Species polygon is filled **black** with semi-transparency (`alpha = 0.5`) so that coastlines and base features remain visible underneath.
- **`coord_sf(expand = FALSE)`** ensures that the map is cropped exactly to the requested bounding box, with no extra padding.  

---

### `plot_anomalies()`
**Description:** Produces a quick base R raster plot of standardized climate anomalies. Useful for fast diagnostics or lightweight PDF export.  
#### Arguments:
- **anom** *(Raster/SpatRaster)*: Standardized anomaly raster to plot.  
- **exp_name** *(character)*: Exposure factor code (e.g., `"sst"`, `"pp"`), included in the title.  
- **extent** *(character, optional)*: Region label to append to the plot title (e.g., `"Global"`, `"Caribbean Sea"`).  
- **carib_box** *(character, default = `'y'`)*: If `'y'`, draws a purple rectangle highlighting the Caribbean Sea bounding box.  
- **mar** *(numeric length-4, default = `c(2,2,2,2)`)*: Plot margins passed to `par(mar=)`.  
#### Notes:
- Use **`viridisLite::turbo()`** colormap for perceptually uniform anomaly shading.  
- Adds **`maps::map("world")`** for coastlines and country outlines.  
- Optionally overlays a bounding box around the Caribbean Sea (`xlim_carib`, `ylim_carib`).  

---

### `plot_anomalies_gg()`
**Description:** Generates a tidy **ggplot2** map of standardized climate anomalies suitable for multi-panel figures and publication.
#### Arguments:
- **anom** *(SpatRaster/Raster\*)*: Standardized anomaly raster to visualize.  
- **exp_name** *(character)*: Exposure factor code (e.g., `"sst"`, `"pp"`); included in the legend title and plot title.  
- **extent** *(character, optional)*: Region label appended to the title (e.g., `"Global"`, `"W. Atlantic"`, `"Caribbean Sea"`).  
- **carib_box** *(character, default = `'n'`)*: If `'y'`, draws a purple rectangle around the Caribbean Sea bounding box (`xlim_carib`, `ylim_carib`).  
- **uscar_box** *(character, default = `'n'`)*: If `'y'`, draws a magenta rectangle around the U.S. Caribbean bounding box (`xlim_uscar`, `ylim_uscar`).  
#### Notes:
- Converts the raster to a tidy data frame (`as.data.frame(anom, xy = TRUE, na.rm = TRUE)`) and maps values with **`geom_raster()`**.  
- Adds land with **`rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")`**.  
- Uses **`viridisLite::turbo()`** via `scale_fill_gradientn()` with legend title `Standardized\nanomaly\n<exp_name>`.  
- Keeps geographic axes correct with **`coord_sf(..., default_crs = sf::st_crs(4326))`** and crops to the data’s bounding box.  
- Formats longitude/latitude tick labels as degrees with W/E and S/N suffixes; enables `guide_axis(check.overlap = TRUE)` to reduce clutter.  
- Optional bounding boxes: Caribbean (purple) and U.S. Caribbean (magenta) are added with **`annotate("rect", ...)`**.  
- For consistent comparisons across panels, one can fix this with a common color range with `scale_fill_gradientn(..., limits = c(-L, L))` when composing multi-plot figures.

---

### `overlap_range()`
**Description:** Creates a ggplot-based map showing the spatial overlap between a species’ distribution range and standardized climate anomalies. Designed for side-by-side comparison with anomaly histograms and categorical summaries.  

#### Arguments:
- **df** *(data.frame)*: Data frame with gridded anomaly values clipped to the species range; must include `x`, `y`, and `anomaly` columns.  
- **species_name** *(character)*: Species name for labeling.  
- **exp_name** *(character)*: Exposure factor code (e.g., `"sst"`, `"pp"`).  
- **domain** *(character)*: Geographic domain label (e.g., `"Caribbean Sea"`, `"W. Atlantic"`).  
- **anom** *(SpatRaster/Raster\*)*: Original anomaly raster (not directly plotted but included for context).  
- **xlim, ylim** *(numeric length-2)*: Map window bounds (West–East, South–North).  
- **main_title** *(character, default = `"all"`)*: Controls the plot title. Options:  
  - `"all"` → includes species name, exposure factor, domain, and cell count.  
  - `"domain"` → shows only domain and cell count.  
- **legend_title** *(character, optional)*: Title for the anomaly legend. Defaults to `"Standardized\nanomaly\n<exp_name>"`.  
- **world** *(sf, optional)*: Pre-supplied world coastlines; if `NULL`, uses **`rnaturalearth`** (medium scale).  
- **world_scale** *(character, default = `"medium"`)*: Resolution passed to `rnaturalearth`.  
- **base_size** *(numeric, default = 12)*: Base text size for theme.  
- **fill_limits** *(numeric length-2, optional)*: Fixes anomaly color scale (useful for cross-plot consistency).  
- **tick_n** *(integer, default = 4)*: Approximate number of axis tick marks.  
- **minor_by** *(integer, default = 1)*: Spacing for minor graticule lines (degrees).  

#### Notes:
- Uses **`geom_tile()`** to display anomaly values within the species range grid.  
- Draws **minor and major graticule lines** with `sf::st_graticule`, clipped to the domain.  
- Coastlines are added with **`geom_sf()`** from `rnaturalearth`.  
- Color scale uses **`viridisLite::turbo()`** with an optional fixed `limits`.  
- Longitude/latitude ticks are formatted with degree symbols and cardinal directions (W/E, S/N).  
- Titles can flexibly display species name, exposure factor, domain, and grid cell count.  

---

### `anom_histogram_gg()`
**Description:** Builds a **ggplot2** histogram of standardized anomaly values (typically from a species‐masked raster) for quick distribution diagnostics and side-by-side comparison with maps.

#### Arguments:
- **r** *(SpatRaster/Raster\*)*: Anomaly raster whose finite values will be binned. The `anom` file that's passed in is already masked to the species footprint and domain.  
- **fill_limits** *(numeric length-2, optional)*: Fixes the x-axis and legend range (e.g., `c(-2, 2)`) to keep scales consistent across panels.

#### Notes:
- Extracts finite raster values via `terra::values()` and bins into **50** bins.  
- Colors bars by bin center using **`scale_fill_gradientn(colors = turbo(bin_num))`**; legend title set to `"Standardized \nanomaly"`.  
- Enforces a common x-range with **`scale_x_continuous(limits = fill_limits)`** so histograms align with map color scales.  
- Adds a dashed vertical line at **0**.

---

### `anom_summary_bars()`
**Description:** Produces a categorical bar chart summarizing the distribution of standardized anomalies within the species/domain footprint. Bars show counts per anomaly bin with percent labels for quick interpretation.

#### Arguments:
- **r** *(SpatRaster/Raster\*)*: Anomaly raster (typically masked to the species range and domain).  
- **species_name** *(character)*: Species label used in the plot title.  
- **exp_name** *(character)*: Exposure factor code (e.g., `"sst"`, `"pp"`) used in the title.  
- **domain** *(character)*: Domain label (e.g., `"Caribbean Sea"`, `"U.S. Caribbean"`) used in the title.  

#### Notes:
- Extracts finite raster values with `terra::values()` and bins them into five categories using:  
  `< -1.5`, `-1.5 to -0.5`, `-0.5 to +0.5`, `0.5 to 1.5`, `> 1.5`.  
- Colors bars with a fixed five-color palette derived from **`viridisLite::turbo()`**, keyed to bin midpoints (`-2, -1, 0, 1, 2`) for stable category coloring.  
- Displays **counts** as bar heights and **percentages** as labels above bars.  
- Uses `scale_y_continuous(expand = expansion(mult = c(0, 0.05)))` to add headroom for labels and `coord_cartesian(clip = "off")` to avoid clipping.  
- Rotates x-axis labels slightly for readability and removes major x-grid lines to reduce clutter.  
- Title combines species, exposure code, and domain as:  
  `"<species_name> | <exp_name> | <domain>"`.  
- Styled with **`theme_minimal()`** and dark axis text/lines for publication-ready legibility.  

---

## Process species
**Description:** Outer-loop driver that processes a single species shapefile: reads and validates geometry, standardizes the species name, prepares output folders, computes/clips the plotting extent, and opens two multi-page PDF devices for downstream figures (distribution+anomalies and exposure+overlap). This function sets up the per-species workspace; the **inner loop** (runs per exposure factor) adds pages to these PDFs.

### `process_species()`
- Requires **sp_file** *(character)*: File path to a species distribution shapefile (`.shp`). Must be lon/lat (EPSG:4326) or convertible.
- **Dependencies & assumptions**
  - Requires project-level globals to be defined beforehand: `out_dir`, `xlim_nwa`, `ylim_nwa`.
  - Assumes the **Utility Functions** are available in scope: `name_from_shp()`, `species_name_clean()`, `bbox_with_pad()`, and `safe_open_pdf()`.
  - Coordinate reference system should be geographic (EPSG:4326). If inputs differ, reproject prior to calling or adapt the loader accordingly.

**Geometry & naming**
  - Reads polygon(s) with `sf::st_read(sp_file, quiet = TRUE)` and repairs issues via `st_make_valid()`.
  - Derives a human-readable species name with `name_from_shp()` and a file-safe slug using `species_name_clean()` (e.g., `"Atlantic Herring"` → `"atlantic-herring"`).
  - **Console logging**: Prints a clear species header and basic status to the console to trace progress in batch runs.

**Output structure**
  - Creates a per-species subfolder: `file.path(out_dir, <slug>)`.
  - Prepares two multi-page PDFs (opened now; pages added later by the inner loop):
    - `<slug>_Distribution-Anomalies.pdf` – species distribution + anomaly maps.
    - `<slug>_Exposure-Overlap.pdf` – overlap map, histogram, and categorical summary.

**PDF handling**
  - Uses `safe_open_pdf()` (prefers `cairo_pdf`, falls back to base `pdf()`).
  - Immediately records the current device numbers and registers paired `on.exit()` blocks that safely close each device if still open.

**Extent computation & clipping**
  - Computes a tight plotting extent from the species geometry with `bbox_with_pad(sp, pad = 0.05)` for a small visual margin.
  - If the species’ bbox exceeds the broader **W. Atlantic** domain (`xlim_nwa`, `ylim_nwa`), it clips the working bbox to those limits. This avoids the computational expense of global plotting and keeps maps focused/reproducible.

**Side effects and return value**
  - Creates folders and opens two PDF devices (left open for the inner loop to append pages).
  - Note that page generation happens in the inner loop and invisibly returns nothing; its primary purpose is setup and device lifecycle management.

## Inner loop: Process Exposure Factor

### Write out plots

## References
