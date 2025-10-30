# Caribbean CVA: Exposure Factor Analyses
The overall goal of the exposure factor analyses is to compare projected future ocean conditions against their past. To do so, we utilize spatially explicit ocean model projections (both historical and future ocean projections) to create a standardized anomaly map (i.e., a static comparison) and provide accompanying analyses. The values of the standardized anomalies are expressed in standard deviation units, which allows the result to be comparable across variables. For each exposure factor anomaly, we mapped the gridded cells of the standardized exposure anomaly and overlapped a given species distribution. Values from the overlapped cells were then compiled into frequency distributions (histograms) and categorized (bar plots). 

This script synthesizes spatial biological and oceanographic information to help evaluate climate exposures for federally managed and ecologically important species in the U.S. Caribbean. The first set of maps produced shows a species' distribution at two scales and exposure factors at four scales. These provide context for the **exposure overlap figures**, which display the spatial overlap between the species distribution and the anomalies, along with a histogram of this data and a categorical bar summary. For the Caribbean CVA, the workflow to produce these maps involves looping over 25 species and 13 exposure factors, resulting in a total of 650 pages. These are grouped into PDFs by species and exposure factor for review by the CVA experts. 

The exposure overlap analyses below are presented within three geographic scopes: a given species distribution within the Western Atlantic Ocean (5°S–72°N, 99°W–40°W), the wider Caribbean region (6°N–28°N, 92°W–57°W), and the U.S Caribbean (16°N–20°N, 69°W–63°W). 

## Final products
**A. Exposure–Overlap Panel (3 × 4 grid; panels A–L)**  
For the **W. Atlantic range**, **Caribbean Sea**, and **U.S. Caribbean** domains, the code builds:
1. **Overlap map** (anomaly tiles within the species footprint).
2. **Categorical bar summary** of anomalies (<-1.5, -1.5 to -0.5, -0.5 to +0.5, 0.5 to 1.5, > 1.5)
3. **Histogram** of standardized anomalies (mask-consistent and scale-aligned) grouped into 0.25 standard deviation bins.
4. **LMHV summary** of absolute values for anomalies. The histogram and summary plots are coded as follows: Low (|σ| < 0.5), Moderate (0.5 ≤ |σ| < 1.5), High (1.5 ≤ |σ| < 2.0), Very High (|σ| ≥ 2.0)

The panels are combined into a figure with a **single shared legend**, then saved to:
- A **multi-page PDF** (one page per exposure factor) named `"<species>_Exposure-Overlap.pdf"`.
- A **PNG** per exposure factor under `.../Exposure-Overlap/`.

**A. Distribution + Anomaly Panel (2 × 3 grid)**  
These provide reference plots for the exposure-overlap analyses. For each species–exposure combination, a six-panel figure is created that shows:
1. **Species distribution (W. Atlantic bbox)** — clean polygon map with coastline context.  
2. **Species distribution (Caribbean bbox)** — zoomed map for regional context.  
3. **Global anomaly map** for the selected exposure factor.  
4. **W. Atlantic anomaly map** (subset of global).  
5. **Caribbean Sea anomaly map** (subset of global; U.S. Caribbean box optionally shown).  
6. **U.S. Caribbean anomaly map**.  

---

### Spatial domains (extents)
All figures are rendered for three consistent spatial windows:
- **Northwest Atlantic (W. Atlantic) range** — a padded species-derived bbox clipped to a project-wide W. Atlantic frame.  
- **Caribbean Sea** — fixed `xlim_carib`, `ylim_carib`.  
- **U.S. Caribbean** — fixed `xlim_uscar`, `ylim_uscar`.

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
Input directory structure:
```r
project/
├─ data/
│  ├─ cmip6/                              # NetCDF exposure files (*.nc)
│  └─ species-distribution-shapefiles/    # Species polygons (*.shp + sidecars)
└─ outputs/
   └─ exposure-overlap/                   # PDFs & PNGs will be written here
```

Output directory structure:
```r
project/
└─ outputs/
└─ exposure-overlap-12panel/
└─ <species-slug>/
├─ <species-slug>_Distribution-Anomalies.pdf # 1 page per exposure factor
├─ <species-slug>_Exposure-Overlap-12panel.pdf # 1 page per exposure factor
├─ Distribution-Anomalies/
│ └─ <species-slug>_Distribution-Anomalies<exp>.png
└─ Exposure-Overlap-12panel/
└─ <species-slug>_Exposure-Overlap<exp>.png
```

---

# Workflow and functions

#### 1) Outer loop: per species shapefile
**Function:** `process_species(sp_file)`

**What it does**
- **Read + validate geometry:** `sf::st_read()` then `st_make_valid()`.
- **Standardize names:** Get a display name with `name_from_shp()` and a file-safe slug via `species_name_clean()`.
- **Compute plotting bbox:** `bbox_with_pad(sp, pad = 0.05)`; clip to the W. Atlantic project window if needed.
- **Prepare outputs:** Create `.../<species-slug>/` and open two multi-page PDFs using `safe_open_pdf()`:
  - `.../<species-slug>_Distribution-Anomalies.pdf`
  - `.../<species-slug>_Exposure-Overlap.pdf`
- **Hold devices open** for the inner loop to append one page per exposure factor.


---

#### 2) Inner loop: per exposure factor (NetCDF anomaly file)
**Driver:** a `for` loop over `nc_files` inside `process_species()`.

**Steps (per exposure)**
1. **Identify exposure:** `exp_name <- exp_name_from_nc(nc_path)`.  
2. **Load anomaly raster:**  
   - `anom <- terra::rast(nc_path, sub = "anomaly")`  
   - Replace fill values with `NA`, `rotate()` from 0–360° to −180–180°.  
3. **Crop by domains:**  
   - `anom_range` (species bbox, W. Atlantic), `anom_carib`, `anom_uscar`.  
4. **Build Figure A (2 × 3):**  
   - Species maps: `plot_distribution()` (W. Atlantic & Caribbean).  
   - Anomaly maps: `plot_anomalies_gg()` for Global, W. Atlantic, Caribbean, U.S. Caribbean.  
   - Assemble with `patchwork`, enlarge first row, add panel tags A–F.  
5. **Mask by species footprint:**  
   - Reproject species to raster CRS; `rasterize(..., cover=TRUE)` to include edge cells.  
   - `anom_masked` for W. Atlantic; then crop to Caribbean and U.S. Caribbean.  
6. **Shared color limits:**  
   - Compute `fill_limits` from masked W. Atlantic values to align legends and x-axes across panels.  
7. **Build Figure B (3 × 3):**  
   - Overlap maps: `overlap_range()` (W. Atlantic, Caribbean, U.S. Caribbean).  
   - Histograms: `anom_histogram_gg()` (all three domains).  
   - Categorical summaries: `anom_summary_bars()` (all three domains).  
   - Remove duplicate legends; extract a single legend with **cowplot**; assemble A–I.  
8. **Write outputs:**  
   - Append **one page** to each open PDF device (A → Distribution–Anomalies; B → Exposure–Overlap).  
   - Save **PNGs** for both figures under the species directory.

### These use five plotting functions
#### Distribution and anomalies maps:
- **`plot_distribution(sp, xlim, ylim, title)`** — species range maps with coastlines.  
- **`plot_anomalies_gg(anom, exp_name, extent, carib_box, uscar_box)`** — tidy anomaly maps (global & subdomains).
#### Exposure overlap plots:
- **`overlap_range(df, species_name, exp_name, domain, xlim, ylim, ...)`** — anomaly tiles within species mask + graticules.  
- **`anom_histogram_gg(r, fill_limits)`** — distribution of masked anomaly values (bin-colored, 0-line).  
- **`anom_summary_bars(r, species_name, exp_name, domain)`** — 5-bin categorical summary with counts and % labels.

### Reproducibility and comparability notes
- **Color limits fixed per exposure:** Use `fill_limits` so all panels for a given exposure share the same scale (improves cross-domain comparison). Consider project-wide fixed limits if you plan to compare across species, too.  
- **Validate CRS early:** Ensure species polygons and rasters align (`st_transform()` to `crs(anom)` before rasterize/mask).  
- **Clip early:** Cropping to domains before masking keeps rasters small and speeds up plotting.  

---

# Code description

## Utility Functions
Helper functions used throughout the exposure–overlap analysis. They make it easier to extract clean names, handle bounding boxes, and safely open PDF devices.

### `name_from_shp`
Extracts a **clean species name** from a shapefile path.  
- Removes file extensions (`.shp`)  
- Inserts spaces between CamelCase (e.g., `AtlanticHerring` → `Atlantic Herring`)  
- Replaces underscores, dashes, and dots with spaces. Trims and normalizes spacing  

### `exp_name_from_nc()`
Extracts the exposure factor code from a NetCDF filename. Everything before the first underscore (_) is returned.

### `bbox_with_pad()`
Calculates a bounding box (xlim, ylim) for an sf object (species polygon). Adds a proportional padding (default = 5%) so that plots have margins around the species range

### `safe_open_pdf()`
Opens a PDF graphics device safely.
- Uses cairo_pdf if available (better text rendering, avoids font issues).
- Falls back to base pdf() if Cairo is unavailable.
- Needed for writing multi-page PDFs of figures.

### `species_name_clean()`
Formats species names for filenames by replacing spaces with hyphens.

---

## Plotting Functions
Documents the plotting utilities used to visualize species ranges, climate anomalies, and their overlap. Each function is designed to be composable, so figures can be combined with **patchwork** or embedded into multi-panel PDFs.

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

### `lmhv_histogram_base()`
**Description:**  Computes and (optionally) draws an HMS-style histogram of standardized anomalies.  
Classifies values into Low (L), Moderate (M), High (H), and Very High (V) categories and returns a summary object that can be reused for categorical bar plots or ggplot wrappers.  

#### Arguments:
- **anom_masked** *(SpatRaster or Raster\**)*: Species-masked anomaly raster for a given domain.  
- **species_name, exp_name, domain** *(character, optional)*: Labels used for external tracking; not added to the base plot.  
- **mar, mgp** *(numeric vectors)*: Graphical parameters for margins and axis spacing.  
- **cex_axis, cex_lab, cex_main** *(numeric)*: Text expansion for axes, labels, and titles.  
- **do_plot** *(logical)*: If `TRUE`, produces a histogram plot; if `FALSE`, only computes and returns summary statistics.  

#### Notes:
- Breaks are set at **0.25 σ intervals** from min to max.  
- Colors: green (L), yellow (M), orange (H), red (V).  
- Returns an object of class `"lmhv_hist_summary"`, containing proportions (`Lp`, `Mp`, `Hp`, `Vp`) and a weighted exposure mean (1–4 scale).  
- Stops with an error if no finite values are found in the raster.  

---

### `lmhv_barplot_base()`
**Description:**  Draws a categorical LMHV bar plot (L, M, H, V) from a histogram summary, showing relative proportions and annotating the mean exposure score when available.  

#### Arguments:
- **lmhv_summary** *(lmhv_hist_summary)*: Object returned by `lmhv_histogram_base()`.  
- **mar, mgp** *(numeric vectors)*: Plot margins and axis spacing.  
- **cex_axis, cex_lab, cex_main** *(numeric)*: Text expansion factors.  
- **do_plot** *(logical)*: If `TRUE`, produces the barplot; otherwise returns the summary invisibly.  

#### Notes:
- Fixed bar colors: green (L), yellow (M), orange (H), red (V).  
- Bars are scaled to 0–1 (proportional).  
- If `exp_mean` is finite, it is shown as a numeric label near the top left.  

---

### `lmhv_histogram_as_gg()` and `lmhv_barplot_as_gg()`
**Description:**  Wrapper functions that produce a ggplot-compatible grob of the LMHV histogram or LMHV barplot for composition with **cowplot** or **patchwork**.  

#### Arguments:
- **x** *(lmhv_hist_summary or Raster\**)*: Either a precomputed summary or a masked raster (in which case the summary is generated internally).  
- **species_name, exp_name, domain** *(character, optional)*: Metadata used only if a raster is supplied.  
- **mar, mgp, cex_axis, cex_lab, cex_main** *(numeric)*: Aesthetic controls, as above.  

#### Notes:
- Converts the base histogram into a ggplot object via **`ggplotify::as.ggplot()`**.  
- Adds outer padding to prevent clipping of axis labels in patchwork layouts.  
- Best practice: precompute with `lmhv_histogram_base(do_plot = FALSE)` and pass the summary.  

---

### `lmhv_hist_draw_base()` and `lmhv_bar_draw_base()`
**Description:**  Low-level helpers that directly draw histograms and barplots from an `lmhv_hist_summary` without recomputation or side effects.  
#### Arguments:
- **summary** *(lmhv_hist_summary)*: Input summary object.  
- **mar, mgp, cex_axis, cex_lab, cex_main** *(numeric)*: Plot aesthetics.  
#### Notes:
- These functions are used internally by the `*_as_gg()` wrappers.  
- You would rarely call them directly, unless embedding base plots manually in a graphics device.  

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

---

## Inner loop: Process Exposure Factor
**Description:** Iterates over CMIP6 NetCDF anomaly files for a single species, deriving per-exposure figures. For each exposure factor, it: (1) reads and prepares the anomaly raster, (2) crops to multiple domains, (3) builds a 2×3 distribution/anomaly map panel, (4) computes species-footprint masks and overlap products for three domains, (5) assembles a 3×3 overlap panel with a single shared legend, and (6) appends one page to each open PDF while also exporting standalone PNGs.

**Notes about performance & reproducibility**
  - Cropping early (`crop()`) keeps rasters small; masking with `cover=TRUE` ensures partial cells are included along polygon edges.
  - Fixing `fill_limits` per exposure allows a **like-for-like** comparison across domains.
  - Ensures `sf` geometries are valid upstream (`st_make_valid`) and that CRS mismatches are resolved before masking.

**Exposure metadata & anomaly read**
  - Derives exposure factor code from filename via `exp_name_from_nc(nc_path)`.
  - Reads the **anomaly** layer with `terra::rast(nc_path, sub = "anomaly")`.
  - Replaces large fill values (`> 1e19`) with `NA` and calls `rotate()` to convert 0–360° longitudes to −180–180° for consistent cropping.

**Domain cropping (three views)**
  - **W. Atlantic species range:** `anom_range <- crop(anom, extent(c(lims$xlim, lims$ylim)))`
  - **Caribbean Sea:** `anom_carib <- crop(anom, ext(xlim_carib, ylim_carib))`
  - **U.S. Caribbean:** `anom_uscar <- crop(anom, ext(xlim_uscar, ylim_uscar))`

---

### **Distribution + anomaly maps (2×3)**
  - `p1, p2`: Species polygons over **world land** via `plot_distribution()` for the species bbox and the Caribbean bbox (larger title on `p1`).
  - `p3–p6`: Anomaly maps via `plot_anomalies_gg()` for Global, W. Atlantic, Caribbean, and U.S. Caribbean; optional domain boxes drawn.
  - Arrange as `(p1|p2)/(p3|p4)/(p5|p6)` with **`patchwork`**; first row height boosted (`plot_layout(heights = c(1.5, 1, 1))`) and panel tags `A–F`.
  
**Vectorization & masking (species overlap)**
  - Align CRS: `sp <- st_transform(sp, crs(anom_range))`, then `spv <- terra::vect(sp)`.
  - Build **cover** mask with `rasterize(spv, anom_range, field = 1, background = NA, cover = TRUE)`.
  - Apply mask: `anom_masked <- mask(anom_range, mask_cover)`.
  - Crop masked raster to subdomains: `anom_masked_carib <- crop(..., carib_ext)`, `anom_masked_uscar <- crop(..., uscar_ext)`.

**Shared color limits**
  - Convert masked range to data frame and compute global min/max:
    - `df_mask <- as.data.frame(anom_masked, xy = TRUE, na.rm = TRUE)`
    - `fill_limits <- c(min(df_mask$anomaly), max(df_mask$anomaly))`
  - Reuse `fill_limits` across maps/histograms to keep palettes and x-axes consistent.

---

### Overlap maps (3×3 with shared legend)
  - **Maps:** `overlap_range()` for W. Atlantic, Caribbean, and U.S. Caribbean (`main_title = "domain"`).
  - **Histograms:** `anom_histogram_gg()` for masked range and subdomains using `fill_limits`.
  - **Categorical summaries:** `anom_summary_bars()` for the same three domains.
  - Remove legends on all but one, extract a single legend via `cowplot::get_legend()`.
  - Assemble rows A–I with `cowplot::plot_grid()` and combine with the legend (ncol = 2, narrow legend column).

---

### Write out plots
- **Outputs (per exposure)**
  - **Append PDF pages:**  
    - `dev_over` ← print the 3×3 overlap panel.  
    - `dev_dist` ← print the 2×3 distribution/anomaly panel.
  - **Export PNGs:**  
    - `species_dir/Distribution-Anomalies/<slug>_Distribution-Anomalies_<exp>.png` (≈10×10 in, 200 dpi).  
    - `species_dir/Exposure-Overlap/<slug>_Exposure-Overlap_<exp>.png` (≈11×11 in, 300 dpi).
      
# References
Code and methods adapted in part from Loughren et al. HMS CVA:

Loughran, C. E., Hazen, E. L., Brodie, S., Jacox, M. G., Whitney, F. A., Payne, M. R., et al. (2025). A climate vulnerability assessment of highly migratory species in the Northwest Atlantic Ocean. *PLOS Climate*, 4(8), e0000530. https://doi.org/10.1371/journal.pclm.0000530  

### Additional references
- Crozier, L. G., Siegel, J., & Finnegan, W. B., et al. (2019). Climate vulnerability assessment for Pacific salmon and steelhead in the California Current. *PLOS ONE*, 14(7), e0217711. https://doi.org/10.1371/journal.pone.0217711  
- Frawley, T., Provost, M., Bellquist, L., Ben-Aderet, N., Blondin, H., Brodie, S., et al. (2025). A collaborative climate vulnerability assessment of California marine fishery species. *PLOS Climate*, 4(2), e0000574. https://doi.org/10.1371/journal.pclm.0000574  
- Hare, J. A., Morrison, W. E., Nelson, M. W., Stachura, M. M., Teeters, E. J., Griffis, R. B., et al. (2016). A vulnerability assessment of fish and invertebrates to climate change on the Northeast U.S. Continental Shelf. *PLOS ONE*, 11(2), e0146756. https://doi.org/10.1371/journal.pone.0146756  
- McClure, M. M., Haltuch, M. A., Berger, A. M., Berger, C., Branch, T. A., Brodziak, J., et al. (2023). Climate vulnerability assessment of species in the California Current Large Marine Ecosystem. *Frontiers in Marine Science*, 10, 1111111. https://doi.org/10.3389/fmars.2023.1111111  
- Morrison, W. E., Nelson, M. W., Howard, J. F., Teeters, E. J., Hare, J. A., Griffis, R. B., & Pugliese, R. (2015). Methodology for assessing the vulnerability of marine fish and shellfish species to a changing climate. *NOAA Technical Memorandum NMFS-OSF-3*. U.S. Department of Commerce. https://doi.org/10.7289/V5TM782J
- Craig, J. K., et al. (2025). A climate vulnerability assessment for the U.S. South Atlantic. https://doi.org/10.1371/journal.pclm.0000543
