# Caribbean CVA – Generating Species Distribution Maps

This script loops through all species distribution shapefiles in `./data/species-distribution-shapefiles/` and produces standardized PNG maps showing each species' range within a Caribbean bounding box. The loop automates map creation for any number of shapefiles in the directory.

- All maps share the same bounding box, symbology, and formatting, making them visually comparable.
- Utility functions make it easy to adapt the script to other projects or naming conventions.
- By saving maps as PNGs in a standard folder, outputs can be tracked and displayed directly on GitHub.

---

## Defining the geographic extent

The user sets a **Caribbean bounding box** based on project needs:

- **West:** −92° (Yucatán Channel buffer)  
- **East:** −57° (just east of Barbados)  
- **South:** 6°N (Panama/Colombia coast buffer)  
- **North:** 27.8°N (covers Bahamas)

This bounding box is stored both as:

- Numeric `xlim_carib` / `ylim_carib` vectors for ggplot cropping.
- An `sf` polygon (`carib_bbox_sf`) for potential masking operations.

---

## Loading basemap data and utility functions

The code uses `rnaturalearth` to load a medium-scale global land polygon (`world`), which is used as the basemap for all plots.

Two helper functions streamline file and label handling:

1. **`name_from_shp(path)`**  
   - Extracts the base filename from a shapefile path and cleans it up by inserting spaces in CamelCase, replacing underscores/hyphens, and trimming whitespace.  
   - Example: `"AtlanticHerring.shp"` → `"Atlantic Herring"`.

2. **`slugify(name)`**  
   - Converts a name into a lowercase, hyphen-separated string safe for filenames.  
   - Example: `"Atlantic Herring"` → `"atlantic-herring"`.

---

## Looping through species

The code searches the `./data/species-distribution-shapefiles/` directory for all `.shp` files.  

The loop conducts the following for each shapefile:

1. **Identify the species name** using `name_from_shp()`.
2. **Read the shapefile** into an `sf` object, ensure it’s valid, and reproject to WGS84 (`EPSG:4326`).
3. **Build the plot** (`pA`) with:
   - Land polygons (`world`) as a gray basemap.
   - The species polygon in semi-transparent black.
   - Coordinate cropping to the Caribbean bounding box.
   - A title `(A) Spatial distribution: [Species Name]`.
   - Minimal theme with black axis labels and visible ticks.
4. **Save as PNG**:
   - Output folder: `./outputs/disbribution-maps/`
   - Filename: slugified species name, e.g., `atlantic-herring.png`.
   - Resolution: 400 dpi, size: 8 × 5 inches.

---

## 7. Output

After running, the script produces one PNG file per shapefile in `./outputs/disbribution-maps/`.

All of the generated maps are here:  
[outputs/disbribution-maps on GitHub](https:/)

--- 
## Example Maps
![Atlantic Herring](https://github.com/holden-harris/Caribbean-CVA/blob/main/outputs/disbribution-maps/atlantic-herring.png?raw=true)
![King Mackerel](https://github.com/holden-harris/Caribbean-CVA/blob/main/outputs/disbribution-maps/king-mackerel.png)
![Misty Grouper](https://github.com/holden-harris/Caribbean-CVA/blob/main/outputs/disbribution-maps/misty-grouper.png)


