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

### Outputs

After running, the script produces one PNG file per shapefile in `./outputs/disbribution-maps/`.

PNG fles of the generated maps are available here: [outputs/disbribution-maps](https:/)

# Species Distribution Maps

<table>
<tr>
<td><img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/outputs/disbribution-maps/atlantic-herring.png?raw=true" alt="Atlantic Herring" width="400"/></td>
<td><img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/outputs/disbribution-maps/ballyhoo.png?raw=true" alt="Ballyhoo" width="400"/></td>
</tr>
<tr>
<td><img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/outputs/disbribution-maps/blue-runner.png?raw=true" alt="Blue Runner" width="400"/></td>
<td><img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/outputs/disbribution-maps/diadema.png?raw=true" alt="Diadema" width="400"/></td>
</tr>
<tr>
<td><img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/outputs/disbribution-maps/dolphinfish.png?raw=true" alt="Dolphinfish" width="400"/></td>
<td><img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/outputs/disbribution-maps/gray-angelfish.png?raw=true" alt="Gray Angelfish" width="400"/></td>
</tr>
<tr>
<td><img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/outputs/disbribution-maps/hogfish.png?raw=true" alt="Hogfish" width="400"/></td>
<td><img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/outputs/disbribution-maps/king-mackerel.png?raw=true" alt="King Mackerel" width="400"/></td>
</tr>
<tr>
<td><img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/outputs/disbribution-maps/lane-snapper.png?raw=true" alt="Lane Snapper" width="400"/></td>
<td><img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/outputs/disbribution-maps/misty-grouper.png?raw=true" alt="Misty Grouper" width="400"/></td>
</tr>
<tr>
<td><img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/outputs/disbribution-maps/mutton-snapper.png?raw=true" alt="Mutton Snapper" width="400"/></td>
<td><img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/outputs/disbribution-maps/nassau-grouper.png?raw=true" alt="Nassau Grouper" width="400"/></td>
</tr>
<tr>
<td><img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/outputs/disbribution-maps/queen-conch.png?raw=true" alt="Queen Conch" width="400"/></td>
<td><img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/outputs/disbribution-maps/queen-snapper.png?raw=true" alt="Queen Snapper" width="400"/></td>
</tr>
<tr>
<td><img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/outputs/disbribution-maps/queen-triggerfish.png?raw=true" alt="Queen Triggerfish" width="400"/></td>
<td><img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/outputs/disbribution-maps/rainbow-parrotfish.png?raw=true" alt="Rainbow Parrotfish" width="400"/></td>
</tr>
<tr>
<td><img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/outputs/disbribution-maps/red-grouper.png?raw=true" alt="Red Grouper" width="400"/></td>
<td><img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/outputs/disbribution-maps/redhind.png?raw=true" alt="Redhind" width="400"/></td>
</tr>
<tr>
<td><img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/outputs/disbribution-maps/sea-cucumber.png?raw=true" alt="Sea Cucumber" width="400"/></td>
<td><img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/outputs/disbribution-maps/silk-snapper.png?raw=true" alt="Silk Snapper" width="400"/></td>
</tr>
<tr>
<td><img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/outputs/disbribution-maps/spiny-lobster.png?raw=true" alt="Spiny Lobster" width="400"/></td>
<td><img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/outputs/disbribution-maps/stoplight-parrotfish.png?raw=true" alt="Stoplight Parrotfish" width="400"/></td>
</tr>
<tr>
<td><img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/outputs/disbribution-maps/white-mullet.png?raw=true" alt="White Mullet" width="400"/></td>
<td><img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/outputs/disbribution-maps/yellowfin-grouper.png?raw=true" alt="Yellowfin Grouper" width="400"/></td>
</tr>
<tr>
<td><img src="https://github.com/holden-harris/Caribbean-CVA/blob/main/outputs/disbribution-maps/yellowtail-snapper.png?raw=true" alt="Yellowtail Snapper" width="400"/></td>
<td></td>
</tr>
</table>

