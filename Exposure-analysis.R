##------------------------------------------------------------------------------
# Clear environment and load packages
rm(list = ls()); gc()
library(dplyr)
library(sf)
library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(marmap)
library(raster)

## -----------------------------------------------------------------------------
## Configuration
species_list <- read.csv("./data/master-lists/species-list.csv")


## -----------------------------------------------------------------------------
## Set geographic extent (bounding box)
xlim_carib <- c(-92, -57)
ylim_carib <- c(  6,  29)

## West (Yucatán Channel buffer): −92°
## East (just east of Barbados): −57°
## South (Panama/Colombia coast buffer): 6°N
## North (covers all Bahamas): 29°N
## If you want a tighter box (less open Atlantic/Gulf water) use xlim = c(-90, -58) and ylim = c(7, 28).

carib_bbox_sf <- sf::st_as_sfc( ## sf polygon for masking/overlays
  sf::st_bbox(c(xmin = xlim_carib[1], ymin = ylim_carib[1],
                xmax = xlim_carib[2], ymax = ylim_carib[2]),
              crs = 4326)
)

## Pull marmap bathymetry (takes a min to pull the API, depending on size and resolution)
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
bathy <- marmap::getNOAA.bathy(lon1 = xlim_carib[1], lon2 = xlim_carib[2],
                               lat1 = ylim_carib[1], lat2 = ylim_carib[2], resolution = 1)
bathy_df <- as.data.frame(rasterToPoints(marmap::as.raster(bathy)))
names(bathy_df) <- c("x","y","z")

##------------------------------------------------------------------------------
## Load species and basemap

## species polygon
species <- sf::st_read("./data/species-distribution-shapefiles/AtlanticHerring.shp", quiet = TRUE) |>
  sf::st_make_valid() |>
  sf::st_transform(4326)


##------------------------------------------------------------------------------
## Make maps for species distribution

# Panel A: species distribution
pA <- ggplot() +
  geom_contour(data = bathy_df, aes(x = x, y = y, z = z),
               breaks = c(-200, -1000), color = "gray70",
               size = 0.3, linetype = "dashed") +
  geom_sf(data = world, fill = "gray90", color = "gray40", linewidth = 0.2) +
  geom_sf(data = species, fill = "gray20", color = NA, alpha = 0.95) +
  coord_sf(xlim = xlim_carib, ylim = ylim_carib, expand = FALSE) +
  labs(title = "(A) Species or Stock Distribution", x = "Longitude", y = "Latitude") +
  theme_minimal(); pA







# basemap & bathy



world <- ne_countries(scale = "medium", returnclass = "sf")
bathy <- getNOAA.bathy(lon1 = -100, lon2 = -45, lat1 = -10, lat2 = 35, resolution = 1)
bathy_raster <- marmap::as.raster(bathy)
crop_extent <- extent(-90, -55, -5, 35)
bathy_cropped <- crop(bathy_raster, crop_extent)
bathy_df <- as.data.frame(rasterToPoints(bathy_cropped))
colnames(bathy_df) <- c("x", "y", "z")

## Plot
ggplot() +
  geom_contour(data = bathy_df, aes(x = x, y = y, z = z),
               breaks = c(-200, -1000), color = "gray70", size = 0.3, linetype = "dashed") +
  #geom_sf(data = world, fill = "gray90", color = "gray40") +  # use full map
  geom_sf(data = species_shp, color = "red", size = 0.5, alpha = 0.7) +
  coord_sf(xlim = c(-100, -45), ylim = c(-5, 35), expand = FALSE) +
  labs(title = "Distribution of Atlantic Herring",
       x = "Longitude", y = "Latitude") +
  theme_minimal()
