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

# Assumes these objects already exist:
#   species     : sf polygon of species range (EPSG:4326)
#   world       : sf world countries for land context
#   bathy_df    : data.frame with columns x, y, z (bathymetry), e.g., from marmap::fortify.bathy()
#   xlim_carib  : c(lon_min, lon_max) numeric vector for Caribbean bbox
#   ylim_carib  : c(lat_min, lat_max) numeric vector for Caribbean bbox

pA <- ggplot() +  
#  geom_contour(   ## draw bathymetry as dashed contour lines
#    data = bathy_df,
#    aes(x = x, y = y, z = z),
#    breaks = c(-400, -800, -1600),  ## isobaths
#    color = "blue", alpha = 0.5,                ## light gray lines
#    linewidth = 0.3,                  ## thin lines (use linewidth for geoms in ggplot >=3.4)
#  ) +
  geom_sf(        ## draw land polygons for geographic context
    data = world,
    fill = "gray80",  # fill color for land
    color = "gray30", # separates countries
    linewidth = 0.5
  ) +
  geom_sf(        # draw the species distribution polygon
    data = species,
    fill = "black",                   # solid black fill...
    color = "black",                       # ...with no border stroke
    alpha = 0.5                       # ...but semi-transparent so bathy/coastlines show through
  ) +
  coord_sf(       # set the map window (crop) using Caribbean bbox
    xlim = xlim_carib, ylim = ylim_carib, expand = FALSE  # no padding around the bbox
  ) +
  labs(           # titles and axis labels
    title = "(A) Species or Stock Distribution",
  ) +
  theme_minimal() +                    # clean base theme
  theme(                               
    axis.text  = element_text(color = "black"),  # axis tick labels black
    axis.title = element_text(color = "black"),  # axis title black (if not blank)
    axis.ticks   = element_line(color = "gray90") # show axis ticks (theme_minimal hides them by default)
  ); pA  # print the plot
    