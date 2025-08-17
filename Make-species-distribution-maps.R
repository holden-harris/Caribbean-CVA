##------------------------------------------------------------------------------
## Org:     Harris Analytics and Research LLC | Isla Mar
## Project: Caribbean CVA
## Contact: Holden Earl Harris | holden.earl.harris@gmail.com
## Code:    Make Distribution Maps
##------------------------------------------------------------------------------

##------------------------------------------------------------------------------
# Clear environment and load packages
rm(list = ls()); gc()
library(dplyr)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(raster)

## -----------------------------------------------------------------------------
## Configuration
species_list <- read.csv("./data/master-lists/species-list.csv")

## -----------------------------------------------------------------------------
## Set geographic extent (bounding box)
xlim_carib <- c(-92, -57)
ylim_carib <- c(  6,  27.8)

## West (Yucatán Channel buffer): −92°
## East (just east of Barbados): −57°
## South (Panama/Colombia coast buffer): 6°N
## North (covers all Bahamas): 27.8°N
## If you want a tighter box (less open Atlantic/Gulf water) use xlim = c(-90, -58) and ylim = c(7, 28).

carib_bbox_sf <- sf::st_as_sfc( ## sf polygon for masking/overlays
  sf::st_bbox(c(xmin = xlim_carib[1], ymin = ylim_carib[1],
                xmax = xlim_carib[2], ymax = ylim_carib[2]),
              crs = 4326)
)

## Pull marmap bathymetry 
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

## -----------------------------------------------------------------------------
## Functions to pull names from files
name_from_shp <- function(path) { # Base path -> "Nice Name"
  base <- tools::file_path_sans_ext(basename(path))  # insert space between lower->Upper (CamelCase), e.g., "AtlanticHerring" -> "Atlantic Herring"
  base <- gsub("(?<=[a-z])(?=[A-Z])", " ", base, perl = TRUE)  # replace _, -, . with spaces
  base <- gsub("[_\\.\\-]+", " ", base)  # squeeze multiple spaces, trim
  base <- gsub("\\s+", " ", trimws(base))   # Title Case (keeps Genus species looking right)
  base
}

slugify <- function(name) { ## "Nice Name" -> "nice-name" (for filenames/ids)
  out <- tolower(name)
  out <- gsub("[^a-z0-9]+", "-", out)
  gsub("(^-|-$)", "", out)
}

##------------------------------------------------------------------------------
## Loop through shape files

## Pull file names from folder "species-distribtion-shapefiles"
shp_dir <- "./data/species-distribution-shapefiles/"
shp_files <- list.files(shp_dir, pattern = "\\.shp$", full.names = TRUE, recursive = FALSE)
if (length(shp_files) == 0L) stop("No .shp files found in: ", shp_dir)

## Run loop
for (i in 1:length(shp_files)){
# for (i in 1:7){ ## Use for testing
  ## Set species
  shp <- shp_files[i]
  species_name <- name_from_shp(shp)
  print(paste("Processsing", species_name))
  
  ## Species polygon
  species <- sf::st_read(shp, quiet = TRUE) |>
    sf::st_make_valid() |>
    sf::st_transform(4326)
  
 ## Make maps for species distribution ----------------------------------------
 pA <- ggplot() +  
    geom_sf(        ## draw land polygons for geographic context
      data = world,
      fill = "gray80",  ## fill color for land
      color = "gray30", ## separates countries
      linewidth = 0.5
    ) +
    geom_sf(        ## draw the species distribution polygon
      data = species,
      fill = "black", color = "black",  ## solid black fill  
      alpha = 0.5                       ## but semi-transparent so bathy/coastlines show through
    ) +
    coord_sf(       ## set the map window (crop) using Caribbean bbox
      xlim = xlim_carib, ylim = ylim_carib, expand = FALSE  ## no padding around the bbox
    ) +
    labs(           ## titles and axis labels
      title = paste("(A) Spatial distribution:", species_name),
    ) +
    theme_minimal() + ## clean base theme
    theme(                               
      axis.text  = element_text(color = "black"),   ## axis tick labels black
      axis.title = element_text(color = "black"),   ## axis title black (if not blank)
      axis.ticks   = element_line(color = "gray90") ## show axis ticks (theme_minimal hides them by default)
    )
  
  ## Write out plot -----------------------------------------------------------
  out_dir <- "./outputs/disbribution-maps/"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_name <- paste0(out_dir, slugify(species_name), ".png")
  png(out_name, width = 8, height = 5, units = "in", res = 400)
  plot(pA)
  dev.off()
  }
