##------------------------------------------------------------------------------
## Org:     Harris Analytics and Research LLC | Isla Mar
## Project: Caribbean CVA
## Contact: Holden Earl Harris | holden.earl.harris@gmail.com
## Code:    Exposure analyses
##------------------------------------------------------------------------------

##------------------------------------------------------------------------------
## Clear environment and load packages
rm(list = ls()); gc()
library(dplyr)
library(ncdf4)
library(raster)    
library(terra)
library(sp)
library(sf)
library(grDevices)
library(viridisLite) 
library(ggplot2)

##------------------------------------------------------------------------------
## Set directory paths
exp_dir     <- "./data/cmip6/"
spp_dir     <- "./data/species-distribution-shapefiles/"
out_dir     <- "./outputs/"

## -----------------------------------------------------------------------------
## Set geographic extent (bounding box)
xlim_carib <- c(-92, -57)
ylim_carib <- c(  6,  27.8)

## -----------------------------------------------------------------------------
## Read in exposure anomoly map
f <- file.path(exp_dir, "o200_1985-2014_2020-2049.nc") ## set path
anom    <- rast(f, sub = "anomaly") ## one layer
anom[anom > 1e19] <- NA ## Fix fill values
anom        <- rotate(anom) ## NetCDF longitudes are 0–360 so need to rotate to match our extent, which is -180-180
anom_carib <- crop(anom, ext(xlim_carib, ylim_carib)) ## Crop to Caribbean extent

## Plotcheck -------------------------------------------------------------------
## Define HMS extent (note: western longitudes are negative)
xlim_nwa <- c(-99, 40)
ylim_nwa <- c(-20, 72)
anom_nwa <- crop(anom, extent(c(xlim_nwa, ylim_nwa)))

## -----------------------------------------------------------------------------
## Plotcheck regions
par(mfrow = c(2,2))

# A) Global anomalies
plot(anom,
     main = "A) CMIP anomalies",
     xlab = "", ylab = "",
     col = turbo(100))
map("world", add = TRUE, col = "grey20", lwd = 0.6)
rect(xlim_carib[1], ylim_carib[1], xlim_carib[2], ylim_carib[2], border = "purple", lwd = 2) ## Box Caribbean

# C) HMS domain crop
plot(anom_nwa,
     main = "C) W. Atlantic",
     xlab = "", ylab = "",
     col = turbo(100))
map("world", add = TRUE, col = "grey20", lwd = 0.6)
rect(xlim_carib[1], ylim_carib[1], xlim_carib[2], ylim_carib[2], border = "purple", lwd = 2) ## Box Caribbean

# B) Caribbean crop
plot(anom_carib,
     main = "C) Caribbean",
     xlab = "", ylab = "",
     col = turbo(100))
map("world", add = TRUE, col = "grey20", lwd = 0.6)
rect(xlim_carib[1], ylim_carib[1], xlim_carib[2], ylim_carib[2], border = "purple", lwd = 2) ## Box Caribbean
par(mfrow = c(1,1)) ## Reset plotting layout

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



## -----------------------------------------------------------------------------
## Map overlap
shp_files <- list.files(spp_dir, pattern = "\\.shp$", full.names = TRUE, recursive = FALSE)
if (length(shp_files) == 0L) stop("No .shp files found in: ", spp_dir)

## Run loop
#for (i in 1:length(shp_files)){
# for (i in 1:1){ ## Use for testing
  i = 1
  
  ## Set species
  sp_file <- shp_files[i]
  species_name <- name_from_shp(sp_file)
  print(paste("Processsing", species_name))
  
  sp      <- st_read(sp_file, quiet = TRUE) |> st_make_valid()
  sp      <- st_transform(sp, crs(anom_carib))  ## match raster CRS
  spv     <- vect(sp)

  ## Rasterize polygon to the anomaly grid (include partial cells)
  ## Note: cover=TRUE returns fraction of each cell covered by the polygon
  mask_cover <- rasterize(spv, anom_carib, field = 1, background = NA, cover = TRUE)
  mask01  <- mask_cover; mask01[mask01 > 0] <- 1 ## convert to 1/NA like HMS
  
  ## Mask anomalies by species range
  anom_masked <- mask(anom_carib, mask01)
  
  ##
  ## anom_masked is a SpatRaster/Raster* (masked by species footprint). Convert it to df for ggplot. 
  df_mask <- as.data.frame(anom_masked, xy = TRUE, na.rm = TRUE)
  names(df_mask)[3] <- "anom"   # third column is the values
  
  pC <- ggplot() +
    geom_raster(data = df_mask, aes(x = x, y = y, fill = anom)) +
    geom_sf(data = world, fill = "gray80", color = "gray30", linewidth = 0.5) +
    scale_fill_gradientn(colors = turbo(100), name = "SST st anom") +
    coord_sf(xlim = xlim_carib, ylim = ylim_carib, expand = FALSE) +
    labs(title = paste("(C) Distribution / Exposure Factor Overlap:", species_name), x = "", y = "") +
    theme_minimal() +
    theme(                               
      axis.text  = element_text(color = "black"),   ## axis tick labels black
      axis.title = element_text(color = "black"),   ## axis title black (if not blank)
      axis.ticks   = element_line(color = "gray90") ## show axis ticks (theme_minimal hides them by default)
    ); pC
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  plot(anom_masked,
       main = "",
       xlab = "", ylab = "",
       col = turbo(100))
  map("world", add = TRUE, col = "grey20", lwd = 0.6)
  
  
  
  
  
  plot(anom_masked)
  
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

  ## Make maps for species distribution ----------------------------------------
  pC <- ggplot() +  
    geom_sf(        ## draw land polygons for geographic context
      data = world,
      fill = "gray80",  ## fill color for land
      color = "gray30", ## separates countries
      linewidth = 0.5
    ) +
    geom_raster(        ## draw the species distribution polygon
      data = anom_masked,
    ) +
    coord_sf(       ## set the map window (crop) using Caribbean bbox
      xlim = xlim_carib, ylim = ylim_carib, expand = FALSE  ## no padding around the bbox
    ) +
    labs(           ## titles and axis labels
      title = paste("(C) Distribution / Exposure factor", species_name),
    ) +
    theme_minimal() + ## clean base theme
    theme(                               
      axis.text  = element_text(color = "black"),   ## axis tick labels black
      axis.title = element_text(color = "black"),   ## axis title black (if not blank)
      axis.ticks   = element_line(color = "gray90") ## show axis ticks (theme_minimal hides them by default)
    ); pC  
