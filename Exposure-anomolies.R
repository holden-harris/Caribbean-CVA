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

##------------------------------------------------------------------------------
## Set directory paths
nc_path     <- "./data/cmip6/"
species_shp <- "./data/species-distribution-shapefiles/"
out_dir     <- "./outputs/"

## -----------------------------------------------------------------------------
## Set geographic extent (bounding box)
xlim_carib <- c(-92, -57)
ylim_carib <- c(  6,  27.8)

## -----------------------------------------------------------------------------
f <- file.path(nc_path, "o200_1985-2014_2020-2049.nc") ## set path
anom              <- rast(f, sub = "anomaly") ## one layer
anom[anom > 1e19] <- NA ## Fix fill values
anom        <- rotate(anom) ## NetCDF longitudes are 0–360 so need to rotate to match our extent, which is -180-180
anom_carib        <- crop(anom, ext(xlim_carib, ylim_carib)) ## Crop to Caribbean extent
plot(anom_carib)


## Make color pallete that matches HMS -----------------------------------------
vals_for_bins <- values(anom_carib)  # or values(env_an_masked) if you've created that
breaks <- seq(
  floor(min(vals_for_bins, na.rm = TRUE)),
  ceiling(max(vals_for_bins, na.rm = TRUE)),
  by = 0.25
)

my_colors <- rep("red", length(breaks))  # HMS does length(breaks) here
my_colors[breaks >= -0.5 & breaks <=  0.5] <- "green"
my_colors[breaks <  -0.5 & breaks >= -1.5] <- "yellow"
my_colors[breaks >   0.5 & breaks <=  1.5] <- "yellow"
my_colors[breaks <  -1.5 & breaks >= -2.0] <- "orange"
my_colors[breaks >   1.5 & breaks <=  2.0] <- "orange"

## For raster plotting, col must be one shorter than breaks
cols_for_raster <- if (length(breaks) > 1) my_colors[-1] else my_color 


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


