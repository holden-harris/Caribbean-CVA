## -----------------------------------------------------------------------------
## Org:     Harris Analytics and Research LLC | Isla Mar
## Project: Caribbean CVA
## Contact: Holden Earl Harris | holden@harris-analytics.com
## Code:    Extract quantitative exposure attribute scores
## Notes:   This script reuses the LMHV exposure scoring logic from the exposure
##          anomaly workflow, but writes a table instead of figures.
## -----------------------------------------------------------------------------

################################################################################
##
## Setup

## -----------------------------------------------------------------------------
## Clear environment and load packages
rm(list = ls()); gc()

library(dplyr)
library(terra)
library(sf)

## -----------------------------------------------------------------------------
## Set directory paths
exp_dir     <- "./data/cmip6/"
spp_dir     <- "./data/species-distribution-shapefiles/"
out_dir     <- "./outputs/final-scores-compiled/"

## -----------------------------------------------------------------------------
## Species shape files
shp_files <- list.files(spp_dir, pattern = "\\.shp$", full.names = TRUE, recursive = FALSE)
if (length(shp_files) == 0L) stop("No .shp files found in: ", spp_dir)
print(shp_files)

## -----------------------------------------------------------------------------
## CMIP exposure files
nc_files <- list.files(exp_dir, pattern = "\\.nc$", full.names = TRUE, recursive = FALSE)
if (length(nc_files) == 0L) stop("No .nc files found in: ", exp_dir)
print(nc_files)

## -----------------------------------------------------------------------------
## Set geographic extents (bounding boxes)

## Caribbean Sea
xlim_carib <- c(-92, -57)
ylim_carib <- c(6, 28)
carib_ext  <- c(xlim_carib, ylim_carib)

## U.S. Caribbean extent
xlim_uscar <- c(-69, -63.0)
ylim_uscar <- c(16, 20.0)
uscar_ext  <- c(xlim_uscar, ylim_uscar)

## Western Atlantic Ocean
xlim_nwa <- c(-99, -40)
ylim_nwa <- c(-5, 72)

exposure_factor_key <- c(
  "bs", "bt", "chl", "mld", "msstg", "o200", 
  "ph", "pp", "pr", "sso", "sss", "sst", "swsm")
  

################################################################################
##
## Functions

## -----------------------------------------------------------------------------
## Extract species name from shapefile name
## Example: "AtlanticHerring.shp" -> "Atlantic Herring"
name_from_shp <- function(path) {
  base <- tools::file_path_sans_ext(basename(path))
  base <- gsub("(?<=[a-z])(?=[A-Z])", " ", base, perl = TRUE)
  base <- gsub("[_\\.\\-]+", " ", base)
  base <- gsub("\\s+", " ", trimws(base))
  base
}

## -----------------------------------------------------------------------------
## Extract exposure factor code from nc filename
## Example: "o200_1985-2014_2020-2049.nc" -> "o200"
exp_name_from_nc <- function(path) {
  sub("_.*$", "", basename(path))
}

## -----------------------------------------------------------------------------
## Get bounding box from species distribution, with small padding
bbox_with_pad <- function(sf_obj, pad = 0.05) {
  bb <- st_bbox(sf_obj)
  dx <- as.numeric(bb["xmax"] - bb["xmin"])
  dy <- as.numeric(bb["ymax"] - bb["ymin"])
  
  list(
    xlim = c(bb["xmin"] - pad * dx, bb["xmax"] + pad * dx),
    ylim = c(bb["ymin"] - pad * dy, bb["ymax"] + pad * dy)
  )
}

## -----------------------------------------------------------------------------
## HMS-style LMHV summary function
## This is the same core logic previously used to compute exp_mean for plots
lmhv_histogram_base <- function(anom_masked,
                                species_name,
                                exp_name = "",
                                domain = "",
                                do_plot = FALSE) {
  
  ## Extract raster values
  vals <- terra::values(anom_masked, mat = FALSE)
  vals <- vals[is.finite(vals)]
  
  ## Return NA summary if no valid cells are present
  if (!length(vals)) {
    out <- list(
      Lp = NA_real_, Mp = NA_real_, Hp = NA_real_, Vp = NA_real_,
      exp_mean = NA_real_,
      tally_L = NA_integer_, tally_M = NA_integer_,
      tally_H = NA_integer_, tally_VH = NA_integer_
    )
    class(out) <- c("lmhv_hist_summary", "list")
    return(out)
  }
  
  ## Histogram breaks
  breaks <- seq(floor(min(vals)), ceiling(max(vals)), by = 0.25)
  
  ## Compute histogram without plotting
  h <- hist(vals, breaks = breaks, plot = FALSE)
  
  ## Midpoints and counts
  mids <- h$mids
  cnts <- h$counts
  
  ## LMHV bins
  L <- sum(cnts[mids >= -0.5 & mids <=  0.5], na.rm = TRUE)
  M <- sum(cnts[(mids < -0.5 & mids >= -1.5) | (mids > 0.5 & mids <= 1.5)], na.rm = TRUE)
  H <- sum(cnts[(mids < -1.5 & mids >= -2.0) | (mids > 1.5 & mids <= 2.0)], na.rm = TRUE)
  V <- sum(cnts[mids < -2.0 | mids > 2.0], na.rm = TRUE)
  
  tot <- L + M + H + V
  
  out <- list(
    Lp       = if (tot > 0) L / tot else NA_real_,
    Mp       = if (tot > 0) M / tot else NA_real_,
    Hp       = if (tot > 0) H / tot else NA_real_,
    Vp       = if (tot > 0) V / tot else NA_real_,
    exp_mean = if (tot > 0) ((L * 1) + (M * 2) + (H * 3) + (V * 4)) / tot else NA_real_,
    tally_L  = L,
    tally_M  = M,
    tally_H  = H,
    tally_VH = V
  )
  
  class(out) <- c("lmhv_hist_summary", "list")
  out
}


################################################################################
##
## Core processing function

## -----------------------------------------------------------------------------
## Process one species and return a long table of attribute scores

process_species_scores <- function(sp_file) {
  
  ## Read species polygon
  sp <- sf::st_read(sp_file, quiet = TRUE) |> sf::st_make_valid()
  
  ## Species name
  species_name <- name_from_shp(sp_file)
  
  cat("\n------------------------------------------------------------\n")
  cat("Processing species:", species_name, "\n")
  
  ## Get bounding box and clip to Western Atlantic analysis domain
  lims <- bbox_with_pad(sp, pad = 0.05)
  
  lims$xlim <- c(
    max(lims$xlim[1], xlim_nwa[1]),
    min(lims$xlim[2], xlim_nwa[2])
  )
  
  lims$ylim <- c(
    max(lims$ylim[1], ylim_nwa[1]),
    min(lims$ylim[2], ylim_nwa[2])
  )
  
  ## Store results from each exposure factor
  score_list <- vector("list", length(nc_files))
  
  ## Loop through each exposure nc file
  for (i in seq_along(nc_files)) {
    
    nc_path  <- nc_files[i]
    exp_name <- exp_name_from_nc(nc_path)
    
    cat("  - Exposure", i, "of", length(nc_files), ":", exp_name, "\n")
    
    ## Read anomaly layer
    anom <- rast(nc_path, sub = "anomaly")
    
    ## Fix fill values
    anom[anom > 1e19] <- NA
    
    ## Rotate from 0-360 to -180 to 180 if needed
    anom <- rotate(anom)
    
    ## Crop to species regional extent within Western Atlantic
    anom_range <- crop(anom, ext(c(lims$xlim, lims$ylim)))
    
    ## Transform species polygon to raster CRS
    sp_proj <- st_transform(sp, crs(anom_range))
    sp_vect <- vect(sp_proj)
    
    ## Rasterize species footprint
    mask_cover <- rasterize(sp_vect, anom_range, field = 1, background = NA, cover = TRUE)
    
    ## Mask anomalies to species distribution
    anom_masked <- mask(anom_range, mask_cover)
    
    ## Crop masked anomalies to each reporting extent
    anom_masked_watl  <- anom_masked
    anom_masked_carib <- crop(anom_masked, ext(xlim_carib, ylim_carib))
    anom_masked_uscar <- crop(anom_masked, ext(xlim_uscar, ylim_uscar))
    
    ## Calculate LMHV summaries
    sum_watl  <- lmhv_histogram_base(anom_masked_watl,  species_name, exp_name, "Western Atlantic", do_plot = FALSE)
    sum_carib <- lmhv_histogram_base(anom_masked_carib, species_name, exp_name, "Caribbean Sea", do_plot = FALSE)
    sum_uscar <- lmhv_histogram_base(anom_masked_uscar, species_name, exp_name, "U.S. Caribbean", do_plot = FALSE)
    
    ## Save long-format output table
    score_list[[i]] <- tibble(
      stock_name                   = species_name,
      quantitative_exposure_factor = exp_name,
      spatial_extent  = c("Western Atlantic", "Caribbean Sea", "U.S. Caribbean"),
      attribute_score = c(sum_watl$exp_mean,  sum_carib$exp_mean,  sum_uscar$exp_mean),
      tally_L  = c(sum_watl$tally_L,  sum_carib$tally_L,  sum_uscar$tally_L),
      tally_M  = c(sum_watl$tally_M,  sum_carib$tally_M,  sum_uscar$tally_M),
      tally_H  = c(sum_watl$tally_H,  sum_carib$tally_H,  sum_uscar$tally_H),
      tally_VH = c(sum_watl$tally_VH, sum_carib$tally_VH, sum_uscar$tally_VH)
    )
  }
  
  ## Combine exposure rows for this species
  bind_rows(score_list)
}

################################################################################
##
## Run loop

## -----------------------------------------------------------------------------
## Test one species first
 i <- 1
 attribute_score_table_test <- process_species_scores(shp_files[i])
 print(attribute_score_table_test)

## -----------------------------------------------------------------------------
## Run all species with error handling
attribute_score_list <- vector("list", length(shp_files))

for (i in seq_along(shp_files)) {
  attribute_score_list[[i]] <- tryCatch(
    process_species_scores(shp_files[i]),
    error = function(e) {
      message("[ERROR] process_species_scores failed for: ", shp_files[i])
      message("        ", conditionMessage(e))
      NULL
    }
  )
}

## -----------------------------------------------------------------------------
## Combine all species into one table
attribute_score_table <- bind_rows(attribute_score_list)

## -----------------------------------------------------------------------------
## Add full exposure factor names
exposure_factor_key <- c(
  "bs"    = "Bottom salinity",
  "bt"    = "Bottom temperature",
  "chl"   = "Chlorophyll-a concentration",
  "mld"   = "Mixed layer depth",
  "msstg" = "Mean sea surface temperature gradient",
  "o200"  = "Oxygen at 200m",
  "ph"    = "Surface pH",
  "pp"    = "Primary production",
  "precip"    = "Precipitation",
  "sso"   = "Sea surface oxygen",
  "sss"   = "Sea surface salinity",
  "sst"   = "Sea surface temperature",
  "swsm"   = "Surface wind speed magnitude"
)

attribute_score_table <- attribute_score_table %>%
  mutate(
    full_names = exposure_factor_key[quantitative_exposure_factor]
  )

## Reorder
attribute_score_table <- attribute_score_table %>%
  select(
    stock_name,
    quantitative_exposure_factor,
    full_names,
    spatial_extent,
    attribute_score,
    tally_L, tally_M, tally_H, tally_VH
  )

## -----------------------------------------------------------------------------
## Add sum column
attribute_score_table <- attribute_score_table %>%
  mutate(
    attribute_score = round(attribute_score, 3),
    n_tallies = tally_L + tally_M + tally_H + tally_VH
  )

## -----------------------------------------------------------------------------
## Inspect final table
print(attribute_score_table, n = 40)

## -----------------------------------------------------------------------------
## Write output csv
out_file <- file.path(out_dir, "quantitative-exposure-attribute-scores-all.csv")
write.csv(attribute_score_table, out_file, row.names = FALSE)
cat("\nDone. Output written to:\n", out_file, "\n")
