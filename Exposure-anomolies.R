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
#library(ncdf4)
library(raster)    
library(terra)
library(sp)
library(sf)
#library(grDevices)
library(viridisLite) 
library(ggplot2)
library(rnaturalearth)

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

exp_name <- "o200"

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
maps::map("world", add = TRUE, col = "grey20", lwd = 0.6)
rect(xlim_carib[1], ylim_carib[1], xlim_carib[2], ylim_carib[2], border = "purple", lwd = 2) ## Box Caribbean

# C) HMS domain crop
plot(anom_nwa,
     main = "C) W. Atlantic",
     xlab = "", ylab = "",
     col = turbo(100))
maps::map("world", add = TRUE, col = "grey20", lwd = 0.6)
rect(xlim_carib[1], ylim_carib[1], xlim_carib[2], ylim_carib[2], border = "purple", lwd = 2) ## Box Caribbean

# B) Caribbean crop
plot(anom_carib,
     main = "C) Caribbean",
     xlab = "", ylab = "",
     col = turbo(100))
maps::map("world", add = TRUE, col = "grey20", lwd = 0.6)
rect(xlim_carib[1], ylim_carib[1], xlim_carib[2], ylim_carib[2], border = "purple", lwd = 2) ## Box Caribbean
par(mfrow = c(1,1)) ## Reset plotting layout

## -----------------------------------------------------------------------------
## Functions 
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

## Function to get bounding box from spatial distribution
bbox_with_pad <- function(sf_obj, pad = 0.05){
  bb <- st_bbox(sf_obj)
  dx <- as.numeric(bb["xmax"] - bb["xmin"])
  dy <- as.numeric(bb["ymax"] - bb["ymin"])
  list(xlim = c(bb["xmin"] - pad*dx, bb["xmax"] + pad*dx),
       ylim = c(bb["ymin"] - pad*dy, bb["ymax"] + pad*dy))
}

## -----------------------------------------------------------------------------
## Map overlap

## -----------------------------------------------------------------------------
## Set geographic extent (bounding box)
xlim_carib <- c(-92, -57)
ylim_carib <- c(  6,  27.8)
carib_ext  <- c(xlim_carib, ylim_carib)

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
  
  ## Read in and vectorize species
  sp      <- st_read(sp_file, quiet = TRUE) |> st_make_valid()
  lims <- bbox_with_pad(sp, pad = 0.05)   ## Crop exposure anomolies to species range
  anom_range <- crop(anom, extent(c(lims$xlim, lims$ylim)))
  sp      <- st_transform(sp, crs(anom_range))  ## match raster CRS
  spv     <- vect(sp)

  ## Make and apply cover mask
  ## Rasterize polygon to the anomaly grid (include partial cells)
  ## Note that cover=TRUE returns fraction of each cell covered by the polygon
  
  ## Mask cover species range (or W. Atlantic)
  mask_cover  <- rasterize(spv, anom_range, field = 1, background = NA, cover = TRUE)
  anom_masked <- mask(anom_range, mask_cover)
  
  ## Mask cover Caribbean
  #mask_cover_carib  <- rasterize(spv, anom_carib, field = 1, background = NA, cover = TRUE)
  #carib_poly <- as.polygons(carib_ext)     # SpatVector polygon from bbox
  #anom_masked_carib <- mask(anom_carib, carib_poly)
  anom_masked_carib <- crop(anom_masked, carib_ext)
  
  ## Make dataframe for masked values
  ## anom_masked is a SpatRaster/Raster* (masked by species footprint). Convert it to df for ggplot. 
  df_mask <- as.data.frame(anom_masked, xy = TRUE, na.rm = TRUE)
  names(df_mask)[3] <- "anom"   ## third column is the anomolies values
  
  df_mask_carib <- as.data.frame(anom_masked_carib, xy = TRUE, na.rm = TRUE)
  names(df_mask_carib)[3] <- "anom"   # third column is the values
  
  ## Plots ---------------------------------------------------------------------
  ## Pull worldmap
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") 
  
  ## Plot overlap -- Entire species range
  plot_overlap_range <- ggplot() +
    geom_tile(data = df_mask, aes(x = x, y = y, fill = anom)) +
    geom_sf(data = world, fill = "gray80", color = "gray30", linewidth = 0.5) +
    scale_fill_gradientn(colors = turbo(100), name = exp_name) +
    coord_sf(xlim = lims$xlim, ylim = lims$ylim, expand = FALSE) +
    labs(title = paste0(exp_name, " | ", species_name, " | W. Atlantic Range"), x = "", y = "") +
    theme_minimal() +
    theme(                               
      axis.text  = element_text(color = "black"),   ## axis tick labels black
      axis.title = element_text(color = "black"),   ## axis title black (if not blank)
      axis.ticks   = element_line(color = "gray90") ## show axis ticks (theme_minimal hides them by default)
    ); plot_overlap_range
  
  ## Plot overlap -- Caribbean Only
  plot_overlap_carib <- ggplot() +
    geom_tile(data = df_mask_carib, aes(x = x, y = y, fill = anom)) +
    geom_sf(data = world, fill = "gray80", color = "gray30", linewidth = 0.5) +
    scale_fill_gradientn(colors = turbo(100), name = exp_name) +
    coord_sf(xlim = xlim_carib, ylim = ylim_carib, expand = FALSE) +
    labs(title = paste0(exp_name, " | ", species_name, " | Caribbean"), x = "", y = "") +
    theme_minimal() +
    theme(                               
      axis.text  = element_text(color = "black"),   ## axis tick labels black
      axis.title = element_text(color = "black"),   ## axis title black (if not blank)
      axis.ticks   = element_line(color = "gray90") ## show axis ticks (theme_minimal hides them by default)
    ); plot_overlap_carib
  
  ## Plot together
  p1 <- plot_overlap_range  + theme(legend.position = "none")
  p2 <- plot_overlap_carib  # keep legend here
  comb_overlaps_plots <- cowplot::plot_grid(p1, p2, ncol = 2, labels = c("", "")); comb_overlaps_plots
  
  ## -----------------------------------------------------------------------------
  ## Grid cell count histogram and barplot
  
  ## ---- HMS-style histogram of anomalies (panel D) ------------------------------
  ## anom_masked: SpatRaster/Raster* of anomalies already masked to species range
  ## species_name: string for title
  ## exp_name: label for x-axis (e.g., "SST st anom" or "o200m st anom")
  
  anom_histogram <- function(anom_masked, species_name, exp_name = "", domain = "") {
    
    ## Pull raster values
    vals <- if (inherits(anom_masked, "SpatRaster")) {
      terra::values(anom_masked, mat = FALSE)
    } else {
      raster::values(anom_masked)
    }
    vals <- vals[is.finite(vals)]
    
    ## Set breaks at 0.25 increments spanning the masked data and use same colors as HMS
    breaks <- seq(floor(min(vals)), ceiling(max(vals)), by = 0.25)
    my_colors <- rep("red", length(breaks))
    my_colors[breaks >= -0.5 & breaks <=  0.5] <- "green"
    my_colors[breaks <  -0.5 & breaks >= -1.5] <- "yellow"
    my_colors[breaks >   0.5 & breaks <=  1.5] <- "yellow"
    my_colors[breaks <  -1.5 & breaks >= -2.0] <- "orange"
    my_colors[breaks >   1.5 & breaks <=  2.0] <- "orange"
    
    ## Make the histogram (freq=FALSE so y-axis is proportion/density)
    h <- hist(vals,
              breaks = breaks,
              freq   = FALSE,
              col    = my_colors,
              xlab   = paste(exp_name, "anomalies"),
              ylab   = "Percent",
              main   = paste0(exp_name, " | ", species_name, " | ", domain))
  }

  par(mfrow=c(1,2))
  anom_histogram(anom_masked, species_name, exp_name = "0200", domain = "Entire range")  
  anom_histogram(anom_masked_carib, species_name, exp_name = "0200", domain = "Caribbean")  
  par(mfrow=c(1,1))
    
  # install.packages("cowplot") # if needed
  library(ggplotify)
  library(cowplot)
  
  p1 <- plot_overlap_range  + theme(legend.position = "none")
  p2 <- plot_overlap_carib
  p3 <- as.ggplot(~ anom_histogram(anom_masked,       species_name, exp_name = "O200", domain = "Entire range"))
  p4 <- as.ggplot(~ anom_histogram(anom_masked_carib, species_name, exp_name = "O200", domain = "Caribbean"))
  
  cowplot::plot_grid(
    cowplot::plot_grid(p1, p2, ncol = 2, labels = c("", "")),
    cowplot::plot_grid(p3, p4, ncol = 2, labels = c("", "")),
    ncol = 1,
    rel_heights = c(1, 1)
  )
  
  ## -----------------------------------------------------------------------------
  ## Make histogram of anomalies 
  
  anom_histogram_gg <- function(r, species_name, exp_name, domain) {
    vals <- terra::values(r)
    vals <- vals[is.finite(vals)]
    df <- data.frame(anom = vals)
    bin_num <- 50
    
    ggplot(df, aes(x = anom, fill = ..x..)) +
      geom_histogram(bins = bin_num, color = "gray25") +
      scale_fill_gradientn(colors = turbo(bin_num), name = "Standardized \nanomaly") +
      scale_x_continuous(expand = c(0,0)) +   # no padding on x-axis
      scale_y_continuous(expand = c(0,0)) +   # optional: also removes gap at y=0
      labs(
        title = paste0(species_name, " | ", exp_name, " | ", domain),
        x = "",
        y = "Count"
      ) +
      geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
      theme_minimal(base_size = 12) + 
      theme(
        axis.text   = element_text(color = "black"),
        axis.title  = element_text(color = "black"),
        axis.ticks  = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),   # black x-axis
        axis.line.y = element_line(color = "black"),   # add y-axis for balance
        legend.position = c(0.95, 0.95),
        legend.justification = c("right","top")
      )
  }
  
  ## Make plots
  p_hist_range <- anom_histogram_gg(anom_masked, species_name, "O200", "Entire range"); p_hist_range
  p_hist_carib <- anom_histogram_gg(anom_masked_carib, species_name, "O200", "Caribbean"); p_hist_carib
  
  ## -----------------------------------------------------------------------------
  ## Make summary of anomalies 
  
  anom_summary_bars <- function(r, species_name, exp_name, domain) {
    # 1. get raster values
    vals <- terra::values(r)
    vals <- vals[is.finite(vals)]
    df <- data.frame(anom = vals)
    
    # 2. cut into categories
    df$cat <- cut(
      df$anom,
      breaks = c(-Inf, -1.5, -0.5, 0.5, 1.5, Inf),
      labels = c("< -1.5", "-1.5 to -0.5", "-0.5 to +0.5", "0.5 to 1.5", "> 1.5"),
      right = TRUE
    )
    
    # 3. assign representative values for coloring (bin midpoints)
    cat_vals <- c(-2, -1, 0, 1, 2)  
    cols <- setNames(
      turbo(5)[as.integer(scales::rescale(cat_vals, to = c(1,5)))],
      levels(df$cat)
    )
    
    # 4. summarize counts and percents
    df_sum <- as.data.frame(table(df$cat))
    colnames(df_sum) <- c("Category", "Count")
    df_sum$Category <- factor(df_sum$Category, levels = levels(df$cat))
    df_sum$Percent <- 100 * df_sum$Count / sum(df_sum$Count)
    
    # 5. plot with labels
    ggplot(df_sum, aes(x = Category, y = Count, fill = Category)) +
      geom_col(color = "gray25") +
      geom_text(aes(label = paste0(round(Percent,1), "%")),
                vjust = -0.5, color = "black", size = 3.5) +
      scale_fill_manual(values = cols, guide = "none") +
      scale_y_continuous(
        expand = expansion(mult = c(0, 0.05))  # 0 bottom, 5% top padding
      ) +
      coord_cartesian(clip = "off") +          # don't clip the top labels
      labs(title = paste0(species_name, " | ", exp_name, " | ", domain),
           x = "Anomaly category", y = "Count") +
      theme_minimal(base_size = 12) +
      theme(axis.text=element_text(color="black"),
            axis.title=element_text(color="black"),
            axis.ticks=element_line(color="black"),
            axis.line.x=element_line(color="black"),
            axis.line.y=element_line(color="black"),
            panel.grid.major.x = element_blank(),
            plot.margin = margin(5, 10, 5, 10))
  }
  
  p_summary_range <- anom_summary_bars(anom_masked, species_name, "O200", "Entire range"); p_summary_range
  p_summary_carib <- anom_summary_bars(anom_masked_carib, species_name, "O200", "Caribbean"); p_summary_carib
  
  
  
  
  
  
  
  
  
    
    
    # summarized barplot (L/M/H/V) + mean score text
    
    # 5) L/M/H/V counts & percents (same logic as HMS)
    L <- sum(h$counts[h$breaks >= -0.5 & h$breaks <=  0.5], na.rm = TRUE)
    1M <- sum(h$counts[(h$breaks < -0.5 & h$breaks >= -1.5) |
                        (h$breaks >  0.5 & h$breaks <=  1.5)], na.rm = TRUE)
    H <- sum(h$counts[(h$breaks < -1.5 & h$breaks >= -2.0) |
                        (h$breaks >  1.5 & h$breaks <=  2.0)], na.rm = TRUE)
    V <- sum(h$counts[h$breaks < -2.0 | h$breaks > 2.0], na.rm = TRUE)
    
    tot <- sum(h$counts)
    Lp <- L / tot; Mp <- M / tot; Hp <- H / tot; Vp <- V / tot
    
    # 6) weighted mean exposure score (1,2,3,4 for L,M,H,V)
    exp_fact_mean <- ((L*1) + (M*2) + (H*3) + (V*4)) / (L + M + H + V)
    
    invisible(list(
      hist_object   = h,
      breaks        = breaks,
      Lp = Lp, Mp = Mp, Hp = Hp, Vp = Vp,
      exp_fact_mean = exp_fact_mean
    ))
  }
  
  barplot(height = c(res_D$Lp, res_D$Mp, res_D$Hp, res_D$Vp),
          names.arg = c("L","M","H","V"),
          col = c("green","yellow","orange","red"),
          ylab = "Percent", ylim = c(0,1))
  abline(h = 0)
  text(x = 1, y = 0.85, labels = round(res_D$exp_fact_mean, 1))
  
