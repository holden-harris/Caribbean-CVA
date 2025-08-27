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
library(raster)    
library(terra)
library(sp)
library(sf)
library(viridisLite) 
library(ggplot2)
library(rnaturalearth)
library(scales)
library(ggplotify) 
library(patchwork) 

##------------------------------------------------------------------------------
## Set directory paths
exp_dir     <- "./data/cmip6/"
spp_dir     <- "./data/species-distribution-shapefiles/"
out_dir     <- "./outputs/"

## -----------------------------------------------------------------------------
## Set geographic extents (bounding box)

## Caribbean Sea
xlim_carib <- c(-92, -57)
ylim_carib <- c(  6,  28)
carib_ext  <- c(xlim_carib, ylim_carib)

## U.S. Caribbean Extent
xlim_uscar <- c(-69, -63.0)
ylim_uscar <- c(16, 20.0)
uscar_ext  <- c(xlim_uscar, ylim_uscar)

## W. Atlantic Ocean
xlim_nwa <- c(-99, -40) # W Atlantic: 99°W to 40°W
ylim_nwa <- c(-5, 72)

## -----------------------------------------------------------------------------
## Utility functions 

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
## 
## Plotting Functions

## -----------------------------------------------------------------------------
## Make maps for species distribution 
plot_distribution <- function (sp, xlim, ylim, title = NULL, domain = NULL) {
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") ## Pullworldmap
  p_distribution <- ggplot() +  
    geom_sf(        ## draw land polygons for geographic context
      data = world,
      fill = "gray80",  ## fill color for land
      color = "gray30", ## separates countries
      linewidth = 0.5
    ) +
    geom_sf(        ## draw the species distribution polygon
      data = sp,
      fill = "black", color = "black",  ## solid black fill  
      alpha = 0.5                       ## but semi-transparent so bathy/coastlines show through
    ) +
    coord_sf(       ## set the map window (crop) using Caribbean bbox
      xlim = xlim, ylim = ylim, expand = FALSE  ## no padding around the bbox
    ) +
    labs(           ## titles and axis labels
      title = paste(title),
    ) +
    theme_minimal() + ## clean base theme
    theme(                               
      axis.text  = element_text(color = "black"),   ## axis tick labels black
      axis.title = element_text(color = "black", size = 16),   ## axis title black (if not blank)
      axis.ticks   = element_line(color = "gray90") ## show axis ticks (theme_minimal hides them by default)
    )
  p_distribution
}

## -----------------------------------------------------------------------------
## Plot anomalies (base R)

plot_anomalies <- function (anom, extent = "", carib_box = 'y', mar = c(2,2,2,2)){
  op <- par(no.readonly = TRUE); on.exit(par(op))
  par(mar = mar)
  plot(anom,
       main = paste0 (exp_name, " anomalies: ", extent),
       xlab = "", ylab = "",
       col = turbo(100))
  maps::map("world", add = TRUE, col = "grey20", lwd = 0.6)
  if (carib_box == 'y'){
    rect(xlim_carib[1], ylim_carib[1], xlim_carib[2], ylim_carib[2],
         border = "purple", lwd = 2)
  }
}


## -----------------------------------------------------------------------------
## Plot anomalies with GG-Plots
## Makes compiling with distribution maps easioer

plot_anomalies_gg <- function(anom, extent = "", carib_box = 'n', uscar_box = 'n') {
  df <- as.data.frame(anom, xy = TRUE, na.rm = TRUE)
  names(df)[3] <- "anom"
  
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  xlim_in <- range(df$x); ylim_in <- range(df$y)
  
  lab_lon <- function(x) sprintf("%g°%s", abs(x), ifelse(x < 0, "W", "E"))
  lab_lat <- function(y) sprintf("%g°%s", abs(y), ifelse(y < 0, "S", "N"))
  
  p <- ggplot() +
    geom_raster(data = df, aes(x = x, y = y, fill = anom)) +
    geom_sf(data = world, fill = "gray80", color = "gray30", linewidth = 0.5) +
    scale_fill_gradientn(colors = turbo(100), 
                         name =  paste0(exp_name)) +
    coord_sf(xlim = xlim_in, ylim = ylim_in, expand = FALSE, default_crs = sf::st_crs(4326), clip = "on") +
    scale_x_continuous(labels = lab_lon, guide = guide_axis(check.overlap = TRUE)) +
    scale_y_continuous(labels = lab_lat, guide = guide_axis(check.overlap = TRUE)) +
    labs(title = paste0("Anomalies: ", extent), x = "", y = "") +
    theme_minimal() +
    theme(axis.text = element_text(color = "black"),
          axis.title = element_text(color = "black"),
          axis.ticks = element_line(color = "gray90"))
  
  ## Option: Add rectangle for Caribbean Sea
  if (carib_box == 'y') {
    p <- p + annotate("rect",
                      xmin = xlim_carib[1], xmax = xlim_carib[2],
                      ymin = ylim_carib[1], ymax = ylim_carib[2],
                      fill = NA, color = "purple", linewidth = 1
    )
  }
  ## Option: Add rectangle for U.S. Caribbean
  if (uscar_box == 'y') {
    p <- p + annotate("rect",
                      xmin = xlim_uscar[1], xmax = xlim_uscar[2],
                      ymin = ylim_uscar[1], ymax = ylim_uscar[2],
                      fill = NA, color = "magenta", linewidth = 1
    )
  }
  p
}

  ## -----------------------------------------------------------------------------
  ## Overlap plot 1
  ## Map overlap of spatial range and anomalies

  ## helper: keep only interior breaks, excluding the min/max limits
  .clip_interior <- function(b, lims, eps = 1e-8) {
    b[b > (min(lims) + eps) & b < (max(lims) - eps)]
  }
  
  overlap_range <- function(df, species_name, exp_name, domain, anom,
                            xlim, ylim, main_title = "all", 
                            legend_title = NULL, world = NULL, 
                            world_scale = "medium", base_size = 12, 
                            fill_limits = NULL, tick_n = 4,
                            minor_by = 1) {
    
    if (is.null(world)) {
      world <- rnaturalearth::ne_countries(scale = world_scale, returnclass = "sf")
    }
    if (is.null(legend_title)) legend_title <- paste0("Standardized \nanomaly\n", exp_name)
    if (main_title == "all"){
      main_title <- paste0("Spe: ", species_name, "\nExp: ", exp_name,
                           "\n\n", domain, ", ", nrow(df), " cells")
    } else if (main_title == "domain"){
      main_title <- paste0(domain, ", ", nrow(df), " cells")
    }
    
    # major ticks (rounded) then drop endpoints
    bx_raw <- scales::breaks_pretty(n = tick_n)(range(xlim))
    by_raw <- scales::breaks_pretty(n = tick_n)(range(ylim))
    bx <- .clip_interior(unique(round(bx_raw)), xlim)
    by <- .clip_interior(unique(round(by_raw)), ylim)
    
    # minor graticule
    lon_minor <- seq(floor(min(xlim)), ceiling(max(xlim)), by = minor_by)
    lat_minor <- seq(floor(min(ylim)), ceiling(max(ylim)), by = minor_by)
    
    lab_lon <- function(x) sprintf("%g°%s", abs(x), ifelse(x < 0, "W", "E"))
    lab_lat <- function(y) sprintf("%g°%s", abs(y), ifelse(y < 0, "S", "N"))
    
    grat_minor <- try(sf::st_graticule(lon = lon_minor, lat = lat_minor, crs = 4326), silent = TRUE)
    grat_major <- try(sf::st_graticule(lon = bx,        lat = by,        crs = 4326), silent = TRUE)
    
    ggplot() +
      geom_tile(data = df, aes(x = x, y = y, fill = anomaly)) +
      { if (!inherits(grat_minor, "try-error")) 
        geom_sf(data = grat_minor, color = "gray90", linewidth = 0.05, alpha = 0.75, inherit.aes = FALSE) } +
      { if (!inherits(grat_major, "try-error")) 
        geom_sf(data = grat_major, color = "gray90", linewidth = 0.05, alpha = 0.75, inherit.aes = FALSE) } +
      geom_sf(data = world, fill = "gray80", color = "gray30", linewidth = 0.5) +
      scale_fill_gradientn(colors = turbo(100), name = legend_title, limits = fill_limits) +
      coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, default_crs = sf::st_crs(4326)) +
      scale_x_continuous(breaks = bx, labels = lab_lon, guide = guide_axis(check.overlap = TRUE)) +
      scale_y_continuous(breaks = by, labels = lab_lat, guide = guide_axis(check.overlap = TRUE)) +
      labs(title = main_title, x = "", y = "") +
      theme_minimal(base_size = base_size) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text  = element_text(color = "black"),
            axis.title = element_text(color = "black"),
            axis.ticks = element_line(color = "gray90"))
  }

  ## -----------------------------------------------------------------------------
  ## Overlap plot 2
  ## Make histogram of anomalies 
  
  anom_histogram_gg <- function(r, fill_limits = NULL) {
    vals <- terra::values(r)
    vals <- vals[is.finite(vals)]
    df <- data.frame(anom = vals)
    bin_num <- 50
    
    ggplot(df, aes(x = anom, fill = ..x..)) +
      geom_histogram(bins = bin_num, color = "gray25") +
      scale_fill_gradientn(
        colors = turbo(bin_num), 
        name   = "Standardized \nanomaly",
        limits = fill_limits           # keep legend consistent
      ) +
      scale_x_continuous(
        limits = fill_limits,          # <-- unify x-axis across plots
        expand = c(0,0)
      ) +
      scale_y_continuous(expand = c(0,0)) +
      labs(x = "Standardized anomaly", y = "Count") +
      geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
      theme_minimal(base_size = 12) + 
      theme(
        axis.text   = element_text(color = "black"),
        axis.title  = element_text(color = "black"),
        axis.ticks  = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        legend.position       = c(0.95, 0.95),
        legend.justification  = c("right","top")
      )
  }
  
  ## -----------------------------------------------------------------------------
  ## Overlap plot 3
  ## Summarize overlapped anomalies 
  
  anom_summary_bars <- function(r, species_name, exp_name, domain) {
    ## First get raster values
    vals <- terra::values(r)
    vals <- vals[is.finite(vals)]
    df <- data.frame(anom = vals)
    
    ## Cut into categories
    df$cat <- cut(
      df$anom,
      breaks = c(-Inf, -1.5, -0.5, 0.5, 1.5, Inf),
      labels = c("< -1.5", "-1.5 to -0.5", "-0.5 to +0.5", "0.5 to 1.5", "> 1.5"),
      right = TRUE
    )
    
    ## Assign representative values for coloring (bin midpoints)
    cat_vals <- c(-2, -1, 0, 1, 2)  
    cols <- setNames(
      turbo(5)[as.integer(scales::rescale(cat_vals, to = c(1,5)))],
      levels(df$cat)
    )
    
    ## Summarize counts and percents
    df_sum <- as.data.frame(table(df$cat))
    colnames(df_sum) <- c("Category", "Count")
    df_sum$Category <- factor(df_sum$Category, levels = levels(df$cat))
    df_sum$Percent <- 100 * df_sum$Count / sum(df_sum$Count)
    
    ## Make plot with labels
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
      theme(axis.text = element_text(color="black"),
            axis.text.x = element_text(color = "black", angle = 25, vjust = 0.5, hjust = 0.5),  # rotate labels
            axis.title=element_text(color="black"),
            axis.ticks=element_line(color="black"),
            axis.line.x=element_line(color="black"),
            axis.line.y=element_line(color="black"),
            panel.grid.major.x = element_blank(),
            plot.margin = margin(5, 10, 5, 10))
  }
  
  

## -----------------------------------------------------------------------------
## Read in shape files
shp_files <- list.files(spp_dir, pattern = "\\.shp$", full.names = TRUE, recursive = FALSE)
if (length(shp_files) == 0L) stop("No .shp files found in: ", spp_dir)


## Run loop
#for (i in 1:length(shp_files)){
# for (i in 1:1){ ## Use for testing
i = 2

## Read in species distribution ------------------------------------------------
sp_file <- shp_files[i]
species_name <- name_from_shp(sp_file)
print(paste("Processsing", species_name))
sp  <- sf::st_read(sp_file, quiet = TRUE) |> st_make_valid()

## Conditionally clip lims to the NWA box if it extends beyond it
lims    <- bbox_with_pad(sp, pad = 0.05)   ## Crop exposure anomolies to species range
needs_clip_x <- lims$xlim[1] < xlim_nwa[1] || lims$xlim[2] > xlim_nwa[2]
needs_clip_y <- lims$ylim[1] < ylim_nwa[1] || lims$ylim[2] > ylim_nwa[2]
if (needs_clip_x || needs_clip_y) {
  lims$xlim <- c(max(lims$xlim[1], xlim_nwa[1]),
                 min(lims$xlim[2], xlim_nwa[2]))
  lims$ylim <- c(max(lims$ylim[1], ylim_nwa[1]),
                 min(lims$ylim[2], ylim_nwa[2]))
}



## -----------------------------------------------------------------------------
## Read in exposure anomoly map

exp_nc_file_name <- "o200_1985-2014_2020-2049.nc"
exp_name <- sub("_.*", "", exp_nc_file_name)

f <- file.path(exp_dir, exp_nc_file_name) ## set path
anom <- rast(f, sub = "anomaly") ## one layer
anom[anom > 1e19] <- NA ## Fix fill values
anom  <- rotate(anom) ## NetCDF longitudes are 0–360 so need to rotate to match our extent, which is -180-180

## Crop to ranges: Entire W. Atlantic range, Caribbean Sea, and U.S. Caribbean
anom_range <- crop(anom, extent(c(lims$xlim, lims$ylim)))
anom_carib <- crop(anom, ext(xlim_carib, ylim_carib)) ## Crop to Caribbean bb
anom_uscar <- crop(anom, ext(xlim_uscar, ylim_uscar)) ## Crop to U.S. Carib


  

## -----------------------------------------------------------------------------
## Plot distribution and exposure factors


# Make the six plots
p1 <- plot_distribution(sp, lims$xlim, lims$ylim, title = species_name) +
  theme(plot.title = element_text(size = 16))
p2 <- plot_distribution(sp, xlim_carib, ylim_carib)
p3 <- plot_anomalies_gg(anom,       extent = "Global",      carib_box = 'y', uscar_box = 'y')
p4 <- plot_anomalies_gg(anom_range, extent = "W. Atlantic", carib_box = 'y', uscar_box = 'y')
p5 <- plot_anomalies_gg(anom_carib, extent = "Caribbean Sea", uscar_box = 'y')
p6 <- plot_anomalies_gg(anom_uscar, extent = "U.S. Caribbean")

## Arrange into 2 columns × 3 rows
six_panel <- (
  (p1 | p2) / (p3 | p4) / (p5 | p6)
  ) + 
  plot_layout(heights = c(1.5, 1, 1)) + ## First row 1.5× larger
  plot_annotation(tag_levels = 'A') 
six_panel

  ## -----------------------------------------------------------------------------
  ## Map overlap

  ## Vectorize species
  sp      <- st_transform(sp, crs(anom_range))  ## match raster CRS
  spv     <- vect(sp)

  ## Make and apply cover mask
  ## Rasterize polygon to the anomaly grid (include partial cells)
  ## Note that cover=TRUE returns fraction of each cell covered by the polygon
  
  ## Mask cover species range (or W. Atlantic)
  mask_cover  <- rasterize(spv, anom_range, field = 1, background = NA, cover = TRUE)
  anom_masked <- mask(anom_range, mask_cover)
  
  ## Mask cover Caribbean
  anom_masked_carib <- crop(anom_masked, carib_ext)
  anom_masked_uscar <- crop(anom_masked, uscar_ext)
  
  ## Make dataframe for masked values
  ## anom_masked is a SpatRaster/Raster* (masked by species footprint). Convert it to df for ggplot. 
  df_mask       <- as.data.frame(anom_masked, xy = TRUE, na.rm = TRUE)
  df_mask_carib <- as.data.frame(anom_masked_carib, xy = TRUE, na.rm = TRUE)
  df_mask_uscar <- as.data.frame(anom_masked_uscar, xy = TRUE, na.rm = TRUE)

  ## Compute global min/max
  min_anom_val  <- min(df_mask$anomaly, na.rm = TRUE)
  max_anom_val  <- max(df_mask$anomaly, na.rm = TRUE)
  fill_limits   <- c(min_anom_val, max_anom_val)
  
  ## Make plots ----------------------------------------------------------------
  ## Make overlap plots
  p_overlap_range <- overlap_range(df_mask, species_name, exp_name,
    domain       = "W. Atlantic range",
    xlim         = lims$xlim,
    ylim         = lims$ylim,
    main_title   = "domain", 
  ); p_overlap_range
  
  p_overlap_carib <- overlap_range(df_mask_carib, species_name, exp_name,
    domain       = "Caribbean Sea",
    xlim         = xlim_carib,
    ylim         = ylim_carib,
    main_title   = "domain", 
    fill_limits  = fill_limits
  ); p_overlap_carib
  
  p_overlap_uscar <- overlap_range(df_mask_uscar, species_name, exp_name,
    domain       = "U.S. Caribbean",
    xlim         = xlim_uscar,
    ylim         = ylim_uscar,
 #   main_title   = "domain", 
    fill_limits  = fill_limits
  ); p_overlap_uscar
  
  ## Make histogram plots
  p_hist_range <- anom_histogram_gg(anom_masked, fill_limits); p_hist_range
  p_hist_carib <- anom_histogram_gg(anom_masked_carib, fill_limits); p_hist_carib
  p_hist_uscar <- anom_histogram_gg(anom_masked_uscar, fill_limits = fill_limits); p_hist_uscar
  
  ## Make summary plots
  p_summary_range <- anom_summary_bars(anom_masked, species_name, exp_name, "Entire range"); p_summary_range
  p_summary_carib <- anom_summary_bars(anom_masked_carib, species_name, exp_name, "Caribbean"); p_summary_carib
  p_summary_uscar <- anom_summary_bars(anom_masked_uscar, species_name, exp_name, "U.S. Caribbean"); p_summary_uscar
  
  ## ---------------------------------------------------------------------------
  ## Combine overlap plots
  ## Plot overlap map, histogram, and summary bar chart together
  
  ## Keep one legend (from the Caribbean map, say)
  leg <- cowplot::get_legend(p_overlap_carib + theme(legend.position = "right"))
  
  ## Remove legends from the others
  a <- p_overlap_range  + theme(legend.position = "none")                 
  b <- p_overlap_carib  + theme(legend.position = "none")
  c <- p_overlap_uscar  + theme(legend.position = "none")
    
  d <- p_hist_range     + theme(legend.position = "none", plot.title = element_blank())             
  e <- p_hist_carib     + theme(legend.position = "none", plot.title = element_blank()) 
  f <- p_hist_uscar     + theme(legend.position = "none", plot.title = element_blank()) 

  g <- p_summary_range  + theme(plot.title = element_blank())             
  h <- p_summary_carib  + theme(plot.title = element_blank())           
  i <- p_summary_uscar  + theme(plot.title = element_blank())           
  
  ## Rows (Caribbean left, Range right), with labels
  row1 <- cowplot::plot_grid(a, b, c, ncol = 3, labels = c("A","B","C"), 
                             label_size = 12, label_fontface = "bold")
  row2 <- cowplot::plot_grid(d, e, f, ncol = 3, labels = c("D","E","F"),
                             label_size = 12, label_fontface = "bold")
  row3 <- cowplot::plot_grid(g, h, i, ncol = 3, labels = c("G","H","I"),
                             label_size = 12, label_fontface = "bold")
  
  main <- cowplot::plot_grid(row1, row2, row3, ncol = 1,
                             rel_heights = c(1, 1, 1), align = "hv")
  final_9panel <- cowplot::plot_grid(main, leg, ncol = 2, rel_widths = c(1, 0.08))
  
  
  ## Add the legend on the right
  final_9panel <- cowplot::plot_grid(main, leg, ncol = 2, rel_widths = c(1, 0.12))
  final_9panel
  
  ## Save plot
  out_name <- paste0(out_dir, species_name, "_", exp_name, "_9panel.png")
  ggsave(file.path(out_name), final_9panel, 
         width = 11, height = 11, dpi = 500, bg = "white")
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
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
  
