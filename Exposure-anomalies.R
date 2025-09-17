##------------------------------------------------------------------------------
## Org:     Harris Analytics and Research LLC | Isla Mar
## Project: Caribbean CVA
## Contact: Holden Earl Harris | holden.earl.harris@gmail.com
## Code:    Exposure analyses
##------------------------------------------------------------------------------

################################################################################
##
## Setup

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
out_dir     <- "./outputs/exposure-overlap-12panel/"

## Species shape files
shp_files <- list.files(spp_dir, pattern = "\\.shp$", full.names = TRUE, recursive = FALSE)
if (length(shp_files) == 0L) stop("No .shp files found in: ", spp_dir)
print(shp_files)

## CMIP exposure files
nc_files  <- list.files(exp_dir, pattern = "\\.nc$",  full.names = TRUE, recursive = FALSE)
if (length(nc_files) == 0L) stop("No .nc files found in: ", exp_dir)
print(nc_files)

## -----------------------------------------------------------------------------
## Set geographic extents (bounding boxes)

## Caribbean Sea
xlim_carib <- c(-92, -57)
ylim_carib <- c(  6,  28)
carib_ext  <- c(xlim_carib, ylim_carib)

## U.S. Caribbean Extent
xlim_uscar <- c(-69, -63.0)
ylim_uscar <- c(16, 20.0)
uscar_ext  <- c(xlim_uscar, ylim_uscar)

## W. Atlantic Ocean
xlim_nwa <- c(-99, -40) ## W Atlantic: 99°W to 40°W
ylim_nwa <- c(-5, 72)


################################################################################
##
## Functions

## -----------------------------------------------------------------------------
## Utility functions 

## Extract species name from shp filename
name_from_shp <- function(path) { # Base path -> "Nice Name"
  base <- tools::file_path_sans_ext(basename(path))  # insert space between lower->Upper (CamelCase), e.g., "AtlanticHerring" -> "Atlantic Herring"
  base <- gsub("(?<=[a-z])(?=[A-Z])", " ", base, perl = TRUE)  # replace _, -, . with spaces
  base <- gsub("[_\\.\\-]+", " ", base)  # squeeze multiple spaces, trim
  base <- gsub("\\s+", " ", trimws(base))   # Title Case (keeps Genus species looking right)
  base
}

## Extract exposure factor (exp_name) from nc filename, e.g., "o200_1985-2014_2020-2049.nc" -> "o200"
exp_name_from_nc <- function(path) {
  sub("_.*$", "", basename(path))
}

## Function to get bounding box from spatial distribution
bbox_with_pad <- function(sf_obj, pad = 0.05){
  bb <- st_bbox(sf_obj)
  dx <- as.numeric(bb["xmax"] - bb["xmin"])
  dy <- as.numeric(bb["ymax"] - bb["ymin"])
  list(xlim = c(bb["xmin"] - pad*dx, bb["xmax"] + pad*dx),
       ylim = c(bb["ymin"] - pad*dy, bb["ymax"] + pad*dy))
}

## Open PDF graphics maker ---------------------------------------------------
safe_open_pdf <- function(path, width = 11, height = 11, onefile = TRUE) {
  has_cairo <- tryCatch(isTRUE(base::capabilities("cairo")), error = function(e) FALSE)
  if (has_cairo) {
    ok <- tryCatch({
      grDevices::cairo_pdf(path, width = width, height = height, onefile = onefile)
      TRUE
    }, error = function(e) FALSE)
    if (ok) return(invisible(NULL))
  }
  grDevices::pdf(path, width = width, height = height, onefile = onefile, useDingbats = FALSE)
}

## Make clean species name (e.g., Atlantic Herring to Atlantic-Herring)  
species_name_clean <- function(x) gsub(" ", "-", x)

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

plot_anomalies <- function (anom, exp_name, extent = "", carib_box = 'y', mar = c(2,2,2,2)){
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
## Makes compiling with distribution maps easier

plot_anomalies_gg <- function(anom, exp_name, extent = "", carib_box = 'n', uscar_box = 'n') {
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
                         name =  paste0("Stnd.\nanom\n", exp_name)) +
    coord_sf(xlim = xlim_in, ylim = ylim_in, expand = FALSE, default_crs = sf::st_crs(4326), clip = "on") +
    scale_x_continuous(labels = lab_lon, guide = guide_axis(check.overlap = TRUE)) +
    scale_y_continuous(labels = lab_lat, guide = guide_axis(check.overlap = TRUE)) +
    labs(title = paste0(extent, ": ", exp_name), x = "", y = "") +
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

  ## ---------------------------------------------------------------------------
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
    if (is.null(legend_title)) legend_title <- paste0("Standardized\nanomaly\n", exp_name)
    if (main_title == "all"){
      main_title <- paste0("Spe: ", species_name, "\nExp: ", exp_name,
                           "\n\n", domain, ", ", nrow(df), " cells")
    } else if (main_title == "domain"){
      main_title <- paste0(domain, ", ", nrow(df), " cells")
    }
    
    ## Major ticks (rounded) then drop endpoints
    bx_raw <- scales::breaks_pretty(n = tick_n)(range(xlim))
    by_raw <- scales::breaks_pretty(n = tick_n)(range(ylim))
    bx <- .clip_interior(unique(round(bx_raw)), xlim)
    by <- .clip_interior(unique(round(by_raw)), ylim)
    
    ## Minor graticule
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
  ## HMS-style LMHV histogram (base graphics), wrapped for cowplot/patchwork use
  ## Returns both a plot expression and the LMHV summary you can reuse.
  
  ## --- Tweakable defaults for row 4–5 base plots ---
  .lmhv_mar <- c(5.2, 4.8, 2.2, 1.4) + 0.1  # bottom, left, top, right (lines)
  .lmhv_mgp <- c(2.4, 0.7, 0)               # axis title, tick labels, tick line
  .lmhv_axis_cex <- 0.8
  .lmhv_lab_cex  <- 0.9
  .lmhv_main_cex <- 0.9
  
  lmhv_histogram_base <- function(anom_masked, species_name, exp_name = "", domain = "",
                                  mar = .lmhv_mar, mgp = .lmhv_mgp,
                                  cex_axis = .lmhv_axis_cex, cex_lab = .lmhv_lab_cex,
                                  cex_main = .lmhv_main_cex,
                                  do_plot = TRUE) {   # new flag
    # values
    vals <- if (inherits(anom_masked, "SpatRaster")) terra::values(anom_masked, mat = FALSE)
    else raster::values(anom_masked)
    vals <- vals[is.finite(vals)]
    if (!length(vals)) stop("No finite values in masked anomalies for: ", domain)
    
    # bins & colors
    breaks <- seq(floor(min(vals)), ceiling(max(vals)), by = 0.25)
    cols   <- rep("red", length(breaks))
    cols[breaks >= -0.5 & breaks <=  0.5] <- "green"
    cols[breaks <  -0.5 & breaks >= -1.5] <- "yellow"
    cols[breaks >   0.5 & breaks <=  1.5] <- "yellow"
    cols[breaks <  -1.5 & breaks >= -2.0] <- "orange"
    cols[breaks >   1.5 & breaks <=  2.0] <- "orange"
    
    h <- NULL
    if (do_plot) {
      op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
      par(mar = mar, mgp = mgp, tcl = -0.25, las = 1, xaxs = "i", yaxs = "i",
          xpd = NA, cex.axis = cex_axis, cex.lab = cex_lab, cex.main = cex_main)
      
      h <- hist(vals, breaks = breaks, freq = FALSE, col = cols,
                xlab = "Standardized anomaly",
                ylab = "Percent",
                main = "")
    } else {
      # still compute histogram counts invisibly
      h <- hist(vals, breaks = breaks, plot = FALSE)
    }
    
    mids <- h$mids; cnts <- h$counts
    L <- sum(cnts[mids >= -0.5 & mids <=  0.5], na.rm = TRUE)
    M <- sum(cnts[(mids < -0.5 & mids >= -1.5) | (mids > 0.5 & mids <= 1.5)], na.rm = TRUE)
    H <- sum(cnts[(mids < -1.5 & mids >= -2.0) | (mids > 1.5 & mids <= 2.0)], na.rm = TRUE)
    V <- sum(cnts[mids < -2.0 | mids > 2.0], na.rm = TRUE)
    tot <- L + M + H + V
    
    out <- list(hist_object = h, breaks = breaks,
                Lp = if (tot > 0) L/tot else 0,
                Mp = if (tot > 0) M/tot else 0,
                Hp = if (tot > 0) H/tot else 0,
                Vp = if (tot > 0) V/tot else 0,
                exp_mean = if (tot > 0) ((L*1) + (M*2) + (H*3) + (V*4))/tot else NA_real_)
    class(out) <- c("lmhv_hist_summary","list")
    out
  }
  
  
  lmhv_barplot_base <- function(lmhv_summary,
                                mar = .lmhv_mar, mgp = .lmhv_mgp,
                                cex_axis = .lmhv_axis_cex, cex_lab = .lmhv_lab_cex,
                                cex_main = .lmhv_main_cex,
                                do_plot = TRUE) {
    stopifnot(inherits(lmhv_summary, "lmhv_hist_summary"))
    pcts <- c(lmhv_summary$Lp, lmhv_summary$Mp, lmhv_summary$Hp, lmhv_summary$Vp)
    
    if (do_plot) {
      op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
      par(mar = mar, mgp = mgp, tcl = -0.25, las = 1, xaxs = "i", yaxs = "i",
          xpd = NA, cex.axis = cex_axis, cex.lab = cex_lab, cex.main = cex_main)
      
      bp <- barplot(height = pcts,
                    names.arg = c("L","M","H","V"),
                    col = c("green","yellow","orange","red"),
                    ylim = c(0, 1),
                    xlab = "Anomaly category",
                    ylab = "Percent",
                    main = "")
      abline(h = 0)
      if (is.finite(lmhv_summary$exp_mean)) {
        text(x = bp[1], y = 0.92, labels = round(lmhv_summary$exp_mean, 1), xpd = NA)
      }
    }
    
    invisible(lmhv_summary)  # return the summary regardless, don’t auto-plot if do_plot=FALSE
  }
  
  
  ## slightly larger inner margins for the base plots
  .lmhv_mar <- c(7.2, 6.6, 1.2, 2.0) + 0.1   # bottom, left, top, right
  .lmhv_mgp <- c(2.8, 0.7, 0)
  
  ## wrappers that add an outer padding via ggdraw/draw_plot
  lmhv_histogram_as_gg <- function(x, species_name = NULL, exp_name = NULL, domain = NULL,
                                   mar = .lmhv_mar, mgp = .lmhv_mgp,
                                   cex_axis = .lmhv_axis_cex, cex_lab = .lmhv_lab_cex,
                                   cex_main = .lmhv_main_cex) {
    # compute summary if needed
    summary <- if (inherits(x, "lmhv_hist_summary")) {
      x
    } else {
      lmhv_histogram_base(x, species_name, exp_name, domain,
                          mar = mar, mgp = mgp,
                          cex_axis = cex_axis, cex_lab = cex_lab, cex_main = cex_main,
                          do_plot = FALSE)
    }
    
    # build a ggplot-like grob by DRAWING from the summary
    base_gg <- ggplotify::as.ggplot(function() lmhv_hist_draw_base(summary, mar, mgp, cex_axis, cex_lab, cex_main))
    
    # for older cowplot: set clip at ggdraw, not draw_plot
    cowplot::ggdraw(clip = "off") +
      cowplot::draw_plot(base_gg, x = 0.06, y = 0.14, width = 0.90, height = 0.86)
  }
  
  lmhv_barplot_as_gg <- function(x,
                                 mar = .lmhv_mar, mgp = .lmhv_mgp,
                                 cex_axis = .lmhv_axis_cex, cex_lab = .lmhv_lab_cex,
                                 cex_main = .lmhv_main_cex) {
    # Require a summary; if user passed a raster by mistake, fail clearly
    stopifnot("lmhv_hist_summary" %in% class(x))
    summary <- x
    
    base_gg <- ggplotify::as.ggplot(function() lmhv_bar_draw_base(summary, mar, mgp, cex_axis, cex_lab, cex_main))
    cowplot::ggdraw(clip = "off") +
      cowplot::draw_plot(base_gg, x = 0.06, y = 0.14, width = 0.90, height = 0.86)
  }
  
  ## ---------- helpers that actually draw from a summary (no auto-plot side effects)
  lmhv_hist_draw_base <- function(summary,
                                  mar = .lmhv_mar, mgp = .lmhv_mgp,
                                  cex_axis = .lmhv_axis_cex, cex_lab = .lmhv_lab_cex,
                                  cex_main = .lmhv_main_cex) {
    stopifnot(inherits(summary, "lmhv_hist_summary"))
    h <- summary$hist_object
    
    # colors by bin midpoints
    mids <- h$mids
    cols <- rep("red", length(mids))
    cols[mids >= -0.5 & mids <=  0.5] <- "green"
    cols[(mids < -0.5 & mids >= -1.5) | (mids > 0.5 & mids <= 1.5)] <- "yellow"
    cols[(mids < -1.5 & mids >= -2.0) | (mids > 1.5 & mids <= 2.0)] <- "orange"
    
    op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
    par(mar = mar, mgp = mgp, tcl = -0.25, las = 1, xaxs = "i", yaxs = "i",
        xpd = NA, cex.axis = cex_axis, cex.lab = cex_lab, cex.main = cex_main)
    
    plot(h, freq = FALSE, col = cols,
         xlab = "Standardized anomaly",
         ylab = "Percent",
         main = "")
  }
  
  lmhv_bar_draw_base <- function(summary,
                                 mar = .lmhv_mar, mgp = .lmhv_mgp,
                                 cex_axis = .lmhv_axis_cex, cex_lab = .lmhv_lab_cex,
                                 cex_main = .lmhv_main_cex) {
    stopifnot(inherits(summary, "lmhv_hist_summary"))
    pcts <- c(summary$Lp, summary$Mp, summary$Hp, summary$Vp)
    
    op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
    par(mar = mar, mgp = mgp, tcl = -0.25, las = 1, xaxs = "i", yaxs = "i",
        xpd = NA, cex.axis = cex_axis, cex.lab = cex_lab, cex.main = cex_main)
    
    bp <- barplot(height = pcts,
                  names.arg = c("L","M","H","V"),
                  col = c("green","yellow","orange","red"),
                  ylim = c(0, 1),
                  xlab = "Anomaly category",
                  ylab = "Percent",
                  main = "")
    abline(h = 0)
    if (is.finite(summary$exp_mean)) {
      text(x = bp[1], y = 0.92, labels = round(summary$exp_mean, 1), xpd = NA)
    }
  }

################################################################################
##
## Run loops
  
## Clean graphics devices (optional, but helps if devices leaked earlier)
while (!is.null(grDevices::dev.list())) grDevices::dev.off()
  
## -----------------------------------------------------------------------------
## Process species shape files

## -----------------------------------------------------------------------------
## Core function: Process species shape files and make maps

## For testing
#i = 1; sp_file = shp_files[i]; sp  <- sf::st_read(sp_file, quiet = TRUE) |> st_make_valid()

process_species <- function(sp_file) {
    ## Setup species
    sp  <- sf::st_read(sp_file, quiet = TRUE) |> st_make_valid()
    species_name <- name_from_shp(sp_file)
    sname <- species_name_clean(species_name)
    
    species_dir <- paste0(out_dir, sname)
    if (!dir.exists(species_dir)) dir.create(species_dir, recursive = TRUE)
    
    cat("\n------------------------------------------------------------\n")
    cat(sprintf("Processing species: %s\n", species_name))
    
    ## BBox/clip to NWA. This speeds up processing if the dataset is global. 
    lims    <- bbox_with_pad(sp, pad = 0.05)
    needs_clip_x <- lims$xlim[1] < xlim_nwa[1] || lims$xlim[2] > xlim_nwa[2]
    needs_clip_y <- lims$ylim[1] < ylim_nwa[1] || lims$ylim[2] > ylim_nwa[2]
    if (needs_clip_x || needs_clip_y) {
      lims$xlim <- c(max(lims$xlim[1], xlim_nwa[1]), min(lims$xlim[2], xlim_nwa[2]))
      lims$ylim <- c(max(lims$ylim[1], ylim_nwa[1]), min(lims$ylim[2], ylim_nwa[2]))
    }
    
    ## Open two PDF devices ---------------------------------------------------
    ## Set directory paths
    dist_pdf_path <- file.path(species_dir, paste0(sname, "_Distribution-Anomalies.pdf"))
    over_pdf_path <- file.path(species_dir, paste0(sname, "_Exposure-Overlap-12panel.pdf"))
    
    ## Safe-open PDF graphics device
    safe_open_pdf(dist_pdf_path, width = 11, height = 11, onefile = TRUE)
    dev_dist <- grDevices::dev.cur()
    on.exit({
      if (!is.null(grDevices::dev.list()) && dev_dist %in% unlist(grDevices::dev.list())) {
        grDevices::dev.set(dev_dist); grDevices::dev.off()
      }
    }, add = TRUE)
    
    safe_open_pdf(over_pdf_path, width = 13, height = 17, onefile = TRUE)
    dev_over <- grDevices::dev.cur()
    on.exit({
      if (!is.null(grDevices::dev.list()) && dev_over %in% unlist(grDevices::dev.list())) {
        grDevices::dev.set(dev_over); grDevices::dev.off()
      }
    }, add = TRUE)
    
    ## -------------------------------------------------------------------------
    ## Inner loop over exposures 
    for (seq in seq_along(nc_files)) {
    #for (seq in 1:2) { ## Use for testing first two exposure factors
      
      ## Set file path and pull exposure factor name
      nc_path  <- nc_files[seq]
      exp_name <- exp_name_from_nc(nc_path)
      cat(sprintf("  - Exposure (%d/%d): %s\n", seq, length(nc_files), exp_name)) ## Loop printout
      
      ## Read in anomaly map
      ## Assumptions:
      ## - NetCDF contains a layer named "anomaly" (standardized anomalies)
      ## - Longitudes may be 0..360 and require rotate()
      ## - Fill values may be > 1e19 and should be set to NA
      
      anom <- rast(nc_path, sub = "anomaly") ## one layer
      anom[anom > 1e19] <- NA ## Fix fill values
      anom  <- rotate(anom) ## NetCDF longitudes are 0–360 so need to rotate to match our extent, which is -180-180
      
      ## Crop to ranges: Entire W. Atlantic range, Caribbean Sea, and U.S. Caribbean
      anom_range <- crop(anom, extent(c(lims$xlim, lims$ylim)))
      anom_carib <- crop(anom, ext(xlim_carib, ylim_carib)) ## Crop to Caribbean bb
      anom_uscar <- crop(anom, ext(xlim_uscar, ylim_uscar)) ## Crop to U.S. Carib
      
      ## -----------------------------------------------------------------------------
      ## Figure I: Make distribution maps
      
      p1 <- plot_distribution(sp, lims$xlim, lims$ylim, title = species_name) +
        theme(plot.title = element_text(size = 16))
      p2 <- plot_distribution(sp, xlim_carib, ylim_carib)
      p3 <- plot_anomalies_gg(anom,       exp_name, extent = "Global",      carib_box = 'y', uscar_box = 'y')
      p4 <- plot_anomalies_gg(anom_range, exp_name, extent = "W. Atlantic", carib_box = 'y', uscar_box = 'y')
      p5 <- plot_anomalies_gg(anom_carib, exp_name, extent = "Caribbean Sea", uscar_box = 'y')
      p6 <- plot_anomalies_gg(anom_uscar, exp_name, extent = "U.S. Caribbean")
      
      ## Arrange into 2 columns × 3 rows
      six_panel <- (
        (p1 | p2) / (p3 | p4) / (p5 | p6)
      ) + 
        plot_layout(heights = c(1.5, 1, 1)) + ## First row 1.5× larger
        plot_annotation(tag_levels = 'A') 

      ## -----------------------------------------------------------------------------
      ## Figure II: Make masks and map overlaps
      
      ## Vectorize species
      sp      <- st_transform(sp, crs(anom_range))  ## match raster CRS
      spv     <- vect(sp)
      
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
      
      ## Compute global min/max to use for all the plots
      min_anom_val  <- min(df_mask$anomaly, na.rm = TRUE)
      max_anom_val  <- max(df_mask$anomaly, na.rm = TRUE)
      fill_limits   <- c(min_anom_val, max_anom_val)
      
      ## ---------------------------------------------------------------------------
      ## Call plot functions to make plots
      
      ## Make overlap plots
      p_overlap_range <- overlap_range(df_mask, species_name, exp_name,
                                       domain       = "W. Atlantic range",
                                       xlim         = lims$xlim,
                                       ylim         = lims$ylim,
                                       main_title   = "domain", 
      )
      p_overlap_carib <- overlap_range(df_mask_carib, species_name, exp_name,
                                       domain       = "Caribbean Sea",
                                       xlim         = xlim_carib,
                                       ylim         = ylim_carib,
                                       main_title   = "domain", 
                                       fill_limits  = fill_limits
      )
      
      p_overlap_uscar <- overlap_range(df_mask_uscar, species_name, exp_name,
                                       domain       = "U.S. Caribbean",
                                       xlim         = xlim_uscar,
                                       ylim         = ylim_uscar,
                                       fill_limits  = fill_limits
      )
      
      ## Make histogram plots (now not used for 3x4 plots)
      p_hist_range <- anom_histogram_gg(anom_masked, fill_limits)
      p_hist_carib <- anom_histogram_gg(anom_masked_carib, fill_limits)
      p_hist_uscar <- anom_histogram_gg(anom_masked_uscar, fill_limits = fill_limits)
      
      ## Make summary plots
      p_summary_range <- anom_summary_bars(anom_masked, species_name, exp_name, "Entire range")
      p_summary_carib <- anom_summary_bars(anom_masked_carib, species_name, exp_name, "Caribbean")
      p_summary_uscar <- anom_summary_bars(anom_masked_uscar, species_name, exp_name, "U.S. Caribbean")
      
      ## Compute HMS-style LMHV hist summaries + plots for each domain ----
      sum_range <- lmhv_histogram_base(anom_masked,       species_name, exp_name, "Entire range", do_plot = FALSE)
      sum_carib <- lmhv_histogram_base(anom_masked_carib, species_name, exp_name, "Caribbean", do_plot = FALSE)
      sum_uscar <- lmhv_histogram_base(anom_masked_uscar, species_name, exp_name, "U.S. Caribbean", do_plot = FALSE)
      
      ## ---------------------------------------------------------------------------
      ## Combine overlap plots
      ## Plot overlap map, histogram, and summary bar chart together
      
      ## Make panels for 3x4 exposure overlap figure
      a <- p_overlap_range  + theme(legend.position = "none")                 
      b <- p_overlap_carib  + theme(legend.position = "none")
      c <- p_overlap_uscar  + theme(legend.position = "none")
      d <- p_summary_range  + theme(plot.title = element_blank())             
      e <- p_summary_carib  + theme(plot.title = element_blank())           
      f <- p_summary_uscar  + theme(plot.title = element_blank())           
      g <- lmhv_histogram_as_gg(sum_range)
      h <- lmhv_histogram_as_gg(sum_carib)
      i <- lmhv_histogram_as_gg(sum_uscar)
      j <- lmhv_barplot_as_gg(sum_range)
      k <- lmhv_barplot_as_gg(sum_carib)
      l <- lmhv_barplot_as_gg(sum_uscar)
  
      pad <- ggplot2::theme(plot.margin = ggplot2::margin(6, 6, 6, 6))
      a <- a + pad; b <- b + pad; c <- c + pad
      d <- d + pad; e <- e + pad; f <- f + pad
      g <- g + pad; h <- h + pad; i <- i + pad
      j <- j + pad; k <- k + pad; l <- l + pad

      ## Compile plots into rows (A-L) 
      row1 <- cowplot::plot_grid(a, b, c, ncol = 3, labels = c("A","B","C"),
                                 label_size = 12, label_fontface = "bold")
      row2 <- cowplot::plot_grid(d, e, f, ncol = 3, labels = c("D","E","F"),
                                 label_size = 12, label_fontface = "bold")
      row3 <- cowplot::plot_grid(g, h, i, ncol = 3, labels = c("G","H","I"),
                                 label_size = 12, label_fontface = "bold")
      row4 <- cowplot::plot_grid(j, k, l, ncol = 3, labels = c("J","K","L"),
                                 label_size = 12, label_fontface = "bold")
      
      main_12 <- cowplot::plot_grid(row1, row2, row3, row4, ncol = 1, ## Put plots together
                                    rel_heights = c(1.4,0.85,1.2,1.2), align = "hv")
      leg <- cowplot::get_legend(p_overlap_carib + theme(legend.position = "right", legend.justification = "top")) ## Extract legend from one of the plots
      
      final_12panel <- cowplot::ggdraw() + ## Overlay legend on top-right of the main grid
        cowplot::draw_plot(main_12, 0, 0, 1, 1) + ## Fill the canvas
        cowplot::draw_plot(leg, 0.21, 0.73, 0.25, 0.25) ## x,y,width,height in [0,1]
      
      ## Test --> Make single PDF (Don't use in loop)
      grDevices::cairo_pdf(file.path(species_dir, "test_3x4.pdf"), width = 13, height = 17);  print(final_12panel); dev.off()
      
      ## Write one page to each PDF --------------------------------------------
      grDevices::dev.set(dev_over); print(final_12panel)  # page appended to Exposure-Overlap.pdf
      grDevices::dev.set(dev_dist); print(six_panel)     # page appended to Distribution-Anomalies.pdf
      
      ## -----------------------------------------------------------------------
      ## Also save individual PNG plots
      
      ## Distribution and anomalies maps
      out_dist_dir <- file.path(species_dir, "Distribution-Anomalies")
      if (!dir.exists(out_dist_dir)) dir.create(out_dist_dir, recursive = TRUE)
      out_name_dist_plot <- file.path(out_dist_dir,
        paste0(sname, "_Distribution-Anomalies_", exp_name, ".png"))
      ggsave(
        filename = out_name_dist_plot,
        plot     = six_panel,
        width    = 10, height = 10, dpi = 150, bg = "white"
      )
      
      ## Exposure factor overlap maps
      out_exp_dir <- file.path(species_dir, "Exposure-Overlap-12panel")
      if (!dir.exists(out_exp_dir)) dir.create(out_exp_dir, recursive = TRUE)
      out_name_overlap_plot <- file.path(out_exp_dir,
        paste0(sname, "_Exposure-Overlap_", exp_name, ".png"))
      ggsave(
        filename = out_name_overlap_plot,
        plot     = final_12panel,
        width    = 13, height = 17, dpi = 250, bg = "white"
      )
    }
    
    ## Optional: force a GC to release file handles on Windows immediately
    invisible(gc())
  }


## Test single species
i = 1;  process_species(shp_files[i])
  
## Run loop for all species, with error isolation per species ------------------
## Include visible errors per species 
  for (i in seq_along(shp_files)) {
    tryCatch(
      process_species(shp_files[i]),
      error = function(e) {
        message("[ERROR] process_species failed for: ", shp_files[i])
        message("        ", conditionMessage(e))
      }
    )
  }
  
## Final graphics cleanup: close any straggling devices
while (!is.null(grDevices::dev.list())) grDevices::dev.off()  

  