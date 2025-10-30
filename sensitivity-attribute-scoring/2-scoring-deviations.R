##------------------------------------------------------------------------------
## Caribbean CVA – Preliminary Sensitivity Attribute SD plots (HMS-style)
## Inputs:
##   - score_table_all: tidy table built from the 16 reviewer workbooks with cols:
##       SourceFile, Scorer, stock_name, row_idx,
##       Attribute_name, Attribute_type,
##       Data_quality, Scoring_rank_1..4
## Output:
##   - ./outputs/prework/plots/<species>.pdf     (stacked LMHV bars by attribute)
##   - ./outputs/prework/scorer_colors.csv       (Scorer → color key; legend hidden in plots)
##   - ./outputs/prework/preliminary_sd.csv      (per-species mean SD across attributes)

rm(list = ls()); gc()
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(forcats)
  library(purrr)
})

##------------------------------------------------------------------------------
## Read the compiled scoring table

in_csv  <- "./data/preliminary-scores/score_table_all.csv"   
out_dir <- "./outputs/prework"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE, recursive = TRUE)

score_table_all <- readr::read_csv(in_csv, show_col_types = FALSE)

## Sanity
stopifnot(all(c("Scorer","stock_name","Attribute_name",
                "Scoring_rank_1","Scoring_rank_2","Scoring_rank_3","Scoring_rank_4") %in% names(score_table_all)))

##------------------------------------------------------------------------------
## Shorten/standardize attribute names (match HMS figure labels)

## Short label map (add the "or" variant)
short_map <- c(
  "Adult mobility"                                   = "Adult Mobility",
  "Habitat specificity"                              = "Habitat Specificity",
  "Mobility and dispersal of early life stages"      = "Mobility Dispersal of ELS",
  "Mobility and dispersal or early life stages"      = "Mobility Dispersal of ELS",
  "Predation and competition dynamics"               = "Pred & Comp Dynamics",
  "Prey specificity"                                 = "Prey Specificity",
  "Complexity in reproductive strategy"              = "Repro Strat Sensitivity",
  "Spawning characteristics"                         = "Reproductive Cycle",
  "Specificity in early life history requirements"   = "Specificity EL Hist REQs",
  "Stock size/status"                                = "Stock Size Status",
  "Tolerance to ocean acidification"                 = "Sensitivity to OA",
  "Population growth rate"                           = "Pop Growth Rate",
  "Sensitivity to temperature"                       = "Sensitivity to Temp",
  "Species range"                                    = "Species Range",
  "Other stressors"                                  = "Other Stressors",
  "Genetic diversity"                                = "Genetic Diversity"
)

## Create short labels from Attribute_name
score_table_all <- score_table_all %>%
  mutate(
    Attribute_name  = str_squish(Attribute_name),
    Attribute_short = recode(Attribute_name, !!!short_map, .default = Attribute_name)
  )

## Optional: enforce attribute display order (edit as desired)
#att_order <- c(
#  "Adult Mobility","Habitat Specificity","Mobility_Dispersal of ELS",
#  "Other Stressors","Pop Growth Rate","Prey Specificity","Reproductive Cycle",
#  "Repro Strat Sensitivity","Sensitivity to OA","Sensitivity to Temp",
#  "Site Fidelity","Specificity EL Hist REQs","Stock Size_Status",
#  "Species Range","Genetic Diversity","Pred & Comp Dynamics"
#)
#score_table_all <- score_table_all |>
#  mutate(Attribute_short = factor(Attribute_short, levels = att_order))

##------------------------------------------------------------------------------
##  Prepare long format for plotting and HMS calculations
##   - one row per Scorer × Species × Attribute × Category (L/M/H/V)
##   - num_experts = number of scorers who picked any category for that attribute

## Tallies long-form: one row per Scorer × Species × Attribute × category, val ∈ [0..5]
lmhv_long <- score_table_all %>%
  dplyr::select(Scorer, stock_name, Attribute_name, Attribute_short,
                Scoring_rank_1, Scoring_rank_2, Scoring_rank_3, Scoring_rank_4) %>%
  tidyr::pivot_longer(cols = dplyr::starts_with("Scoring_rank_"),
                      names_to = "category", values_to = "val") %>%
  dplyr::mutate(
    category = dplyr::recode(category,
                             "Scoring_rank_1" = "L",
                             "Scoring_rank_2" = "M",
                             "Scoring_rank_3" = "H",
                             "Scoring_rank_4" = "V"),
    val = suppressWarnings(as.numeric(val))
  ) %>%
  tidyr::replace_na(list(val = 0)) %>%
  dplyr::mutate(category = factor(category, levels = c("L","M","H","V")))

## num_experts = distinct scorers who put ANY tallies (>0) on that attribute
experts_per_attr <- lmhv_long %>%
  dplyr::group_by(stock_name, Attribute_name, Scorer) %>%
  dplyr::summarise(tot = sum(val, na.rm = TRUE), .groups = "drop") %>%
  dplyr::group_by(stock_name, Attribute_name) %>%
  dplyr::summarise(num_experts = dplyr::n_distinct(Scorer[tot > 0]), .groups = "drop")


##------------------------------------------------------------------------------
## Compute per-species × attribute HMS mean and SD
##   - L,M,H,V = sums across scorers
##   - Weighted mean (HMS): ((L*1 + M*2 + H*3 + V*4) / (num_experts * 5))
##   - SD: sd(rep(1,L), rep(2,M), rep(3,H), rep(4,V))


## Per-stock × attribute Data Quality mean/sd across scorers
dq_stats <- score_table_all %>%
  group_by(stock_name, Attribute_name, Attribute_short) %>%
  summarise(
    dq_mean = round(mean(Data_quality, na.rm = TRUE),3),
    dq_sd   = round(sd(Data_quality,   na.rm = TRUE),3),
#    n_dq    = sum(!is.na(Data_quality)), ## number of scoreres
    .groups = "drop"
  )


## Sum tallies across scorers by category
## Ensure Data_quality is numeric once
score_table_all <- score_table_all %>%
  dplyr::mutate(Data_quality = suppressWarnings(as.numeric(Data_quality)))

## Sum tallies across scorers by category + add num_experts, att_mean/att_sd, dq_mean/dq_sd
lmhv_sums <- lmhv_long %>%
  dplyr::group_by(stock_name, Attribute_name, Attribute_short, category) %>%
  dplyr::summarise(n = sum(val, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = category, values_from = n, values_fill = 0) %>%
  ## ensure L/M/H/V exist even if a category was absent
  { for (nm in c("L","M","H","V")) if (!nm %in% names(.)) .[[nm]] <- 0L; . } %>%
  ## num_experts = distinct scorers who placed any tallies on this attribute
  dplyr::left_join(
    lmhv_long %>%
      dplyr::group_by(stock_name, Attribute_name) %>%
      dplyr::summarise(num_experts = dplyr::n_distinct(Scorer[val > 0]), .groups = "drop"),
    by = c("stock_name","Attribute_name")
  ) %>%
  dplyr::mutate(
    lmhv_mean = dplyr::if_else(num_experts > 0,
                              ((L*1) + (M*2) + (H*3) + (V*4)) / (num_experts * 5),
                              NA_real_)
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    lmhv_sd = {
      vec <- c(rep(1, L), rep(2, M), rep(3, H), rep(4, V))
      if (length(vec) >= 2) stats::sd(vec) else NA_real_
    }
  ) %>%
  dplyr::ungroup() %>%
  ## Data Quality mean/sd per stock × attribute (inline)
  dplyr::left_join(
    score_table_all %>%
      dplyr::group_by(stock_name, Attribute_name) %>%
      dplyr::summarise(
        dq_mean = mean(Data_quality, na.rm = TRUE),
        dq_sd   = stats::sd(Data_quality, na.rm = TRUE),
        .groups = "drop"
      ),
    by = c("stock_name","Attribute_name")
  )

write.csv(lmhv_sums, file.path(out_dir, "sp-x-att_score_summaries.csv"), row.names = F)

##------------------------------------------------------------------------------
## Build Scorer → color key (saved to CSV). Legend stays hidden in plots.

library(pals)
scorers <- sort(unique(score_table_all$Scorer))
pal <- pals::glasbey(length(scorers))   # returns exactly N distinct colors

scorer_colors <- tibble::tibble(
  Scorer = scorers,
  color  = pal
)

## Order scorers A–Z and lay them out in a fixed grid
scorer_colors_plot <- scorer_colors %>%
  arrange(Scorer) %>%
  mutate(
    ## choose number of columns; tweak if you have many scorers
    ncol = 4L,
    idx  = row_number(),
    row  = (idx - 1L) %/% ncol + 1L,
    col  = (idx - 1L) %%  ncol + 1L
  )

## Use geom_tile to draw swatches + text labels next to them
p_key <- ggplot(scorer_colors_plot, aes(x = col, y = -row)) +
  ## color tile
  geom_tile(aes(fill = I(color)), width = 0.25, height = 0.8, position = position_nudge(x = -0.35)) +
  ## label
  geom_text(aes(label = Scorer), hjust = 0, nudge_x = -0.05, size = 4.0, color = "black") +
  ## title + styling
  labs(
    title = "Scorer Color Key",
  ) +
  ## remove axes, fix coord so spacing looks nice
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.10))) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text  = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 11)
  )
p_key

##------------------------------------------------------------------------------
## Plot (into make per-species panel 
##   - stacked bars L/M/H/V by scorer (colors), one facet per attribute
##   - annotate facet labels with m= and sd= (rounded to 2)
##   - legend hidden; separate scorer_colors.csv is saved already

## --- Build df_counts for plotting and give it Attribute_short
df_counts <- lmhv_long %>%
  filter(stock_name == stock_name) %>%
  mutate(
    Attribute_short = recode(Attribute_name, !!!short_map, .default = Attribute_name),
    Attribute_short = forcats::fct_inorder(Attribute_short),
    category        = factor(category, levels = c("L","M","H","V"))
  )

## --- Make df_annot no Attribute_short here, just label text
df_annot <- lmhv_sums %>%
  filter(stock_name == stock_name) %>%
  mutate(lab = sprintf("m= %.1f  sd= %.2f", round(lmhv_mean, 1), round(lmhv_sd, 2))) %>%
  distinct(stock_name, Attribute_name, .keep_all = TRUE) %>%  # 1 row per attribute
  select(stock_name, Attribute_name, lab)

## --- Join back on BOTH keys; Attribute_short remains from df_counts
df_counts <- df_counts %>%
  left_join(df_annot, by = c("stock_name","Attribute_name"))


## --- Plotting function ---
make_species_plot <- function(plot_stock) {
  
  #plot_stock <- "Atlantic thread herring" ## for testing
  
  ## Filter to one stock + enforce L-M-H-V order
  df_counts_plot <- df_counts %>%
    dplyr::filter(stock_name == plot_stock) %>%
    dplyr::mutate(category = factor(category, levels = c("L","M","H","V")))
  stopifnot(length(unique(df_counts_plot$stock_name)) == 1)

  ## Label positions per facet (max stacked tallies + margin) — use *filtered* data
  lab_pos <- df_counts_plot %>%
    dplyr::group_by(Attribute_short, category) %>%
    dplyr::summarise(n = sum(val, na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(Attribute_short) %>%
    dplyr::summarise(ymax = max(n, na.rm = TRUE), .groups = "drop") %>%
    dplyr::left_join(df_counts_plot %>% dplyr::distinct(Attribute_short, lab),
                     by = "Attribute_short") %>%
    dplyr::mutate(x = 2.5, y = ymax + 0.5)

  ## Colors per scorer (universal)
  cols <- scorer_colors$color
  names(cols) <- scorer_colors$Scorer

  ## Facet strip labels: "<Attribute>\n<m and sd>"
  strip_labs <- df_counts_plot %>%
    dplyr::distinct(Attribute_short, lab) %>%
    dplyr::mutate(strip = paste0(as.character(Attribute_short), "\n", lab))
  lab_map <- setNames(strip_labs$strip, strip_labs$Attribute_short)

  ## Aggregate to one row per Scorer × category × facet (clean borders)
  df_counts_agg <- df_counts_plot %>%
    dplyr::group_by(stock_name, Attribute_short, Scorer, category, lab) %>%
    dplyr::summarise(val = sum(val, na.rm = TRUE), .groups = "drop")

  ## Compute a single, universal y-limit across all facets for this stock
  y_max <- df_counts_plot %>%
    dplyr::group_by(Attribute_short, category) %>%
    dplyr::summarise(n = sum(val, na.rm = TRUE), .groups = "drop") %>%
    dplyr::summarise(max_n = max(n, na.rm = TRUE)) %>%
    dplyr::pull(max_n)
  y_max <- ifelse(is.finite(y_max), y_max, 0)
  
  ## Make page title
  n_scorers <- length(unique(df_counts_plot$Scorer))
  page_title <- paste0(plot_stock, " (", n_scorers, " scorers)")
  
  ##  --- Plot ---
  p <- ggplot(df_counts_agg, aes(x = category, y = val, fill = Scorer)) +
    geom_col(width = 0.85, position = "stack", color = "black", linewidth = 0.3) +
    facet_wrap(~ Attribute_short, ncol = 3, labeller = as_labeller(lab_map)) +
    labs(title = page_title, x = NULL, y = NULL) +
    scale_x_discrete(limits = c("L","M","H","V")) +
    scale_fill_manual(values = cols, guide = "none") +
    
    ## Universal y-axis + whole-number ticks + small headroom
    scale_y_continuous(
      limits = c(0, y_max),
      expand = expansion(mult = c(0, 0.06)),
      breaks = scales::breaks_pretty(n = 6),
      labels = scales::number_format(accuracy = 1)  ## labels shown as whole numbers
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.text.x = element_text(size = 10, face = "bold", lineheight = 0.95),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(size = 9, color = "black"),
      axis.text.y = element_text(size = 9, color = "black"),
      plot.title  = element_text(size = 14, face = "bold")
    )
  p
}

## Test 
first_species <- lmhv_long %>% dplyr::distinct(stock_name) %>% dplyr::arrange(stock_name) %>% dplyr::pull(stock_name) %>% .[1]
print(first_species)
p <- make_species_plot(first_species)
print(p)


## ------------------------------------------------------------
## Compile plots into a (multi-page) PDF

dir.create(dirname(out_dir), showWarnings = FALSE, recursive = TRUE) ## Make directory if it doesnt exist
out_pdf <- file.path(out_dir, "prework_all_species.pdf")
species_vec <- unique(score_table_all$stock_name) ## Set up species vector

## Set up PDF writer
grDevices::pdf(out_pdf, width = 8.5, height = 11)
on.exit(grDevices::dev.off(), add = TRUE)

## Loop through stock names  
for (sp in species_vec) {
  p <- tryCatch(make_species_plot(sp), error = function(e) NULL)
  if (!is.null(p)) {
    print(p)
  } else {
    message("Skipped (no plot or error): ", sp)
  }
}
dev.off()
message("Wrote PDF: ", out_pdf)





##------------------------------------------------------------------------------
##Export per-species preliminary S 
prelim_sd <- lmhv_sums %>%
  group_by(stock_name) %>%
  summarise(
    prelim_sd = round(mean(HMS_sd, na.rm = TRUE), 2),
    n_attributes = sum(!is.na(HMS_sd)),
    .groups = "drop"
  ) %>%
  arrange(stock_name)

readr::write_csv(prelim_sd, file.path(out_dir, "preliminary_sd.csv"))
message("Wrote: ", file.path(out_dir, "preliminary_sd.csv"))
message("Wrote: ", file.path(out_dir, "scorer_colors.csv"))
