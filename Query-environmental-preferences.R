##------------------------------------------------------------------------------
##
## Org:     Isla Mar, Harris Analytics and Research LLC
## Project: Caribbean CVA
## Contact: Holden Earl Harris | holden.earl.harris@gmail.com
## Code:    Query Fishbase Database
##          Find environmental preference

rm(list=ls()); gc()
library(dplyr)
library(dplyr)
library(rfishbase)

# Load species list
sp_list <- read.csv("./species-list.csv")  

# Query environmental and trait data from fishbase database (takes a minute)
species_fb <- rfishbase::species(sp_list$scientific_name)
growth_fb  <- popgrowth(sp_list$scientific_name)
ecology_fb <- ecology(sp_list$scientific_name)
repro_fb  <- reproduction(sp_list$scientific_name)


# Just pull selected columns to select from species_db
selected_cols_species <- c(
  "Species", 
  "FBname", "BodyShapeI",
  "DepthRangeShallow", "DepthRangeDeep", 
  "DepthRangeComShallow", "DepthRangeComDeep", 
  "LongevityWild", "LongevityCaptive", 
  "Vulnerability", "VulnerabilityClimate", 
  "Length", "CommonLength", "Weight"
)

# Left join the selected columns from species_db to species_list
sp_attritbutes <- 
  sp_list %>%
  left_join(
    select(species_fb, all_of(selected_cols_species)),
    by = c("scientific_name" = "Species")
  )

# Summarize growth rates, just get Linf and and K values
growth_summary <- growth_fb %>%
  group_by(Species) %>%
  summarise(
    n_growth_studies = n(),
    K_avg = mean(K, na.rm = TRUE),
    Loo_avg = mean(Loo, na.rm = TRUE),
    .groups = "drop"
  )

sp_attritbutes <- # Join into attributes data frame
  sp_attritbutes %>%
  left_join(
    growth_summary,
    by = c("scientific_name" = "Species")
  )

# Pull selected columns to select from repro_db
selected_cols_repro <- c(
  "Species", 
  "ReproMode", "Spawning", 
  "RepGuild1", "RepGuild2", "AddInfos"
)

# Left join the selected columns from repro_db to attributes dataframe
sp_attritbutes <- 
  sp_attritbutes %>%
  left_join(
    select(repro_fb, all_of(selected_cols_repro)),
    by = c("scientific_name" = "Species")
  )

# Pull selected columns to select from ecology_db
selected_cols_ecology <- c(
  "Species", 
  "DietTroph", "DietSeTroph", 
  "DietRemark", "AddRems"
)

# Left join the selected columns from repro_db to attributes dataframe
sp_attritbutes <- 
  sp_attritbutes %>%
  left_join(
    select(ecology_fb, all_of(selected_cols_ecology)),
    by = c("scientific_name" = "Species")
  )



