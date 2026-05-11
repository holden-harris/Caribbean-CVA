##------------------------------------------------------------------------------
## Org:     Isla Mar, Harris Analytics and Research LLC
## Project: Caribbean CVA
## Contact: Holden Earl Harris | holden.earl.harris@gmail.com
## Code:    Query FishBase database for environmental and trait data
##------------------------------------------------------------------------------

# Clear environment and load packages
rm(list = ls()); gc()
library(dplyr)
library(rfishbase)

#------------------------------------------------------------------------------
# Set folder directories
dir_in  <- "./data/master-lists/"
dir_out <- "./outputs"

#------------------------------------------------------------------------------
# Load species list
sp_list <- read.csv(paste0(dir_in, "./species-list.csv"))  
# A species list will need a column with 'scientific_name' to match in FishBase

#------------------------------------------------------------------------------
# Query FishBase tables
# These API calls may take ~1 minute depending on list size
#------------------------------------------------------------------------------
spp_raw     <- species(sp_list$scientific_name)
growth_raw  <- popgrowth(sp_list$scientific_name)
ecology_raw <- ecology(sp_list$scientific_name)
repro_raw   <- reproduction(sp_list$scientific_name)

#------------------------------------------------------------------------------
# Select and join species-level attributes
selected_cols_species <- c(
  "Species", "FBname", "BodyShapeI",
  "DepthRangeShallow", "DepthRangeDeep", 
  "DepthRangeComShallow", "DepthRangeComDeep", 
  "LongevityWild", "LongevityCaptive", 
  "Vulnerability", "VulnerabilityClimate", 
  "Length", "CommonLength", "Weight"
)

sp_attributes <- sp_list %>% # join selected columns into sp_list
  left_join(select(spp_raw, all_of(selected_cols_species)),
            by = c("scientific_name" = "Species"))

#------------------------------------------------------------------------------
# Summarize growth data (K, Loo) by species
growth_summary <- growth_raw %>%
  group_by(Species) %>%
  summarise(
    n_growth_studies = n(),
    K_avg   = mean(K, na.rm = TRUE),
    Loo_avg = mean(Loo, na.rm = TRUE),
    .groups = "drop"
  )

sp_attributes <- sp_attributes %>% # left-join into species attributes dataframe
  left_join(growth_summary, by = c("scientific_name" = "Species"))

#------------------------------------------------------------------------------
# Select and join reproduction fields
selected_cols_repro <- c(
  "Species", "ReproMode", "Spawning", 
  "RepGuild1", "RepGuild2", "AddInfos"
)

sp_attributes <- sp_attributes %>%
  left_join(select(repro_raw, all_of(selected_cols_repro)),
            by = c("scientific_name" = "Species"))

#------------------------------------------------------------------------------
# Select and join ecological traits (trophic level, remarks)
selected_cols_ecology <- c(
  "Species", "DietTroph", "DietSeTroph", 
  "DietRemark", "AddRems"
)

sp_attributes <- sp_attributes %>%
  left_join(select(ecology_raw, all_of(selected_cols_ecology)),
            by = c("scientific_name" = "Species"))

#------------------------------------------------------------------------------
# Write out table
write.csv(sp_attributes, 
          paste0(dir_out, "fishbase_species_attributes.csv"), row.names = FALSE)
