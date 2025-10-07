rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('functions/custom_functions.R', local=util)
setwd(wd)

library(rnaturalearth)
library(ggplot2)
library(tidyterra)
library(patchwork)
world <- vect(ne_countries(scale = "medium", returnclass = "sf"))
westam <- vect(ne_states(country = c("Canada", "Mexico", "United States of America"), returnclass = "sf"))

data <- readRDS(file.path(wd, 'output/model', 'data_24sept2025.rds'))
datasets <- readRDS(file.path(wd, 'output/model', 'datasets_24sept2025.rds'))

# Climate data
clim_pred <- readRDS(file.path(wd, "output", "climate",  "climpredictors_24sept2025.rds"))
clim_pred$soilmoist_mjj <- clim_pred$soilmoist_2m_mjj*100 # in percentage
clim_pred$gdd <- clim_pred$gdd/100 # in x100 degC
clim_pred$gdd_amjjas <- clim_pred$gdd_amjjas/100 # in x100 degC
clim_pred$vpd_mjj <- clim_pred$vpd_mjj # in hPa?
clim_pred <- na.omit(clim_pred)

# Get shapefile name
temp <- stringr::str_split_fixed(datasets$species_name, " ", 3)
datasets$species_name_simplified <- paste(temp[,1], temp[,2], sep=" ")
species_treerings <- data.frame(species = unique(datasets$species_name_simplified))
species_treerings <- data.frame(species = species_treerings[species_treerings != 'NA ',])
species_little <- read.csv(file.path(wd, 'data', 'ustreeatlas', 'Little_datatable.csv'))
species_little <- species_little[,c('Latin.Name', 'SHP..')]
names(species_little) <- c('species', 'shp')
species_little[species_little$species == 'Chamaecyparis nootkatensis', 'species'] <- 'Callitropsis nootkatensis'
species_little[species_little$species == 'Libocedrus decurrens', 'species'] <- 'Calocedrus decurrens'
species_treerings <- merge(species_treerings, species_little,  all.x = TRUE)

# GDD experience across species range
species_gddrange <- readRDS(file.path(wd, 'output', 'climate', 'species_gddrange_experienced.rds'))

for(i in 1:nrow(species_treerings)){
  
  sp <- species_treerings[i, 'species']
  shpname <- species_treerings[i, 'shp']
  if(is.na(shpname)){next}
  
  rangemap <- aggregate(vect(file.path(wd, 'data', 'ustreeatlas', 'shp', shpname, paste0(shpname, '.shp'))))
  crs(rangemap) <- 'EPSG:4267'
  records <- vect(unique(datasets[datasets$species_name_simplified == sp, c('north_lat', 'east_lon')]), 
                  geom = c('east_lon', 'north_lat'))
  crs(records) <- 'EPSG:4267'
  map <- ggplot() +
    geom_spatvector(data = westam, color = 'white', fill = 'grey90') +
    geom_spatvector(data = rangemap, color = 'white', fill = '#51ba69') +
    geom_spatvector(data = records, color = 'white', fill = '#BA6951', shape = 21, stroke = 1.1, size = 1.5) +
    ggplot2::theme_minimal() +
    theme(axis.title = ggplot2::element_blank()) +
    ggplot2::coord_sf(xlim = c(-124.925, -99.025), ylim = c(25.065,52.925), expand = FALSE)
  
  tree_species <- which(data$species_idxs == which(data$uniq_species_ids == unique(datasets[datasets$species_name_simplified %in% sp, 'species_code'])))
  obs_species <- unlist(sapply(tree_species, function(t){data$tree_start_idxs[t]:data$tree_end_idxs[t]}))
  obs_gdd <- data$gdd_amjjas_obs[obs_species]*100
  
  exp_gdd <- species_gddrange[[sp]]
  
  gdd_plot <- ggplot() +
    geom_histogram(data = data.frame(gdd = exp_gdd), aes(x = gdd, y=after_stat(density)), fill="#51ba69", alpha=0.6, color = 'white', bins = 30) +
    geom_histogram(data = data.frame(gdd = obs_gdd), aes(x = gdd, y=after_stat(density)), fill="#BA6951", alpha=0.6, color = 'white', bins = 30) +
    theme_classic() +
    labs(x = 'GDD Apr.-Sept. (°C)', y = 'Density')
  
  comb <- map + gdd_plot + plot_layout(ncol = 2, width = c(0.5,0.4))
  
  ggsave(comb, file = file.path(wd, 'figures', 'representativeness', paste0( tolower(gsub(" ", "_", sp)), '_withmap.pdf')),
         height = 4, width = 8)
  
}

# And all species?
obs_gdd <- data$gdd_amjjas_obs*100

exp_gdd <- species_gddrange[[sp]]

ggplot() +
  geom_histogram(data = data.frame(gdd = exp_gdd), aes(x = gdd, y=after_stat(density)), fill="#51ba69", alpha=0.6, color = 'white', bins = 30) +
  geom_histogram(data = data.frame(gdd = obs_gdd), aes(x = gdd, y=after_stat(density)), fill="#BA6951", alpha=0.6, color = 'white', bins = 30) +
  theme_classic() +
  labs(x = 'GDD Apr.-Sept. (°C)', y = 'Density')


