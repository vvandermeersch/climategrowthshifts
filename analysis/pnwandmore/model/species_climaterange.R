rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('functions/custom_functions.R', local=util)
setwd(wd)

library(terra)

#### Observed range!

# Model data
data <- readRDS(file.path(wd, 'output/model', 'data_24sept2025.rds'))
datasets <- readRDS(file.path(wd, 'output/model', 'datasets_24sept2025.rds'))

# Climate data
clim_pred <- readRDS(file.path(wd, "output", "climate",  "climpredictors_24sept2025.rds"))
clim_pred$soilmoist_mjj <- clim_pred$soilmoist_2m_mjj*100 # in percentage
clim_pred$gdd <- clim_pred$gdd/100 # in x100 degC
clim_pred$gdd_amjjas <- clim_pred$gdd_amjjas/100 # in x100 degC
clim_pred$vpd_mjj <- clim_pred$vpd_mjj # in hPa?
clim_pred <- na.omit(clim_pred)
climate_bystand <- aggregate(cbind(gdd, soilmoist_mjj) ~ dataset, clim_pred,  function (x) quantile(x, c(0.05, 0.5, 0.95)))
climate_bystand <- do.call(data.frame, climate_bystand)
names(climate_bystand) <- c('dataset', paste0('gdd_q', c(5,50,95)), paste0('sm_q', c(5,50,95)))

# Species x stands x climate
species_bystand <- datasets[, c('dataset', 'species_code')]
species_climate_bystand <- merge(species_bystand, clim_pred)
species_climate_obs <- aggregate(cbind(gdd, soilmoist_mjj, vpd_mjj) ~ species_code, species_climate_bystand,  function (x) quantile(x, c(0.05, 0.5, 0.95)))
species_climate_obs <- do.call(data.frame, species_climate_obs)
names(species_climate_obs) <- c('species', paste0('gdd_obs_q', c(5,50,95)), paste0('sm_obs_q', c(5,50,95)), paste0('vpd_obs_q', c(5,50,95)))

no_stands_byspecies <- aggregate(dataset ~ species_code, datasets, function(x){length(unique(x))})
names(no_stands_byspecies) <- c("species", 'no_datasets')



#### "True" range!
years <- 1980:2023
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

species_gddrange <- list()
for(i in 1:nrow(species_treerings)){ 
  sp <- species_treerings[i, 'species']
  cat(paste0(sp,'\n'))
  shpname <- species_treerings[i, 'shp']
  if(is.na(shpname)){next}
  rangemap <- aggregate(vect(file.path(wd, 'data', 'ustreeatlas', 'shp', shpname, paste0(shpname, '.shp'))))
  crs(rangemap) <- 'EPSG:4267'
  
  gdd_range <- c()
  for(y in years){
    gdd <- rast(file.path(wd, 'output', 'climate', 'rasters', paste0('gdd_amjjas_',y,'.nc')))
    gdd_range_y <- mask(gdd, rangemap)
    gdd_range <- c(gdd_range, values(gdd_range_y, na.rm = T))
  }
  gdd_range <- list(gdd_range)
  names(gdd_range) <- sp
  species_gddrange <- c(species_gddrange, gdd_range)
}
saveRDS(species_gddrange, file.path(wd, 'output', 'climate', 'species_gddrange_experienced.rds'))

# All species
allspecies_rangemap <- c()
for(i in 1:nrow(species_treerings)){ 
  sp <- species_treerings[i, 'species']
  cat(paste0(sp,'\n'))
  shpname <- species_treerings[i, 'shp']
  if(is.na(shpname)){next}
  rangemap <- aggregate(vect(file.path(wd, 'data', 'ustreeatlas', 'shp', shpname, paste0(shpname, '.shp'))))
  crs(rangemap) <- 'EPSG:4267'
  
  rangemap <- crop(rangemap, ext(-124.925, -89.025, 25.065,52.925))
  if(sum(dim(rangemap)==0)==2){next}
  
  if(i == 1){
    allspecies_rangemap <- rangemap
  }else{
    allspecies_rangemap <- union(allspecies_rangemap, rangemap)
  }
}
saveRDS(allspecies_rangemap, file.path(wd, 'output', 'allspecies_rangemap.rds'))
 
gdd_range <- c()
for(y in years){
  gdd <- rast(file.path(wd, 'output', 'climate', 'rasters', paste0('gdd_amjjas_',y,'.nc')))
  gdd_range_y <- mask(gdd, allspecies_rangemap)
  gdd_range <- c(gdd_range, values(gdd_range_y, na.rm = T))
}
gdd_range <- list(gdd_range)
names(gdd_range) <- 'allspecies'
saveRDS(gdd_range, file.path(wd, 'output', 'climate', 'allspecies_gddrange_experienced.rds'))
