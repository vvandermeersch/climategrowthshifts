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
library(terra)
library(tidyterra)
library(patchwork)
library(colorspace)
world <- vect(ne_countries(scale = "medium", returnclass = "sf"))
westam <- vect(ne_states(country = c("Canada", "Mexico", "United States of America"), returnclass = "sf"))

# Model data
data <- readRDS(file.path(wd, 'output/model', 'data_24sept2025.rds'))
datasets <- readRDS(file.path(wd, 'output/model', 'datasets_24sept2025.rds'))
group_keys <- interaction(
  round(datasets$north_lat,2),
  round(datasets$south_lat,2),
  round(datasets$east_lon,2),
  round(datasets$west_lon,2),
  drop = TRUE
)
datasets$grouped_stand <- paste0("S", as.integer(factor(group_keys)))
base_samples <- readRDS(file.path(wd, "output/model/samples_24sept2025_partialpooling_2clades_centered_extended.rds"))

# Climate data
clim_pred <- readRDS(file.path(wd, "output", "climate",  "climpredictors_07oct2025.rds"))
clim_pred$soilmoist_mjj <- clim_pred$soilmoist_2m_mjj*100 # in percentage
clim_pred$gdd <- clim_pred$gdd/100 # in x100 degC
clim_pred$gdd_amjjas <- clim_pred$gdd_amjjas/100 # in x100 degC
clim_pred$vpd_mjj <- clim_pred$vpd_mjj # in hPa?
clim_pred <- na.omit(clim_pred)
clim_pred <- merge(clim_pred, datasets[,c('dataset', 'grouped_stand', 'east_lon', 'north_lat')])

gdd0 = 10;
sm0 = 25;
vpd0 = 8;

# Species
species <- unique(datasets[,c('species_code', 'species_name')])
species[species$species_code == 'PISA', 'species_name'] <- "Pinus sabiniana"
temp <- stringr::str_split_fixed(species$species_name, " ", 3)
species$species_name_simplified <- paste(temp[,1], temp[,2], sep=" ")
species <- species[species$species_code %in% data$uniq_species_ids,]
species$subclade <- NA
species[grepl('Pinus', species$species_name_simplified), 'subclade'] <- 'Pinus'
species[grepl('Larix|Pseudotsuga', species$species_name_simplified), 'subclade'] <- 'Laricoideae'
species[grepl('Picea', species$species_name_simplified), 'subclade'] <- 'Picea'
species[grepl('Abies|Tsuga', species$species_name_simplified), 'subclade'] <- 'Abietoideae'
species[grepl('Juniperus|Callitropsis|Calocedrus|Chamaecyparis|Thuja|Taxodium|Sequoiadendron', species$species_name_simplified), 'subclade'] <- 'Cupressaceae'
species[grepl('Quercus', species$species_name_simplified), 'subclade'] <- 'Quercus'
species[grepl('Populus', species$species_name_simplified), 'subclade'] <- 'Populus'
species[grepl('Platanus', species$species_name_simplified), 'subclade'] <- 'Platanus'

# Climate slopes
betas_climate <- data.frame(
  species = data$uniq_species_ids,
  clade = data$clade_idxs,
  alpha_q5 = util$ensemble_mcmc_quantile_est(base_samples[['alpha']], 0.05),
  alpha_q50 = util$ensemble_mcmc_quantile_est(base_samples[['alpha']], 0.50),
  alpha_q95 = util$ensemble_mcmc_quantile_est(base_samples[['alpha']], 0.95),
  betagdd_q5 = sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('beta_gdd[', sp, ']')]], 0.05)),
  betagdd_q50 = sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('beta_gdd[', sp, ']')]], 0.5)),
  betagdd_q95= sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('beta_gdd[', sp, ']')]], 0.95)),
  betasm_q5 = sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('beta_sm[', sp, ']')]], 0.05)),
  betasm_q50 = sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('beta_sm[', sp, ']')]], 0.5)),
  betasm_q95= sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('beta_sm[', sp, ']')]], 0.95)),
  betavpd_q5 = sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('beta_vpd[', sp, ']')]], 0.05)),
  betavpd_q50 = sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('beta_vpd[', sp, ']')]], 0.5)),
  betavpd_q95= sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('beta_vpd[', sp, ']')]], 0.95))
)
betas_climate <- merge(betas_climate, species, by.x = 'species', by.y = 'species_code')

# Look at one species first
species_little <- read.csv(file.path(wd, 'data', 'ustreeatlas', 'Little_datatable.csv'))
species_little <- species_little[,c('Latin.Name', 'SHP..')]
names(species_little) <- c('species_name_simplified', 'shp')
species_little[species_little$species_name_simplified == 'Chamaecyparis nootkatensis', 'species_name_simplified'] <- 'Callitropsis nootkatensis'
species_little[species_little$species_name_simplified == 'Libocedrus decurrens', 'species_name_simplified'] <- 'Calocedrus decurrens'
species <- merge(species, species_little,  all.x = TRUE)

# Absolute
for(i in 1:nrow(species)){
  
  speciesid <- species[i, 'species_code']
  shpname <- species[i, 'shp']
  if(is.na(shpname)){next}
  
  focustrees <- which(data$species_idxs == which(data$uniq_species_ids == speciesid))
  focusstands <-  data$uniq_stand_ids[unique(data$stand_idxs[focustrees])]
  growth_trends <- unique(clim_pred[clim_pred$grouped_stand %in% focusstands, 
                                    c('grouped_stand', 'east_lon', 'north_lat', 'year',
                                      'gdd', 'soilmoist_mjj', 'vpd_mjj')])
  growth_trends$log_rw_onlyclimate <- betas_climate[betas_climate$species == speciesid, 'alpha_q50'] +
    betas_climate[betas_climate$species == speciesid, 'betagdd_q50'] * (growth_trends$gdd-gdd0) +
    betas_climate[betas_climate$species == speciesid, 'betasm_q50'] * (growth_trends$soilmoist_mjj-sm0) +
    betas_climate[betas_climate$species == speciesid, 'betavpd_q50'] * (growth_trends$vpd_mjj-vpd0)
  growth_trends$rw_onlyclimate <- exp(growth_trends$log_rw_onlyclimate)
  growth_trends_before2000 <- aggregate(rw_onlyclimate ~ grouped_stand + east_lon + north_lat, 
                                        data = growth_trends[growth_trends$year < 2000,], FUN = mean)
  names(growth_trends_before2000)[4] <- 'growth_before2000'
  growth_trends_after2000 <- aggregate(rw_onlyclimate ~ grouped_stand + east_lon + north_lat, 
                                       data = growth_trends[growth_trends$year >= 2000,], FUN = mean)
  names(growth_trends_after2000)[4] <- 'growth_after2000'
  growth_trends_summ <- merge(growth_trends_before2000, growth_trends_after2000)
  growth_trends_summ$trends <- growth_trends_summ$growth_after2000 - growth_trends_summ$growth_before2000
  
  print(range(growth_trends_summ$trends))
  
  rangemap <- aggregate(vect(file.path(wd, 'data', 'ustreeatlas', 'shp', shpname, paste0(shpname, '.shp'))))
  crs(rangemap) <- 'EPSG:4267'
  records <- vect(growth_trends_summ, 
                  geom = c('east_lon', 'north_lat'))
  crs(records) <- 'EPSG:4267'
  trend_points <- ggplot() +
    geom_spatvector(data = westam, color = 'white', fill = 'grey95', linewidth = 0.7) +
    geom_spatvector(data = rangemap, color = 'white', fill = '#51ba69', alpha = 0.2) +
    geom_spatvector(data = records, aes(fill = trends), color = 'black', shape = 21, stroke = 0.5, size = 4) +
    scale_fill_continuous_diverging("Blue-Red 3", rev = TRUE,
                                    name = 'Growth trends, climate only\n(after 2000 - before 2000, mm)') +
    ggplot2::theme_void() +
    theme(axis.title = ggplot2::element_blank(), legend.position = 'inside',
          legend.position.inside = c(0.25, 0.15), legend.direction = 'horizontal',
          legend.title.position = 'top', legend.key.width = unit(1, 'cm'),
          legend.key.height = unit(0.4, 'cm'), legend.ticks = element_blank(),
          legend.frame = element_rect(color = 'black'), panel.border = element_rect(color = 'black')) +
    ggplot2::coord_sf(xlim = c(-126.925, -99.025), ylim = c(22.065,52.925), expand = FALSE)
  
  growth_trends_summ$trends <- ifelse(growth_trends_summ$trends > 0.2, 0.2, ifelse(growth_trends_summ$trends < -0.2, -0.2, growth_trends_summ$trends))
  records <- vect(growth_trends_summ, 
                  geom = c('east_lon', 'north_lat'))
  crs(records) <- 'EPSG:4267'
  trend_hex <- ggplot() +
    geom_spatvector(data = terra::aggregate(westam), color = 'grey50', fill = 'grey97', linewidth = 0.2) +
    geom_spatvector(data = rangemap, color = NA, fill = '#51ba69', alpha = 0.2) +
    stat_summary_hex(aes(x = east_lon, y = north_lat, z = trends), 
                     data = growth_trends_summ, 
                     fun = mean, binwidth = 2, alpha = 1, color = 'white')+
    geom_spatvector(data = rangemap, color = 'grey50', fill = NA) +
    geom_spatvector(data = records, aes(fill = trends), color = 'black', shape = 21, stroke = 0.5, size = 2) +
    scale_fill_continuous_diverging("Blue-Red 3", rev = TRUE, 
                                    name = 'Growth trends, climate only\n(after 2000 - before 2000, mm)',
                                    limits = c(-0.2,0.2), breaks = seq(-0.2,0.2,0.1), labels = c('<-0.2', '-0.1', '0', '0.1', '>0.2')) +
    ggplot2::theme_void() +
    theme(axis.title = ggplot2::element_blank(), legend.position = 'inside',
          legend.position.inside = c(0.25, 0.15), legend.direction = 'horizontal',
          legend.title.position = 'top', legend.key.width = unit(1, 'cm'),
          legend.key.height = unit(0.4, 'cm'), legend.ticks = element_blank(),
          legend.frame = element_rect(color = 'black'), panel.border = element_rect(color = 'black', fill = NA)) +
    ggplot2::coord_sf(xlim = c(-126.925, -99.025), ylim = c(22.065,52.925), expand = FALSE)
  
  comb <- trend_points + plot_spacer() + trend_hex + plot_layout(width = c(1,0.05,1))
  ggsave(comb, file = file.path(wd, 'figures', 'newyork2025', 'climate_growth_trends', 'per_species', 
                                paste0(species[species$species_code == speciesid, 'species_name_simplified'], '.pdf')), width = 12, height = 8)
  
  
}

# With min and max time series
for(i in 1:nrow(species)){
  
  speciesid <- species[i, 'species_code']
  
  focustrees <- which(data$species_idxs == which(data$uniq_species_ids == speciesid))
  focusstands <-  data$uniq_stand_ids[unique(data$stand_idxs[focustrees])]
  growth_trends <- unique(clim_pred[clim_pred$grouped_stand %in% focusstands, 
                                    c('grouped_stand', 'east_lon', 'north_lat', 'year',
                                      'gdd', 'soilmoist_mjj', 'vpd_mjj')])
  growth_trends$log_rw_onlyclimate <- betas_climate[betas_climate$species == speciesid, 'alpha_q50'] +
    betas_climate[betas_climate$species == speciesid, 'betagdd_q50'] * (growth_trends$gdd-gdd0) +
    betas_climate[betas_climate$species == speciesid, 'betasm_q50'] * (growth_trends$soilmoist_mjj-sm0) +
    betas_climate[betas_climate$species == speciesid, 'betavpd_q50'] * (growth_trends$vpd_mjj-vpd0)
  growth_trends$rw_onlyclimate <- exp(growth_trends$log_rw_onlyclimate)
  growth_trends_before2000 <- aggregate(rw_onlyclimate ~ grouped_stand + east_lon + north_lat, 
                                        data = growth_trends[growth_trends$year < 2000,], FUN = mean)
  names(growth_trends_before2000)[4] <- 'growth_before2000'
  growth_trends_after2000 <- aggregate(rw_onlyclimate ~ grouped_stand + east_lon + north_lat, 
                                       data = growth_trends[growth_trends$year >= 2000,], FUN = mean)
  names(growth_trends_after2000)[4] <- 'growth_after2000'
  growth_trends_summ <- merge(growth_trends_before2000, growth_trends_after2000)
  growth_trends_summ$trends <- growth_trends_summ$growth_after2000 - growth_trends_summ$growth_before2000
  
  max_trend_stand <- growth_trends_summ[growth_trends_summ$trends == max(growth_trends_summ$trends), 'grouped_stand']
  min_trend_stand <- growth_trends_summ[growth_trends_summ$trends == min(growth_trends_summ$trends), 'grouped_stand']
  
  growth_trends_summ$trends <- ifelse(growth_trends_summ$trends > 0.2, 0.2, ifelse(growth_trends_summ$trends < -0.2, -0.2, growth_trends_summ$trends))
  records <- vect(growth_trends_summ, 
                  geom = c('east_lon', 'north_lat'))
  crs(records) <- 'EPSG:4267'
  trend_hex <- ggplot() +
    geom_spatvector(data = terra::aggregate(westam), color = 'grey50', fill = 'grey97', linewidth = 0.2) +
    stat_summary_hex(aes(x = east_lon, y = north_lat, z = trends), 
                     data = growth_trends_summ, 
                     fun = mean, binwidth = 2, alpha = 1, color = 'white')+
    geom_spatvector(data = records[records$grouped_stand %in% c(max_trend_stand, min_trend_stand)], aes(fill = trends), color = 'black', shape = 21, stroke = 0.5, size = 4) +
    scale_fill_continuous_diverging("Blue-Red 3", rev = TRUE, 
                                    name = 'Growth trends, climate only\n(after 2000 - before 2000, mm)',
                                    limits = c(-0.2,0.2), breaks = seq(-0.2,0.2,0.1), labels = c('<-0.2', '-0.1', '0', '0.1', '>0.2')) +
    ggplot2::theme_void() +
    theme(axis.title = ggplot2::element_blank(), legend.position = 'inside',
          legend.position.inside = c(0.25, 0.12), legend.direction = 'horizontal',
          legend.title.position = 'top', legend.key.width = unit(1, 'cm'),
          legend.key.height = unit(0.4, 'cm'), legend.ticks = element_blank(),
          legend.frame = element_rect(color = 'black'), panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.8)) +
    ggplot2::coord_sf(xlim = c(-126.925, -99.025), ylim = c(22.065,52.925), expand = FALSE)
  
  ylimits <- c(floor(min(growth_trends[growth_trends$grouped_stand %in% c(max_trend_stand, min_trend_stand), 'rw_onlyclimate'])*10)/10,
               floor(max(growth_trends[growth_trends$grouped_stand %in% c(max_trend_stand, min_trend_stand), 'rw_onlyclimate'])*10+1)/10)
  
  maxtrend_ts <- ggplot(data = growth_trends[growth_trends$grouped_stand == max_trend_stand, ]) +
    geom_vline(aes(xintercept = 2000), linetype = 'dotted', color = 'grey70') +
    geom_smooth(aes(x = year, y = rw_onlyclimate), method ='loess', se = FALSE,
                linetype = "dashed", color = "#9baee2") +
    geom_line(aes(x = year, y = rw_onlyclimate)) +
    geom_point(aes(x = year, y = rw_onlyclimate)) +
    scale_y_continuous(position = 'right', limits = ylimits) + 
    theme_classic() +
    labs(x = '', y = 'Ring width (mm)')
  
  mintrend_ts <-ggplot(data = growth_trends[growth_trends$grouped_stand == min_trend_stand, ]) +
    geom_vline(aes(xintercept = 2000), linetype = 'dotted', color = 'grey70') +
    geom_smooth(aes(x = year, y = rw_onlyclimate), method ='loess', se = FALSE,
                linetype = "dashed", color = "#db8e8e") +
    geom_line(aes(x = year, y = rw_onlyclimate)) +
    geom_point(aes(x = year, y = rw_onlyclimate)) +
    scale_y_continuous(position = 'right', limits = ylimits) + 
    theme_classic() +
    labs(x = '', y = 'Ring width (mm)')
  
  design <- 
    "AB
   AC"
  comb <- trend_hex + mintrend_ts + maxtrend_ts + plot_layout(design = design, width = c(0.6,0.5))
  
  ggsave(comb, file = file.path(wd, 'figures', 'newyork2025', 'climate_growth_trends', 'per_species', 'with_2timeseries',
                                paste0(species[species$species_code == speciesid, 'species_name_simplified'], '.pdf')), width = 12, height = 8)
  
  
}


# Relative
for(i in 1:nrow(species)){
  
  speciesid <- species[i, 'species_code']
  shpname <- species[i, 'shp']
  if(is.na(shpname)){next}
  
  focustrees <- which(data$species_idxs == which(data$uniq_species_ids == speciesid))
  focusstands <-  data$uniq_stand_ids[unique(data$stand_idxs[focustrees])]
  growth_trends <- unique(clim_pred[clim_pred$grouped_stand %in% focusstands, 
                                    c('grouped_stand', 'east_lon', 'north_lat', 'year',
                                      'gdd', 'soilmoist_mjj', 'vpd_mjj')])
  growth_trends$log_rw_onlyclimate <- betas_climate[betas_climate$species == speciesid, 'alpha_q50'] +
    betas_climate[betas_climate$species == speciesid, 'betagdd_q50'] * (growth_trends$gdd-gdd0) +
    betas_climate[betas_climate$species == speciesid, 'betasm_q50'] * (growth_trends$soilmoist_mjj-sm0) +
    betas_climate[betas_climate$species == speciesid, 'betavpd_q50'] * (growth_trends$vpd_mjj-vpd0)
  growth_trends$rw_onlyclimate <- exp(growth_trends$log_rw_onlyclimate)
  growth_trends_before2000 <- aggregate(rw_onlyclimate ~ grouped_stand + east_lon + north_lat, 
                                        data = growth_trends[growth_trends$year < 2000,], FUN = mean)
  names(growth_trends_before2000)[4] <- 'growth_before2000'
  growth_trends_after2000 <- aggregate(rw_onlyclimate ~ grouped_stand + east_lon + north_lat, 
                                       data = growth_trends[growth_trends$year >= 2000,], FUN = mean)
  names(growth_trends_after2000)[4] <- 'growth_after2000'
  growth_trends_summ <- merge(growth_trends_before2000, growth_trends_after2000)
  growth_trends_summ$trends <- ((growth_trends_summ$growth_after2000 - growth_trends_summ$growth_before2000)/growth_trends_summ$growth_before2000)*100
  
  print(range(growth_trends_summ$trends))
  
  rangemap <- aggregate(vect(file.path(wd, 'data', 'ustreeatlas', 'shp', shpname, paste0(shpname, '.shp'))))
  crs(rangemap) <- 'EPSG:4267'
  records <- vect(growth_trends_summ, 
                  geom = c('east_lon', 'north_lat'))
  crs(records) <- 'EPSG:4267'
  trend_points <- ggplot() +
    geom_spatvector(data = westam, color = 'white', fill = 'grey95', linewidth = 0.7) +
    geom_spatvector(data = rangemap, color = 'white', fill = '#51ba69', alpha = 0.2) +
    geom_spatvector(data = records, aes(fill = trends), color = 'black', shape = 21, stroke = 0.5, size = 4) +
    scale_fill_continuous_diverging("Blue-Red 3", rev = TRUE,
                                    name = 'Growth trends, climate only\n(after 2000 - before 2000, mm)') +
    ggplot2::theme_void() +
    theme(axis.title = ggplot2::element_blank(), legend.position = 'inside',
          legend.position.inside = c(0.25, 0.15), legend.direction = 'horizontal',
          legend.title.position = 'top', legend.key.width = unit(1, 'cm'),
          legend.key.height = unit(0.4, 'cm'), legend.ticks = element_blank(),
          legend.frame = element_rect(color = 'black'), panel.border = element_rect(color = 'black', fill = NA)) +
    ggplot2::coord_sf(xlim = c(-126.925, -99.025), ylim = c(22.065,52.925), expand = FALSE)
  
  growth_trends_summ$trends <- ifelse(growth_trends_summ$trends > 10, 10, ifelse(growth_trends_summ$trends < -10, -10, growth_trends_summ$trends))
  records <- vect(growth_trends_summ, 
                  geom = c('east_lon', 'north_lat'))
  crs(records) <- 'EPSG:4267'
  trend_hex <- ggplot() +
    geom_spatvector(data = terra::aggregate(westam), color = 'grey50', fill = 'grey97', linewidth = 0.2) +
    geom_spatvector(data = rangemap, color = NA, fill = '#51ba69', alpha = 0.2) +
    stat_summary_hex(aes(x = east_lon, y = north_lat, z = trends), 
                     data = growth_trends_summ, 
                     fun = mean, binwidth = 2, alpha = 1, color = 'white')+
    geom_spatvector(data = rangemap, color = 'grey50', fill = NA) +
    geom_spatvector(data = records, aes(fill = trends), color = 'black', shape = 21, stroke = 0.5, size = 2) +
    scale_fill_continuous_diverging("Blue-Red 3", rev = TRUE, 
                                    name = 'Growth trends, climate only\n(relative change, %)',
                                    limits = c(-10,10), breaks = seq(-10,10,5), labels = c('<-10', '-5', '0', '5', '>10')) +
    ggplot2::theme_void() +
    theme(axis.title = ggplot2::element_blank(), legend.position = 'inside',
          legend.position.inside = c(0.25, 0.15), legend.direction = 'horizontal',
          legend.title.position = 'top', legend.key.width = unit(1, 'cm'),
          legend.key.height = unit(0.4, 'cm'), legend.ticks = element_blank(),
          legend.frame = element_rect(color = 'black'), panel.border = element_rect(color = 'black', fill = NA)) +
    ggplot2::coord_sf(xlim = c(-126.925, -99.025), ylim = c(22.065,52.925), expand = FALSE)
  
  comb <- trend_points + plot_spacer() + trend_hex + plot_layout(width = c(1,0.05,1))
  ggsave(comb, file = file.path(wd, 'figures', 'newyork2025', 'climate_growth_trends', 'per_species', 'relative_change',
                                paste0(species[species$species_code == speciesid, 'species_name_simplified'], '.pdf')), width = 12, height = 8)
  
  
}

# Growth trends across all species
growth_trends_all <- data.frame()
for(i in 1:nrow(species)){
  
  speciesid <- species[i, 'species_code']
  shpname <- species[i, 'shp']
  if(is.na(shpname)){next}
  
  focustrees <- which(data$species_idxs == which(data$uniq_species_ids == speciesid))
  focusstands <-  data$uniq_stand_ids[unique(data$stand_idxs[focustrees])]
  growth_trends <- unique(clim_pred[clim_pred$grouped_stand %in% focusstands, 
                                    c('grouped_stand', 'east_lon', 'north_lat', 'year',
                                      'gdd', 'soilmoist_mjj', 'vpd_mjj')])
  growth_trends$log_rw_onlyclimate <- betas_climate[betas_climate$species == speciesid, 'alpha_q50'] +
    betas_climate[betas_climate$species == speciesid, 'betagdd_q50'] * (growth_trends$gdd-gdd0) +
    betas_climate[betas_climate$species == speciesid, 'betasm_q50'] * (growth_trends$soilmoist_mjj-sm0) +
    betas_climate[betas_climate$species == speciesid, 'betavpd_q50'] * (growth_trends$vpd_mjj-vpd0)
  growth_trends$rw_onlyclimate <- exp(growth_trends$log_rw_onlyclimate)
  growth_trends_before2000 <- aggregate(rw_onlyclimate ~ grouped_stand + east_lon + north_lat, 
                                        data = growth_trends[growth_trends$year < 2000,], FUN = mean)
  names(growth_trends_before2000)[4] <- 'growth_before2000'
  growth_trends_after2000 <- aggregate(rw_onlyclimate ~ grouped_stand + east_lon + north_lat, 
                                       data = growth_trends[growth_trends$year >= 2000,], FUN = mean)
  names(growth_trends_after2000)[4] <- 'growth_after2000'
  growth_trends_summ <- merge(growth_trends_before2000, growth_trends_after2000)
  growth_trends_summ$trends <- growth_trends_summ$growth_after2000 - growth_trends_summ$growth_before2000
  growth_trends_summ$rel_trends <- (growth_trends_summ$growth_after2000 - growth_trends_summ$growth_before2000)/growth_trends_summ$growth_before2000*100
  growth_trends_summ$species <- speciesid
  
  growth_trends_all <- rbind(growth_trends_all, growth_trends_summ)
}
records <- vect(growth_trends_all, 
                geom = c('east_lon', 'north_lat'))
crs(records) <- 'EPSG:4267'
growth_trends <- ggplot() +
  geom_spatvector(data = terra::aggregate(westam), color = 'grey50', fill = 'grey97', linewidth = 0.2) +
  stat_summary_hex(aes(x = east_lon, y = north_lat, z = trends), 
                   data = growth_trends_all, 
                   fun = mean, binwidth = 2, alpha = 1, color = 'white')+
  geom_spatvector(data = records, aes(fill = trends), color = 'black', shape = 21, stroke = 0.5, size = 2) +
  scale_fill_continuous_diverging("Blue-Red 3", rev = TRUE, 
                                  name = 'Growth trends, climate only\n(after 2000 - before 2000, mm)') +
  ggplot2::theme_void() +
  theme(axis.title = ggplot2::element_blank(), legend.position = 'inside',
        legend.position.inside = c(0.25, 0.12), legend.direction = 'horizontal',
        legend.title.position = 'top', legend.key.width = unit(0.8, 'cm'),
        legend.key.height = unit(0.3, 'cm'), legend.ticks = element_blank(),
        legend.frame = element_rect(color = 'black'), panel.border = element_rect(color = 'black', fill = NA)) +
  ggplot2::coord_sf(xlim = c(-126.925, -99.025), ylim = c(22.065,52.925), expand = FALSE)

growth_rel_trends <- ggplot() +
  geom_spatvector(data = terra::aggregate(westam), color = 'grey50', fill = 'grey97', linewidth = 0.2) +
  stat_summary_hex(aes(x = east_lon, y = north_lat, z = rel_trends), 
                   data = growth_trends_all, 
                   fun = mean, binwidth = 2, alpha = 1, color = 'white')+
  scale_fill_continuous_diverging("Blue-Red 3", rev = TRUE, 
                                  name = 'Growth trends, climate only\n(relative, %)') +
  ggplot2::theme_void() +
  theme(axis.title = ggplot2::element_blank(), legend.position = 'inside',
        legend.position.inside = c(0.25, 0.12), legend.direction = 'horizontal',
        legend.title.position = 'top', legend.key.width = unit(0.8, 'cm'),
        legend.key.height = unit(0.3, 'cm'), legend.ticks = element_blank(),
        legend.frame = element_rect(color = 'black'), panel.border = element_rect(color = 'black', fill = NA)) +
  ggplot2::coord_sf(xlim = c(-126.925, -99.025), ylim = c(22.065,52.925), expand = FALSE)

no_species <- ggplot() +
  geom_spatvector(data = terra::aggregate(westam), color = 'grey50', fill = 'grey97', linewidth = 0.2) +
  stat_summary_hex(aes(x = east_lon, y = north_lat, z = as.numeric(as.factor(species))), 
                   data = growth_trends_all, 
                   fun = function(x) length(unique(x)), binwidth = 2, alpha = 1, color = 'white')+
  stat_summary_hex(aes(x = east_lon, y = north_lat, z = as.numeric(as.factor(species)), label = after_stat(value)), 
                   data = growth_trends_all, 
                   fun = function(x) length(unique(x)), binwidth = 2, alpha = 1, color = 'white', geom = 'text')+
  scale_fill_continuous_sequential("BluGrn", rev = TRUE, 
                                   name = 'Number of species') +
  ggplot2::theme_void() +
  theme(axis.title = ggplot2::element_blank(), legend.position = 'inside',
        legend.position.inside = c(0.25, 0.12), legend.direction = 'horizontal',
        legend.title.position = 'top', legend.key.width = unit(0.8, 'cm'),
        legend.key.height = unit(0.3, 'cm'), legend.ticks = element_blank(),
        legend.frame = element_rect(color = 'black'), panel.border = element_rect(color = 'black', fill = NA)) +
  ggplot2::coord_sf(xlim = c(-126.925, -99.025), ylim = c(22.065,52.925), expand = FALSE)

no_stands <- ggplot() +
  geom_spatvector(data = terra::aggregate(westam), color = 'grey50', fill = 'grey97', linewidth = 0.2) +
  stat_summary_hex(aes(x = east_lon, y = north_lat, z = as.numeric(as.factor(grouped_stand))), 
                   data = growth_trends_all, 
                   fun = function(x) length(unique(x)), binwidth = 2, alpha = 1, color = 'white')+
  stat_summary_hex(aes(x = east_lon, y = north_lat, z = as.numeric(as.factor(grouped_stand)), label = after_stat(value)), 
                   data = growth_trends_all, 
                   fun = function(x) length(unique(x)), binwidth = 2, alpha = 1, color = 'white', geom = 'text')+
  scale_fill_continuous_sequential("Teal", rev = TRUE, 
                                   name = 'Number of stands') +
  ggplot2::theme_void() +
  theme(axis.title = ggplot2::element_blank(), legend.position = 'inside',
        legend.position.inside = c(0.25, 0.12), legend.direction = 'horizontal',
        legend.title.position = 'top', legend.key.width = unit(0.8, 'cm'),
        legend.key.height = unit(0.3, 'cm'), legend.ticks = element_blank(),
        legend.frame = element_rect(color = 'black'), panel.border = element_rect(color = 'black', fill = NA)) +
  ggplot2::coord_sf(xlim = c(-126.925, -99.025), ylim = c(22.065,52.925), expand = FALSE)

comb <- no_species + plot_spacer() + no_stands + plot_spacer() + growth_trends + plot_spacer() + growth_rel_trends + plot_layout(width = c(1,0.05,1, 0.05, 1,  0.05, 1))
ggsave(comb, file = file.path(wd, 'figures', 'newyork2025', 'climate_growth_trends', 'growth_trends.pdf'), width = 22, height = 7.5)

# Climate trends
climate_trends <- unique(clim_pred[, 
                                  c('grouped_stand', 'east_lon', 'north_lat', 'year',
                                    'gdd', 'soilmoist_mjj', 'vpd_mjj')])

climate_trends_before2000 <- aggregate(gdd ~ grouped_stand + east_lon + north_lat, 
                                      data = climate_trends[climate_trends$year < 2000,], FUN = mean)
names(climate_trends_before2000)[4] <- 'before2000'
climate_trends_after2000 <- aggregate(gdd ~ grouped_stand + east_lon + north_lat, 
                                       data = climate_trends[climate_trends$year >= 2000,], FUN = mean)
names(climate_trends_after2000)[4] <- 'after2000'
climate_trends_summ <- merge(climate_trends_before2000, climate_trends_after2000)
climate_trends_summ$trends <- climate_trends_summ$after2000 - climate_trends_summ$before2000

records <- vect(climate_trends_summ, 
                geom = c('east_lon', 'north_lat'))
crs(records) <- 'EPSG:4267'
gdd_trends <- ggplot() +
  geom_spatvector(data = terra::aggregate(westam), color = 'grey50', fill = 'grey97', linewidth = 0.2) +
  stat_summary_hex(aes(x = east_lon, y = north_lat, z = trends), 
                   data = climate_trends_summ, 
                   fun = mean, binwidth = 2, alpha = 1, color = 'white')+
  geom_spatvector(data = records, aes(fill = trends), color = 'black', shape = 21, stroke = 0.5, size = 2) +
  scale_fill_continuous_diverging("Blue-Red 3", rev = TRUE, 
                                  name = 'GDD trends (x100 degC)') +
  ggplot2::theme_void() +
  theme(axis.title = ggplot2::element_blank(), legend.position = 'inside',
        legend.position.inside = c(0.25, 0.12), legend.direction = 'horizontal',
        legend.title.position = 'top', legend.key.width = unit(0.8, 'cm'),
        legend.key.height = unit(0.3, 'cm'), legend.ticks = element_blank(),
        legend.frame = element_rect(color = 'black'), panel.border = element_rect(color = 'black', fill = NA)) +
  ggplot2::coord_sf(xlim = c(-126.925, -99.025), ylim = c(22.065,52.925), expand = FALSE)

climate_trends <- unique(clim_pred[, 
                                   c('grouped_stand', 'east_lon', 'north_lat', 'year',
                                     'gdd', 'soilmoist_mjj', 'vpd_mjj')])

climate_trends_before2000 <- aggregate(soilmoist_mjj ~ grouped_stand + east_lon + north_lat, 
                                       data = climate_trends[climate_trends$year < 2000,], FUN = mean)
names(climate_trends_before2000)[4] <- 'before2000'
climate_trends_after2000 <- aggregate(soilmoist_mjj ~ grouped_stand + east_lon + north_lat, 
                                      data = climate_trends[climate_trends$year >= 2000,], FUN = mean)
names(climate_trends_after2000)[4] <- 'after2000'
climate_trends_summ <- merge(climate_trends_before2000, climate_trends_after2000)
climate_trends_summ$trends <- climate_trends_summ$after2000 - climate_trends_summ$before2000

records <- vect(climate_trends_summ, 
                geom = c('east_lon', 'north_lat'))
crs(records) <- 'EPSG:4267'
sm_trends <- ggplot() +
  geom_spatvector(data = terra::aggregate(westam), color = 'grey50', fill = 'grey97', linewidth = 0.2) +
  stat_summary_hex(aes(x = east_lon, y = north_lat, z = trends), 
                   data = climate_trends_summ, 
                   fun = mean, binwidth = 2, alpha = 1, color = 'white')+
  geom_spatvector(data = records, aes(fill = trends), color = 'black', shape = 21, stroke = 0.5, size = 2) +
  scale_fill_continuous_diverging("Blue-Red 3", rev = TRUE, 
                                  name = 'Soil moisture trends (%)') +
  ggplot2::theme_void() +
  theme(axis.title = ggplot2::element_blank(), legend.position = 'inside',
        legend.position.inside = c(0.25, 0.12), legend.direction = 'horizontal',
        legend.title.position = 'top', legend.key.width = unit(0.8, 'cm'),
        legend.key.height = unit(0.3, 'cm'), legend.ticks = element_blank(),
        legend.frame = element_rect(color = 'black'), panel.border = element_rect(color = 'black', fill = NA)) +
  ggplot2::coord_sf(xlim = c(-126.925, -99.025), ylim = c(22.065,52.925), expand = FALSE)

climate_trends_before2000 <- aggregate(vpd_mjj ~ grouped_stand + east_lon + north_lat, 
                                       data = climate_trends[climate_trends$year < 2000,], FUN = mean)
names(climate_trends_before2000)[4] <- 'before2000'
climate_trends_after2000 <- aggregate(vpd_mjj ~ grouped_stand + east_lon + north_lat, 
                                      data = climate_trends[climate_trends$year >= 2000,], FUN = mean)
names(climate_trends_after2000)[4] <- 'after2000'
climate_trends_summ <- merge(climate_trends_before2000, climate_trends_after2000)
climate_trends_summ$trends <- climate_trends_summ$after2000 - climate_trends_summ$before2000

records <- vect(climate_trends_summ, 
                geom = c('east_lon', 'north_lat'))
crs(records) <- 'EPSG:4267'
vpd_trends <- ggplot() +
  geom_spatvector(data = terra::aggregate(westam), color = 'grey50', fill = 'grey97', linewidth = 0.2) +
  stat_summary_hex(aes(x = east_lon, y = north_lat, z = trends), 
                   data = climate_trends_summ, 
                   fun = mean, binwidth = 2, alpha = 1, color = 'white')+
  geom_spatvector(data = records, aes(fill = trends), color = 'black', shape = 21, stroke = 0.5, size = 2) +
  scale_fill_continuous_diverging("Blue-Red 3", rev = FALSE, 
                                  name = 'VPD trends (hPa)') +
  ggplot2::theme_void() +
  theme(axis.title = ggplot2::element_blank(), legend.position = 'inside',
        legend.position.inside = c(0.25, 0.12), legend.direction = 'horizontal',
        legend.title.position = 'top', legend.key.width = unit(0.8, 'cm'),
        legend.key.height = unit(0.3, 'cm'), legend.ticks = element_blank(),
        legend.frame = element_rect(color = 'black'), panel.border = element_rect(color = 'black', fill = NA)) +
  ggplot2::coord_sf(xlim = c(-126.925, -99.025), ylim = c(22.065,52.925), expand = FALSE)

comb <- gdd_trends + plot_spacer() + sm_trends + plot_spacer() + vpd_trends + plot_layout(width = c(1,0.05,1, 0.05, 1))
ggsave(comb, file = file.path(wd, 'figures', 'newyork2025', 'climate_growth_trends', 'climate_trends.pdf'), width = 14, height = 6)
