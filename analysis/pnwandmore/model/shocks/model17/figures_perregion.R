rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(cmdstanr)
library(rnaturalearth)
library(terra)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)

folder <- '/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model/'
gq_samples <- readRDS(file.path(folder, 'model17/gq_samples_model17_HSGP_PIPOexcept1.rds'))
data <- readRDS(file.path(folder,'data_26mar2026_PIPOexcept1_19202024.rds'))
datasets <- readRDS(file.path(folder,'datasets_26mar2026_PIPOexcept1_19202024.rds'))

world <- vect(ne_countries(scale = "medium", returnclass = "sf"))
westam <- vect(ne_states(country = c("Canada", "Mexico", "United States of America"), returnclass = "sf"))
par(mfrow = c(1,1), cex.main = 0.9, mar = c(4,4,1,1))
hill <- rast(file.path(wd, 'data/srtm', 'hill.tif'))
mtn <- crop(vect(file.path(wd, 'data/usgs/physio_shp', 'physio.shp')),
            ext(c(-125, -104.2, 30, 50)))


plot(crop(hill,  ext(c(-125, -104.2, 30, 49))),
     col = grey(c(0:100)/100), legend = F, maxcell = Inf,
     xlim = c(-125, -104.2), ylim = c(30, 49), box = FALSE,
     clip = TRUE)
plot(westam, xlim = c(-130, -104.6), ylim = c(30, 49),
     box = FALSE, clip = TRUE, col = '#d9edf090', border = 'white', add = T)
plot(aggregate(buffer(mtn[mtn$PROVINCE == 'COLORADO PLATEAUS'], width = 1)), col = '#ffe1bd50', border = 'white', add = TRUE, lwd = 2)
plot(aggregate(buffer(mtn[mtn$PROVINCE == 'SOUTHERN ROCKY MOUNTAINS'], width = 1)), col = '#ffe1bd50', border = 'white', add = TRUE, lwd = 2)
plot(aggregate(buffer(mtn[mtn$PROVINCE == 'NORTHERN ROCKY MOUNTAINS'], width = 1)), col = '#ffe1bd50', border = 'white', add = TRUE, lwd = 2)
plot(aggregate(buffer(mtn[mtn$SECTION == 'SIERRA NEVADA'], width = 1)), col = '#ffe1bd50', border = 'white', add = TRUE, lwd = 2)
plot(aggregate(buffer(mtn[mtn$SECTION %in% c('NORTHERN CASCADE MOUNTAINS', 'MIDDLE CASCADE MOUNTAINS','SOUTHERN CASCADE MOUNTAINS')], width = 1)), 
     col = '#ffe1bd50', border = 'white', add = TRUE, lwd = 2)
plot(aggregate(buffer(mtn[mtn$SECTION %in% c('SONORAN DESERT', 'MEXICAN HIGHLAND','SACRAMENTO')], width = 1)), 
     col = '#ffe1bd50', border = 'white', add = TRUE, lwd = 2)
points(x = data$uniq_stand_lon, y = data$uniq_stand_lat,
       pch = 19, col = 'white', cex = 1.3)
points(x = data$uniq_stand_lon, y = data$uniq_stand_lat,
       pch = 19, col = '#ffb861', cex = 0.7)


stands <- vect(data.frame(x = as.numeric(data$uniq_stand_lon), y = as.numeric(data$uniq_stand_lat), id = data$uniq_stand_ids))


# ------------------------------- #
# PROPORTION OF SHOCKS PER REGION #
# ------------------------------- #

par(mfrow = c(5,1))
# Desert and highlands
extent <- aggregate(buffer(mtn[mtn$SECTION %in% c('SONORAN DESERT', 'MEXICAN HIGHLAND','SACRAMENTO')], width = 1))
region <- 'desert_highlands'
region_name <- 'Southern Desert and Highlands'
stands_region <- crop(stands, extent)
for(i in 1:data$N_all_years){
  gq_samples[[paste0('num_sck_perregion[',region,',',i,']')]] <- 0
}
stands_region <- which(data$uniq_stand_ids %in% stands_region$id)
ntrees <- sum(data$stand_idxs %in% stands_region)
for(s in stands_region){
  idxs <- seq(1 + data$N_all_years * (s-1), data$N_all_years * s, 1)
  for(i in 1:length(idxs)){
    gq_samples[[paste0('num_sck_perregion[',region,',',i,']')]] <- gq_samples[[paste0('num_sck_perregion[',region,',',i,']')]] +
      gq_samples[[paste0('num_conc_sck_stdlvl[',idxs[i],']')]]/length(stands_region)*100
  }
}
idxs <- seq(1, data$N_all_years, 1)
names <- paste0('num_sck_perregion[',region, ',', idxs,']')
util$plot_disc_pushforward_quantiles(gq_samples, names, 
                                     ylab = 'Average proportion of shocks\n(% of trees within a stand)', 
                                     display_ylim = c(0, 50),
                                     xticklabs = data$all_years,
                                     main = paste0(region_name, ' (',ntrees, ' trees, ', length(stands_region), ' stands)'))

# Colorado Plateau
extent <- aggregate(buffer(mtn[mtn$PROVINCE == 'COLORADO PLATEAUS'], width = 1))
region <- 'colorado_plateau'
region_name <- 'Colorado Plateau'
stands_region <- crop(stands, extent)
for(i in 1:data$N_all_years){
  gq_samples[[paste0('num_sck_perregion[',region,',',i,']')]] <- 0
}
stands_region <- which(data$uniq_stand_ids %in% stands_region$id)
ntrees <- sum(data$stand_idxs %in% stands_region)
for(s in stands_region){
  idxs <- seq(1 + data$N_all_years * (s-1), data$N_all_years * s, 1)
  for(i in 1:length(idxs)){
    gq_samples[[paste0('num_sck_perregion[',region,',',i,']')]] <- gq_samples[[paste0('num_sck_perregion[',region,',',i,']')]] +
      gq_samples[[paste0('num_conc_sck_stdlvl[',idxs[i],']')]]/length(stands_region)*100
  }
}
idxs <- seq(1, data$N_all_years, 1)
names <- paste0('num_sck_perregion[',region, ',', idxs,']')
util$plot_disc_pushforward_quantiles(gq_samples, names, 
                                     ylab = 'Average proportion of shocks\n(% of trees within a stand)', 
                                     display_ylim = c(0, 50),
                                     xticklabs = data$all_years,
                                     main = paste0(region_name, ' (',ntrees, ' trees, ', length(stands_region), ' stands)'))

# Southern Rockies
extent <- aggregate(buffer(mtn[mtn$PROVINCE == 'SOUTHERN ROCKY MOUNTAINS'], width = 1))
region <- 'southern_rockies'
region_name <- 'Southern Rocky Mountains'
stands_region <- crop(stands, extent)
for(i in 1:data$N_all_years){
  gq_samples[[paste0('num_sck_perregion[',region,',',i,']')]] <- 0
}
stands_region <- which(data$uniq_stand_ids %in% stands_region$id)
ntrees <- sum(data$stand_idxs %in% stands_region)
for(s in stands_region){
  idxs <- seq(1 + data$N_all_years * (s-1), data$N_all_years * s, 1)
  for(i in 1:length(idxs)){
    gq_samples[[paste0('num_sck_perregion[',region,',',i,']')]] <- gq_samples[[paste0('num_sck_perregion[',region,',',i,']')]] +
      gq_samples[[paste0('num_conc_sck_stdlvl[',idxs[i],']')]]/length(stands_region)*100
  }
}
idxs <- seq(1, data$N_all_years, 1)
names <- paste0('num_sck_perregion[',region, ',', idxs,']')
util$plot_disc_pushforward_quantiles(gq_samples, names, 
                                     ylab = 'Average proportion of shocks\n(% of trees within a stand)', 
                                     display_ylim = c(0, 50),
                                     xticklabs = data$all_years,
                                     main = paste0(region_name, ' (',ntrees, ' trees, ', length(stands_region), ' stands)'))

# Sierra Nevada
extent <- aggregate(buffer(mtn[mtn$SECTION == 'SIERRA NEVADA'], width = 1))
region <- 'sierra_nevada'
region_name <- 'Sierra Nevada'
stands_region <- crop(stands, extent)
for(i in 1:data$N_all_years){
  gq_samples[[paste0('num_sck_perregion[',region,',',i,']')]] <- 0
}
stands_region <- which(data$uniq_stand_ids %in% stands_region$id)
ntrees <- sum(data$stand_idxs %in% stands_region)
for(s in stands_region){
  idxs <- seq(1 + data$N_all_years * (s-1), data$N_all_years * s, 1)
  for(i in 1:length(idxs)){
    gq_samples[[paste0('num_sck_perregion[',region,',',i,']')]] <- gq_samples[[paste0('num_sck_perregion[',region,',',i,']')]] +
      gq_samples[[paste0('num_conc_sck_stdlvl[',idxs[i],']')]]/length(stands_region)*100
  }
}
idxs <- seq(1, data$N_all_years, 1)
names <- paste0('num_sck_perregion[',region, ',', idxs,']')
util$plot_disc_pushforward_quantiles(gq_samples, names, 
                                     ylab = 'Average proportion of shocks\n(% of trees within a stand)', 
                                     display_ylim = c(0, 50),
                                     xticklabs = data$all_years,
                                     main = paste0(region_name, ' (',ntrees, ' trees, ', length(stands_region), ' stands)'))

# Cascade range
extent <- aggregate(buffer(mtn[mtn$SECTION %in% c('NORTHERN CASCADE MOUNTAINS', 'MIDDLE CASCADE MOUNTAINS','SOUTHERN CASCADE MOUNTAINS')], width = 1))
region <- 'cascade_range'
region_name <- 'Cascade Range'
stands_region <- crop(stands, extent)
for(i in 1:data$N_all_years){
  gq_samples[[paste0('num_sck_perregion[',region,',',i,']')]] <- 0
}
stands_region <- which(data$uniq_stand_ids %in% stands_region$id)
ntrees <- sum(data$stand_idxs %in% stands_region)
for(s in stands_region){
  idxs <- seq(1 + data$N_all_years * (s-1), data$N_all_years * s, 1)
  for(i in 1:length(idxs)){
    gq_samples[[paste0('num_sck_perregion[',region,',',i,']')]] <- gq_samples[[paste0('num_sck_perregion[',region,',',i,']')]] +
      gq_samples[[paste0('num_conc_sck_stdlvl[',idxs[i],']')]]/length(stands_region)*100
  }
}
idxs <- seq(1, data$N_all_years, 1)
names <- paste0('num_sck_perregion[',region, ',', idxs,']')
util$plot_disc_pushforward_quantiles(gq_samples, names, 
                                     ylab = 'Average proportion of shocks\n(% of trees within a stand)', 
                                     display_ylim = c(0, 50),
                                     xticklabs = data$all_years,
                                     main = paste0(region_name, ' (',ntrees, ' trees, ', length(stands_region), ' stands)'))



# ---------------------------- #
# IMPACTS OF SHOCKS PER REGION #
# ---------------------------- #

par(mfcol = c(5,2))
# Desert and highlands
extent <- aggregate(buffer(mtn[mtn$SECTION %in% c('SONORAN DESERT', 'MEXICAN HIGHLAND','SACRAMENTO')], width = 1))
region <- 'desert_highlands'
region_name <- 'Southern Desert and Highlands'
stands_region <- crop(stands, extent)
for(i in 1:data$N_all_years){
  gq_samples[[paste0('relchg_sck_perregion[',region,',',i,']')]] <- 0
}
stands_region <- which(data$uniq_stand_ids %in% stands_region$id)
ntrees <- sum(data$stand_idxs %in% stands_region)
first_obsyear <- min(sapply(which(data$stand_idxs %in% stands_region), function(t){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  return(min(data$years[idxs]))
}))
last_obsyear <- max(sapply(which(data$stand_idxs %in% stands_region), function(t){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  return(max(data$years[idxs]))
}))
year_idxs <- which(data$all_years %in% first_obsyear:last_obsyear)
for(s in stands_region){
  idxs <- seq(1 + data$N_all_years * (s-1), data$N_all_years * s, 1)
  for(i in 1:length(idxs)){
    gq_samples[[paste0('relchg_sck_perregion[',region,',',i,']')]] <- gq_samples[[paste0('relchg_sck_perregion[',region,',',i,']')]] +
      gq_samples[[paste0('avg_exp_m1_delta_sck_stdlvl[',idxs[i],']')]]/length(stands_region)*100
  }
}
idxs <- seq(1, data$N_all_years, 1)
names <- paste0('relchg_sck_perregion[',region, ',', idxs,']')
# non-observed years will not appear
for(i in idxs[!(idxs %in% year_idxs)]){
  gq_samples[[paste0('relchg_sck_perregion[',region,',',i,']')]] <- gq_samples[[paste0('relchg_sck_perregion[',region,',',i,']')]]+1e6 
}
util$plot_disc_pushforward_quantiles(gq_samples, names, 
                                     ylab = 'Average relative change in\ngrowth due to shocks (%)', 
                                     display_ylim = c(-50, 50),
                                     xticklabs = data$all_years,
                                     main = paste0(region_name, ' (',ntrees, ' trees, ', length(stands_region), ' stands)'))

# Colorado Plateau
extent <- aggregate(buffer(mtn[mtn$PROVINCE == 'COLORADO PLATEAUS'], width = 1))
region <- 'colorado_plateau'
region_name <- 'Colorado Plateau'
stands_region <- crop(stands, extent)
for(i in 1:data$N_all_years){
  gq_samples[[paste0('relchg_sck_perregion[',region,',',i,']')]] <- 0
}
stands_region <- which(data$uniq_stand_ids %in% stands_region$id)
ntrees <- sum(data$stand_idxs %in% stands_region)
first_obsyear <- min(sapply(which(data$stand_idxs %in% stands_region), function(t){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  return(min(data$years[idxs]))
}))
last_obsyear <- max(sapply(which(data$stand_idxs %in% stands_region), function(t){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  return(max(data$years[idxs]))
}))
year_idxs <- which(data$all_years %in% first_obsyear:last_obsyear)
for(s in stands_region){
  idxs <- seq(1 + data$N_all_years * (s-1), data$N_all_years * s, 1)
  for(i in 1:length(idxs)){
    gq_samples[[paste0('relchg_sck_perregion[',region,',',i,']')]] <- gq_samples[[paste0('relchg_sck_perregion[',region,',',i,']')]] +
      gq_samples[[paste0('avg_exp_m1_delta_sck_stdlvl[',idxs[i],']')]]/length(stands_region)*100
  }
}
idxs <- seq(1, data$N_all_years, 1)
names <- paste0('relchg_sck_perregion[',region, ',', idxs,']')
# non-observed years will not appear
for(i in idxs[!(idxs %in% year_idxs)]){
  gq_samples[[paste0('relchg_sck_perregion[',region,',',i,']')]] <- gq_samples[[paste0('relchg_sck_perregion[',region,',',i,']')]]+1e6 
}
util$plot_disc_pushforward_quantiles(gq_samples, names, 
                                     ylab = 'Average relative change in\ngrowth due to shocks (%)', 
                                     display_ylim = c(-50, 50),
                                     xticklabs = data$all_years,
                                     main = paste0(region_name, ' (',ntrees, ' trees, ', length(stands_region), ' stands)'))

# Southern Rockies
extent <- aggregate(buffer(mtn[mtn$PROVINCE == 'SOUTHERN ROCKY MOUNTAINS'], width = 1))
region <- 'southern_rockies'
region_name <- 'Southern Rocky Mountains'
stands_region <- crop(stands, extent)
for(i in 1:data$N_all_years){
  gq_samples[[paste0('relchg_sck_perregion[',region,',',i,']')]] <- 0
}
stands_region <- which(data$uniq_stand_ids %in% stands_region$id)
ntrees <- sum(data$stand_idxs %in% stands_region)
first_obsyear <- min(sapply(which(data$stand_idxs %in% stands_region), function(t){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  return(min(data$years[idxs]))
}))
last_obsyear <- max(sapply(which(data$stand_idxs %in% stands_region), function(t){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  return(max(data$years[idxs]))
}))
year_idxs <- which(data$all_years %in% first_obsyear:last_obsyear)
for(s in stands_region){
  idxs <- seq(1 + data$N_all_years * (s-1), data$N_all_years * s, 1)
  for(i in 1:length(idxs)){
    gq_samples[[paste0('relchg_sck_perregion[',region,',',i,']')]] <- gq_samples[[paste0('relchg_sck_perregion[',region,',',i,']')]] +
      gq_samples[[paste0('avg_exp_m1_delta_sck_stdlvl[',idxs[i],']')]]/length(stands_region)*100
  }
}
idxs <- seq(1, data$N_all_years, 1)
names <- paste0('relchg_sck_perregion[',region, ',', idxs,']')
# non-observed years will not appear
for(i in idxs[!(idxs %in% year_idxs)]){
  gq_samples[[paste0('relchg_sck_perregion[',region,',',i,']')]] <- gq_samples[[paste0('relchg_sck_perregion[',region,',',i,']')]]+1e6 
}
util$plot_disc_pushforward_quantiles(gq_samples, names, 
                                     ylab = 'Average relative change in\ngrowth due to shocks (%)', 
                                     display_ylim = c(-50, 50),
                                     xticklabs = data$all_years,
                                     main = paste0(region_name, ' (',ntrees, ' trees, ', length(stands_region), ' stands)'))

# Sierra Nevada
extent <- aggregate(buffer(mtn[mtn$SECTION == 'SIERRA NEVADA'], width = 1))
region <- 'sierra_nevada'
region_name <- 'Sierra Nevada'
stands_region <- crop(stands, extent)
for(i in 1:data$N_all_years){
  gq_samples[[paste0('relchg_sck_perregion[',region,',',i,']')]] <- 0
}
stands_region <- which(data$uniq_stand_ids %in% stands_region$id)
ntrees <- sum(data$stand_idxs %in% stands_region)
first_obsyear <- min(sapply(which(data$stand_idxs %in% stands_region), function(t){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  return(min(data$years[idxs]))
}))
last_obsyear <- max(sapply(which(data$stand_idxs %in% stands_region), function(t){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  return(max(data$years[idxs]))
}))
year_idxs <- which(data$all_years %in% first_obsyear:last_obsyear)
for(s in stands_region){
  idxs <- seq(1 + data$N_all_years * (s-1), data$N_all_years * s, 1)
  for(i in 1:length(idxs)){
    gq_samples[[paste0('relchg_sck_perregion[',region,',',i,']')]] <- gq_samples[[paste0('relchg_sck_perregion[',region,',',i,']')]] +
      gq_samples[[paste0('avg_exp_m1_delta_sck_stdlvl[',idxs[i],']')]]/length(stands_region)*100
  }
}
idxs <- seq(1, data$N_all_years, 1)
names <- paste0('relchg_sck_perregion[',region, ',', idxs,']')
# non-observed years will not appear
for(i in idxs[!(idxs %in% year_idxs)]){
  gq_samples[[paste0('relchg_sck_perregion[',region,',',i,']')]] <- gq_samples[[paste0('relchg_sck_perregion[',region,',',i,']')]]+1e6 
}
util$plot_disc_pushforward_quantiles(gq_samples, names, 
                                     ylab = 'Average relative change in\ngrowth due to shocks (%)', 
                                     display_ylim = c(-50, 50),
                                     xticklabs = data$all_years,
                                     main = paste0(region_name, ' (',ntrees, ' trees, ', length(stands_region), ' stands)'))

# Cascade range
extent <- aggregate(buffer(mtn[mtn$SECTION %in% c('NORTHERN CASCADE MOUNTAINS', 'MIDDLE CASCADE MOUNTAINS','SOUTHERN CASCADE MOUNTAINS')], width = 1))
region <- 'cascade_range'
region_name <- 'Cascade Range'
stands_region <- crop(stands, extent)
for(i in 1:data$N_all_years){
  gq_samples[[paste0('relchg_sck_perregion[',region,',',i,']')]] <- 0
}
stands_region <- which(data$uniq_stand_ids %in% stands_region$id)
ntrees <- sum(data$stand_idxs %in% stands_region)
first_obsyear <- min(sapply(which(data$stand_idxs %in% stands_region), function(t){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  return(min(data$years[idxs]))
}))
last_obsyear <- max(sapply(which(data$stand_idxs %in% stands_region), function(t){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  return(max(data$years[idxs]))
}))
year_idxs <- which(data$all_years %in% first_obsyear:last_obsyear)
for(s in stands_region){
  idxs <- seq(1 + data$N_all_years * (s-1), data$N_all_years * s, 1)
  for(i in 1:length(idxs)){
    gq_samples[[paste0('relchg_sck_perregion[',region,',',i,']')]] <- gq_samples[[paste0('relchg_sck_perregion[',region,',',i,']')]] +
      gq_samples[[paste0('avg_exp_m1_delta_sck_stdlvl[',idxs[i],']')]]/length(stands_region)*100
  }
}
idxs <- seq(1, data$N_all_years, 1)
names <- paste0('relchg_sck_perregion[',region, ',', idxs,']')
# non-observed years will not appear
for(i in idxs[!(idxs %in% year_idxs)]){
  gq_samples[[paste0('relchg_sck_perregion[',region,',',i,']')]] <- gq_samples[[paste0('relchg_sck_perregion[',region,',',i,']')]]+1e6 
}
util$plot_disc_pushforward_quantiles(gq_samples, names, 
                                     ylab = 'Average relative change in\ngrowth due to shocks (%)', 
                                     display_ylim = c(-50, 50),
                                     xticklabs = data$all_years,
                                     main = paste0(region_name, ' (',ntrees, ' trees, ', length(stands_region), ' stands)'))



# -------------------------------------- #
# IMPACTS OF CLIM. PREDICTORS PER REGION #
# -------------------------------------- #

# Desert and highlands
extent <- aggregate(buffer(mtn[mtn$SECTION %in% c('SONORAN DESERT', 'MEXICAN HIGHLAND','SACRAMENTO')], width = 1))
region <- 'desert_highlands'
region_name <- 'Southern Desert and Highlands'
stands_region <- crop(stands, extent)
for(i in 1:data$N_all_years){
  gq_samples[[paste0('relchg_clim_perregion[',region,',',i,']')]] <- 0
}
stands_region <- which(data$uniq_stand_ids %in% stands_region$id)
ntrees <- sum(data$stand_idxs %in% stands_region)
first_obsyear <- min(sapply(which(data$stand_idxs %in% stands_region), function(t){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  return(min(data$years[idxs]))
}))
last_obsyear <- max(sapply(which(data$stand_idxs %in% stands_region), function(t){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  return(max(data$years[idxs]))
}))
year_idxs <- which(data$all_years %in% first_obsyear:last_obsyear)
for(s in stands_region){
  idxs <- seq(1 + data$N_all_years * (s-1), data$N_all_years * s, 1)
  for(i in 1:length(idxs)){
    gq_samples[[paste0('relchg_clim_perregion[',region,',',i,']')]] <- gq_samples[[paste0('relchg_clim_perregion[',region,',',i,']')]] +
      gq_samples[[paste0('avg_exp_m1_clim_stdlvl[',idxs[i],']')]]/length(stands_region)*100
  }
}
idxs <- seq(1, data$N_all_years, 1)
names <- paste0('relchg_clim_perregion[',region, ',', idxs,']')
# non-observed years will not appear
for(i in idxs[!(idxs %in% year_idxs)]){
  gq_samples[[paste0('relchg_clim_perregion[',region,',',i,']')]] <- gq_samples[[paste0('relchg_clim_perregion[',region,',',i,']')]]+1e6 
}
util$plot_disc_pushforward_quantiles(gq_samples, names, 
                                     ylab = 'Average relative change in\ngrowth due to climate (%)', 
                                     display_ylim = c(-50, 50),
                                     xticklabs = data$all_years,
                                     main = paste0(region_name, ' (',ntrees, ' trees, ', length(stands_region), ' stands)'))

# Colorado Plateau
extent <- aggregate(buffer(mtn[mtn$PROVINCE == 'COLORADO PLATEAUS'], width = 1))
region <- 'colorado_plateau'
region_name <- 'Colorado Plateau'
stands_region <- crop(stands, extent)
for(i in 1:data$N_all_years){
  gq_samples[[paste0('relchg_clim_perregion[',region,',',i,']')]] <- 0
}
stands_region <- which(data$uniq_stand_ids %in% stands_region$id)
ntrees <- sum(data$stand_idxs %in% stands_region)
first_obsyear <- min(sapply(which(data$stand_idxs %in% stands_region), function(t){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  return(min(data$years[idxs]))
}))
last_obsyear <- max(sapply(which(data$stand_idxs %in% stands_region), function(t){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  return(max(data$years[idxs]))
}))
year_idxs <- which(data$all_years %in% first_obsyear:last_obsyear)
for(s in stands_region){
  idxs <- seq(1 + data$N_all_years * (s-1), data$N_all_years * s, 1)
  for(i in 1:length(idxs)){
    gq_samples[[paste0('relchg_clim_perregion[',region,',',i,']')]] <- gq_samples[[paste0('relchg_clim_perregion[',region,',',i,']')]] +
      gq_samples[[paste0('avg_exp_m1_clim_stdlvl[',idxs[i],']')]]/length(stands_region)*100
  }
}
idxs <- seq(1, data$N_all_years, 1)
names <- paste0('relchg_clim_perregion[',region, ',', idxs,']')
# non-observed years will not appear
for(i in idxs[!(idxs %in% year_idxs)]){
  gq_samples[[paste0('relchg_clim_perregion[',region,',',i,']')]] <- gq_samples[[paste0('relchg_clim_perregion[',region,',',i,']')]]+1e6 
}
util$plot_disc_pushforward_quantiles(gq_samples, names, 
                                     ylab = 'Average relative change in\ngrowth due to climate (%)', 
                                     display_ylim = c(-50, 50),
                                     xticklabs = data$all_years,
                                     main = paste0(region_name, ' (',ntrees, ' trees, ', length(stands_region), ' stands)'))

# Southern Rockies
extent <- aggregate(buffer(mtn[mtn$PROVINCE == 'SOUTHERN ROCKY MOUNTAINS'], width = 1))
region <- 'southern_rockies'
region_name <- 'Southern Rocky Mountains'
stands_region <- crop(stands, extent)
for(i in 1:data$N_all_years){
  gq_samples[[paste0('relchg_clim_perregion[',region,',',i,']')]] <- 0
}
stands_region <- which(data$uniq_stand_ids %in% stands_region$id)
ntrees <- sum(data$stand_idxs %in% stands_region)
first_obsyear <- min(sapply(which(data$stand_idxs %in% stands_region), function(t){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  return(min(data$years[idxs]))
}))
last_obsyear <- max(sapply(which(data$stand_idxs %in% stands_region), function(t){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  return(max(data$years[idxs]))
}))
year_idxs <- which(data$all_years %in% first_obsyear:last_obsyear)
for(s in stands_region){
  idxs <- seq(1 + data$N_all_years * (s-1), data$N_all_years * s, 1)
  for(i in 1:length(idxs)){
    gq_samples[[paste0('relchg_clim_perregion[',region,',',i,']')]] <- gq_samples[[paste0('relchg_clim_perregion[',region,',',i,']')]] +
      gq_samples[[paste0('avg_exp_m1_clim_stdlvl[',idxs[i],']')]]/length(stands_region)*100
  }
}
idxs <- seq(1, data$N_all_years, 1)
names <- paste0('relchg_clim_perregion[',region, ',', idxs,']')
# non-observed years will not appear
for(i in idxs[!(idxs %in% year_idxs)]){
  gq_samples[[paste0('relchg_clim_perregion[',region,',',i,']')]] <- gq_samples[[paste0('relchg_clim_perregion[',region,',',i,']')]]+1e6 
}
util$plot_disc_pushforward_quantiles(gq_samples, names, 
                                     ylab = 'Average relative change in\ngrowth due to climate (%)', 
                                     display_ylim = c(-50, 50),
                                     xticklabs = data$all_years,
                                     main = paste0(region_name, ' (',ntrees, ' trees, ', length(stands_region), ' stands)'))

# Sierra Nevada
extent <- aggregate(buffer(mtn[mtn$SECTION == 'SIERRA NEVADA'], width = 1))
region <- 'sierra_nevada'
region_name <- 'Sierra Nevada'
stands_region <- crop(stands, extent)
for(i in 1:data$N_all_years){
  gq_samples[[paste0('relchg_clim_perregion[',region,',',i,']')]] <- 0
}
stands_region <- which(data$uniq_stand_ids %in% stands_region$id)
ntrees <- sum(data$stand_idxs %in% stands_region)
first_obsyear <- min(sapply(which(data$stand_idxs %in% stands_region), function(t){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  return(min(data$years[idxs]))
}))
last_obsyear <- max(sapply(which(data$stand_idxs %in% stands_region), function(t){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  return(max(data$years[idxs]))
}))
year_idxs <- which(data$all_years %in% first_obsyear:last_obsyear)
for(s in stands_region){
  idxs <- seq(1 + data$N_all_years * (s-1), data$N_all_years * s, 1)
  for(i in 1:length(idxs)){
    gq_samples[[paste0('relchg_clim_perregion[',region,',',i,']')]] <- gq_samples[[paste0('relchg_clim_perregion[',region,',',i,']')]] +
      gq_samples[[paste0('avg_exp_m1_clim_stdlvl[',idxs[i],']')]]/length(stands_region)*100
  }
}
idxs <- seq(1, data$N_all_years, 1)
names <- paste0('relchg_clim_perregion[',region, ',', idxs,']')
# non-observed years will not appear
for(i in idxs[!(idxs %in% year_idxs)]){
  gq_samples[[paste0('relchg_clim_perregion[',region,',',i,']')]] <- gq_samples[[paste0('relchg_clim_perregion[',region,',',i,']')]]+1e6 
}
util$plot_disc_pushforward_quantiles(gq_samples, names, 
                                     ylab = 'Average relative change in\ngrowth due to climate (%)', 
                                     display_ylim = c(-50, 50),
                                     xticklabs = data$all_years,
                                     main = paste0(region_name, ' (',ntrees, ' trees, ', length(stands_region), ' stands)'))
# Cascade range
extent <- aggregate(buffer(mtn[mtn$SECTION %in% c('NORTHERN CASCADE MOUNTAINS', 'MIDDLE CASCADE MOUNTAINS','SOUTHERN CASCADE MOUNTAINS')], width = 1))
region <- 'cascade_range'
region_name <- 'Cascade Range'
stands_region <- crop(stands, extent)
for(i in 1:data$N_all_years){
  gq_samples[[paste0('relchg_clim_perregion[',region,',',i,']')]] <- 0
}
stands_region <- which(data$uniq_stand_ids %in% stands_region$id)
ntrees <- sum(data$stand_idxs %in% stands_region)
first_obsyear <- min(sapply(which(data$stand_idxs %in% stands_region), function(t){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  return(min(data$years[idxs]))
}))
last_obsyear <- max(sapply(which(data$stand_idxs %in% stands_region), function(t){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  return(max(data$years[idxs]))
}))
year_idxs <- which(data$all_years %in% first_obsyear:last_obsyear)
for(s in stands_region){
  idxs <- seq(1 + data$N_all_years * (s-1), data$N_all_years * s, 1)
  for(i in 1:length(idxs)){
    gq_samples[[paste0('relchg_clim_perregion[',region,',',i,']')]] <- gq_samples[[paste0('relchg_clim_perregion[',region,',',i,']')]] +
      gq_samples[[paste0('avg_exp_m1_clim_stdlvl[',idxs[i],']')]]/length(stands_region)*100
  }
}
idxs <- seq(1, data$N_all_years, 1)
names <- paste0('relchg_clim_perregion[',region, ',', idxs,']')
# non-observed years will not appear
for(i in idxs[!(idxs %in% year_idxs)]){
  gq_samples[[paste0('relchg_clim_perregion[',region,',',i,']')]] <- gq_samples[[paste0('relchg_clim_perregion[',region,',',i,']')]]+1e6 
}
util$plot_disc_pushforward_quantiles(gq_samples, names, 
                                     ylab = 'Average relative change in\ngrowth due to climate (%)', 
                                     display_ylim = c(-50, 50),
                                     xticklabs = data$all_years,
                                     main = paste0(region_name, ' (',ntrees, ' trees, ', length(stands_region), ' stands)'))


