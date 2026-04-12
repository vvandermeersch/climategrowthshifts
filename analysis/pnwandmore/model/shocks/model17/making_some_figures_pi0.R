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


plot(mtn)
text(mtn, labels = mtn$PROVINCE)

plot(aggregate(buffer(mtn[mtn$PROVINCE == ''], width = 1)))
plot(aggregate(buffer(mtn[mtn$PROVINCE == ''], width = 1)))
plot(aggregate(buffer(mtn[mtn$PROVINCE == 'CASCADE-SIERRA MOUNTAINS'], width = 1)))

# ACROSS ALL STANDS
par(mfrow = c(3,1), cex.main = 0.9, mar = c(4,6,1,1), cex.axis = 1, cex.lab = 1.3)
# Number of shocks across all stands
for(i in 1:data$N_all_years){
  gq_samples[[paste0('num_sck[',i,']')]]<- 0
}

for(s in 1:data$N_stands){
  idxs <- seq(1 + data$N_all_years * (s-1), data$N_all_years * s, 1)
  for(i in 1:length(idxs)){
    gq_samples[[paste0('num_sck[',i,']')]] <- gq_samples[[paste0('num_sck[',i,']')]] +
      gq_samples[[paste0('num_conc_sck_stdlvl[',idxs[i],']')]]/data$N_stands
  }
}
for(i in 1:data$N_all_years){
  gq_samples[[paste0('num_sck[',i,']')]] <- gq_samples[[paste0('num_sck[',i,']')]]*100
}

names <- paste0('num_sck[',1:data$N_all_years,']')
util$plot_disc_pushforward_quantiles(gq_samples, names, 
                                     ylab = 'Average proportion of shocks\n(% of trees within a stand)', 
                                     display_ylim = c(0, 25),
                                     xticklabs = data$all_years)
# abline(h = 8, lty = 3)
selected_years <- c(1925, 1954, 1956, 1963, 2002)
abline(v = which(data$all_years %in% selected_years), lty = 3, col = 'grey80')

# Relative reduction of growth due to shocks
for(i in 1:data$N_all_years){
  gq_samples[[paste0('avg_exp_m1_delta_sck[',i,']')]]<- 0
}

for(s in 1:data$N_stands){
  idxs <- seq(1 + data$N_all_years * (s-1), data$N_all_years * s, 1)
  for(i in 1:length(idxs)){
    gq_samples[[paste0('avg_exp_m1_delta_sck[',i,']')]] <- gq_samples[[paste0('avg_exp_m1_delta_sck[',i,']')]] +
      gq_samples[[paste0('avg_exp_m1_delta_sck_stdlvl[',idxs[i],']')]]/data$N_stands
  }
}

for(i in 1:data$N_all_years){
  gq_samples[[paste0('avg_exp_m1_delta_sck[',i,']')]] <- gq_samples[[paste0('avg_exp_m1_delta_sck[',i,']')]]*100
}

names <- paste0('avg_exp_m1_delta_sck[',1:data$N_all_years,']')
util$plot_disc_pushforward_quantiles(gq_samples, names, 
                                     ylab = 'Average relative change in growth\ndue to shocks (%)', 
                                     display_ylim = c(-30, 25),
                                     xticklabs = data$all_years)
abline(h = 0, lty = 3)
abline(v = which(data$all_years %in% selected_years), lty = 3, col = 'grey80')

# Relative reduction of growth due to climate predictors
for(i in 1:data$N_all_years){
  gq_samples[[paste0('avg_exp_m1_clim[',i,']')]]<- 0
}

for(s in 1:data$N_stands){
  idxs <- seq(1 + data$N_all_years * (s-1), data$N_all_years * s, 1)
  for(i in 1:length(idxs)){
    gq_samples[[paste0('avg_exp_m1_clim[',i,']')]] <- gq_samples[[paste0('avg_exp_m1_clim[',i,']')]] +
      gq_samples[[paste0('avg_exp_m1_clim_stdlvl[',idxs[i],']')]]/data$N_stands
  }
}

for(i in 1:data$N_all_years){
  gq_samples[[paste0('avg_exp_m1_clim[',i,']')]] <- gq_samples[[paste0('avg_exp_m1_clim[',i,']')]]*100
}

names <- paste0('avg_exp_m1_clim[',1:data$N_all_years,']')
util$plot_disc_pushforward_quantiles(gq_samples, names, 
                                     ylab = 'Average relative change in growth\ndue to climate predictors (%)', 
                                     display_ylim = c(-30, 20),
                                     xticklabs = data$all_years)
abline(h = 0, lty = 3)
abline(v = which(data$all_years %in% selected_years), lty = 3, col = 'grey80')


for(i in which(data$all_years %in% selected_years)){
  gq_samples[[paste0('prop_expl_shock[',i,']')]] <- 0
  gq_samples[[paste0('prop_expl_shock[',i,']')]] <- 
    gq_samples[[paste0('avg_exp_m1_delta_sck[',i,']')]]/(gq_samples[[paste0('avg_exp_m1_clim[',i,']')]] + gq_samples[[paste0('avg_exp_m1_delta_sck[',i,']')]])*100
}

par(mfrow = c(1,1), cex.main = 0.9, mar = c(4,6,1,1), cex.axis = 0.9, cex.lab = 0.9)
names <- paste0('prop_expl_shock[',which(data$all_years %in% selected_years),']')
util$plot_disc_pushforward_quantiles(gq_samples, names, 
                                     ylab = 'Proportion of change\nexplained by shocks (%)', 
                                     display_ylim = c(0, 300),
                                     xticklabs = selected_years)


for(s in 1:data$N_stands){
  gq_samples[[paste0('num_sck_perstand[',s,']')]]<- 0
}

for(s in 1:data$N_stands){
  idxs <- seq(1 + data$N_all_years * (s-1), data$N_all_years * s, 1)
  for(i in 1:length(idxs)){
    gq_samples[[paste0('num_sck_perstand[',s,']')]] <- gq_samples[[paste0('num_sck_perstand[',s,']')]] +
      gq_samples[[paste0('num_conc_sck_stdlvl[',idxs[i],']')]]/data$N_stand_years[s]
  }
}
for(s in 1:data$N_stands){
  gq_samples[[paste0('num_sck_perstand[',s,']')]] <- gq_samples[[paste0('num_sck_perstand[',s,']')]]*100
}

names <- paste0('num_sck_perstand[',1:data$N_stands,']')
util$plot_disc_pushforward_quantiles(gq_samples, names, 
                                     ylab = 'Average proportion of shocks\n(% of years)', 
                                     display_ylim = c(0, 30),
                                     xticklabs = 1:data$N_stands)



idxs <- which(rep(data$all_years, data$N_stands)==1963)
names <- paste0('num_conc_sck_stdlvl[',idxs,']')
util$plot_disc_pushforward_quantiles(gq_samples, names, 
                                     ylab = 'Average proportion of shocks\n(% of years)', 
                                     display_ylim = c(0, 1),
                                     xticklabs = 1:data$N_stands)


states <- unique(datasets[, c('state', 'grouped_stand')])
for(s in unique(states$state)){
  for(i in 1:data$N_all_years){
    gq_samples[[paste0('num_sck_perstate[',s,',',i,']')]] <- 0
  }
}

for(state in unique(states$state)){
  stands <- which(data$uniq_stand_ids %in% states[states$state == state, 'grouped_stand'])
  for(s in stands){
    idxs <- seq(1 + data$N_all_years * (s-1), data$N_all_years * s, 1)
    for(i in 1:length(idxs)){
      gq_samples[[paste0('num_sck_perstate[',state,',',i,']')]] <- gq_samples[[paste0('num_sck_perstate[',state,',',i,']')]] +
        gq_samples[[paste0('num_conc_sck_stdlvl[',idxs[i],']')]]/length(stands)
    }
  }
}

for(s in unique(states$state)){
  for(i in 1:data$N_all_years){
    gq_samples[[paste0('num_sck_perstate[',s,',',i,']')]] <- gq_samples[[paste0('num_sck_perstate[',s,',',i,']')]]*100
  }
}

par(mfrow = c(3,2), cex.main = 0.9, mar = c(4,6,1,1), cex.axis = 0.8, cex.lab = 1)
selected_years <- data$all_years
for(y in selected_years){
  idx <- which(data$all_years==y)
  names <- paste0('num_sck_perstate[',unique(states$state), ',', idx,']')
  util$plot_disc_pushforward_quantiles(gq_samples, names, 
                                       ylab = 'Average proportion of shocks\n(% of trees within a stand)', 
                                       display_ylim = c(0, 100),
                                       xticklabs = unique(states$state),
                                       main = y)
}









