rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(cmdstanr)
library(terra)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('mcmc_custom_functions.R', local = util)
setwd(wd)

folder <- '/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model'

fit <- readRDS(file.path(wd, 'output/model/model21bis', 'fit_model21bis_HGSP_24species_365stands_withinit_4threads_climateshocks.rds'))
params <- c(
  "mu_alpha", "sigma_alpha", "alpha",
  "sigma_alpha_stand", "alpha_stand",
  "mu_beta_gdd", "sigma_beta_gdd", "beta_gdd",
  "mu_beta_pre", "sigma_beta_pre", "beta_pre",
  "mu_beta_vpd", "sigma_beta_vpd", "beta_vpd",
  "delta_clim",
  "tau_clim", # "kappa_clim_free", "kappa_clim",
  'log_kappa_clim', 'kappa_clim',
  "f_tilde_ind",
  "mu_log_rho", "sigma_log_rho", "rho_merged",
  "mu_log_gamma", "sigma_log_gamma", "log_gamma_merged", 'gamma_merged',
  "mu_log_tau_conc", "sigma_log_tau_conc", "tau_conc",
  "mu_phi", "sigma_phi", "logit_phi_sck", "phi_sck", "beta_phi_vpd", "beta_phi_pre",
  "mu_omega_conc", "sigma_omega_conc", "logit_omega_conc_sck", "omega_conc_sck",
  "mu_omega_shutdown", "sigma_omega_shutdown", "logit_omega_shutdown", "omega_shutdown",
  "thetas_idio", "tau_idio",
  "sigma"
)
base_samples <- fit$draws(params)
names <- dimnames(base_samples)$variable
base_samples <- lapply(1:dim(base_samples)[3],
                       function(k){t(matrix(base_samples[1:dim(base_samples)[1],1:dim(base_samples)[2],k],
                                            nrow = dim(base_samples)[1], ncol = dim(base_samples)[2]))})
names(base_samples) <- names

pre0 <- 5
vpd0 <- 23

data <- readRDS(file.path(folder,'data_10july2026_24species_365stands_19502022.rds'))
sites <- vect(data.frame(x = as.numeric(data$uniq_stand_lon), y = as.numeric(data$uniq_stand_lat)))
sites$id <- 1:nrow(sites)

states <- sapply(1:data$N_stands, function(s){
  unique(substr(data$uniq_tree_ids[data$stand_idxs == s], 1, 2))
})

northamerica <- geodata::gadm(country = geodata::country_codes("North America")$ISO3, level = 1, resolution = 2, path = file.path(wd, 'data/gadm'))

extent <- ext(c(-126, -102, 30,49))
us <- crop(geodata::gadm(country = "USA",level = 1, resolution = 1, path = file.path(wd, 'data')), extent)
canmex <- crop(geodata::gadm(country = c('MEX', 'CAN'),level =0, resolution = 1, path = file.path(wd, 'data')), extent)
mtn <- crop(vect(file.path(wd, 'data/usgs/physio_shp', 'physio.shp')), extent)

par(mfrow = c(1,1))
plot(us, ext = extent, clip = FALSE, lwd = 0.5, 
     box = FALSE, buffer = FALSE, axes = FALSE,
     mar = c(0,2,0.5,0), col = 'grey90', border = 'white')
plot(canmex, ext = extent, clip = FALSE, lwd = 0.5, 
     box = FALSE, buffer = FALSE, axes = FALSE,
     mar = c(0,2,0.5,0), col = 'grey90', border = 'white', add = T)
lines(aggregate(rbind(us, canmex)), col = 'white', lwd = 3)

colplat <- aggregate(buffer(mtn[mtn$PROVINCE %in% c('COLORADO PLATEAUS', 'SOUTHERN ROCKY MOUNTAINS')], width = 1))
plot(colplat, col = '#DB474360', border = 'white', add = TRUE, lwd = 2.5)
normnt <- aggregate(buffer(mtn[mtn$PROVINCE %in% c('NORTHERN ROCKY MOUNTAINS') |
                                 mtn$SECTION %in% c('NORTHERN CASCADE MOUNTAINS', 
                                                    'MIDDLE CASCADE MOUNTAINS','SOUTHERN CASCADE MOUNTAINS')], width = 1))
plot(normnt,  col = '#7C873E60', border = 'white', add = TRUE, lwd = 2.5)
midmnt <- aggregate(buffer(mtn[mtn$PROVINCE %in% c('MIDDLE ROCKY MOUNTAINS')], width = 1))
plot(midmnt,  col = '#8064A280', border = 'white', add = TRUE, lwd = 2.5)
sienev <- aggregate(buffer(mtn[mtn$SECTION == 'SIERRA NEVADA'], width = 1))
plot(sienev, col = "#5495CF60", border = 'white', add = TRUE, lwd = 2.5) 
deserts <- aggregate(buffer(mtn[mtn$SECTION %in% c('SONORAN DESERT', 'MEXICAN HIGHLAND','SACRAMENTO')], width = 1))
plot(deserts, col = '#F5AF4D60', border = 'white', add = TRUE, lwd = 2.5)
# lines(canmex, lwd = 0.5)
points(sites, pch = 20, cex = 1.4, col = 'white')
points(sites, pch = 20, cex = 0.7, col = 'grey20')

legend(
  "bottomleft", inset = c(-0.05, 0.1),
  legend = c('Northern Rockies & Cascades', 'Middle Rockies', 'Sierra Nevada', 'Southern Rockies & Colorado Plateau', 'Southwest Deserts & Highlands'),
  col = c('#7C873E80', "#8064A280", "#5495CF80", '#DB474380', '#F5AF4D80'),
  pch = 15, cex = 0.75, bty = 'n', xpd = NA, pt.cex = 1.5
)

# Colorado plateau and more
par(mfrow = c(1,1))
sites_here <- terra::intersect(sites, colplat)$id
base_samples[['mean_phi']] <- 0
for(s in sites_here){
  years <- seq(1+data$N_all_years*(s-1), data$N_all_years*s, 1)
  for(i in 1:length(years)){
    
    logit_phi <- base_samples[[paste0('logit_phi_sck[', s, ']')]]+
      (data$vpd_obs[years[i]]-vpd0)*base_samples[['beta_phi_vpd']]+
      (data$pre_obs[years[i]]-pre0)*base_samples[['beta_phi_pre']]
    phi <- boot::inv.logit(logit_phi)
    
    base_samples[['mean_phi']]  <- base_samples[['mean_phi']]  + phi
  }
}
base_samples[['mean_phi']]  <- base_samples[['mean_phi']] /(data$N_all_years*length(sites_here))
util$plot_expectand_pushforward(base_samples[['mean_phi']] , 30, flim = c(0, 0.3), col = '#DB4743',
                                border = '#DB474360',  ylim = c(0,100), display_name = bquote('Average'~phi))


# Sierra Nevada
sites_here <- terra::intersect(sites, sienev)$id
base_samples[['mean_phi']] <- 0
for(s in sites_here){
  years <- seq(1+data$N_all_years*(s-1), data$N_all_years*s, 1)
  for(i in 1:length(years)){
    
    logit_phi <- base_samples[[paste0('logit_phi_sck[', s, ']')]]+
      (data$vpd_obs[years[i]]-vpd0)*base_samples[['beta_phi_vpd']]+
      (data$pre_obs[years[i]]-pre0)*base_samples[['beta_phi_pre']]
    phi <- boot::inv.logit(logit_phi)
    
    base_samples[['mean_phi']]  <- base_samples[['mean_phi']]  + phi
  }
}
base_samples[['mean_phi']]  <- base_samples[['mean_phi']] /(data$N_all_years*length(sites_here))
util$plot_expectand_pushforward(base_samples[['mean_phi']] , 30, flim = c(0, 0.3), col = '#5495CF',
                                border = '#5495CF60', add = T)

# Middle Rockies
sites_here <- terra::intersect(sites, midmnt)$id
base_samples[['mean_phi']] <- 0
for(s in sites_here){
  years <- seq(1+data$N_all_years*(s-1), data$N_all_years*s, 1)
  for(i in 1:length(years)){
    
    logit_phi <- base_samples[[paste0('logit_phi_sck[', s, ']')]]+
      (data$vpd_obs[years[i]]-vpd0)*base_samples[['beta_phi_vpd']]+
      (data$pre_obs[years[i]]-pre0)*base_samples[['beta_phi_pre']]
    phi <- boot::inv.logit(logit_phi)
    
    base_samples[['mean_phi']]  <- base_samples[['mean_phi']]  + phi
  }
}
base_samples[['mean_phi']]  <- base_samples[['mean_phi']] /(data$N_all_years*length(sites_here))
util$plot_expectand_pushforward(base_samples[['mean_phi']] , 30, flim = c(0, 0.3), col = '#8064A2',
                                border = '#8064A260', add = T)

# Northern Mountains
sites_here <- terra::intersect(sites, normnt)$id
base_samples[['mean_phi']] <- 0
for(s in sites_here){
  years <- seq(1+data$N_all_years*(s-1), data$N_all_years*s, 1)
  for(i in 1:length(years)){
    
    logit_phi <- base_samples[[paste0('logit_phi_sck[', s, ']')]]+
      (data$vpd_obs[years[i]]-vpd0)*base_samples[['beta_phi_vpd']]+
      (data$pre_obs[years[i]]-pre0)*base_samples[['beta_phi_pre']]
    phi <- boot::inv.logit(logit_phi)
    
    base_samples[['mean_phi']]  <- base_samples[['mean_phi']]  + phi
  }
}
base_samples[['mean_phi']]  <- base_samples[['mean_phi']] /(data$N_all_years*length(sites_here))
util$plot_expectand_pushforward(base_samples[['mean_phi']], 30, flim = c(0, 0.3), col = '#7C873E',
                                border = '#7C873E60', add = T)

# Deserts
sites_here <- terra::intersect(sites, deserts)$id
base_samples[['mean_phi']] <- 0
for(s in sites_here){
  years <- seq(1+data$N_all_years*(s-1), data$N_all_years*s, 1)
  for(i in 1:length(years)){
    
    logit_phi <- base_samples[[paste0('logit_phi_sck[', s, ']')]]+
      (data$vpd_obs[years[i]]-vpd0)*base_samples[['beta_phi_vpd']]+
      (data$pre_obs[years[i]]-pre0)*base_samples[['beta_phi_pre']]
    phi <- boot::inv.logit(logit_phi)
    
    base_samples[['mean_phi']]  <- base_samples[['mean_phi']]  + phi
  }
}
base_samples[['mean_phi']]  <- base_samples[['mean_phi']] /(data$N_all_years*length(sites_here))
util$plot_expectand_pushforward(base_samples[['mean_phi']], 30, flim = c(0, 0.3), col = '#F5AF4D',
                                border = '#F5AF4D60', add = T)

legend(
  "bottomleft", inset = c(0.6, 0.7),
  legend = c('Northern Rockies & Cascades', 'Middle Rockies', 'Sierra Nevada', 'Southern Rockies & Colorado Plateau', 'Southwest Deserts & Highlands'),
  col = c('#7C873E80', "#8064A280", "#5495CF80", '#DB474380', '#F5AF4D80'),
  pch = 15, cex = 0.75, bty = 'n', xpd = NA, pt.cex = 1.5
)


sites_here <- terra::intersect(sites, normnt)$id
idxs <- c(sapply(sites_here, function(s) seq(1+data$N_all_years*(s-1), data$N_all_years*s, 1)))
hist(data$vpd_obs[idxs], breaks = seq(4,50, 1), col = '#7C873E60', border = 'white', 
     ylim = c(0,300), xlim = c(5,45), xlab = 'VPD (hPa)', main ='')
sites_here <- terra::intersect(sites, sienev)$id
idxs <- c(sapply(sites_here, function(s) seq(1+data$N_all_years*(s-1), data$N_all_years*s, 1)))
hist(data$vpd_obs[idxs], breaks = seq(4,50, 1), col = '#5495CF60', border = 'white', add = T)
sites_here <- terra::intersect(sites, colplat)$id
idxs <- c(sapply(sites_here, function(s) seq(1+data$N_all_years*(s-1), data$N_all_years*s, 1)))
hist(data$vpd_obs[idxs], breaks = seq(4,50, 1), col = '#DB474360', border = 'white', add = T) 

