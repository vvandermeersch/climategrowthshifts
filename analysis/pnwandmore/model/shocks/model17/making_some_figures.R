rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(cmdstanr)
library(rnaturalearth)
library(terra)
world <- vect(ne_countries(scale = "medium", returnclass = "sf"))
westam <- vect(ne_states(country = c("Canada", "Mexico", "United States of America"), returnclass = "sf"))

setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)


folder <- '/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model/'
datasets <- readRDS(file.path(folder,'datasets_26mar2026_TSME_standclimate_18962024.rds'))
data <- readRDS(file.path(folder,'data_26mar2026_TSME_standclimate_18962024.rds'))
fit <- readRDS(file.path(folder, 'model17/fit_HSGP_26mar2026_TSME_model17.rds'))

par(mfrow = c(1,1), cex.main = 0.9, mar = c(4,4,1,1))
plot(westam, xlim = c(-130, -104), ylim = c(30, 50),
     box = FALSE, clip = TRUE, col = 'grey97')
points(x = data$uniq_stand_lon, y = data$uniq_stand_lat,
       pch = 19, col = 'white', cex = 2)
points(x = data$uniq_stand_lon, y = data$uniq_stand_lat,
       pch = 19, col = util$c_light, cex = 1)



# Generate predictions!
mod_gq <- cmdstan_model(file.path(wd, 'model/stan/model17/hsgp', 'model17_HSGP_indGP_only1species_withGQ.stan'))
data_gq <- data
data_gq$grainsize <-  ceiling(data_gq$N_stands/4)
data_gq$uniq_tree_ids <- NULL
data_gq$uniq_species_ids <- NULL
data_gq$uniq_stand_ids <- NULL
data_gq$uniq_stand_lat <- NULL
data_gq$uniq_stand_lon <- NULL
fit_gq <- mod_gq$generate_quantities(fit, data = data_gq, seed = 5838293, parallel_chains = 4)
gc()
gq_samples <- fit_gq$draws(c('log_rw_pred', 'f', 'f_ind', 'delta_sck',
                             'num_conc_sck_stdlvl', 'num_idsc_sck_stdlvl', 'num_dbl_sck_stdlvl',
                             'num_sck_stdlvl', 'avg_exp_m1_delta_sck_stdlvl', 
                             'avg_exp_m1_clim_stdlvl'))
gc()
names <- dimnames(gq_samples)$variable
gq_samples <- lapply(1:dim(gq_samples)[3], 
                     function(k){t(matrix(gq_samples[1:dim(gq_samples)[1],1:dim(gq_samples)[2],k], 
                                          nrow = dim(gq_samples)[1], ncol = dim(gq_samples)[2]))})
names(gq_samples) <- names

par(mfrow = c(3,1), cex.main = 0.9, mar = c(4,6,1,1), cex.axis = 1, cex.lab = 1.3)
# Number of shocks across all stands
for(i in 1:data$N_all_years){
  gq_samples[[paste0('num_sck[',i,']')]]<- 0
}

for(s in 1:data$N_stands){
  idxs <- seq(1 + data$N_all_years * (s-1), data$N_all_years * s, 1)
  for(i in 1:length(idxs)){
    gq_samples[[paste0('num_sck[',i,']')]] <- gq_samples[[paste0('num_sck[',i,']')]] +
      gq_samples[[paste0('num_sck_stdlvl[',idxs[i],']')]]/data$N_stands
  }
}
for(i in 1:data$N_all_years){
  gq_samples[[paste0('num_sck[',i,']')]] <- gq_samples[[paste0('num_sck[',i,']')]]*100
}

names <- paste0('num_sck[',1:data$N_all_years,']')
util$plot_disc_pushforward_quantiles(gq_samples, names, 
                                     ylab = 'Average proportion of shocks\n(% of trees within a stand)', 
                                     display_ylim = c(0, 20),
                                     xticklabs = data$all_years)
abline(v = which(data$all_years %in% c(1899, 1916, 1956, 1991, 2011)), lty = 3)

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
                                     display_ylim = c(-20, 25),
                                     xticklabs = data$all_years)
abline(h = 0, col = 'grey70', lty = 2, lwd = 2)
abline(v = which(data$all_years %in% c(1899, 1916, 1956, 1991, 2011)), lty = 3)

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
                                     display_ylim = c(-20, 25),
                                     xticklabs = data$all_years)
abline(h = 0, col = 'grey70', lty = 2, lwd = 2)
abline(v = which(data$all_years %in% c(1899, 1916, 1956, 1991, 2011)), lty = 3)


for(i in which(data$all_years %in% c(1899, 1916, 1956, 1991, 2011))){
  gq_samples[[paste0('prop_expl_shock[',i,']')]] <- 0
  gq_samples[[paste0('prop_expl_shock[',i,']')]] <- 
    gq_samples[[paste0('avg_exp_m1_delta_sck[',i,']')]]/(gq_samples[[paste0('avg_exp_m1_clim[',i,']')]] + gq_samples[[paste0('avg_exp_m1_delta_sck[',i,']')]])*100
}

par(mfrow = c(1,1), cex.main = 0.9, mar = c(4,6,1,1), cex.axis = 0.9, cex.lab = 0.9)
names <- paste0('prop_expl_shock[',which(data$all_years %in% c(1899, 1916, 1956, 1991, 2011)),']')
util$plot_disc_pushforward_quantiles(gq_samples, names, 
                                     ylab = 'Proportion of change\nexplained by shocks (%)', 
                                     display_ylim = c(0, 100),
                                     xticklabs = c(1899, 1916, 1956, 1991, 2011))

