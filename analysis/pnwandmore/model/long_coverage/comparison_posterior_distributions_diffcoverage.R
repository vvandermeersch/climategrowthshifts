rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('functions/custom_functions.R', local=util)
setwd(wd)

data <- readRDS(file.path(wd, 'output/model/climatena', 'data_24sept2025_long_pinus.rds'))
samples <- readRDS(file.path(wd, "output/model/climatena", "samples_24sept2025_differentunits_partialpooling_byclade_centered_extended_constrainedrhos_updatedpriors_long_pinus.rds"))

short_samples <- readRDS(file.path(wd, 'output/model/climatena', 'samples_15sept2025_2predictors_partialpooling_2clades_centered_extended_short.rds'))
short_data <- readRDS(file.path(wd, 'output/model/climatena', 'data_11august2025_short.rds'))
short_species <- which(short_data$uniq_species_ids %in% data$uniq_species_ids)
short_rhos <- sapply(short_species, function(sp){
  util$ensemble_mcmc_quantile_est(short_samples[[paste0('rho_sp[', sp, ']')]], 0.5)
})
long_rhos <- sapply(1:data$N_species, function(sp){
  util$ensemble_mcmc_quantile_est(samples[[paste0('rho_sp[', sp, ']')]], 0.5)
})
short_betas_gdd <- sapply(short_species, function(sp){
  util$ensemble_mcmc_quantile_est(short_samples[[paste0('beta_gdd[', sp, ']')]], 0.5)
})
long_betas_gdd <- sapply(1:data$N_species, function(sp){
  util$ensemble_mcmc_quantile_est(samples[[paste0('beta_gdd[', sp, ']')]], 0.5)
})
short_betas_pre <- sapply(short_species, function(sp){
  util$ensemble_mcmc_quantile_est(short_samples[[paste0('beta_winterprec[', sp, ']')]], 0.5)
})
long_betas_pre <- sapply(1:data$N_species, function(sp){
  util$ensemble_mcmc_quantile_est(samples[[paste0('beta_winterprec[', sp, ']')]], 0.5)
})

# Rhos
par(mfrow = c(2,1), mar = c(3,4.5,2,2), cex.lab = 1.2)
names <- sapply(1:data$N_species,
                function(sp) paste0('rho_sp[', sp, ']'))
util$plot_disc_pushforward_quantiles(samples, names, baseline_values = short_rhos, main = 'Rhos, long-coverage model (black = short)')
names <- sapply(short_species,
                function(sp) paste0('rho_sp[', sp, ']'))
util$plot_disc_pushforward_quantiles(short_samples, names, baseline_values = long_rhos, main = 'Rhos, hort-coverage model (black = long)')

# GDD
par(mfrow = c(2,1), mar = c(3,4.5,2,2), cex.lab = 1.2)
names <- sapply(1:data$N_species,
                function(sp) paste0('beta_gdd[', sp, ']'))
util$plot_disc_pushforward_quantiles(samples, names, baseline_values = short_betas_gdd, main = 'GDD slopes, long-coverage model (black = short)')
names <- sapply(short_species,
                function(sp) paste0('beta_gdd[', sp, ']'))
util$plot_disc_pushforward_quantiles(short_samples, names, baseline_values = long_betas_gdd, main = 'GDD, slopes short-coverage model (black = long)')

# Winter prec.
par(mfrow = c(2,1), mar = c(3,4.5,2,2), cex.lab = 1.2)
names <- sapply(1:data$N_species,
                function(sp) paste0('beta_winterprec[', sp, ']'))
util$plot_disc_pushforward_quantiles(samples, names, baseline_values = short_betas_pre, main = 'Winter prec. slopes, long-coverage model (black = short)')
names <- sapply(short_species,
                function(sp) paste0('beta_winterprec[', sp, ']'))
util$plot_disc_pushforward_quantiles(short_samples, names, baseline_values = long_betas_pre, main = 'Winter prec. slopes short-coverage model (black = long)')

#
par(mfrow = c(2,1), mar = c(3,4.5,2,2), cex.lab = 1.2)
stand <- 3
names <- sapply(1:data$N_all_years,
                function(sp) paste0('f_tilde_sh[', stand, ',', sp, ']'))
util$plot_realizations(samples, names, plot_xs = data$all_years, N_plots = 10,
                       main = "Long-coverage model")
# abline(v=2002, lwd=1, lty=2, col="grey50")

short_stand <- unique(short_data$stand_idx[which(short_data$uniq_tree_ids %in% data$uniq_tree_ids[data$stand_idxs == stand])])
names <- sapply(1:short_data$N_all_years,
                function(sp) paste0('f_tilde_sh[', short_stand, ',', sp, ']'))
util$plot_realizations(short_samples, names, plot_xs = short_data$all_years, N_plots = 10,
                       display_xlim = range(data$all_years), main = "Short-coverage model")
# abline(v=2002, lwd=1, lty=2, col="grey50")


data$uniq_tree_ids[data$stand_idxs == stand]
data$species_idxs[data$stand_idxs == stand]
short_data$uniq_tree_ids[short_data$stand_idxs == short_stand]
short_data$species_idxs[short_data$stand_idxs == short_stand]
