rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('functions/custom_functions.R', local=util)
setwd(wd)

data <- readRDS(file.path(wd, 'output/model', 'data_11mar2026_oregonTSME_standclimate_18962024.rds'))
fit <- readRDS(file.path(wd, 'output/model/model16', "fit_11mar2026_oregonTSME_18962024_model16_indGP.rds")) 

samples <- fit$draws(c('rho_ind', 'gamma_ind', 
                       'rho_sp', 'gamma_sp'))


names <- dimnames(samples)$variable
samples <- lapply(1:dim(samples)[3], 
                  function(k){t(matrix(samples[1:dim(samples)[1],1:dim(samples)[2],k], 
                                       nrow = dim(samples)[1], ncol = dim(samples)[2]))})
names(samples) <- names
util$check_all_expectand_diagnostics(samples)

util$plot_expectand_pushforward(samples[['rho_sp[1]']], 30)
util$plot_expectand_pushforward(samples[['rho_ind']], 30)
util$plot_expectand_pushforward(samples[['gamma_sp[1]']], 30)
util$plot_expectand_pushforward(samples[['gamma_ind']], 30)



samples <- fit$draws(
  c("phi_sck0", "tau_phi_sck", "alpha_phi_sck","phi_sck", 
    "omega_conc_sck0", "tau_omega_conc_sck", "alpha_omega_conc_sck", "omega_conc_sck",
    "mu_logdelta_omega_nonconc_sck", "tau_logdelta_omega_nonconc_sck", "logdelta_omega_nonconc_sck", "omega_nonconc_sck"))

names <- dimnames(samples)$variable
samples <- lapply(1:dim(samples)[3], 
                  function(k){t(matrix(samples[1:dim(samples)[1],1:dim(samples)[2],k], 
                                       nrow = dim(samples)[1], ncol = dim(samples)[2]))})
names(samples) <- names
util$check_all_expectand_diagnostics(samples)


names <- paste0('phi_sck[', 1:data$N_stands, ']')
util$plot_disc_pushforward_quantiles(samples, names)
names <- paste0('omega_nonconc_sck[', 1:data$N_stand_species, ']')
util$plot_disc_pushforward_quantiles(samples, names)
names <- paste0('omega_conc_sck[', 1:data$N_stand_species, ']')
util$plot_disc_pushforward_quantiles(samples, names)
