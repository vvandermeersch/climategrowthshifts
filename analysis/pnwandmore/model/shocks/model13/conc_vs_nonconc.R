rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(terra)
library(rnaturalearth)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)


data <-readRDS(file.path(wd, 'output/model', 'data_30jan2025_gymnosperms_standclimate_19502024_7species.rds'))
fit_conconly <- readRDS(file.path(wd, 'output/model/model13', 'fit_16feb2026_gymnosperms_standclimate_19502024_7species_model13_updatedpriors_onlyconcordance.rds'))
fit_nonconc <- readRDS(file.path(wd, 'output/model/model13', 'fit_16feb2026_gymnosperms_standclimate_19502024_7species_model13_updatedpriors.rds'))

# Parameter samples
param_samples <- fit_conconly$draws(variables = c("mu_alpha","mu_gdd", "mu_vpd", "mu_pre", 
                                              "tau_alpha","tau_gdd", "tau_vpd", "tau_pre",
                                              "alpha","beta_gdd", "beta_vpd", "beta_pre",
                                              
                                              "mu_tau_sck", "tau_tau_sck", "tau_sck",
                                              "phi_sck", "omega_conc_sck",
                                              
                                              "mu_rho", "tau_rho", "rho_sp",
                                              "mu_gamma", "tau_gamma", "gamma_sp",
                                              
                                              "rho_sh",
                                              "mu_kappa", "tau_kappa", "kappa_sh",
                                              
                                              "phi_sck0", "tau_phi_sck", "alpha_phi_sck",
                                              "omega_conc_sck0", "mu_omega_conc_sck", "tau_omega_conc_sck", "alpha_omega_conc_sck", "omega_conc_sck",
                                              # "mu_logdelta_omega_nonconc_sck", "tau_logdelta_omega_nonconc_sck", "logdelta_omega_nonconc_sck", "omega_nonconc_sck",
                                              
                                              "sigma"))
names <- dimnames(param_samples)$variable
param_samples_fit_conconly <- lapply(1:dim(param_samples)[3], function(k){t(matrix(param_samples[1:dim(param_samples)[1],1:dim(param_samples)[2],k], 
                                                                      nrow = dim(param_samples)[1], ncol = dim(param_samples)[2]))})
names(param_samples_fit_conconly) <- names


param_samples <- fit_nonconc$draws(variables = c("mu_alpha","mu_gdd", "mu_vpd", "mu_pre", 
                                                  "tau_alpha","tau_gdd", "tau_vpd", "tau_pre",
                                                  "alpha","beta_gdd", "beta_vpd", "beta_pre",
                                                  
                                                  "mu_tau_sck", "tau_tau_sck", "tau_sck",
                                                  "phi_sck", "omega_conc_sck",
                                                  
                                                  "mu_rho", "tau_rho", "rho_sp",
                                                  "mu_gamma", "tau_gamma", "gamma_sp",
                                                  
                                                  "rho_sh",
                                                  "mu_kappa", "tau_kappa", "kappa_sh",
                                                  
                                                  "phi_sck0", "tau_phi_sck", "alpha_phi_sck",
                                                  "omega_conc_sck0", "mu_omega_conc_sck", "tau_omega_conc_sck", "alpha_omega_conc_sck", "omega_conc_sck",
                                                  "mu_logdelta_omega_nonconc_sck", "tau_logdelta_omega_nonconc_sck", "logdelta_omega_nonconc_sck", "omega_nonconc_sck",
                                                  
                                                  "sigma"))
names <- dimnames(param_samples)$variable
param_samples_fit_nonconc <- lapply(1:dim(param_samples)[3], function(k){t(matrix(param_samples[1:dim(param_samples)[1],1:dim(param_samples)[2],k], 
                                                                                   nrow = dim(param_samples)[1], ncol = dim(param_samples)[2]))})
names(param_samples_fit_nonconc) <- names




par(mfrow = c(3,2), mar = c(1.5,5,1,1), cex.main = 1, cex.lab = 1.2)
names <- paste0("phi_sck[",1:data$N_stands,"]")
util$plot_disc_pushforward_quantiles(param_samples_fit_conconly, names, 
                                     xlab="Year", ylab=bquote(phi[shock ~ (stand)]^concordant),
                                     display_ylim=c(-0.05, 1.05), main = 'Only concordant shocks')
util$plot_disc_pushforward_quantiles(param_samples_fit_nonconc, names, 
                                     xlab="Year", ylab=bquote(phi[shock ~ (stand)]^concordant),
                                     display_ylim=c(-0.05, 1.05), main = 'Concordant & non-concordant shocks')

names <- paste0("omega_conc_sck[",1:data$N_stand_species,"]")
util$plot_disc_pushforward_quantiles(param_samples_fit_conconly, names, 
                                     xlab="Year", ylab=bquote(omega[shock ~ (tree)]^concordant),
                                     display_ylim=c(-0.05, 1.05), main = '')
util$plot_disc_pushforward_quantiles(param_samples_fit_nonconc, names, 
                                     xlab="Year", ylab=bquote(omega[shock ~ (tree)]^concordant),
                                     display_ylim=c(-0.05, 1.05), main = '')

names <- paste0("omega_nonconc_sck[",1:data$N_stand_species,"]")
plot.new()
util$plot_disc_pushforward_quantiles(param_samples_fit_nonconc, names, 
                                     xlab="Year", ylab=bquote(omega[shock ~ (tree)]^{non-concordant}),
                                     display_ylim=c(-0.05, 1.05), main = '')

par(mfrow = c(1,1), mar = c(5,5,1,1))
plot(x = NULL, y = NULL,
     xlim = c(0,0.2), ylim = c(0,0.2),
     xlab = bquote(phi[shock ~ (stand)]^concordant *' \u00D7 '* omega[shock ~ (stand)]^concordant),
     ylab = bquote((1-phi[shock ~ (stand)]^concordant) *' \u00D7 '* omega[shock ~ (stand)]^{non-concordant}),
     asp = 1)
saveme <- c()
for(i in 1:data$N_stand_species){
  
  sp <- unique(data$species_idxs[data$stand_species_idxs == i])
  phi <- paste0("phi_sck[",sp,"]")
  omega_conc <- paste0("omega_conc_sck[",i,"]")
  omega_nonconc <- paste0("omega_nonconc_sck[",i,"]")
  
  param_samples_fit_nonconc[['phi_omega_nonconc']] <- 
    (1 - param_samples_fit_nonconc[[phi]])*param_samples_fit_nonconc[[omega_nonconc]]
  
  param_samples_fit_nonconc[['phi_omega_conc']] <- 
    (param_samples_fit_nonconc[[phi]])*param_samples_fit_nonconc[[omega_conc]]
  
  conc <- util$ensemble_mcmc_quantile_est(param_samples_fit_nonconc[['phi_omega_conc']], c(0.05,0.5, 0.95))
  nonconc <- util$ensemble_mcmc_quantile_est(param_samples_fit_nonconc[['phi_omega_nonconc']], c(0.05,0.5, 0.95))
  
  if(nonconc['50%'] > 0.085){saveme <- c(saveme, i)}
  
  segments(x0 = conc['50%'], y0 = nonconc['5%'], y1 = nonconc['95%'])
  segments(x0 = conc['5%'], x1 = conc['95%'], y0 = nonconc['50%'])
  points(x = conc['50%'], y = nonconc['50%'], pch = 20)
}
abline(a = 0, b = 1, lty = 2, col = 'grey80')




library(cmdstanr)

# Generate quantities
mod_gq <- cmdstan_model(file.path(wd, 'model/stan', 'model13_updatedpriors_wGQ.stan'))
data_gq <- data
data_gq$grainsize <-  ceiling(data_gq$N_stands/8)
data_gq$uniq_tree_ids <- NULL
data_gq$uniq_species_ids <- NULL
data_gq$uniq_stand_ids <- NULL
data_gq$N_clades <- 1
data_gq$uniq_stand_lat <- NULL
data_gq$uniq_stand_lon <- NULL
fit_gq <- mod_gq$generate_quantities(fit_nonconc, data = data_gq, seed = 5838293, parallel_chains = 4)
gq_samples <- fit_gq$draws()
names <- dimnames(gq_samples)$variable
gq_samples <- lapply(1:dim(gq_samples)[3], function(k){t(matrix(gq_samples[1:dim(gq_samples)[1],1:dim(gq_samples)[2],k], 
                                                                nrow = dim(gq_samples)[1], ncol = dim(gq_samples)[2]))})
names(gq_samples) <- names


# Look at the substand 81
tree_idxs <- which(data$stand_species_idxs == 46)
par(mfrow = c(9,3), mar = c(2,3,0,1), mgp = c(2,0.5,0),
    cex.lab = 0.8, cex.axis = 0.7)
for(t in tree_idxs){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  stand <- data$stand_idxs[t]
  names <- paste0("log_rw_pred[",idxs,"]")

  util$plot_conn_pushforward_quantiles(gq_samples, names, data$years[idxs]-1949,
                                       xlab="Year", ylab="Log ring width (per mm)",
                                       display_ylim=c(-8, 2), display_xlim = range(data$years[idxs])-1949)
  points(data$years[idxs]-1949, log(data$rw_obs[idxs]), pch=16, cex=1, col="white")
  points(data$years[idxs]-1949, log(data$rw_obs[idxs]), pch=16, cex=0.5, col="black")
  text(x = .5, y = -6.5, labels = data$uniq_tree_ids[t], cex = 0.8, adj = 0)
  text(x = .5, y = -7.7, labels = paste0('Species:', data$uniq_species_ids[data$species_idxs[t]]), cex = 0.8, adj = 0)

  # start <- 1+(stand-1)*data$N_all_years
  # end <- stand*data$N_all_years
  # lines(1:data$N_all_years, data$pre_obs[start:end]/10, col = 'darkblue')
  # lines(1:data$N_all_years, data$vpd_obs[start:end]/20, col = 'darkred')
  
  names <- paste0("sck_state[",idxs,"]")
  util$plot_disc_pushforward_quantiles(gq_samples, names, data$years[idxs],
                                       xlab="Year", ylab="Shock state",
                                       display_ylim=c(-0.05, 1.05))
  
  names <- paste0("delta_sck[",idxs,"]")
  util$plot_conn_pushforward_quantiles(gq_samples, names, data$years[idxs]-1949,
                                       xlab="Year", ylab="Shock amplitude",
                                       display_ylim=c(-5, 1), display_xlim = range(data$years[idxs])-1949)
}

