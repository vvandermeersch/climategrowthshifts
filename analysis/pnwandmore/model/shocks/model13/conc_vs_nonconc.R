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
fit_conconly <- readRDS(file.path(wd, 'output/model/model13', 'fit_30jan2025_gymnosperms_standclimate_19502024_7species_model13_updatedpriors_onlyconcordance.rds'))
fit_nonconc <- readRDS(file.path(wd, 'output/model/model13', 'fit_30jan2025_gymnosperms_standclimate_19502024_7species_model13_updatedpriors.rds'))

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
                                              "omega_conc_sck0", "mu_omega_conc_sck", "tau_omega_conc_sck", "alpha_omega_conc_sck",
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
                                                  "omega_conc_sck0", "mu_omega_conc_sck", "tau_omega_conc_sck", "alpha_omega_conc_sck",
                                                  "mu_logdelta_omega_nonconc_sck", "tau_logdelta_omega_nonconc_sck", "logdelta_omega_nonconc_sck", "omega_nonconc_sck",
                                                  
                                                  "sigma"))
names <- dimnames(param_samples)$variable
param_samples_fit_nonconc <- lapply(1:dim(param_samples)[3], function(k){t(matrix(param_samples[1:dim(param_samples)[1],1:dim(param_samples)[2],k], 
                                                                                   nrow = dim(param_samples)[1], ncol = dim(param_samples)[2]))})
names(param_samples_fit_nonconc) <- names




par(mfrow = c(2,1), mar = c(2.5,4,1,1), cex.main = 0.7)
names <- paste0("phi_sck[",1:data$N_stands,"]")
util$plot_disc_pushforward_quantiles(param_samples_fit_conconly, names, 
                                     xlab="Year", ylab=bquote(phi[shock]),
                                     display_ylim=c(-0.05, 0.85), main = 'Only concordant shocks')
util$plot_disc_pushforward_quantiles(param_samples_fit_nonconc, names, 
                                     xlab="Year", ylab=bquote(phi[shock]),
                                     display_ylim=c(-0.05, 0.85), main = 'Concordant & non-concordant shocks')

