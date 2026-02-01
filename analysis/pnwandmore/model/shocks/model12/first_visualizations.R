rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(terra)
library(rnaturalearth)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)


data <-readRDS(file.path(wd, 'output/model', 'data_29jan2025_gymnosperms_standclimate_19801990.rds'))
fit <- readRDS(file.path(wd, 'output/model', 'fit_29jan2025_gymnosperms_standclimate_19801990.rds'))

param_samples <- fit$draws(variables = c("mu_alpha","mu_gdd", "mu_vpd", "mu_pre", 
                                         "tau_alpha","tau_gdd", "tau_vpd", "tau_pre",
                                         "alpha","beta_gdd", "beta_vpd", "beta_pre",
                                         
                                         "mu_tau_sck", "tau_tau_sck", "tau_sck",
                                         "phi_sck", "omega_conc_sck",
                                         
                                         "mu_rho", "tau_rho", "rho_sp",
                                         "mu_gamma", "tau_gamma", "gamma_sp",
                                         
                                         "rho_sh",
                                         "mu_kappa", "tau_kappa", "kappa_sh",
                                         
                                         "sigma"))

# Transform as a Rstan object (to be able to use Mike's functions)
samples <- lapply(1:dim(param_samples)[3], function(k){t(matrix(param_samples[1:dim(param_samples)[1],1:dim(param_samples)[2],k], 
                                                                nrow = dim(param_samples)[1], ncol = dim(param_samples)[2]))})
names(samples) <- dimnames(param_samples)$variable
util$check_all_expectand_diagnostics(samples)

util$plot_pairs_by_chain(samples[['tau_sck[16]']], 'tau_sck[16]',
                         samples[['rho_sp[16]']], 'rho_sp[16]')

util$plot_pairs_by_chain(samples[['mu_tau_sck[1]']], 'mu_tau_sck[1]',
                         log(samples[['tau_tau_sck[1]']]), 'log(tau_tau_sck[1])')

subsamples <- util$filter_expectands(samples, c('tau_sck[1]', 'tau_sck[2]'))
util$check_all_expectand_diagnostics(subsamples)
util$plot_pairs_by_chain(samples[['tau_sck[1]']], 'tau_sck[1]',
                         samples[['tau_sck[2]']], 'tau_sck[2]')

subsamples <- util$filter_expectands(samples, c('mu_tau_sck[1]', 'tau_tau_sck[1]'))
util$check_all_expectand_diagnostics(subsamples)
util$plot_pairs_by_chain(samples[['mu_tau_sck[1]']], 'mu_tau_sck[1]',
                         samples[['tau_tau_sck[1]']], 'tau_tau_sck[1]')

subsamples <- util$filter_expectands(samples, c('rho_sp[1]', 'rho_sp[2]'))
util$check_all_expectand_diagnostics(subsamples)
util$plot_pairs_by_chain(samples[['rho_sp[1]']], 'rho_sp[1]',
                         samples[['rho_sp[2]']], 'rho_sp[2]')


subsamples <- util$filter_expectands(samples, c('rho_sh'))
util$check_all_expectand_diagnostics(subsamples)
