rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(cmdstanr)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('mcmc_custom_functions.R', local = util)
setwd(wd)

folder <- '/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model'

data <- readRDS(file.path(folder,'data_6july2026_21species_115stands_19502022.rds'))

states <- sapply(1:data$N_stands, function(s){
  unique(substr(data$uniq_tree_ids[data$stand_idxs == s], 1, 2))
})

states_stsp <- sapply(1:data$N_stand_species, function(s){
  unique(substr(data$uniq_tree_ids[data$stand_species_idxs == s], 1, 2))
})

fit <- readRDS(file.path(wd, 'output/model/model21bis', 'fit_model21bis_HGSP_21species_115stands_4threads_climateshocks_updated.rds'))
fit$time()
fit$diagnostic_summary()

params <- c(
  "mu_alpha", 
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
  "mu_omega_conc", "sigma_omega_species", "sigma_omega_stand", "alpha_omega_species", "alpha_omega_stand",
  "mu_omega_shutdown", "sigma_omega_shutdown", "logit_omega_shutdown", "omega_shutdown",
  "thetas_idio", "tau_idio",
  # 'thetas_baseline', 'omega_thetas', 'thetas_idio',
  "sigma"
)


base_samples <- fit$draws(params)
names <- dimnames(base_samples)$variable
base_samples <- lapply(1:dim(base_samples)[3],
                       function(k){t(matrix(base_samples[1:dim(base_samples)[1],1:dim(base_samples)[2],k],
                                            nrow = dim(base_samples)[1], ncol = dim(base_samples)[2]))})
names(base_samples) <- names

# base_samples <- readRDS(file.path(wd, 'output/model/model21bis', 'basesamples_model21bis_HGSP_11species_99stands_withinit.rds'))

params <- c(
  "mu_alpha", 
  "sigma_alpha_stand", "alpha_stand",
  "mu_beta_gdd", "sigma_beta_gdd", "beta_gdd",
  "mu_beta_pre", "sigma_beta_pre", "beta_pre",
  "mu_beta_vpd", "sigma_beta_vpd", "beta_vpd",
  "tau_clim", # "kappa_clim_free", "kappa_clim",
  'log_kappa_clim',
  "mu_log_rho", "sigma_log_rho", "rho_merged",
  "mu_log_gamma", "sigma_log_gamma", "log_gamma_merged",
  "mu_log_tau_conc", "sigma_log_tau_conc", "tau_conc",
  "mu_phi", "sigma_phi", "logit_phi_sck", "beta_phi_vpd", "beta_phi_pre",
  "mu_omega_conc", "sigma_omega_species", "sigma_omega_stand", "alpha_omega_species", "alpha_omega_stand", 
  "mu_omega_shutdown", "sigma_omega_shutdown", "logit_omega_shutdown", 
  "thetas_idio", "tau_idio",
  "sigma"
)

for(p in params){
  message(paste0('\n\n---------\nParameter(s):',p))
  rest <- util$filter_expectands(base_samples, p, 
                                 check_arrays = TRUE)
  util$check_all_expectand_diagnostics(rest, min_ess_hat_per_chain=50)
}
