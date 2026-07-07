rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(cmdstanr)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)

folder <- '/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model'

data <- readRDS(file.path(folder,'data_25june2026_11species_99stands_19502020.rds'))

fit <- readRDS(file.path(wd, 'output/model/model21bis', 'fit_model21bis_HGSP_11species_99stands_withinit.rds'))

params <- c(
  "mu_alpha", "sigma_alpha", "alpha",
  "mu_alpha_stand", "sigma_alpha_stand", 
  # "alpha_stand",
  "mu_beta_gdd", "sigma_beta_gdd", "beta_gdd",
  "mu_beta_pre", "sigma_beta_pre", "beta_pre",
  "mu_beta_vpd", "sigma_beta_vpd", "beta_vpd",
  # "delta_clim", 
  "tau_clim", "kappa_clim_free",
  # "f_tilde_ind",
  "mu_log_rho", "sigma_log_rho", "rho_merged",
  "mu_log_gamma", "sigma_log_gamma", "log_gamma_merged", 'gamma_merged',
  "tau_conc",
  "mu_phi", "sigma_phi", "logit_phi_sck",
  "mu_omega_conc", "sigma_omega_conc", "logit_omega_conc_sck",
  "mu_omega_shutdown", "sigma_omega_shutdown", "logit_omega_shutdown",
  "thetas_idio", "tau_idio",
  "sigma"
)

base_samples <- fit$draws(params)
names <- dimnames(base_samples)$variable
base_samples <- lapply(1:dim(base_samples)[3],
                  function(k){t(matrix(base_samples[1:dim(base_samples)[1],1:dim(base_samples)[2],k],
                                       nrow = dim(base_samples)[1], ncol = dim(base_samples)[2]))})
names(base_samples) <- names


for(p in params){
  message(paste0('\n\n---------\nParameter(s):',p))
  rest <- util$filter_expectands(base_samples, p, 
                                 check_arrays = TRUE)
  util$check_all_expectand_diagnostics(rest)
}

util$plot_pairs_by_chain(base_samples[['mu_log_rho']], 'mu_log_rho',
                         base_samples[['sigma_log_rho']], 'sigma_log_rho')

util$plot_pairs_by_chain(base_samples[['mu_log_gamma']], 'mu_log_gamma',
                         base_samples[['sigma_log_gamma']], 'sigma_log_gamma')

util$plot_pairs_by_chain(base_samples[['tau_clim']], 'tau_clim',
                         base_samples[['tau_idio']], 'tau_idio')

for(sp in 1:data$N_species){
  util$plot_pairs_by_chain(base_samples[[paste0('rho_merged[',sp,']')]], paste0('rho_merged[',sp,']'),
                           base_samples[[paste0('gamma_merged[',sp,']')]], paste0('gamma_merged[',sp,']'))
}

for(sp in 1:data$N_species){
  util$plot_pairs_by_chain(base_samples[[paste0('alpha[',sp,']')]], paste0('alpha[',sp,']'),
                           base_samples[[paste0('beta_vpd[',sp,']')]], paste0('beta_vpd[',sp,']'))
}

par(mfrow = c(2,2), mar = c(4,4,1,1))
util$plot_disc_pushforward_quantiles(base_samples, paste0('phi_sck[',1:data$N_stands,']'), display_ylim = c(0,1),
                                     ylab = "p(stand concordance)")
prior <- rbeta(1e6, 2.3, 6.07)
lines(density(prior)$y, density(prior)$x, col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_disc_pushforward_quantiles(base_samples, paste0('omega_conc_sck[',1:data$N_stands,']'), display_ylim = c(0,1),
                                     ylab = "p(tree shock | stand concordance)")
prior <- rbeta(1e6, 4.24,  2.80)
lines(density(prior)$y, density(prior)$x, col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_disc_pushforward_quantiles(base_samples, paste0('omega_shutdown[',1:data$N_stands,']'), display_ylim = c(0,1),
                                     ylab = "p(growth shutdown | tree shock)")
prior <- rbeta(1e6, 1.66, 6.86)
lines(density(prior)$y, density(prior)$x, col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_disc_pushforward_quantiles(base_samples, paste0('tau_conc[',1:data$N_species,']'), display_ylim = c(0,1.5),
                                     ylab = "tau_conc")
