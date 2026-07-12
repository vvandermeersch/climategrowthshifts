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

fit <- readRDS(file.path(wd, 'output/model/model21bis', 'fit_model21bis_HGSP_21species_115stands_withinit_4threads_geommeanconstraint.rds'))
fit$time()
fit$diagnostic_summary()

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
  "tau_conc",
  "mu_phi", "sigma_phi", "logit_phi_sck", "phi_sck",
  "mu_omega_conc", "sigma_omega_conc", "logit_omega_conc_sck", "omega_conc_sck",
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
  "mu_alpha", "sigma_alpha", "alpha",
  "sigma_alpha_stand", "alpha_stand",
  "mu_beta_gdd", "sigma_beta_gdd", "beta_gdd",
  "mu_beta_pre", "sigma_beta_pre", "beta_pre",
  "mu_beta_vpd", "sigma_beta_vpd", "beta_vpd",
  "tau_clim", # "kappa_clim_free", "kappa_clim",
  'log_kappa_clim', 'kappa_clim',
  "mu_log_rho", "sigma_log_rho", "rho_merged",
  "mu_log_gamma", "sigma_log_gamma", "log_gamma_merged", 'gamma_merged',
  "tau_conc",
  "mu_phi", "sigma_phi", "logit_phi_sck", 
  "mu_omega_conc", "sigma_omega_conc", "logit_omega_conc_sck", 
  "mu_omega_shutdown", "sigma_omega_shutdown", "logit_omega_shutdown", 
  "thetas_idio", "tau_idio",
  # 'thetas_baseline', 
  # 'omega_thetas', 
  'thetas_idio',
  "sigma"
)

params <- c('delta_clim')

for(p in params){
  message(paste0('\n\n---------\nParameter(s):',p))
  rest <- util$filter_expectands(base_samples, p, 
                                 check_arrays = TRUE)
  util$check_all_expectand_diagnostics(rest)
}

util$plot_pairs_by_chain(base_samples[['mu_log_rho']], 'mu_log_rho',
                         log(base_samples[['sigma_log_rho']]), 'log(sigma_log_rho)')

util$plot_pairs_by_chain(base_samples[['mu_log_gamma']], 'mu_log_gamma',
                         log(base_samples[['sigma_log_gamma']]), 'log(sigma_log_gamma)')

util$plot_pairs_by_chain(base_samples[['tau_clim']], 'tau_clim',
                         base_samples[['tau_idio']], 'tau_idio')

util$plot_pairs_by_chain(base_samples[['alpha[1]']], 'alpha[1]',
                         base_samples[['alpha_stand[2]']], 'alpha_stand[2]')

util$plot_pairs_by_chain(base_samples[['mu_alpha']], 'mu_alpha',
                         log(base_samples[['sigma_alpha']]), 'log(sigma_alpha)')


util$plot_pairs_by_chain(base_samples[['omega_conc_sck[6]']], 'omega_conc_sck[6]',
                         base_samples[['phi_sck[6]']], 'phi_sck[6]')

util$plot_pairs_by_chain(base_samples[['mu_phi']], 'mu_phi',
                         log(base_samples[['sigma_phi']]), 'log(sigma_phi)')
util$plot_pairs_by_chain(base_samples[['mu_omega_conc']], 'mu_omega_conc',
                         log(base_samples[['sigma_omega_conc']]), 'log(sigma_omega_conc)')
util$plot_pairs_by_chain(base_samples[['mu_omega_shutdown']], 'mu_omega_shutdown',
                         log(base_samples[['sigma_omega_shutdown']]), 'log(sigma_omega_shutdown)')

for(sp in 1:data$N_species){
  util$plot_pairs_by_chain(base_samples[[paste0('rho_merged[',sp,']')]], paste0('rho_merged[',sp,']'),
                           base_samples[[paste0('gamma_merged[',sp,']')]], paste0('gamma_merged[',sp,']'))
}

for(sp in 1:data$N_species){
  util$plot_pairs_by_chain(base_samples[[paste0('alpha[',sp,']')]], paste0('alpha[',sp,']'),
                           base_samples[[paste0('beta_vpd[',sp,']')]], paste0('beta_vpd[',sp,']'))
}

for(sp in 1:data$N_species){
  util$plot_pairs_by_chain(base_samples[[paste0('alpha[',sp,']')]], paste0('alpha[',sp,']'),
                           base_samples[[paste0('kappa_clim[',sp,']')]], paste0('kappa_clim[',sp,']'))
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

# alpha slope
layout(matrix(c(1,2,3,3), ncol = 2, byrow = T))
util$plot_expectand_pushforward(base_samples[['mu_alpha']], 50,
                                display_name = 'mu_alpha',
                                flim = c(-log(10),log(10)))
prior <- rnorm(1e6, 0, log(10)/2.32)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[['sigma_alpha']], 50,
                                display_name = 'sigma_alpha', 
                                flim = c(0,log(10)))
prior <- rnorm(1e6, 0, log(10)/2.32)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_disc_pushforward_quantiles(base_samples, paste0('alpha[',1:data$N_species,']'), display_ylim = c(-3,3),
                                     ylab = "alpha")

# alpha_stand slope
layout(matrix(c(1,2,3,3), ncol = 2, byrow = T))
plot.new()
# util$plot_expectand_pushforward(base_samples[['mu_alpha_stand']], 50,
#                                 display_name = 'mu_alpha_stand',http://127.0.0.1:32015/graphics/9e462507-8d98-4175-aa79-8c76d2daf6b1.png
#                                 flim = c(-log(10),log(10)))
# prior <- rnorm(1e6, 0, log(10)/2.32)
# lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[['sigma_alpha_stand']], 50,
                                display_name = 'sigma_alpha_stand', 
                                flim = c(0,log(10)))
prior <- rnorm(1e6, 0, log(10)/2.32)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_disc_pushforward_quantiles(base_samples, paste0('alpha_stand[',1:data$N_stands,']'), display_ylim = c(-3,3),
                                     ylab = "alpha_stand")

# Individual GP
par(mar = c(4,4,1,1))
layout(matrix(c(1,2,3,3), ncol = 2, byrow = T))
util$plot_expectand_pushforward(base_samples[['mu_log_rho']], 50,
                                display_name = 'mu_log_rho',
                                flim = c(1,5))
prior <- rnorm(1e6, log(15), 0.5)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[['sigma_log_rho']], 50,
                                display_name = 'sigma_log_rho', 
                                flim = c(0,1))
prior <- rnorm(1e6, 0, 0.5)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_disc_pushforward_quantiles(base_samples, paste0('rho_merged[',1:data$N_species,']'), display_ylim = c(1,20),
                                     ylab = "rho_merged", xlab = 'Species')

# GDD slope
layout(matrix(c(1,2,3,3), ncol = 2, byrow = T))
util$plot_expectand_pushforward(base_samples[['mu_beta_gdd']], 50,
                                display_name = 'mu_beta_gdd',
                                flim = c(-0.25,0.25))
prior <- rnorm(1e6, 0, log(1.3) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[['sigma_beta_gdd']], 50,
                                display_name = 'sigma_beta_gdd', 
                                flim = c(0,0.25))
prior <- rnorm(1e6, 0, log(1.3) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_disc_pushforward_quantiles_sign(base_samples, paste0('beta_gdd[',1:data$N_species,']'), display_ylim = c(-0.15,0.15),
                                     ylab = "beta_gdd", xticklabs = data$uniq_species_ids)

# VPD slope
layout(matrix(c(1,2,3,3), ncol = 2, byrow = T))
util$plot_expectand_pushforward(base_samples[['mu_beta_vpd']], 50,
                                display_name = 'mu_beta_vpd',
                                flim = c(-0.25,0.25))
prior <- rnorm(1e6, 0, log(1.3) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[['sigma_beta_vpd']], 50,
                                display_name = 'sigma_beta_vpd', 
                                flim = c(0,0.25))
prior <- rnorm(1e6, 0, log(1.3) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_disc_pushforward_quantiles_sign(base_samples, paste0('beta_vpd[',1:data$N_species,']'), display_ylim = c(-0.15,0.15),
                                     ylab = "beta_vpd", xticklabs = data$uniq_species_ids)

# Winter precipitation slope
layout(matrix(c(1,2,3,3), ncol = 2, byrow = T))
util$plot_expectand_pushforward(base_samples[['mu_beta_pre']], 50,
                                display_name = 'mu_beta_pre',
                                flim = c(-0.25,0.25))
prior <- rnorm(1e6, 0, log(1.3) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[['sigma_beta_pre']], 50,
                                display_name = 'sigma_beta_pre', 
                                flim = c(0,0.25))
prior <- rnorm(1e6, 0, log(1.3) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_disc_pushforward_quantiles_sign(base_samples, paste0('beta_pre[',1:data$N_species,']'), display_ylim = c(-0.15,0.15),
                                     ylab = "beta_pre", xticklabs = data$uniq_species_ids)


# Missing mesoenvironment
layout(matrix(c(1,2,3,3), ncol = 2, byrow = T))
util$plot_expectand_pushforward(base_samples[['tau_clim']], 50,
                                display_name = 'tau_clim',
                                flim = c(0,log(2)))
prior <- rnorm(1e6, 0, log(2)/2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
plot.new()
util$plot_disc_pushforward_quantiles(base_samples, paste0('kappa_clim[',1:data$N_species,']'), display_ylim = c(0,2),
                                     ylab = "kappa_clim")
# kappa = 2: twice more sensitive than the average species (kappa = 1)

# p(stand concordance)
par(mar = c(4,4,1,1))
layout(matrix(c(1,2,3,3), ncol = 2, byrow = T))
util$plot_expectand_pushforward(base_samples[['mu_phi']], 50,
                                display_name = 'mu_phi',
                                flim = c(-4,-1))
prior <- rnorm(1e6, -2.17, 0.40)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[['sigma_phi']], 50,
                                display_name = 'sigma_phi', 
                                flim = c(0,2))
prior <- rnorm(1e6, 0, 1)
lines(density(prior), col = util$c_mid_teal, lwd = 2, lty = 2)
util$plot_disc_pushforward_quantiles(base_samples, paste0('phi_sck[',order(data$uniq_stand_lat),']'), display_ylim = c(0,1),
                                     ylab = "p(stand concordance)", xlab = 'Stands (order by latitude)')
abline(v = 0.5, lty = 2, col = 'grey30')
for(st in unique(states)){
  abline(v = max(which(states == st))+ 0.5, lty = 2, col = 'grey30') 
}

# p(tree shock | stand concordance)
layout(matrix(c(1,2,3,3), ncol = 2, byrow = T))
util$plot_expectand_pushforward(base_samples[['mu_omega_conc']], 50,
                                display_name = 'mu_omega_conc',
                                flim = c(-2.5,1))
prior <- rnorm(1e6, -0.49, 0.46)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[['sigma_omega_conc']], 50,
                                display_name = 'sigma_omega_conc', 
                                flim = c(0,4))
prior <- rnorm(1e6, 0, 1.5)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_disc_pushforward_quantiles(base_samples, paste0('omega_conc_sck[',1:data$N_stand_species,']'), display_ylim = c(0,1),
                                     ylab = "p(tree shock | stand concordance)", xlab = 'Stands x species')


# p(growth shutdown | tree shock)
layout(matrix(c(1,2,3,3), ncol = 2, byrow = T))
util$plot_expectand_pushforward(base_samples[['mu_omega_shutdown']], 50,
                                display_name = 'mu_omega_shutdown',
                                flim = c(-8,0))
prior <- rnorm(1e6, -3.40, 0.61)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[['sigma_omega_shutdown']], 50,
                                display_name = 'sigma_omega_shutdown', 
                                flim = c(0,4))
prior <- rnorm(1e6, 0, 0.95)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_disc_pushforward_quantiles(base_samples, paste0('omega_shutdown[',1:data$N_stand_species,']'), display_ylim = c(0,1),
                                     ylab = "p(growth shutdown | tree shock)", xlab = 'Stands x species')


# p(tree shock | stand concordance)
par(mfrow=c(1,1))
util$plot_disc_pushforward_quantiles_sign(base_samples, paste0('tau_conc[',phylo_order,']'), display_ylim = c(0,3),
                                          xticklabs = data$uniq_species_ids[phylo_order],
                                          ylab = "Concordant shock amplitude", xlab = '', ignore_sign = T)


# p(idiosyncratic stuff)
thetas_baseline <- c(100, 20, 0.001)
omega_thetas <- 10
alphas = thetas_baseline / omega_thetas + c(1, 1, 1)
thetas_idio <- gtools::rdirichlet(1e6, alphas)
layout(matrix(c(1,2,3,3, 4,4), ncol = 2, byrow = T))
util$plot_expectand_pushforward(base_samples[['tau_idio']], 50,
                                display_name = 'tau_idio',
                                flim = c(0,log(5)))
prior <- rnorm(1e6, 0, log(5) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
plot.new()
util$plot_disc_pushforward_quantiles(base_samples, paste0('thetas_idio[',1:data$N_trees,',2]'), display_ylim = c(0,0.8),
                                     ylab = "p(idiosyncratic shock)")
probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
plot_quantiles <- quantile(thetas_idio[,2], probs)
rect(xleft= -50, xright = -20, ybottom = plot_quantiles[1], ytop = plot_quantiles[9],
     col = paste0(util$c_light_teal,30), border = NA)
rect(xleft= -50, xright = -20, ybottom = plot_quantiles[2], ytop = plot_quantiles[8],
     col = paste0(util$c_light_teal,50), border = NA)
rect(xleft= -50, xright = -20, ybottom = plot_quantiles[3], ytop = plot_quantiles[7],
     col = paste0(util$c_mid_teal,50), border = NA)
rect(xleft= -50, xright = -20, ybottom = plot_quantiles[4], ytop = plot_quantiles[6],
     col = paste0(util$c_mid_teal,50), border = NA)
segments(x0 = -50, x1 =-20, y0 = plot_quantiles[5], col=util$c_dark_teal, lwd=2)
util$plot_disc_pushforward_quantiles(base_samples, paste0('thetas_idio[',1:data$N_trees,',3]'), display_ylim = c(0,0.2),
                                     ylab = "p(idiosyncratic shutdown)", xlab = "Trees")

plot_quantiles <- quantile(thetas_idio[,3], probs)
rect(xleft= -50, xright = -20, ybottom = plot_quantiles[1], ytop = plot_quantiles[9],
     col = paste0(util$c_light_teal,30), border = NA)
rect(xleft= -50, xright = -20, ybottom = plot_quantiles[2], ytop = plot_quantiles[8],
     col = paste0(util$c_light_teal,50), border = NA)
rect(xleft= -50, xright = -20, ybottom = plot_quantiles[3], ytop = plot_quantiles[7],
     col = paste0(util$c_mid_teal,50), border = NA)
rect(xleft= -50, xright = -20, ybottom = plot_quantiles[4], ytop = plot_quantiles[6],
     col = paste0(util$c_mid_teal,50), border = NA)
segments(x0 = -50, x1 =-20, y0 = plot_quantiles[5], col=util$c_dark_teal, lwd=2)
