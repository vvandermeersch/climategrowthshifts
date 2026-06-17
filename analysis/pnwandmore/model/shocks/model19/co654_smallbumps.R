rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(cmdstanr)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)


folder <- '/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model'

data <- readRDS(file.path(folder,'data_18may2026_PIPO_1stand_co654_19502008.rds'))

fit <- readRDS(file.path(folder, 'model19', 'fit_model19_HSGP_PIPO_1stand_co654_19502008.rds'))
samples <- fit$draws(c('alpha',
                       'beta_gdd', 'beta_vpd', 'beta_pre',
                       
                       'rho_sp', 'gamma_sp', 'f_tilde',
                       'rho_sh', 'gamma_sh', 'f_tilde_sh',
                       'rho_ind', 'gamma_ind', 'f_ind_tilde',
                       
                       'mu_conc', 'tau_conc',
                       'phi_sck',
                       'omega_conc_sck',
                       'omega_shutdown',
                       
                       'tau_idio', 
                       'thetas_idio',
                       
                       'sigma'
))
names <- dimnames(samples)$variable
samples <- lapply(1:dim(samples)[3],
                  function(k){t(matrix(samples[1:dim(samples)[1],1:dim(samples)[2],k],
                                       nrow = dim(samples)[1], ncol = dim(samples)[2]))})
names(samples) <- names
util$check_all_expectand_diagnostics(samples)
stand_samples <- fit$draws(c('f_sh'))
names <- dimnames(stand_samples)$variable
stand_samples <- lapply(1:dim(stand_samples)[3],
                  function(k){t(matrix(stand_samples[1:dim(stand_samples)[1],1:dim(stand_samples)[2],k],
                                       nrow = dim(stand_samples)[1], ncol = dim(stand_samples)[2]))})
names(stand_samples) <- names
util$check_all_expectand_diagnostics(stand_samples)
tree_samples <- readRDS(file.path(folder, 'model19', 'treesamples_model19_HSGP_PIPO_1stand_co654_19502008.rds'))

util$plot_pairs_by_chain(base_samples[['tau_conc']], 'tau_conc',
                         base_samples[['mu_conc']], 'mu_conc')

util$plot_pairs_by_chain(base_samples[['tau_conc']], 'tau_conc',
                         base_samples[['sigma']], 'sigma')

par(mfrow = c(1,3))
util$plot_expectand_pushforward(base_samples[['mu_conc']], 30,
                                display_name = 'mu_conc',
                                flim = c(0,3))
prior <- rlnorm(1e6, 0.154, 0.546)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(base_samples[['tau_conc']], 30,
                                display_name = 'tau_conc',
                                flim = c(0,1.5))
prior <- rnorm(1e6, 0, 1.5 / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)

util$plot_expectand_pushforward(base_samples[['tau_idio']], 30,
                                display_name = 'tau_idio',
                                flim = c(0,1.5))
prior <- rnorm(1e6, 0, log(5) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
