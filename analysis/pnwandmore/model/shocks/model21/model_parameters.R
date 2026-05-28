rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(cmdstanr)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)


folder <- '/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model'

data <- readRDS(file.path(folder,'data_18may2026_PIPO_1stand_19502008.rds'))

base_samples <- readRDS(file.path(folder, 'model21', 'basesamples_model21_HGSP_PIPO_1stand_co654_19502008_v3.rds'))
util$check_all_expectand_diagnostics(base_samples)

util$plot_pairs_by_chain(base_samples[['tau_conc']], 'tau_conc',
                         base_samples[['mu_conc']], 'mu_conc')

util$plot_pairs_by_chain(base_samples[['sigma']], 'sigma',
                         base_samples[['tau_clim']], 'tau_clim')

util$plot_pairs_by_chain(base_samples[['rho_sh']], 'rho_sh',
                         base_samples[['gamma_sh']], 'gamma_sh')

par(mfrow = c(2,3))
util$plot_expectand_pushforward(base_samples[['sigma']], 30,
                                display_name = 'sigma',
                                flim = c(0,0.5))
prior <- rnorm(1e6, 0, log(1.3)/2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(base_samples[['tau_clim']], 30,
                                display_name = 'tau_clim',
                                flim = c(0,1))
prior <- rnorm(1e6, 0, log(2)/2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(base_samples[['mu_conc']], 30,
                                display_name = 'mu_conc',
                                flim = c(0,3))
prior <- rlnorm(1e6, -0.79, 0.83)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(base_samples[['tau_conc']], 30,
                                display_name = 'tau_conc',
                                flim = c(0,1.5))
prior <- rnorm(1e6, 0, 3 / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(base_samples[['tau_idio']], 30,
                                display_name = 'tau_idio',
                                flim = c(0,1.5))
prior <- rnorm(1e6, 0, log(5) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)

par(mfrow = c(1,3))
util$plot_expectand_pushforward(base_samples[['phi_sck[1]']], 30,
                                display_name = 'p(stand concordance)',
                                flim = c(0,1))
prior <- rbeta(1e6, 2.3, 6.07)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(base_samples[['omega_conc_sck[1]']], 30,
                                display_name = 'p(tree shock | stand concordance)',
                                flim = c(0,1))
prior <- rbeta(1e6, 5.22, 3.08)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(base_samples[['omega_shutdown[1]']], 30,
                                display_name = 'p(growth shutdown | tree shock)',
                                flim = c(0,1))
prior <- rbeta(1e6, 1.26, 5.97)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)


util$plot_pairs_by_chain(base_samples[['delta_clim[1,7]']], 'delta_clim[1,7]',
                         base_samples[['omega_conc_sck[1]']], 'omega_conc_sck[1]')

util$plot_pairs_by_chain(base_samples[['omega_conc_sck[1]']], 'omega_conc_sck[1]',
                         base_samples[['phi_sck[1]']], 'phi_sck[1]')

par(mfrow = c(3,1))
util$plot_disc_pushforward_quantiles(base_samples, paste0('thetas_idio[',1:data$N_trees,',1]'), display_ylim = c(0,1))
util$plot_disc_pushforward_quantiles(base_samples, paste0('thetas_idio[',1:data$N_trees,',2]'), display_ylim = c(0,1))
util$plot_disc_pushforward_quantiles(base_samples, paste0('thetas_idio[',1:data$N_trees,',3]'), display_ylim = c(0,1))

par(mfrow = c(1,1))
util$plot_disc_pushforward_quantiles(base_samples, paste0('delta_clim_raw[1,',1:data$N_all_years,']'), 
                                     display_ylim = c(-2,2))
