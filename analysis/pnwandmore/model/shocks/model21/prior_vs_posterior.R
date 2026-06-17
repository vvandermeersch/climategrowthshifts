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


par(mfrow = c(2,2))
util$plot_expectand_pushforward(base_samples[['alpha']], 30,
                                display_name = 'alpha',
                                flim = c(-log(5),log(5)))
prior <- rnorm(1e6, 0, log(5)/2.32)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[['beta_gdd']], 30,
                                display_name = 'beta_gdd',
                                flim = c(-log(1.8),log(1.8)))
prior <- rnorm(1e6, 0, log(1.8)/2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[['beta_pre']], 30,
                                display_name = 'beta_pre',
                                flim = c(-log(1.8),log(1.8)))
prior <- rnorm(1e6, 0, log(1.8)/2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[['beta_vpd']], 30,
                                display_name = 'beta_vpd',
                                flim = c(-log(1.8),log(1.8)))
prior <- rnorm(1e6, 0, log(1.8)/2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)

par(mfrow = c(3,2))
util$plot_expectand_pushforward(base_samples[['rho_sp']], 30,
                                display_name = 'rho_sp',
                                flim = c(10,80))
abline(v = data$N_all_years, lty = 2, lwd = 2)
prior <- rlnorm(1e6, 3.7, 0.35)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[['gamma_sp']], 30,
                                display_name = 'gamma_sp',
                                flim = c(0,log(10)))
prior <- rnorm(1e6, 0, log(10)/2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[['rho_sh']], 30,
                                display_name = 'rho_sh',
                                flim = c(1,10))
abline(v = data$N_all_years, lty = 2, lwd = 2)
prior <- rlnorm(1e6, 1.7, 0.26)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[['gamma_sh']], 30,
                                display_name = 'gamma_sh',
                                flim = c(0,log(3)))
prior <- rnorm(1e6, 0, log(3)/2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[['rho_ind']], 30,
                                display_name = 'rho_ind',
                                flim = c(1,10))
abline(v = data$N_all_years, lty = 2, lwd = 2)
prior <- rlnorm(1e6, 1.4, 0.25)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[['gamma_ind']], 30,
                                display_name = 'gamma_ind',
                                flim = c(0,log(3)))
prior <- rnorm(1e6, 0, log(3)/2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)

par(mfrow = c(1,2))
util$plot_expectand_pushforward(base_samples[['tau_clim']], 30,
                                display_name = 'tau_clim',
                                flim = c(0,log(2)))
prior <- rnorm(1e6, 0, log(2)/2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[['sigma']], 30,
                                display_name = 'sigma',
                                flim = c(0,log(1.3)))
prior <- rnorm(1e6, 0, log(1.3)/2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)

par(mfrow = c(3,2))
util$plot_expectand_pushforward(base_samples[['phi_sck[1]']], 30,
                                display_name = 'p(stand concordance)',
                                flim = c(0,1))
prior <- rbeta(1e6, 2.3, 6.07)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[['omega_conc_sck[1]']], 30,
                                display_name = 'p(tree shock | stand concordance)',
                                flim = c(0,1))
prior <- rbeta(1e6, 3.48, 1.86)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[['omega_shutdown[1]']], 30,
                                display_name = 'p(growth shutdown | tree shock)',
                                flim = c(0,1))
prior <- rbeta(1e6, 1.26, 5.97)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
plot.new()
util$plot_expectand_pushforward(base_samples[['mu_conc']], 30,
                                display_name = 'mu_conc',
                                flim = c(0,3))
prior <- rlnorm(1e6, -0.79, 0.83)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[['tau_conc']], 30,
                                display_name = 'tau_conc',
                                flim = c(0,1.5))
prior <- rnorm(1e6, 0, 3 / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)

par(mfrow = c(2,2))
thetas_baseline <- c(100, 5, 0.5)
omega_thetas <- 4
alphas = thetas_baseline / omega_thetas + c(1, 1, 1)
thetas_idio <- gtools::rdirichlet(1e6, alphas)
util$plot_disc_pushforward_quantiles(base_samples, paste0('thetas_idio[',1:data$N_trees,',1]'), display_ylim = c(0,1))
lines(density(thetas_idio[,1])$y, density(thetas_idio[,1])$x, col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_disc_pushforward_quantiles(base_samples, paste0('thetas_idio[',1:data$N_trees,',2]'), display_ylim = c(0,1))
lines(density(thetas_idio[,2])$y, density(thetas_idio[,2])$x, col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_disc_pushforward_quantiles(base_samples, paste0('thetas_idio[',1:data$N_trees,',3]'), display_ylim = c(0,1))
lines(density(thetas_idio[,3])$y, density(thetas_idio[,3])$x, col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[['tau_idio']], 30,
                                display_name = 'tau_idio',
                                flim = c(0,1.5))
prior <- rnorm(1e6, 0, log(5) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)

