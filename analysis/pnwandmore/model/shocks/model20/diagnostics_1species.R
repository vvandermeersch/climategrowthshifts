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

# fit <- readRDS(file.path(folder, 'model19', 'fit_model19_HSGP_PIPO_10stands_19502013_fullidio.rds'))
base_samples <- readRDS(file.path(folder, 'model20', 'basesamples_model20_HSGP_PIPO_1stand_co654_19502008.rds'))
util$check_all_expectand_diagnostics(base_samples)

util$plot_pairs_by_chain(base_samples[['tau_sck']], 'tau_sck',
                         base_samples[['mu_sck']], 'mu_sck')

util$plot_pairs_by_chain(base_samples[['rho_sh']], 'rho_sh',
                         base_samples[['gamma_sh']], 'gamma_sh')

util$plot_pairs_by_chain(base_samples[['thetas_conc[1,2]']], 'thetas_conc[1,2]',
                         base_samples[['thetas_conc[1,3]']], 'thetas_conc[1,3]')

util$plot_pairs_by_chain(base_samples[['thetas_conc[1,2]']], 'thetas_conc[1,2]',
                         base_samples[['thetas_conc[1,1]']], 'thetas_conc[1,1]')

util$plot_pairs_by_chain(base_samples[['omega_shutdown[1]']], 'omega_shutdown[1]',
                         base_samples[['thetas_conc[1,3]']], 'thetas_conc[1,3]')

par(mfrow = c(1,4))
util$plot_expectand_pushforward(base_samples[['tau_small']], 30,
                                display_name = 'tau_small',
                                flim = c(0,1.5))
prior <- rnorm(1e6, 0, log(1.5)/2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(base_samples[['mu_sck']], 30,
                                display_name = 'mu_sck',
                                flim = c(0,3))
prior <- rlnorm(1e6, 0.154, 0.546)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(base_samples[['tau_sck']], 30,
                                display_name = 'tau_sck',
                                flim = c(0,1.5))
prior <- rnorm(1e6, 0, 1.5 / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)

util$plot_expectand_pushforward(base_samples[['tau_idio']], 30,
                                display_name = 'tau_idio',
                                flim = c(0,1.5))
prior <- rnorm(1e6, 0, log(5) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)



par(mfrow = c(2,3))
util$plot_expectand_pushforward(base_samples[['phi_conc[1]']], 30,
                                display_name = 'phi_conc[1]',
                                flim = c(0,1))
prior <- rbeta(1e6, 2, 7)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
n <- 1e5
thetas_baseline <- c(5, 100, 20)
omega_thetas <- 3
alphas = thetas_baseline / omega_thetas + c(1, 1, 1)
thetas_conc <- gtools::rdirichlet(n, alphas)
util$plot_expectand_pushforward(base_samples[['thetas_conc[1,1]']], 30,
                                display_name = 'thetas_conc[1,1]',
                                flim = c(0,1))
lines(density(thetas_conc[,1]), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(base_samples[['thetas_conc[1,2]']], 30,
                                display_name = 'thetas_conc[1,2]',
                                flim = c(0,1))
lines(density(thetas_conc[,2]), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(base_samples[['thetas_conc[1,3]']], 30,
                                display_name = 'thetas_conc[1,3]',
                                flim = c(0,1))
lines(density(thetas_conc[,3]), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(base_samples[['omega_shutdown[1]']], 30,
                                display_name = 'omega_shutdown[1]',
                                flim = c(0,1))
prior <- rbeta(1e6, 1.15, 4.3)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
