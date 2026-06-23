rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(cmdstanr)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)

folder <- '/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model'

data1 <- readRDS(file.path(folder,'data_18may2026_PIPO_30stands_19502022.rds'))
base_samples1 <- readRDS(file.path(wd, 'output/model/model21bis', 'basesamples_model21bis_HGSP_PIPO_30stands.rds'))
fit1 <-  readRDS(file.path(wd, 'output/model/model21bis', 'fit_model21bis_HGSP_PIPO_30stands.rds'))
util$check_all_expectand_diagnostics(base_samples1)

data2 <- readRDS(file.path(folder,'data_18may2026_PIPO_30stands_19502022_sample2.rds'))
base_samples2 <- readRDS(file.path(wd, 'output/model/model21bis', 'basesamples_model21bis_HGSP_PIPO_30stands_sample2.rds'))
fit2 <-  readRDS(file.path(wd, 'output/model/model21bis', 'fit_model21bis_HGSP_PIPO_30stands_sample2.rds'))
util$check_all_expectand_diagnostics(base_samples2)


par(mfrow = c(2,2), mar = c(4,4,1,1))
util$plot_expectand_pushforward(base_samples1[['alpha']], 50,
                                display_name = 'alpha',
                                flim = c(-log(10),log(10)))
util$plot_expectand_pushforward(base_samples2[['alpha']], 50,
                                display_name = 'alpha', col = util$c_light,
                                flim = c(-log(10),log(10)), add = T)
prior <- rnorm(1e6, 0, log(10)/2.32)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples1[['beta_gdd']], 50,
                                display_name = 'beta_gdd',
                                flim = c(-log(1.8),log(1.8)))
util$plot_expectand_pushforward(base_samples2[['beta_gdd']], 50,
                                display_name = 'beta_gdd', col = util$c_light,
                                flim = c(-log(1.8),log(1.8)), add = T)
prior <- rnorm(1e6, 0, log(1.8)/2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples1[['beta_pre']], 50,
                                display_name = 'beta_pre',
                                flim = c(-log(1.8),log(1.8)))
util$plot_expectand_pushforward(base_samples2[['beta_pre']], 50,
                                display_name = 'beta_pre', col = util$c_light,
                                flim = c(-log(1.8),log(1.8)), add = T)
prior <- rnorm(1e6, 0, log(1.8)/2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples1[['beta_vpd']], 50,
                                display_name = 'beta_vpd',
                                flim = c(-log(1.8),log(1.8)))
util$plot_expectand_pushforward(base_samples2[['beta_vpd']], 50,
                                display_name = 'beta_vpd', col = util$c_light,
                                flim = c(-log(1.8),log(1.8)), add = T)
prior <- rnorm(1e6, 0, log(1.8)/2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)


par(mfrow = c(2,2), mar = c(4,4,1,1))
util$plot_expectand_pushforward(base_samples1[['rho_merged']], 50,
                                display_name = 'rho_merged',
                                flim = c(0,30))
util$plot_expectand_pushforward(base_samples2[['rho_merged']], 50,
                                display_name = 'rho_merged', 
                                flim = c(0,30), col = util$c_light, add = T)
abline(v = data1$N_all_years, lty = 2, lwd = 2)
prior <- rlnorm(1e6, log(15), 0.5)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples1[['gamma_merged']], 50,
                                display_name = 'gamma_merged',
                                flim = c(0,log(2.3)))
util$plot_expectand_pushforward(base_samples2[['gamma_merged']], 50,
                                display_name = 'gamma_merged',
                                flim = c(0,log(2.3)), col = util$c_light, add = T)
prior <- rnorm(1e6, 0, log(2.3) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples1[['rho_sh']], 50,
                                display_name = 'rho_sh',
                                flim = c(10,100))
util$plot_expectand_pushforward(base_samples2[['rho_sh']], 50,
                                display_name = 'rho_sh',
                                flim = c(10,100), col = util$c_light, add = T)
abline(v = data1$N_all_years, lty = 2, lwd = 2)
prior <- rlnorm(1e6, 3.0, 0.42)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples1[['gamma_sh']], 30,
                                display_name = 'gamma_sh',
                                flim = c(0,log(3)))
util$plot_expectand_pushforward(base_samples2[['gamma_sh']], 30,
                                display_name = 'gamma_sh',
                                flim = c(0,log(3)), col = util$c_light, add = T)
prior <- rnorm(1e6, 0, log(3)/2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)


par(mfrow = c(2,2))
util$plot_expectand_pushforward(base_samples1[['tau_clim']], 50,
                                display_name = 'tau_clim',
                                flim = c(0,log(2)))
util$plot_expectand_pushforward(base_samples2[['tau_clim']], 50,
                                display_name = 'tau_clim',
                                flim = c(0,log(2)), col = util$c_light, add = T)
prior <- rnorm(1e6, 0, log(2)/2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples1[['sigma']], 50,
                                display_name = 'sigma',
                                flim = c(0,log(1.3)))
util$plot_expectand_pushforward(base_samples2[['sigma']], 50,
                                display_name = 'sigma',
                                flim = c(0,log(1.3)), col = util$c_light, add = T)
prior <- rnorm(1e6, 0, log(1.1) /2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples1[[paste0('tau_conc')]], 60,
                                display_name = 'tau_conc',
                                flim = c(0,3))
util$plot_expectand_pushforward(base_samples2[[paste0('tau_conc')]], 60,
                                display_name = 'tau_conc',
                                flim = c(0,3), col = util$c_light, add = T)
prior <- rnorm(1e6, 0, 3 / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples1[['tau_idio']], 60,
                                display_name = 'tau_idio',
                                flim = c(0,1.5))
util$plot_expectand_pushforward(base_samples2[['tau_idio']], 60,
                                display_name = 'tau_idio',
                                flim = c(0,1.5), col = util$c_light, add = T)
prior <- rnorm(1e6, 0, log(5) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)