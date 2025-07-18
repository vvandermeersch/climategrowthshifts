rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)

# Model data
data <- readRDS(file.path(wd, 'output/model', 'data_11july2025_legacies.rds'))

# Posterior quantification - diagnostics
fit <- readRDS(file.path(wd, 'output/model', 'fit_11july2025_completepooling_legacies.rds')) # run on Margot
diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit)
base_samples <- util$filter_expectands(samples,
                                       c('alpha', 
                                         'beta_gdd', 'beta_sm', 
                                         'beta_vpd', 'rho_beta_vpd',
                                         'rho_sp', 'gamma_sp',
                                         'f_tilde_sh',
                                         'rho_sh', 'kappa_sh',
                                         'sigma'),
                                       check_arrays=TRUE)
util$check_all_expectand_diagnostics(base_samples)


util$plot_pairs_by_chain(log(samples[[paste0('rho_beta_vpd')]]), paste0('log(rho_beta_vpd)'), samples[[paste0('beta_vpd')]], paste0('beta_vpd'))
par(mfrow=c(1, 1), mar = c(4,4,1,1))
names <- sapply(1:4,
                function(sp) paste0('beta_vpd_leg[', sp, ']'))
util$plot_disc_pushforward_quantiles(samples, names,
                                     xlab="Lag",
                                     xticklabs=c(0:3),
                                     ylab="VPD effect")

par(mfrow=c(3, 1), mar = c(4,4,1,1))
util$plot_expectand_pushforward(samples[['beta_gdd']], 200,
                                flim = c(-0.15,0.15),
                                display_name="beta_gdd")
util$plot_expectand_pushforward(samples[['beta_sm']], 200,
                                flim = c(-0.15,0.15),
                                display_name="beta_sm")
util$plot_expectand_pushforward(samples[['beta_vpd']], 200,
                                flim = c(-0.15,0.15),
                                display_name="beta_vpd")


util$plot_pairs_by_chain(samples[[paste0('rho_sp[33]')]], paste0('rho_sp[33]'), samples[[paste0('beta_vpd')]], paste0('beta_vpd'))
