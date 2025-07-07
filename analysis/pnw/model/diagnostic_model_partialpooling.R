rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnw"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)

# Model data
data <- readRDS(file.path(wd, 'output/model', 'data_03july2025.rds'))

# Posterior quantification - diagnostics
fit <- readRDS(file.path(wd, 'output/model', 'fit_03july2025_partialpooling_2clades_centered.rds')) # run on Margot
diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit)
base_samples <- util$filter_expectands(samples,
                                       c('alpha', 
                                         'mu_gdd', 'mu_sm', 'mu_vpd',
                                         'tau_gdd', 'tau_sm', 'tau_vpd',
                                         'beta_gdd', 'beta_sm', 'beta_vpd',
                                         'f_tilde_sh',
                                         'rho_sh', 'kappa_sh',
                                         'sigma'),
                                       check_arrays=TRUE)
util$check_all_expectand_diagnostics(base_samples)

par(mfrow=c(1, 2), mar = c(4,4,1,1))
util$plot_expectand_pushforward(samples[['mu_gdd[1]']], 20,
                                flim = c(-0.15,0.25),
                                display_name="beta_gdd")
util$plot_expectand_pushforward(samples[['mu_gdd[2]']], 20,
                                flim = c(-0.15,0.25),
                                display_name="beta_gdd")

par(mfrow=c(1, 2), mar = c(4,4,1,1))
util$plot_expectand_pushforward(samples[['mu_vpd[1]']],50,
                                flim = c(-0.15,0.25),
                                display_name="beta_vpd")
util$plot_expectand_pushforward(samples[['mu_vpd[2]']], 50,
                                flim = c(-0.15,0.25),
                                display_name="beta_vpd")

par(mfrow=c(1, 2), mar = c(4,4,1,1))
util$plot_expectand_pushforward(samples[['mu_sm[1]']], 50,
                                flim = c(-0.15,0.25),
                                display_name="beta_sm")
util$plot_expectand_pushforward(samples[['mu_sm[2]']], 50,
                                flim = c(-0.15,0.25),
                                display_name="beta_sm")



par(mfrow=c(1, 1), mar = c(4,4,1,1))
names <- sapply(1:data$N_species,
                function(sp) paste0('rho_sp[', sp, ']'))
util$plot_disc_pushforward_quantiles(samples, names,
                                     xlab="Species",
                                     ylab="Rho")


par(mfrow=c(1, 1), mar = c(4,4,1,1))
names <- sapply(1:data$N_species,
                function(sp) paste0('beta_vpd[', sp, ']'))
util$plot_disc_pushforward_quantiles(samples, names,
                                     xlab="Species",
                                     ylab="beta_vpd")

par(mfrow=c(1, 1), mar = c(4,4,1,1))
names <- sapply(1:data$N_species,
                function(sp) paste0('beta_gdd[', sp, ']'))
util$plot_disc_pushforward_quantiles(samples, names,
                                     xlab="Species",
                                     ylab="beta_gdd")

