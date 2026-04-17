rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(cmdstanr)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)


folder <- '/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model/'

data <- readRDS(file.path(folder,'data_26mar2026_WestCoast_19202024.rds'))
fit <- readRDS(file.path(folder, 'model17/fit_HSGP_26mar2026_WestCoast_19202024_pi0.rds'))
gc()
fit$diagnostic_summary()

# Load samples
samples <- fit$draws(c('mu_alpha', 'tau_alpha', 'alpha',
                       
                       'mu_gdd', 'tau_gdd', 'beta_gdd', 
                       'mu_vpd', 'tau_vpd', 'beta_vpd', 
                       'mu_pre', 'tau_pre', 'beta_pre', 
                       
                       'mu_rho', 'tau_rho', 'rho_sp', 
                       'mu_gamma', 'tau_gamma', 'gamma_sp',
                       
                       'rho_sh', # 'gamma_sh',
                       'mu_kappa', 'tau_kappa', 'kappa_sh', 
                       
                       'rho_ind', 'gamma_ind',
                       
                       'tau_sck', 'mu_tau_sck', 'tau_tau_sck',
                       'phi_sck0', 'mu_phi_sck', 'tau_phi_sck', 'alpha_phi_sck',
                       'omega_conc_sck0', 'tau_omega_conc_sck', 'alpha_omega_conc_sck',
                       # 'pi_idsc_sck0', 'tau_pi_idsc_sck', 'alpha_pi_idsc_sck',
                       
                       'sigma'
))
names <- dimnames(samples)$variable
chains <- 1:4
samples <- lapply(1:dim(samples)[3], 
                  function(k){t(matrix(samples[1:dim(samples)[1],chains,k], 
                                       nrow = dim(samples)[1], ncol =  length(chains)))})
names(samples) <- names
util$check_all_expectand_diagnostics(samples)


# Intercept and slopes
util$check_all_expectand_diagnostics(util$filter_expectands(samples, c('mu_alpha', 'tau_alpha', 'alpha'),
                                                            check_arrays = T))
layout(matrix(c(1,2,3,3), ncol = 2, byrow = T))
par(mar = c(1,4,4,1), cex.lab = 1)
util$plot_expectand_pushforward(samples[['mu_alpha']], 50, flim= c(-2,2), main = 'mu_alpha')
prior <- rnorm(1e6, 0, log(5)/2.32)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(samples[['tau_alpha']], 50, flim= c(0,2), main = 'tau_alpha')
prior <- rnorm(1e6, 0, (5^0.25)/2.32)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
par(mar = c(3,4,2,1))
util$plot_disc_pushforward_quantiles(samples, paste0('alpha[',1:data$N_species,']'),
                                     ylab = 'alpha')


util$check_all_expectand_diagnostics(util$filter_expectands(samples, c('mu_gdd', 'tau_gdd', 'beta_gdd'),
                                                            check_arrays = T))
layout(matrix(c(1,2,3,3), ncol = 2, byrow = T))
par(mar = c(1,4,4,1), cex.lab = 1)
util$plot_expectand_pushforward(samples[['mu_gdd']], 50, flim= c(-0.3,0.3), main = 'mu_gdd')
prior <- rnorm(1e6, 0, log(1.8) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(samples[['tau_gdd']], 50, flim= c(0,0.1), main = 'tau_gdd')
prior <- rnorm(1e6, 0, log(1.8^0.25) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
par(mar = c(3,4,2,1))
util$plot_disc_pushforward_quantiles(samples, paste0('beta_gdd[',1:data$N_species,']'),
                                     ylab = 'beta_gdd', display_ylim = c(-0.1, 0.1))

util$check_all_expectand_diagnostics(util$filter_expectands(samples, c('mu_pre', 'tau_pre', 'beta_pre'),
                                                            check_arrays = T))
layout(matrix(c(1,2,3,3), ncol = 2, byrow = T))
par(mar = c(1,4,4,1), cex.lab = 1)
util$plot_expectand_pushforward(samples[['mu_pre']], 50, flim= c(-0.3,0.3), main = 'mu_pre')
prior <- rnorm(1e6, 0, log(1.8) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(samples[['tau_pre']], 50, flim= c(0,0.1), main = 'tau_pre')
prior <- rnorm(1e6, 0, log(1.8^0.25) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
par(mar = c(3,4,2,1))
util$plot_disc_pushforward_quantiles(samples, paste0('beta_pre[',1:data$N_species,']'),
                                     ylab = 'beta_pre', display_ylim = c(-0.1, 0.1))

util$check_all_expectand_diagnostics(util$filter_expectands(samples, c('mu_vpd', 'tau_vpd', 'beta_vpd'),
                                                            check_arrays = T))
layout(matrix(c(1,2,3,3), ncol = 2, byrow = T))
par(mar = c(1,4,4,1), cex.lab = 1)
util$plot_expectand_pushforward(samples[['mu_vpd']], 50, flim= c(-0.3,0.3), main = 'mu_vpd')
prior <- rnorm(1e6, 0, log(1.8) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(samples[['tau_vpd']], 50, flim= c(0,0.1), main = 'tau_vpd')
prior <- rnorm(1e6, 0, log(1.8^0.25) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
par(mar = c(3,4,2,1))
util$plot_disc_pushforward_quantiles(samples, paste0('beta_vpd[',1:data$N_species,']'),
                                     ylab = 'beta_vpd', display_ylim = c(-0.1, 0.1))

# Gaussian processes

layout(matrix(c(1,2,3,3), ncol = 2, byrow = T))
par(mar = c(1,4,4,1), cex.lab = 1)
util$plot_expectand_pushforward(samples[['mu_rho']], 50, flim= c(0,80), main = 'mu_rho')
prior <- rlnorm(1e6, 3.7, 0.35)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(samples[['tau_rho']], 50, flim= c(0,30), main = 'tau_rho')
prior <- rnorm(1e6, 0, 6 / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
prior <- rnorm(1e6, 0, 30 / 2.57)
lines(density(prior), col = util$c_mid_teal, lwd = 1.5, lty = 1)
par(mar = c(3,4,2,1))
util$plot_disc_pushforward_quantiles(samples, paste0('rho_sp[',1:data$N_species,']'),
                                     ylab = 'rho_sp', display_ylim = c(0, 80))

layout(matrix(c(1,2,3,3), ncol = 2, byrow = T))
par(mar = c(1,4,4,1), cex.lab = 1)
util$plot_expectand_pushforward(samples[['mu_gamma']], 50, flim= c(0,2), main = 'mu_gamma')
prior <- rnorm(1e6, 0, log(10) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(samples[['tau_gamma']], 50, flim= c(0,0.5), main = 'tau_gamma')
prior <- rnorm(1e6, 0, 0.23 / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
par(mar = c(3,4,2,1))
util$plot_disc_pushforward_quantiles(samples, paste0('gamma_sp[',1:data$N_species,']'),
                                     ylab = 'gamma_sp', display_ylim = c(0, 0.8))


par(mfrow = c(1,2))
util$plot_expectand_pushforward(samples[['rho_ind']], 50, flim= c(1,8), 'rho_ind')
prior <- rlnorm(1e6, 0.80, 0.40)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
prior <- rlnorm(1e6, 1.4, 0.35)
lines(density(prior), col = util$c_mid_teal, lwd = 2, lty = 1)
util$plot_expectand_pushforward(samples[['gamma_ind']], 50, flim= c(0,1), 'gamma_ind')
prior <- rnorm(1e6, 0, log(3) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)

par(mfrow = c(1,1))
util$plot_expectand_pushforward(samples[['rho_sh']], 50, flim= c(1,8), 'rho_sh')
prior <- rlnorm(1e6, 0.4, 0.3)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)

layout(matrix(c(1,2,3,3), ncol = 2, byrow = T))
par(mar = c(1,4,4,1), cex.lab = 1)
util$plot_expectand_pushforward(samples[['mu_kappa']], 50, flim= c(0,2), main = 'mu_kappa')
prior <- rlnorm(1e6, 0, 0.41 / 2.32)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(samples[['tau_kappa']], 50, flim= c(0,0.5), main = 'tau_kappa')
prior <- rnorm(1e6, 0, 0.2 / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
par(mar = c(3,4,2,1))
util$plot_disc_pushforward_quantiles(samples, paste0('kappa_sh[',1:data$N_species,']'),
                                     ylab = 'kappa_sh', display_ylim = c(0, 1))


# Shocks
layout(matrix(c(1,2,3,3), ncol = 2, byrow = T))
par(mar = c(1,4,4,1), cex.lab = 1)
util$plot_expectand_pushforward(samples[['mu_tau_sck']], 50, flim= c(0, 5), main = 'mu_tau_sck')
prior <- rnorm(1e6, 0, 5 / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(samples[['tau_tau_sck']], 50, flim= c(0,2), main = 'tau_tau_sck')
prior <- rnorm(1e6, 0, 1 / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
par(mar = c(3,4,2,1))
util$plot_disc_pushforward_quantiles(samples, paste0('tau_sck[',1:data$N_species,']'),
                                     ylab = 'tau_sck', display_ylim = c(0, 4))


shocks <- fit$draws(c('phi_sck', 'omega_conc_sck'))
names <- dimnames(shocks)$variable
shocks <- lapply(1:dim(shocks)[3], 
                 function(k){t(matrix(shocks[1:dim(shocks)[1],1:dim(shocks)[2],k], 
                                      nrow = dim(shocks)[1], ncol = dim(shocks)[2]))})
names(shocks) <- names
par(mfrow = c(2,1), mar = c(4,4,1,1))
util$plot_disc_pushforward_quantiles(shocks, paste0('phi_sck[', 1:data$N_stands, ']'), display_ylim = c(0,1),
                                     ylab = 'p(stand conc.)')
util$plot_disc_pushforward_quantiles(shocks, paste0('omega_conc_sck[', 1:data$N_stand_species, ']'), display_ylim = c(0,1),
                                     ylab = 'p(tree sck | stand conc.)')
util$plot_disc_pushforward_quantiles(shocks, paste0('pi_idsc_sck[', 1:data$N_stand_species, ']'), display_ylim = c(0,0.5),
                                     ylab = 'p(tree idiosync. sck)')

layout(matrix(c(1,2,3,3), ncol = 2, byrow = T))  
par(mar = c(1,4,4,1), cex.lab = 1)
util$plot_expectand_pushforward(samples[['phi_sck0']], 50, flim= c(0, 1), main = 'phi_sck0')
prior <- rbeta(1e6, 2, 20)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(samples[['tau_phi_sck']], 50, flim= c(0,5), main = 'tau_phi_sck')
prior <- rnorm(1e6, 0, 1)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
par(mar = c(3,4,2,1))
util$plot_disc_pushforward_quantiles(shocks, paste0('phi_sck[', 1:data$N_stands, ']'), display_ylim = c(0,1),
                                     ylab = 'p(stand conc.)')

layout(matrix(c(1,2,3,3), ncol = 2, byrow = T))
par(mar = c(1,4,4,1), cex.lab = 1)
util$plot_expectand_pushforward(samples[['omega_conc_sck0']], 50, flim= c(0, 1), main = 'omega_conc_sck0')
prior <- rbeta(1e6, 12.28, 34)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(samples[['tau_omega_conc_sck']], 50, flim= c(0,5), main = 'tau_omega_conc_sck')
prior <- rnorm(1e6, 0, 1)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
par(mar = c(3,4,2,1))
util$plot_disc_pushforward_quantiles(shocks, paste0('omega_conc_sck[', 1:data$N_stand_species, ']'), display_ylim = c(0,1),
                                     ylab = 'p(tree sck | stand conc.)')


layout(matrix(c(1,2,3,3), ncol = 2, byrow = T))
par(mar = c(1,4,4,1), cex.lab = 1)
util$plot_expectand_pushforward(samples[['pi_idsc_sck0']], 50, flim= c(0, 1), main = 'pi_idsc_sck0')
prior <- rbeta(1e6, 2, 100)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(samples[['tau_pi_idsc_sck']], 50, flim= c(0,5), main = 'tau_pi_idsc_sck')
prior <- rnorm(1e6, 0, 0.5/2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
par(mar = c(3,4,2,1))
util$plot_disc_pushforward_quantiles(shocks, paste0('pi_idsc_sck[', 1:data$N_stand_species, ']'), display_ylim = c(0,0.5),
                                     ylab = 'p(tree idiosync. sck)')


