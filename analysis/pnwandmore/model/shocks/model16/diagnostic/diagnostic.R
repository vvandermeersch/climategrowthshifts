rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(cmdstanr)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)


folder <- '/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model/'

data <- readRDS(file.path(folder,'data_11mar2026_idahoPIPO_standclimate_18962024.rds'))
fit <- readRDS(file.path(folder, 'model16/fit_11mar2026_idahoPIPO.rds'))
gc()


# samples <- fit$draws(c('mu_alpha', 'tau_alpha', 'alpha',
# 
#                        'mu_gdd', 'tau_gdd', 'beta_gdd',
#                        'mu_vpd', 'tau_vpd', 'beta_vpd',
#                        'mu_pre', 'tau_pre', 'beta_pre',
# 
#                        'mu_rho', 'tau_rho','rho_sp',
#                        'mu_gamma', 'tau_gamma', 'gamma_sp',
# 
#                        'rho_ind', 'gamma_ind',
#                        'rho_sh'
#                        ))

samples <- fit$draws(c('alpha', 'beta_gdd', 'beta_vpd', 'beta_pre',
                       
                       'rho_sp', 'gamma_sp',
                       
                       'rho_ind', 'gamma_ind',
                       'rho_sh'
))


names <- dimnames(samples)$variable
samples <- lapply(1:dim(samples)[3], 
                  function(k){t(matrix(samples[1:dim(samples)[1],1:dim(samples)[2],k], 
                                       nrow = dim(samples)[1], ncol = dim(samples)[2]))})
names(samples) <- names
util$check_all_expectand_diagnostics(samples)

#### Alphas ####
util$plot_pairs_by_chain(samples[['mu_alpha[1]']], 'mu_alpha[1]',
                         log(samples[['tau_alpha[1]']]), 'log(tau_alpha)')
util$plot_pairs_by_chain(samples[['alpha[5]']], 'alpha[5]',
                         log(samples[['tau_alpha[1]']]), 'log(tau_alpha)')

par(mfrow = c(2,1))
util$plot_expectand_pushforward(samples[['mu_alpha[1]']], 50, flim= c(-2,2), 'mu_alpha')
mu <- rnorm(1e6, 0, log(5)/2.32)
lines(density(mu), col = 'darkblue')

util$plot_expectand_pushforward(samples[['tau_alpha[1]']], 50, flim= c(0,2), 'tau_alpha')
tau <- rnorm(1e6, 0, log(5^0.25)/2.32)
lines(density(tau), col = 'darkblue')
tau <- rnorm(1e6, 0, log(5)/2.32)
lines(density(tau), col = 'darkblue', lty = 2)

par(mfrow = c(5,2), mar = c(4,4,1,1))
for(i in 1:10){
  util$plot_expectand_pushforward(samples[[paste0('alpha[',i,']')]], 50, flim= c(-2,2))
}

#### Betas ####
util$plot_pairs_by_chain(samples[['mu_gdd[1]']], 'mu_gdd[1]',
                         log(samples[['tau_gdd[1]']]), 'log(tau_gdd)')
util$plot_pairs_by_chain(samples[['beta_gdd[5]']], 'beta_gdd[5]',
                         log(samples[['tau_gdd[1]']]), 'log(tau_gdd)')

par(mfrow = c(5,2), mar = c(4,4,1,1))
for(i in 1:10){
  util$plot_expectand_pushforward(samples[[paste0('beta_gdd[',i,']')]], 50, flim= c(-0.25,0.25))
}

par(mfrow = c(2,1))
util$plot_expectand_pushforward(samples[['mu_gdd[1]']], 50, flim= c(-0.25,0.25))
mu <- rnorm(1e6, 0, log(1.8) / 2.57)
lines(density(mu), col = 'darkblue')
util$plot_expectand_pushforward(samples[['tau_gdd[1]']], 50, flim= c(0,0.5))
tau <- rnorm(1e6, 0, log(1.8^0.25) / 2.57)
lines(density(tau), col = 'darkblue')

util$plot_pairs_by_chain(samples[['mu_pre[1]']], 'mu_pre[1]',
                         log(samples[['tau_pre[1]']]), 'log(tau_pre)')
util$plot_pairs_by_chain(samples[['beta_pre[6]']], 'beta_pre[6]',
                         log(samples[['tau_pre[1]']]), 'log(tau_pre)')

par(mfrow = c(5,2), mar = c(4,4,1,1))
for(i in 1:10){
  util$plot_expectand_pushforward(samples[[paste0('beta_pre[',i,']')]], 50, flim= c(-0.1,0.15))
}

par(mfrow = c(2,1))
util$plot_expectand_pushforward(samples[['mu_pre[1]']], 50, flim= c(-0.25,0.25))
mu <- rnorm(1e6, 0, log(1.8) / 2.57)
lines(density(mu), col = 'darkblue')
util$plot_expectand_pushforward(samples[['tau_pre[1]']], 50, flim= c(0,0.5))
tau <- rnorm(1e6, 0, log(1.8^0.25) / 2.57)
lines(density(tau), col = 'darkblue')

util$plot_pairs_by_chain(samples[['mu_vpd[1]']], 'mu_vpd[1]',
                         log(samples[['tau_vpd[1]']]), 'log(tau_vpd)')

par(mfrow = c(2,1))
util$plot_expectand_pushforward(samples[['mu_vpd[1]']], 50, flim= c(-0.25,0.25))
mu <- rnorm(1e6, 0, log(1.8) / 2.57)
lines(density(mu), col = 'darkblue')
util$plot_expectand_pushforward(samples[['tau_vpd[1]']], 50, flim= c(0,0.5))
tau <- rnorm(1e6, 0, log(1.8^0.25) / 2.57)
lines(density(tau), col = 'darkblue')



#### Long-term GPs ####
util$plot_pairs_by_chain(samples[['mu_rho[1]']], 'mu_rho[1]',
                         log(samples[['tau_rho[1]']]), 'log(tau_rho)')
util$plot_pairs_by_chain(samples[['rho_sp[8]']], 'rho_sp[8]',
                         log(samples[['tau_rho[1]']]), 'log(tau_rho)')

par(mfrow = c(2,1))
util$plot_expectand_pushforward(samples[['mu_rho[1]']], 50, flim= c(0,100), bquote(mu[rho]))
mu <- rlnorm(1e6, 3.55, 0.24)
lines(density(mu), col = 'darkblue')

util$plot_expectand_pushforward(samples[['tau_rho[1]']], 50, flim= c(0,50), bquote(tau[rho]))
tau <- rnorm(1e6, 0, 20 / 2.57)
lines(density(tau), col = 'darkblue')

par(mfrow = c(5,2), mar = c(4,4,1,1))
for(i in 1:10){
  util$plot_expectand_pushforward(samples[[paste0('rho_sp[',i,']')]], 50, flim= c(0,150))
}

#### Short-term GPs ####
util$plot_pairs_by_chain(samples[['rho_ind']], 'rho_ind',
                         samples[['gamma_ind']], 'gamma_ind')

util$plot_pairs_by_chain(samples[['rho_ind']], 'rho_ind',
                         samples[['mu_rho[1]']], 'mu_rho[1]')

par(mfrow = c(2,1))
util$plot_expectand_pushforward(samples[['rho_ind']], 50, flim= c(0,25), display_name = 'Short-term rho (tree)')
mu <- rlnorm(1e6, 0.80, 0.40)
lines(density(mu), col = 'darkblue')

util$plot_expectand_pushforward(samples[['gamma_ind']], 50, flim= c(0,1), display_name = 'Short-term gamma')
tau <- rnorm(1e6, 0, log(3) / 2.57)
lines(density(tau), col = 'darkblue')

util$plot_expectand_pushforward(samples[['rho_sh']], 50, flim= c(0,25), display_name = 'Short-term rho (stand)')
mu <- rlnorm(1e6, 1.7, 0.26)
lines(density(mu), col = 'darkblue')


#### GP vs GP
util$plot_pairs_by_chain(samples[['rho_ind']], 'rho_ind',
                         samples[['rho_sp[2]']], 'rho_sp[2]')
abline(a = 0, b = 1)


#### THE SHOKS ####

samples <- fit$draws(c('phi_sck0', 'tau_phi_sck', 'alpha_phi_sck',
                       'omega_conc_sck0', 'tau_omega_conc_sck', 'alpha_omega_conc_sck',
                       'mu_logdelta_omega_nonconc_sck', 'tau_logdelta_omega_nonconc_sck', 'logdelta_omega_nonconc_sck'))
names <- dimnames(samples)$variable
samples <- lapply(1:dim(samples)[3], 
                  function(k){t(matrix(samples[1:dim(samples)[1],1:dim(samples)[2],k], 
                                       nrow = dim(samples)[1], ncol = dim(samples)[2]))})
names(samples) <- names
util$check_all_expectand_diagnostics(samples)


util$plot_pairs_by_chain(boot::logit(samples[['phi_sck0']]), 'logit(phi_sck0)',
                         log(samples[['tau_phi_sck']]), 'log(tau_phi_sck)')
util$plot_pairs_by_chain(samples[['phi_sck0']], 'phi_sck0',
                         log(samples[['tau_phi_sck']]), 'log(tau_phi_sck)')
util$plot_pairs_by_chain(boot::logit(samples[['omega_conc_sck0']]), 'logit(omega_conc_sck0)',
                         log(samples[['tau_omega_conc_sck']]), 'log(tau_omega_conc_sck)')
util$plot_pairs_by_chain(samples[['mu_logdelta_omega_nonconc_sck']], 'mu_logdelta_omega_nonconc_sck',
                         log(samples[['tau_logdelta_omega_nonconc_sck']]), 'log(tau_logdelta_omega_nonconc_sck)')


par(mfrow = c(2,1))
util$plot_expectand_pushforward(samples[['phi_sck0']], 50, flim= c(0,1), 'phi_sck0')
mu <- rbeta(1e6, 2, 20)
lines(density(mu), col = 'darkblue')

util$plot_expectand_pushforward(samples[['tau_phi_sck']], 50, flim= c(0,5), 'tau_phi_sck')
tau <- rnorm(1e6, 0, 1)
lines(density(tau), col = 'darkblue')

par(mfrow = c(5,2), mar = c(4,4,1,1))
for(i in 1:44){
  util$plot_expectand_pushforward(boot::inv.logit(samples[[paste0('alpha_phi_sck[',i,']')]]), 30, flim= c(0,1), display_name = paste0('alpha_phi_sck[',i,']'))
}

util$plot_pairs_by_chain(samples[['alpha_phi_sck[26]']], 'alpha_phi_sck[26]',
                         log(samples[['tau_phi_sck']]), 'log(tau_phi_sck)')


par(mfrow = c(2,1))
util$plot_expectand_pushforward(samples[['omega_conc_sck0']], 50, flim= c(0,1), 'omega_conc_sck0')
mu <- rbeta(1e6, 30, 20)
lines(density(mu), col = 'darkblue')

util$plot_expectand_pushforward(samples[['tau_omega_conc_sck']], 50, flim= c(0,5), 'tau_omega_conc_sck')
tau <- rnorm(1e6, 0, 1)
lines(density(tau), col = 'darkblue')

util$plot_pairs_by_chain(samples[['omega_conc_sck0']], 'omega_conc_sck0',
                         samples[['tau_omega_conc_sck']], 'tau_omega_conc_sck')

util$plot_pairs_by_chain(samples[['omega_conc_sck0']], 'omega_conc_sck0',
                         samples[['phi_sck0']], 'phi_sck0')

par(mfrow = c(5,2), mar = c(4,4,1,1))
for(i in 1:52){
  util$plot_expectand_pushforward(boot::inv.logit(samples[[paste0('alpha_omega_conc_sck[',i,']')]]), 30, flim= c(0,1),
                                  display_name = paste0('alpha_omega_conc_sck[',i,']'))
}

util$plot_pairs_by_chain(samples[['omega_conc_sck0']], 'omega_conc_sck0',
                         samples[['tau_omega_conc_sck']], 'tau_omega_conc_sck')

util$plot_pairs_by_chain(boot::inv.logit(samples[['alpha_omega_conc_sck[1]']]), 'omega_conc_sck[1]',
                         log(samples[['tau_omega_conc_sck']]), 'log(tau_omega_conc_sck)')

util$plot_pairs_by_chain(samples[['alpha_omega_conc_sck[52]']], 'alpha_omega_conc_sck[52]',
                         samples[['alpha_phi_sck[44]']], 'alpha_phi_sck[44]')
abline(a = 0, b = 1)



par(mfrow = c(2,1))
util$plot_expectand_pushforward(samples[['mu_logdelta_omega_nonconc_sck']], 50, flim= c(-5,5), 'mu_logdelta_omega_nonconc_sck')
mu <- rnorm(1e6, log(8), log(2)/2.57)
lines(density(mu), col = 'darkblue')
mu <- rnorm(1e6, log(5), log(3)/2.57)
lines(density(mu), col = 'darkblue', lty = 2)

util$plot_expectand_pushforward(samples[['tau_logdelta_omega_nonconc_sck']], 50, flim= c(0,5), 'tau_logdelta_omega_nonconc_sck')
tau <- rnorm(1e6, 0, 0.3/2.57)
lines(density(tau), col = 'darkblue')
tau <- rnorm(1e6, 0, 1/2.57)
lines(density(tau), col = 'darkblue', lty = 2)

util$plot_expectand_pushforward(samples[['logdelta_omega_nonconc_sck[49]']], 50, flim= c(0,5), 'logdelta_omega_nonconc_sck[49]')



samples <- fit$draws(c('omega_nonconc_sck', 'mu_logdelta_omega_nonconc_sck', 'tau_logdelta_omega_nonconc_sck'))
names <- dimnames(samples)$variable
samples <- lapply(1:dim(samples)[3], 
                  function(k){t(matrix(samples[1:dim(samples)[1],1:dim(samples)[2],k], 
                                       nrow = dim(samples)[1], ncol = dim(samples)[2]))})
names(samples) <- names
util$check_all_expectand_diagnostics(samples)
util$plot_disc_pushforward_quantiles(samples, paste0('omega_nonconc_sck[', 1:data$N_stand_species, ']'))

util$plot_pairs_by_chain(samples[['omega_nonconc_sck[33]']], 'omega_nonconc_sck[33]',
                         samples[['tau_logdelta_omega_nonconc_sck']], 'tau_logdelta_omega_nonconc_sck')

util$plot_pairs_by_chain(samples[['mu_logdelta_omega_nonconc_sck']], 'mu_logdelta_omega_nonconc_sck',
                         samples[['tau_logdelta_omega_nonconc_sck']], 'tau_logdelta_omega_nonconc_sck')


samples <- fit$draws(c('mu_logdelta_omega_nonconc_sck', 'tau_logdelta_omega_nonconc_sck', 'gamma_ind', 'phi_sck0'))
names <- dimnames(samples)$variable
samples <- lapply(1:dim(samples)[3], 
                  function(k){t(matrix(samples[1:dim(samples)[1],1:dim(samples)[2],k], 
                                       nrow = dim(samples)[1], ncol = dim(samples)[2]))})
names(samples) <- names

util$plot_pairs_by_chain(samples[['phi_sck0']], 'phi_sck0',
                         samples[['gamma_ind']], 'gamma_ind')




samples <- fit$draws(c('tau_alpha', 'tau_rho', 'tau_gamma',
                       'tau_kappa', 'tau_sck', 'tau_tau_sck', 'tau_phi_sck',
                       'tau_omega_conc_sck', 'tau_logdelta_omega_nonconc_sck', 'sigma'))
names <- dimnames(samples)$variable
samples <- lapply(1:dim(samples)[3], 
                  function(k){t(matrix(samples[1:dim(samples)[1],1:dim(samples)[2],k], 
                                       nrow = dim(samples)[1], ncol = dim(samples)[2]))})
names(samples) <- names

par(mfrow = c(3,3))
util$plot_expectand_pushforward(samples[['tau_alpha[1]']], 50, flim= c(0,3), 'tau_alpha')
tau <- rnorm(1e6, 0, log(5^0.25)/2.32)
lines(density(tau), col = 'darkblue')
tau <- rnorm(1e6, 0, log(5)/2.32) # what I should put
lines(density(tau), col = 'darkgreen', lty = 2)

util$plot_expectand_pushforward(samples[['tau_rho[1]']], 50, flim= c(0,50), 'tau_rho')
tau <- rnorm(1e6, 0, 20 / 2.57)
lines(density(tau), col = 'darkblue')

util$plot_expectand_pushforward(samples[['tau_gamma[1]']], 50, flim= c(0,1), 'tau_gamma')
tau <- rnorm(1e6, 0, 0.23 / 2.57)
lines(density(tau), col = 'darkblue')

util$plot_expectand_pushforward(samples[['tau_kappa[1]']], 50, flim= c(0,1), 'tau_kappa')
tau <- rnorm(1e6, 0, 0.2 / 2.57)
lines(density(tau), col = 'darkblue')

util$plot_expectand_pushforward(samples[['tau_tau_sck[1]']], 50, flim= c(0,2), 'tau_tau_sck')
tau <- rnorm(1e6, 0, log(2) / 2.57)
lines(density(tau), col = 'darkblue')

util$plot_expectand_pushforward(samples[['tau_phi_sck']], 50, flim= c(0,5), 'tau_phi_sck')
tau <- rnorm(1e6, 0, 1)
lines(density(tau), col = 'darkblue')

util$plot_expectand_pushforward(samples[['tau_omega_conc_sck']], 50, flim= c(0,5), 'tau_omega_conc_sck')
tau <- rnorm(1e6, 0, 1)
lines(density(tau), col = 'darkblue')

util$plot_expectand_pushforward(samples[['tau_logdelta_omega_nonconc_sck']], 50, flim= c(0,2), 'tau_logdelta_omega_nonconc_sck')
tau <- rnorm(1e6, 0, 0.3/2.57)
lines(density(tau), col = 'darkblue')
tau <- rnorm(1e6, 0, 1/2.57)
lines(density(tau), col = 'darkgreen', lty = 2)

util$plot_expectand_pushforward(samples[['sigma']], 50, flim= c(0,0.3), 'sigma')
tau <- rnorm(1e6, 0, log(1.1) / 2.57)
lines(density(tau), col = 'darkblue')



util$plot_pairs_by_chain(samples[['tau_phi_sck']], 'tau_phi_sck',
                         samples[['tau_rho[1]']], 'tau_rho[1]')



samples <- fit$draws(c('phi_sck', 'omega_conc_sck', 'omega_nonconc_sck'))
names <- dimnames(samples)$variable
samples <- lapply(1:dim(samples)[3], 
                  function(k){t(matrix(samples[1:dim(samples)[1],1:dim(samples)[2],k], 
                                       nrow = dim(samples)[1], ncol = dim(samples)[2]))})
names(samples) <- names

par(mfrow= c(3,1))
util$plot_disc_pushforward_quantiles(samples, paste0('phi_sck[',1:data$N_stand_species,']'), ylab = 'phi_sck')
util$plot_disc_pushforward_quantiles(samples, paste0('omega_conc_sck[',1:data$N_stand_species,']'), ylab = 'omega_conc_sck')
util$plot_disc_pushforward_quantiles(samples, paste0('omega_nonconc_sck[',1:data$N_stand_species,']'), ylab = 'omega_nonconc_sck')
