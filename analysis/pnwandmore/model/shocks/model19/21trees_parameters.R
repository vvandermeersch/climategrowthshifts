rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(cmdstanr)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)


folder <- '/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model'

data <- readRDS(file.path(folder,'data_18may2026_PIPO_10stands_19502013.rds'))

fit <- readRDS(file.path(folder, 'model19', 'fit_model19_HSGP_PIPO_10stands_19502013_fullidio.rds'))
base_samples <- fit$draws(c('alpha', 
                            'beta_gdd', 'beta_vpd', 'beta_pre',
                            
                            'rho_sp', 'gamma_sp',
                            'rho_sh', 'gamma_sh',
                            'rho_ind', 'gamma_ind',
                            
                            'mu_conc', 'tau_conc',
                            'phi_sck',
                            'omega_conc_sck',
                            'omega_shutdown',
                            
                            'tau_small_idio', 'tau_large_idio',
                            'thetas_idio',
                            
                            'sigma'
))
names <- dimnames(base_samples)$variable
base_samples <- lapply(1:dim(base_samples)[3],
                       function(k){t(matrix(base_samples[1:dim(base_samples)[1],1:dim(base_samples)[2],k],
                                            nrow = dim(base_samples)[1], ncol = dim(base_samples)[2]))})
names(base_samples) <- names
util$check_all_expectand_diagnostics(base_samples)


par(mfrow = c(3,2))
for(s in 1:data$N_stands){
  util$plot_expectand_pushforward(base_samples[[paste0('phi_sck[',s,']')]], 30,
                                  display_name = 'Probability of concordant state',
                                  flim = c(0,1), add = (s != 1))
}
prior <- rbeta(1e6, 2, 13)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
plot.new()
for(s in 1:data$N_stands){
  util$plot_expectand_pushforward(base_samples[[paste0('omega_conc_sck[',s,']')]], 30,
                                  display_name = 'Probability of shock | concordant',
                                  flim = c(0,1), add = (s != 1))
}
prior <- rbeta(1e6, 2, 2)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
for(s in 1:data$N_stands){
  util$plot_expectand_pushforward(base_samples[[paste0('omega_shutdown[',s,']')]], 30,
                                  display_name = 'Probability of shutdown | shock',
                                  flim = c(0,1), add = (s != 1))
}
prior <- rbeta(1e6, 2, 2)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(base_samples[['mu_conc']], 30,
                                display_name = 'mu_conc',
                                flim = c(0,log(30)))
prior <- rnorm(1e6, log(20), log(10)/2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(base_samples[['tau_conc']], 30,
                                display_name = 'tau_conc',
                                flim = c(0,2))
prior <- rnorm(1e6, 0, 2 / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)

par(mfrow = c(1,4))
util$plot_disc_pushforward_quantiles(base_samples, paste0('thetas_idio[',1:data$N_trees,',1]'), display_ylim = c(0,1))
util$plot_disc_pushforward_quantiles(base_samples, paste0('thetas_idio[',1:data$N_trees,',2]'), display_ylim = c(0,1))
util$plot_disc_pushforward_quantiles(base_samples, paste0('thetas_idio[',1:data$N_trees,',3]'), display_ylim = c(0,1))
util$plot_disc_pushforward_quantiles(base_samples, paste0('thetas_idio[',1:data$N_trees,',4]'), display_ylim = c(0,1))

par(mfrow=c(4,1))
util$plot_expectand_pushforward(base_samples[['sigma']], 30,
                                display_name = 'sigma',
                                flim = c(0,1))
util$plot_expectand_pushforward(base_samples[['tau_small_idio']], 30,
                                display_name = 'tau_idio',
                                flim = c(0,1))
prior <- rnorm(1e6, 0, log(1.5) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(base_samples[['tau_large_idio']], 30,
                                display_name = 'tau_idio',
                                flim = c(0,2))
prior <- rnorm(1e6, 0, log(5) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(base_samples[['tau_conc']], 30,
                                display_name = 'tau_conc',
                                flim = c(0,2))
prior <- rnorm(1e6, 0, 2 / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)

par(mfrow = c(3,2))

util$plot_expectand_pushforward(base_samples[['rho_sp']], 30,
                                display_name = 'rho_species',
                                flim = c(1, 100))
util$plot_expectand_pushforward(base_samples[['gamma_sp']], 30,
                                display_name = 'gamma_species',
                                flim = c(0, 1))

util$plot_expectand_pushforward(base_samples[['rho_sh']], 30,
                                display_name = 'rho_short',
                                flim = c(1, 10))
util$plot_expectand_pushforward(base_samples[['gamma_sh']], 30,
                                display_name = 'gamma_short',
                                flim = c(0, 1))

util$plot_expectand_pushforward(base_samples[['rho_ind']], 30,
                                display_name = 'rho_ind',
                                flim = c(1, 10))
util$plot_expectand_pushforward(base_samples[['gamma_ind']], 30,
                                display_name = 'gamma_ind',
                                flim = c(0, 1))

