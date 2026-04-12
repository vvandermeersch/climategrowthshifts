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
fit_with <- readRDS(file.path(folder, 'model17/fit_HSGP_11mar2026_idahoPIPO_model17_newpriors.rds'))
fit_without <- readRDS(file.path(folder, 'model17/fit_HSGP_11mar2026_idahoPIPO_model17_newpriors_pi0.rds'))
gc()


# Load samples
samples <- fit_with$draws(c('alpha', 'beta_gdd', 'beta_vpd', 'beta_pre',
                       
                       'rho_sp', 'gamma_sp',
                       'rho_sh', 'gamma_sh',
                       'rho_ind', 'gamma_ind',
                       
                       'tau_sck',
                       'phi_sck0', 'tau_phi_sck', 'alpha_phi_sck',
                       'omega_conc_sck0', 'tau_omega_conc_sck', 'alpha_omega_conc_sck',
                       
                       'pi_idsc_sck',
                       
                       'sigma'
))
names <- dimnames(samples)$variable
samples <- lapply(1:dim(samples)[3], 
                  function(k){t(matrix(samples[1:dim(samples)[1],1:dim(samples)[2],k], 
                                       nrow = dim(samples)[1], ncol = dim(samples)[2]))})
names(samples) <- names
samples_with <- samples

samples <- fit_without$draws(c('alpha', 'beta_gdd', 'beta_vpd', 'beta_pre',
                            
                            'rho_sp', 'gamma_sp',
                            'rho_sh', 'gamma_sh',
                            'rho_ind', 'gamma_ind',
                            
                            'tau_sck',
                            'phi_sck0', 'tau_phi_sck', 'alpha_phi_sck',
                            'omega_conc_sck0', 'tau_omega_conc_sck', 'alpha_omega_conc_sck',
                            
                            # 'pi_idsc_sck',
                            
                            'sigma'
))
names <- dimnames(samples)$variable
samples <- lapply(1:dim(samples)[3], 
                  function(k){t(matrix(samples[1:dim(samples)[1],1:dim(samples)[2],k], 
                                       nrow = dim(samples)[1], ncol = dim(samples)[2]))})
names(samples) <- names
samples_without <- samples





# Shocks
par(mfrow = c(4,2), mar = c(4,4,1,1))
util$plot_expectand_pushforward(samples_with[['tau_sck']], 50, flim= c(0,6), 'tau_sck')
util$plot_expectand_pushforward(samples_without[['tau_sck']], 50, flim= c(0,6), 'tau_sck',
                                add = TRUE, col = util$c_light)
prior <- rnorm(1e6, 0, 5 / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
plot.new()
util$plot_expectand_pushforward(samples_with[['phi_sck0']], 50, flim= c(0,1), 'phi_sck0')
util$plot_expectand_pushforward(samples_without[['phi_sck0']], 50, flim= c(0,1), 'tau_sck',
                                add = TRUE, col = util$c_light)
prior <- rbeta(1e6, 2, 20)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(samples_with[['tau_phi_sck']], 50, flim= c(0,5), 'tau_phi_sck',
                                ylim = c(0,2))
util$plot_expectand_pushforward(samples_without[['tau_phi_sck']], 50, flim= c(0,5), 'tau_phi_sck',
                                add = TRUE, col = util$c_light)
prior <- rnorm(1e6, 0, 1)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(samples_with[['omega_conc_sck0']], 50, flim= c(0,1), 'omega_conc_sck0',
                                ylim = c(0,30))
util$plot_expectand_pushforward(samples_without[['omega_conc_sck0']], 50, flim= c(0,1), 'omega_conc_sck0',
                                add = TRUE, col = util$c_light)
prior <- rbeta(1e6, 12.28, 34)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(samples_with[['tau_omega_conc_sck']], 50, flim= c(0,5), 'tau_omega_conc_sck')
util$plot_expectand_pushforward(samples_without[['tau_omega_conc_sck']], 50, flim= c(0,5), 'tau_omega_conc_sck',
                                add = TRUE, col = util$c_light)
prior <- rnorm(1e6, 0, 1)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(samples_with[['pi_idsc_sck']], 50, flim= c(0,1), 'pi_idsc_sck')
prior <- rbeta(1e6, 2, 100)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)



shocks <- fit_with$draws(c('phi_sck', 'omega_conc_sck'))
names <- dimnames(shocks)$variable
shocks <- lapply(1:dim(shocks)[3], 
                 function(k){t(matrix(shocks[1:dim(shocks)[1],1:dim(shocks)[2],k], 
                                      nrow = dim(shocks)[1], ncol = dim(shocks)[2]))})
names(shocks) <- names
shocks_with <- shocks
shocks <- fit_without$draws(c('phi_sck', 'omega_conc_sck'))
names <- dimnames(shocks)$variable
shocks <- lapply(1:dim(shocks)[3], 
                 function(k){t(matrix(shocks[1:dim(shocks)[1],1:dim(shocks)[2],k], 
                                      nrow = dim(shocks)[1], ncol = dim(shocks)[2]))})
names(shocks) <- names
shocks_without <- shocks
par(mfcol = c(2,2), mar = c(4,4,1,1), cex.main = 0.9)
util$plot_disc_pushforward_quantiles(shocks_with, paste0('phi_sck[', 1:data$N_stands, ']'), display_ylim = c(0,0.7),
                                     ylab = 'p(stand conc.)', main = 'With idiosyncratic shocks')
util$plot_disc_pushforward_quantiles(shocks_with, paste0('omega_conc_sck[', 1:data$N_stands, ']'), display_ylim = c(0,1),
                                     ylab = 'p(tree sck | stand conc.)')
util$plot_disc_pushforward_quantiles(shocks_without, paste0('phi_sck[', 1:data$N_stands, ']'), display_ylim = c(0,0.7),
                                     ylab = 'p(stand conc.)', main = 'Without idiosyncratic shocks')
util$plot_disc_pushforward_quantiles(shocks_without, paste0('omega_conc_sck[', 1:data$N_stands, ']'), display_ylim = c(0,1),
                                     ylab = 'p(tree sck | stand conc.)')

