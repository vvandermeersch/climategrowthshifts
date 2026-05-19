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

par(mfrow = c(3,7), mar = c(2,3,1,1), mgp = c(2,0.4,0))
for(t in 1:data$N_trees){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  plot(x = data$years[idxs], y = data$rw_obs[idxs], type = 'l',
       bty = 'n', ylab = 'rw (mm)', xlab = '',
       ylim = c(0,4), xlim = range(data$all_years), col = 'grey30')
  sdw <- which(data$rw_obs[idxs] == 0)
  points(x = data$years[idxs][sdw], y = data$rw_obs[idxs][sdw],
         pch = 20, cex = 0.7, col = 'darkred')
  abline(v = data$years[idxs][sdw], lty = 2, col = 'darkred',
         lwd = 0.7)
}

par(mfrow = c(3,7), mar = c(2,3,1,1), mgp = c(2,0.4,0))
for(t in 1:data$N_trees){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  plot(x = data$years[idxs], y = log(data$rw_obs[idxs]), type = 'l',
       bty = 'n', ylab = 'rw (mm)', xlab = '',
       ylim = c(-4,1.5), xlim = range(data$all_years), col = 'grey30')
  sdw <- which(data$rw_obs[idxs] == 0)
  # points(x = data$years[idxs][sdw], y = data$rw_obs[idxs][sdw],
  #        pch = 20, cex = 0.7, col = 'darkred')
  abline(v = c(1956, 2002, 2007), lty = 2, col = 'darkred',
         lwd = 0.7)
}

fit <- readRDS(file.path(folder, 'model18', 'fit_model18_HSGP_PIPO_1stand_19502008_v3.rds'))
base_samples <- fit$draws(c('alpha', 
                       'beta_gdd', 'beta_vpd', 'beta_pre',
                       
                       'rho_sp', 'gamma_sp',
                       'rho_sh', 'gamma_sh',
                       'rho_ind', 'gamma_ind',
                       
                       'mu_sck', 'tau_sck',
                       'phi_sck',
                       'omega_conc_sck',
                       'omega_shutdown',
                       
                       'sigma'
))
names <- dimnames(base_samples)$variable
base_samples <- lapply(1:dim(base_samples)[3],
                  function(k){t(matrix(base_samples[1:dim(base_samples)[1],1:dim(base_samples)[2],k],
                                       nrow = dim(base_samples)[1], ncol = dim(base_samples)[2]))})
names(base_samples) <- names
util$check_all_expectand_diagnostics(base_samples)


par(mfrow = c(3,2))
util$plot_expectand_pushforward(base_samples[['phi_sck[1]']], 30,
                                display_name = 'Probability of concordant state',
                                flim = c(0,1))
prior <- rbeta(1e6, 2, 13)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
plot.new()
util$plot_expectand_pushforward(base_samples[['omega_conc_sck[1]']], 30,
                                display_name = 'Probability of shock | concordant',
                                flim = c(0,1))
prior <- rbeta(1e6, 2, 2)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(base_samples[['omega_shutdown[1]']], 30,
                                display_name = 'Probability of shutdown | shock',
                                flim = c(0,1))
prior <- rbeta(1e6, 1, 10)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(base_samples[['mu_sck']], 30,
                                display_name = 'mu_shock',
                                flim = c(0,3))
prior <- rnorm(1e6, log(20), log(10)/2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(base_samples[['tau_sck']], 30,
                                display_name = 'tau_sck',
                                flim = c(0,2))
prior <- rnorm(1e6, 0, 2 / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)

par(mfrow = c(1,1))
util$plot_expectand_pushforward(base_samples[['sigma']], 30,
                                display_name = 'sigma',
                                flim = c(0,1))
prior <- rnorm(1e6, 0, log(1.5) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)

par(mfrow = c(3,2))
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

util$plot_expectand_pushforward(base_samples[['omega_conc_sck[1]']], 30,
                                display_name = 'Probability of shock | concordant',
                                flim = c(0,1))
prior <- rbeta(1e6, 3, 12)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(base_samples[['omega_shutdown[1]']], 30,
                                display_name = 'Probability of shutdown | shock',
                                flim = c(0,1))
prior <- rbeta(1e6, 1, 25)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(base_samples[['mu_sck']], 30,
                                display_name = 'mu_shock',
                                flim = c(0,3))
prior <- rnorm(1e6, log(11), log(7)/2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(base_samples[['tau_sck']], 30,
                                display_name = 'tau_sck',
                                flim = c(0,2))
prior <- rnorm(1e6, 0, 2 / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)

par(mfrow = c(1,1))
util$plot_expectand_pushforward(base_samples[['sigma']], 30,
                                display_name = 'sigma',
                                flim = c(0,1))
prior <- rnorm(1e6, 0, log(1.5) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)

fsh_samples <- fit$draws(c('f_sh'))
names <- dimnames(fsh_samples)$variable
fsh_samples <- lapply(1:dim(fsh_samples)[3],
                       function(k){t(matrix(fsh_samples[1:dim(fsh_samples)[1],1:dim(fsh_samples)[2],k],
                                            nrow = dim(fsh_samples)[1], ncol = dim(fsh_samples)[2]))})
names(fsh_samples) <- names
par(mfrow = c(1,1))
util$plot_conn_pushforward_quantiles(fsh_samples, paste0('f_sh[1,',1:data$N_all_years,']'), data$all_years)
abline(v = c(1956, 2002, 2007), lty = 2)


