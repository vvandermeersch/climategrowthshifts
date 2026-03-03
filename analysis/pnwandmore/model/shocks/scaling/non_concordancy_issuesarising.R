rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(terra)
library(rnaturalearth)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)

data <- readRDS(file.path(wd, 'output/model', 'data_15jan2025_PIPO_standclimate_19802024.rds'))
fit <- readRDS(file.path(wd, 'output/model', 'fit_15jan2025_PIPO_standclimate_19802024_standpoolingonly_nonconcordant_updatedpriors2.rds'))

param_samples <- fit$draws(variables = c("alpha",
                                         "beta_gdd", "beta_vpd", "beta_pre", 
                                         "tau_sck", "sigma",
                                         "phi_sck0", "mu_phi_sck", "tau_phi_sck", 
                                         "alpha_tilde_phi_sck", "phi_sck",
                                         
                                         "omega_conc_sck0", "mu_omega_conc_sck", "tau_omega_conc_sck", 
                                         "alpha_omega_conc_sck", "omega_conc_sck",
                                         
                                         "omega_nonconc_sck0", "mu_omega_nonconc_sck", "tau_omega_nonconc_sck", 
                                         "alpha_omega_nonconc_sck", "omega_nonconc_sck",
                                         "rho_sp", "rho_sh", "rho_sh", 'sigma'))

samples <- lapply(1:dim(param_samples)[3], function(k){t(matrix(param_samples[1:dim(param_samples)[1],1:dim(param_samples)[2],k], 
                                                                nrow = dim(param_samples)[1], ncol = dim(param_samples)[2]))})
names(samples) <- dimnames(param_samples)$variable
util$check_all_expectand_diagnostics(samples)

util$plot_pairs_by_chain(samples[['alpha_omega_nonconc_sck[30]']], 'alpha_omega_nonconc_sck[30]',
                         log(samples[['tau_omega_nonconc_sck']]), 'tau_omega_nonconc_sck')

util$plot_pairs_by_chain(samples[['tau_sck[1]']], 'tau_sck[1]',
                         samples[['omega_nonconc_sck0']], 'omega_nonconc_sck0')

util$plot_pairs_by_chain(samples[['phi_sck0']], 'phisck0',
                         samples[['omega_conc_sck0']], 'omega_conc_sck0')

util$plot_pairs_by_chain(samples[['tau_omega_nonconc_sck']], 'tau_omega_nonconc_sck',
                         samples[['tau_phi_sck']], 'tau_phi_sck')

util$plot_pairs_by_chain(, 'phisck0',
                         matrix(samples[['tau_phi_sck']][1,1:1000], ncol = 1000), 'tau_phi_sck')

par(mfrow = c(1,1))
util$plot_expectand_pushforward(matrix(samples[['tau_phi_sck']][2:4,1:1000], ncol = 1000), 50,
                                flim = c(0,5), ylim = c(0,5),
                                display_name = bquote(tau[phi[shock]]))
util$plot_expectand_pushforward(matrix(samples[['tau_phi_sck']][1,1:1000], ncol = 1000), 50,
                                flim = c(0,5), ylim = c(0,5),
                                display_name = bquote(tau[phi[shock]]),
                                col = 'orange', add = TRUE)
y <- rnorm(1e6,0,1/2.57)
lines(density(y), col = 'lightblue', lwd = 3)


par(mfrow = c(1,1))
util$plot_expectand_pushforward(matrix(samples[['mu_phi_sck']][2:4,1:1000], ncol = 1000), 50,
                                flim = c(-5,5), ylim = c(0,5),
                                display_name = bquote(mu[phi[shock]]))
util$plot_expectand_pushforward(matrix(samples[['mu_phi_sck']][1,1:1000], ncol = 1000), 50,
                                flim = c(-5,5), ylim = c(0,5),
                                display_name = bquote(mu[phi[shock]]),
                                col = 'orange', add = TRUE)

par(mfrow = c(1,1))
util$plot_expectand_pushforward(matrix(samples[['phi_sck0']][2:4,1:1000], ncol = 1000), 50,
                                flim = c(0,1), ylim = c(0,50),
                                display_name = bquote(mu[phi[shock]]))
util$plot_expectand_pushforward(matrix(samples[['phi_sck0']][1,1:1000], ncol = 1000), 50,
                                flim = c(0,1), ylim = c(0,50),
                                display_name = bquote(mu[phi[shock]]),
                                col = 'orange', add = TRUE)

par(mfrow = c(1,1))
util$plot_expectand_pushforward(matrix(samples[['omega_nonconc_sck0']][2:4,1:1000], ncol = 1000), 50,
                                flim = c(0,0.3), ylim = c(0,100),
                                display_name = bquote(mu[phi[shock]]))
util$plot_expectand_pushforward(matrix(samples[['omega_nonconc_sck0']][1,1:1000], ncol = 1000), 50,
                                flim = c(0,0.3), ylim = c(0,100),
                                display_name = bquote(mu[phi[shock]]),
                                col = 'orange', add = TRUE)

par(mfrow = c(1,1))
util$plot_expectand_pushforward(matrix(samples[['tau_omega_nonconc_sck']][2:4,1:1000], ncol = 1000), 50,
                                flim = c(0,5), ylim = c(0,5),
                                display_name = bquote(mu[phi[shock]]))
util$plot_expectand_pushforward(matrix(samples[['tau_omega_nonconc_sck']][1,1:1000], ncol = 1000), 50,
                                flim = c(0,5), ylim = c(0,5),
                                display_name = bquote(mu[phi[shock]]),
                                col = 'orange', add = TRUE)
y <- rnorm(1e6,0,1/2.57)
lines(density(y), col = 'lightblue', lwd = 3)


par(mfrow = c(1,1))
util$plot_expectand_pushforward(matrix(samples[['tau_sck[1]']][2:4,1:1000], ncol = 1000), 50,
                                flim = c(0,10), ylim = c(0,5),
                                display_name = bquote(mu[phi[shock]]))
util$plot_expectand_pushforward(matrix(samples[['tau_sck[1]']][1,1:1000], ncol = 1000), 50,
                                flim = c(0,10), ylim = c(0,5),
                                display_name = bquote(mu[phi[shock]]),
                                col = 'orange', add = TRUE)
y <- rnorm(1e6,0,1/2.57)
lines(density(y), col = 'lightblue', lwd = 3)





param_samples <- fit$draws(variables = c("alpha",
                                         "beta_gdd", "beta_vpd", "beta_pre", 
                                         "tau_sck", "sigma",
                                         "phi_sck0", "mu_phi_sck", "tau_phi_sck", 
                                         
                                         "omega_conc_sck0", "mu_omega_conc_sck", "tau_omega_conc_sck", 

                                         
                                         "omega_nonconc_sck0", "mu_omega_nonconc_sck", "tau_omega_nonconc_sck", 

                                         "rho_sp", "rho_sh", "rho_sh", 'sigma'))

samples <- lapply(1:dim(param_samples)[3], function(k){t(matrix(param_samples[1:dim(param_samples)[1],1:dim(param_samples)[2],k], 
                                                                nrow = dim(param_samples)[1], ncol = dim(param_samples)[2]))})
names(samples) <- dimnames(param_samples)$variable
util$check_all_expectand_diagnostics(samples)



qg_samples <-  fit$draws(variables = c('log_rw_pred', 'sck_state', 'delta_sck'))
names <- dimnames(qg_samples)$variable
qg_samples_chain1 <- lapply(1:dim(qg_samples)[3], function(k){t(matrix(qg_samples[1:dim(qg_samples)[1],1,k], 
                                                                nrow = dim(qg_samples)[1], ncol = 1))})
names(qg_samples_chain1) <- names

qg_samples_chains234 <- lapply(1:dim(qg_samples)[3], function(k){t(matrix(qg_samples[1:dim(qg_samples)[1],2:4,k], 
                                                                       nrow = dim(qg_samples)[1], ncol = 3))})
names(qg_samples_chains234) <- names
gc()

par(mfrow = c(2,4), cex.lab = 1.2)



t <- 444
idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
names <- paste0("log_rw_pred[",idxs,"]")
util$plot_conn_pushforward_quantiles(qg_samples_chain1, names, data$years[idxs]-1979,
                                     xlab="Year", ylab="Log ring width (per mm)", 
                                     display_ylim=c(-8, 2), display_xlim = range(data$years[idxs])-1979)
points(data$years[idxs]-1979, log(data$rw_obs[idxs]), pch=16, cex=1, col="white")
points(data$years[idxs]-1979, log(data$rw_obs[idxs]), pch=16, cex=0.5, col="black")

names <- paste0("sck_state[",idxs,"]")
util$plot_disc_pushforward_quantiles(qg_samples_chain1, names, data$years[idxs],
                                     xlab="Year", ylab="Shock state", 
                                     display_ylim=c(-0.05, 1.05))

names <- paste0("delta_sck[",idxs,"]")
util$plot_conn_pushforward_quantiles(qg_samples_chain1, names, data$years[idxs]-1979,
                                     xlab="Year", ylab="Shock amplitude", 
                                     display_ylim=c(-7, 1))

idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
names <- paste0("log_rw_pred[",idxs,"]")

util$plot_conn_pushforward_quantiles(qg_samples_chains234, names, data$years[idxs]-1979,
                                     xlab="Year", ylab="Log ring width (per mm)", 
                                     display_ylim=c(-8, 2), display_xlim = range(data$years[idxs])-1979)
points(data$years[idxs]-1979, log(data$rw_obs[idxs]), pch=16, cex=1, col="white")
points(data$years[idxs]-1979, log(data$rw_obs[idxs]), pch=16, cex=0.5, col="black")

names <- paste0("sck_state[",idxs,"]")
util$plot_disc_pushforward_quantiles(qg_samples_chains234, names, data$years[idxs],
                                     xlab="Year", ylab="Shock state", 
                                     display_ylim=c(-0.05, 1.05))

names <- paste0("delta_sck[",idxs,"]")
util$plot_conn_pushforward_quantiles(qg_samples_chains234, names, data$years[idxs]-1979,
                                     xlab="Year", ylab="Shock amplitude", 
                                     display_ylim=c(-7, 1))




quantile(boot::logit(rbeta(1e6, 1, 72)), c(0.025,0.975))
quantile(rnorm(1e6, -5.4,1.25), c(0.025,0.975))

quantile(rbeta(1e6, 1, 72), c(0.025,0.975))
quantile(boot::inv.logit(rnorm(1e6, -5.4,1.25)), c(0.025,0.975))



par(mfrow = c(1,1))
for(s in 1:data$N_stand_species){
  
  util$plot_expectand_pushforward(matrix(samples[['omega_nonconc_sck[1]']][2:4,1:1000], ncol = 1000), 50,
                                  flim = c(0,1), ylim = c(0,100),
                                  display_name = bquote(mu[phi[shock]]),
                                  col = util$c_dark)
  util$plot_expectand_pushforward(matrix(samples[['omega_nonconc_sck[1]']][1,1:1000], ncol = 1000), 50,
                                  flim = c(0,1), ylim = c(0,100),
                                  display_name = bquote(mu[phi[shock]]),
                                  col = util$c_dark_teal, add = TRUE)
  
  util$plot_expectand_pushforward(matrix(samples[['omega_conc_sck[1]']][2:4,1:1000], ncol = 1000), 50,
                                  flim = c(0,1), ylim = c(0,100),
                                  display_name = bquote(mu[phi[shock]]),
                                  col = util$c_light, add = TRUE)
  util$plot_expectand_pushforward(matrix(samples[['omega_conc_sck[1]']][1,1:1000], ncol = 1000), 50,
                                  flim = c(0,1), ylim = c(0,100),
                                  display_name = bquote(mu[phi[shock]]),
                                  col = util$c_light_teal, add = TRUE)
}



util$plot_expectand_pushforward(matrix(samples[['omega_nonconc_sck0']][2:4,1:1000], ncol = 1000), 50,
                                flim = c(0,1), ylim = c(0,100),
                                display_name = bquote(mu[phi[shock]]),
                                col = util$c_dark)
util$plot_expectand_pushforward(matrix(samples[['omega_nonconc_sck0']][1,1:1000], ncol = 1000), 50,
                                flim = c(0,1), ylim = c(0,100),
                                display_name = bquote(mu[phi[shock]]),
                                col = util$c_dark_teal, add = TRUE)

util$plot_expectand_pushforward(matrix(samples[['omega_conc_sck0']][2:4,1:1000], ncol = 1000), 50,
                                flim = c(0,1), ylim = c(0,100),
                                display_name = bquote(mu[phi[shock]]),
                                col = util$c_light, add = TRUE)
util$plot_expectand_pushforward(matrix(samples[['omega_conc_sck0']][1,1:1000], ncol = 1000), 50,
                                flim = c(0,1), ylim = c(0,100),
                                display_name = bquote(mu[phi[shock]]),
                                col = util$c_light_teal, add = TRUE)



par(mfrow = c(1,1))
for(s in 1:data$N_stand_species){
  
  util$plot_expectand_pushforward(matrix(samples[['alpha_omega_nonconc_sck[1]']][2:4,1:1000], ncol = 1000), 50,
                                  flim = c(-30,30), ylim = c(0,1),
                                  display_name = bquote(mu[phi[shock]]),
                                  col = util$c_dark)
  util$plot_expectand_pushforward(matrix(samples[['alpha_omega_nonconc_sck[1]']][1,1:1000], ncol = 1000), 50,
                                  flim = c(-30,30), ylim = c(0,1),
                                  display_name = bquote(mu[phi[shock]]),
                                  col = util$c_dark_teal, add = TRUE)
  
  util$plot_expectand_pushforward(matrix(samples[['alpha_omega_conc_sck[1]']][2:4,1:1000], ncol = 1000), 50,
                                  flim = c(-30,30), ylim = c(0,1),
                                  display_name = bquote(mu[phi[shock]]),
                                  col = util$c_light, add = TRUE)
  util$plot_expectand_pushforward(matrix(samples[['alpha_omega_conc_sck[1]']][1,1:1000], ncol = 1000), 50,
                                  flim = c(-30,30), ylim = c(0,1),
                                  display_name = bquote(mu[phi[shock]]),
                                  col = util$c_light_teal, add = TRUE)
}

