rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(cmdstanr)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)


folder <- '/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model'

data <- readRDS(file.path(folder,'data_22apr2026_PIED_utah_19202024.rds'))
datasets <- readRDS(file.path(folder,'datasets_22apr2026_PSME_except2_19202024.rds'))

# fit <-  readRDS(file.path(folder, 'model19/fit_model19_HSGP_PIED_utah_pi0.rds'))

samples <- readRDS(file.path(folder, 'model19/basesamples_model19_HSGP_PIED_utah.rds'))
util$check_all_expectand_diagnostics(samples)

par(mfrow = c(2,2))
util$plot_expectand_pushforward(samples[['phi_sck0']], 30,
                                display_name = 'Probability of concordant state',
                                flim = c(0,1))
prior <- rbeta(1e6, 2, 13)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
plot.new()
util$plot_expectand_pushforward(samples[['omega_conc_sck0']], 30,
                                display_name = 'Probability of shock | concordant',
                                flim = c(0,1))
prior <- rbeta(1e6, 2, 4)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(samples[['omega_shutdown0']], 30,
                                display_name = 'Probability of shutdown | shock',
                                flim = c(0,1))
prior <- rbeta(1e6, 1, 25)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)

util$plot_expectand_pushforward(samples[['sigma']], 30,
                                display_name = 'sigma',
                                flim = c(0,log(1.3)))

util$plot_expectand_pushforward(samples[['mu_conc']], 30,
                                display_name = 'mu_conc',
                                flim = c(0,3))

util$plot_pairs_by_chain(samples[['phi_sck0']], 'phi_sck0',
                         samples[['omega_conc_sck0']], 'omega_conc_sck0')

util$plot_pairs_by_chain(samples[['tau_conc']], 'tau_conc',
                         samples[['tau_idio_small']], 'tau_small')

util$plot_pairs_by_chain(samples[['tau_conc']], 'tau_conc',
                         samples[['phi_sck0']], 'phi_sck0')

util$plot_pairs_by_chain(samples[['tau_conc']], 'tau_conc',
                         samples[['phi_sck0']], 'phi_sck0')


util$plot_pairs_by_chain(samples[['tau_idio_small']], 'tau_idio_small',
                         samples[['sigma']], 'sigma')


util$plot_disc_pushforward_quantiles(samples, paste0('alpha_omega_shutdown[',1:data$N_stands,']'))


util$plot_pairs_by_chain(samples[[paste0('thetas_idio[',1,',2]')]], paste0('thetas_idio[',1,',2]'),
                         samples[['tau_idio_small']], 'tau_idio_small')


par(mfrow = c(2,2))
util$plot_disc_pushforward_quantiles(samples, paste0('thetas_idio[',1:data$N_trees,',1]'), display_ylim = c(0,1))
util$plot_disc_pushforward_quantiles(samples, paste0('thetas_idio[',1:data$N_trees,',2]'), display_ylim = c(0,1))
util$plot_disc_pushforward_quantiles(samples, paste0('thetas_idio[',1:data$N_trees,',3]'), display_ylim = c(0,1))
util$plot_disc_pushforward_quantiles(samples, paste0('thetas_idio[',1:data$N_trees,',4]'), display_ylim = c(0,1))


fsh_samples <- fit$draws(c('f_sh'))
names <- dimnames(fsh_samples)$variable
fsh_samples <- lapply(1:dim(fsh_samples)[3],
                  function(k){t(matrix(fsh_samples[1:dim(fsh_samples)[1],1:dim(fsh_samples)[2],k],
                                       nrow = dim(fsh_samples)[1], ncol = dim(fsh_samples)[2]))})
names(fsh_samples) <- names


treesamples <- readRDS(file.path(folder, 'model19/treesamples_model19_HSGP_PIED_utah_pi0.rds'))

par(mfcol = c(5,2), cex.lab = 0.8, cex.axis = 0.8, mar = c(1.5,4,2.5,1), mgp = c(1.5, 0.2, 0),tck = -0.03)
t <- 1

chains <- 3
idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
names <- paste0("log_rw_pred[",idxs,"]")
temp <- lapply(treesamples[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Log ring width (per mm)", 
                                     display_ylim=c(-5, 2), display_xlim = range(data$all_years),
                                     main = paste0('Chain ', chains, collapse = ' '))
points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=1, col="white")
points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=0.5, col="black")
text(x = 1920, y = 1.5, labels = data$uniq_tree_ids[t], cex = 0.8, adj = 0, col = "#474747")
abline(v = data$years[idxs[which(data$rw_obs[idxs] == 0)]], col = 'grey80', lty = 2)

names <- paste0("f[",idxs,"]")
temp <- lapply(treesamples[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Short-term GP", 
                                     display_ylim=c(-5, 2), display_xlim = range(data$all_years))

names <- paste0("f_ind[",idxs,"]")
temp <- lapply(treesamples[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Short-term GP", 
                                     display_ylim=c(-5, 2), display_xlim = range(data$all_years))

names <- paste0("delta_sck[",idxs,"]")
temp <- lapply(treesamples[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Shock amplitude", 
                                     display_ylim=c(-5, 2), display_xlim = range(data$all_years))
abline(v = data$years[idxs[which(data$rw_obs[idxs] == 0)]], col = 'grey80', lty = 2)

stand <- data$stand_idxs[t]
stand_idxs <- data$all_years_idxs[idxs]
names <- paste0("f_sh[",stand, ',', stand_idxs,"]")
temp <- lapply(fsh_samples[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Stand-level GP", 
                                     display_ylim=c(-5, 2), display_xlim = range(data$all_years))
abline(v = data$years[idxs[which(data$rw_obs[idxs] == 0)]], col = 'grey80', lty = 2)


chains <- c(1,2,4)
idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
names <- paste0("log_rw_pred[",idxs,"]")
temp <- lapply(treesamples[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Log ring width (per mm)", 
                                     display_ylim=c(-5, 2), display_xlim = range(data$all_years),
                                     main = paste0('Chain ', chains, collapse = ' '))
points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=1, col="white")
points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=0.5, col="black")
text(x = 1920, y = 1.5, labels = data$uniq_tree_ids[t], cex = 0.8, adj = 0, col = "#474747")
abline(v = data$years[idxs[which(data$rw_obs[idxs] == 0)]], col = 'grey80', lty = 2)

names <- paste0("f[",idxs,"]")
temp <- lapply(treesamples[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Short-term GP", 
                                     display_ylim=c(-5, 2), display_xlim = range(data$all_years))

names <- paste0("f_ind[",idxs,"]")
temp <- lapply(treesamples[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Short-term GP", 
                                     display_ylim=c(-5, 2), display_xlim = range(data$all_years))

names <- paste0("delta_sck[",idxs,"]")
temp <- lapply(treesamples[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Shock amplitude", 
                                     display_ylim=c(-5, 2), display_xlim = range(data$all_years))
abline(v = data$years[idxs[which(data$rw_obs[idxs] == 0)]], col = 'grey80', lty = 2)

stand <- data$stand_idxs[t]
stand_idxs <- data$all_years_idxs[idxs]
names <- paste0("f_sh[",stand, ',', stand_idxs,"]")
temp <- lapply(fsh_samples[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Stand-level GP", 
                                     display_ylim=c(-5, 2), display_xlim = range(data$all_years))
abline(v = data$years[idxs[which(data$rw_obs[idxs] == 0)]], col = 'grey80', lty = 2)



