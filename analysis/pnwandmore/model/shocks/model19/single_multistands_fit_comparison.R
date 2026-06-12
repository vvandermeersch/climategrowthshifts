rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(cmdstanr)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)


folder <- '/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model'

data_1s <- readRDS(file.path(folder,'data_18may2026_PIPO_1stand_19502008.rds'))
base_samples_1s <- readRDS(file.path(folder, 'model19', 'basesamples_model19_HSGP_PIPO_1stand_19502008_idio3mods.rds'))
data <- readRDS(file.path(folder,'data_18may2026_PIPO_10stands_19502013.rds'))
base_samples <- readRDS(file.path(folder, 'model19', 'basesamples_model19_HSGP_PIPO_10stands_19502013_idio3mods.rds'))


par(mfrow = c(1,2), mar = c(4,4,1,1))
util$plot_expectand_pushforward(base_samples[['mu_conc']], 30,
                                display_name = 'mu_conc',
                                flim = c(0,2))
util$plot_expectand_pushforward(base_samples_1s[['mu_conc']], 30,
                                display_name = 'mu_conc',
                                flim = c(0,2), col = util$c_mid_teal,
                                add = T)
prior <- rlnorm(1e6, 0.154, 0.546)
lines(density(prior), col = '#95955b', lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(base_samples[['tau_conc']], 30,
                                display_name = 'tau_conc',
                                flim = c(0,1.5))
util$plot_expectand_pushforward(base_samples_1s[['tau_conc']], 30,
                                display_name = 'tau_conc',
                                flim = c(0,1.5), col = util$c_mid_teal,
                                add = T)
prior <- rnorm(1e6, 0, 1.5 / 2.57)
lines(density(prior), col = '#95955b', lwd = 1.5, lty = 2)



par(mfrow = c(1,3))

util$plot_expectand_pushforward(base_samples[['phi_sck[1]']], 30,
                                display_name = 'phi_sck[1]',
                                flim = c(0,1), ylim = c(0,15))
util$plot_expectand_pushforward(base_samples_1s[['phi_sck[1]']], 30,
                                display_name = 'phi_sck[1]',
                                flim = c(0,1), col = util$c_mid_teal,
                                add = T)
prior <- rbeta(1e6, 2, 13)
lines(density(prior), col = '#95955b', lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(base_samples[['omega_conc_sck[1]']], 30,
                                display_name = 'omega_conc_sck[1]',
                                flim = c(0,1))
util$plot_expectand_pushforward(base_samples_1s[['omega_conc_sck[1]']], 30,
                                display_name = 'omega_conc_sck[1]',
                                flim = c(0,1), col = util$c_mid_teal,
                                add = T)
prior <- rbeta(1e6, 4.3, 1.15)
lines(density(prior), col = '#95955b', lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(base_samples[['omega_shutdown[1]']], 30,
                                display_name = 'omega_shutdown[1]',
                                flim = c(0,1))
util$plot_expectand_pushforward(base_samples_1s[['omega_shutdown[1]']], 30,
                                display_name = 'omega_shutdown[1]',
                                flim = c(0,1), col = util$c_mid_teal,
                                add = T)
prior <- rbeta(1e6, 1.5, 4.3)
lines(density(prior), col = '#95955b', lwd = 1.5, lty = 2)


par(mfrow = c(1,3))
util$plot_disc_pushforward_quantiles(base_samples, paste0('thetas_idio[',which(data$stand_idxs == 1),',1]'), display_ylim = c(0,1),
                                     ylab = 'p(no idiosyncratic shock)')
util$plot_disc_pushforward_quantiles(base_samples, paste0('thetas_idio[',which(data$stand_idxs == 1),',2]'), display_ylim = c(0,1),
                                     ylab = 'p(idiosyncratic depres. growth)')
util$plot_disc_pushforward_quantiles(base_samples, paste0('thetas_idio[',which(data$stand_idxs == 1),',3]'), display_ylim = c(0,1),
                                     ylab = 'p(idiosyncratic shutdown)')



util$plot_disc_pushforward_quantiles(base_samples_1s, paste0('thetas_idio[',which(data$stand_idxs == 1),',1]'), display_ylim = c(0,1),
                                     ylab = 'p(no idiosyncratic shock)')
util$plot_disc_pushforward_quantiles(base_samples_1s, paste0('thetas_idio[',which(data$stand_idxs == 1),',2]'), display_ylim = c(0,1),
                                     ylab = 'p(idiosyncratic depres. growth)')
util$plot_disc_pushforward_quantiles(base_samples_1s, paste0('thetas_idio[',which(data$stand_idxs == 1),',3]'), display_ylim = c(0,1),
                                     ylab = 'p(idiosyncratic shutdown)')





par(mfrow = c(3,2))
util$plot_expectand_pushforward(base_samples[['rho_sp']], 30,
                                display_name = 'rho_sp',
                                flim = c(15,80))
util$plot_expectand_pushforward(base_samples_1s[['rho_sp']], 30,
                                display_name = 'rho_sp',
                                flim = c(15,80), col = util$c_mid_teal,
                                add = T)
prior <- rlnorm(1e6, 3.7, 0.35)
lines(density(prior), col = '#95955b', lwd = 1.5, lty = 2)

util$plot_expectand_pushforward(base_samples[['gamma_sp']], 30,
                                display_name = 'gamma_sp',
                                flim = c(0,2))
util$plot_expectand_pushforward(base_samples_1s[['gamma_sp']], 30,
                                display_name = 'gamma_sp',
                                flim = c(0,2), col = util$c_mid_teal,
                                add = T)
prior <- rnorm(1e6, 0, log(10) / 2.57)
lines(density(prior), col = '#95955b', lwd = 1.5, lty = 2)

util$plot_expectand_pushforward(base_samples[['rho_sh']], 30,
                                display_name = 'rho_sh',
                                flim = c(0,5))
util$plot_expectand_pushforward(base_samples_1s[['rho_sh']], 30,
                                display_name = 'rho_sh',
                                flim = c(0,5), col = util$c_mid_teal,
                                add = T)
prior <- rlnorm(1e6, 1.7, 0.26)
lines(density(prior), col = '#95955b', lwd = 1.5, lty = 2)

util$plot_expectand_pushforward(base_samples[['gamma_sh']], 30,
                                display_name = 'gamma_sh',
                                flim = c(0,1))
util$plot_expectand_pushforward(base_samples_1s[['gamma_sh']], 30,
                                display_name = 'gamma_sh',
                                flim = c(0,1), col = util$c_mid_teal,
                                add = T)
prior <- rnorm(1e6, 0, log(3) / 2.57)
lines(density(prior), col = '#95955b', lwd = 1.5, lty = 2)

util$plot_expectand_pushforward(base_samples[['rho_ind']], 30,
                                display_name = 'rho_ind',
                                flim = c(1,10))
util$plot_expectand_pushforward(base_samples_1s[['rho_ind']], 30,
                                display_name = 'rho_ind',
                                flim = c(1,10), col = util$c_mid_teal,
                                add = T)
prior <- rlnorm(1e6, 1.4, 0.25)
lines(density(prior), col = '#95955b', lwd = 1.5, lty = 2)

util$plot_expectand_pushforward(base_samples[['gamma_ind']], 30,
                                display_name = 'gamma_ind',
                                flim = c(0,1))
util$plot_expectand_pushforward(base_samples_1s[['gamma_ind']], 30,
                                display_name = 'gamma_ind',
                                flim = c(0,1), col = util$c_mid_teal,
                                add = T)
prior <- rnorm(1e6, 0, log(3) / 2.57)
lines(density(prior), col = '#95955b', lwd = 1.5, lty = 2)




# Stand-level GPs
fit <- readRDS(file.path(folder, 'model19', 'fit_model19_HSGP_PIPO_10stands_19502013_idio3mods.rds'))
samples <- fit$draws(c('f_sh'))
names <- dimnames(samples)$variable
samples <- lapply(1:dim(samples)[3],
                  function(k){t(matrix(samples[1:dim(samples)[1],1:dim(samples)[2],k],
                                       nrow = dim(samples)[1], ncol = dim(samples)[2]))})
names(samples) <- names
stand_samples <- samples

fit_1s <- readRDS(file.path(folder, 'model19', 'fit_model19_HSGP_PIPO_1stand_19502008_idio3mods.rds'))
samples <- fit_1s$draws(c('f_sh'))
names <- dimnames(samples)$variable
samples <- lapply(1:dim(samples)[3],
                  function(k){t(matrix(samples[1:dim(samples)[1],1:dim(samples)[2],k],
                                       nrow = dim(samples)[1], ncol = dim(samples)[2]))})
names(samples) <- names
stand_samples_1s <- samples


# Trees!
tree_samples <- readRDS(file.path(folder, 'model19', 'treesamples_model19_HSGP_PIPO_10stands_19502013_idio3mods.rds'))
tree_samples_1s <- readRDS(file.path(folder, 'model19', 'treesamples_model19_HSGP_PIPO_1stand_19502008_idio3mods.rds'))

chains <- c(1)

par(mfrow=c(3,3))
t <- 18
idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
names <- paste0("log_rw_pred[",idxs,"]")
temp <- lapply(tree_samples[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Log ring width (per mm)", 
                                     display_ylim=c(-3, 2), display_xlim = range(data$all_years))
points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=1, col="white")
points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=0.5, col="black")
text(x = 1950, y = 1.5, labels = data$uniq_tree_ids[t], cex = 0.8, adj = 0, col = "#474747")
sdw <- which(data$rw_obs[idxs] == 0)
for(p in sdw){
  points(x = data$years[idxs][p], y = log(mean(data$rw_obs[idxs][p-1], data$rw_obs[idxs][p+1])),
         pch = 4, cex = 1, col = 'black')
}

names <- paste0("f[",idxs,"]")
temp <- lapply(tree_samples[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Indivual GP (allometry)", 
                                     display_ylim=c(-3, 2), display_xlim = range(data$all_years),
                                     main = paste0('Chain ', chains))

names <- paste0("f_ind[",idxs,"]")
temp <- lapply(tree_samples[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Indivual GP (canopy)", 
                                     display_ylim=c(-3, 2), display_xlim = range(data$all_years))

names <- paste0("f_sh[",data$stand_idxs[t],",",which(data$all_years %in% data$years[idxs]),"]")
temp <- lapply(stand_samples[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Stand-level GP",
                                     display_ylim=c(-3, 2))

names <- paste0("delta_sck[",idxs,"]")
temp <- lapply(tree_samples[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Shock amplitude", 
                                     display_ylim=c(-3, 2), display_xlim = range(data$all_years))


names <- paste0("shutdown[",idxs,"]")
temp <- lapply(tree_samples[names], function(m) m[chains, , drop = FALSE])
util$plot_disc_pushforward_quantiles(temp, names, 
                                     xlab="Year", ylab="Shutdown state", 
                                     display_ylim=c(0, 1))

names <- paste0("conc_state[",idxs,"]")
temp <- lapply(tree_samples[names], function(m) m[chains, , drop = FALSE])
util$plot_disc_pushforward_quantiles(temp, names, 
                                     xlab="Year", ylab="Concordant state", 
                                     display_ylim=c(0, 1))

names <- paste0("idio_state[",idxs,"]")
temp <- lapply(tree_samples[names], function(m) m[chains, , drop = FALSE])
util$plot_disc_pushforward_quantiles(temp, names, 
                                     xlab="Year", ylab="Idiosyncratic state", 
                                     display_ylim=c(0, 1))

chains <- 1:4

par(mfrow=c(3,3), cex.axis = 1, cex.lab = 2)
t <- 18
idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
names <- paste0("log_rw_pred[",idxs,"]")
temp <- lapply(tree_samples_1s[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Log ring width (per mm)", 
                                     display_ylim=c(-3, 2), display_xlim = range(data$all_years))
points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=1, col="white")
points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=0.5, col="black")
text(x = 1950, y = 1.5, labels = data$uniq_tree_ids[t], cex = 0.8, adj = 0, col = "#474747")
sdw <- which(data$rw_obs[idxs] == 0)
for(p in sdw){
  points(x = data$years[idxs][p], y = log(mean(data$rw_obs[idxs][p-1], data$rw_obs[idxs][p+1])),
         pch = 4, cex = 1, col = 'black')
}

names <- paste0("f[",idxs,"]")
temp <- lapply(tree_samples_1s[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Indivual GP (allometry)", 
                                     display_ylim=c(-3, 2), display_xlim = range(data$all_years),
                                     main = paste('Chain', chains, collapse = ' - '))

names <- paste0("f_ind[",idxs,"]")
temp <- lapply(tree_samples_1s[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Indivual GP (canopy)", 
                                     display_ylim=c(-3, 2), display_xlim = range(data$all_years))

names <- paste0("f_sh[",data$stand_idxs[t],",",which(data$all_years %in% data$years[idxs]),"]")
temp <- lapply(stand_samples_1s[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Stand-level GP",
                                     display_ylim=c(-3, 2))

names <- paste0("delta_sck[",idxs,"]")
temp <- lapply(tree_samples_1s[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Shock amplitude", 
                                     display_ylim=c(-3, 2), display_xlim = range(data$all_years))


names <- paste0("shutdown[",idxs,"]")
temp <- lapply(tree_samples_1s[names], function(m) m[chains, , drop = FALSE])
util$plot_disc_pushforward_quantiles(temp, names, 
                                     xlab="Year", ylab="Shutdown state", 
                                     display_ylim=c(0, 1))

names <- paste0("conc_state[",idxs,"]")
temp <- lapply(tree_samples_1s[names], function(m) m[chains, , drop = FALSE])
util$plot_disc_pushforward_quantiles(temp, names, 
                                     xlab="Year", ylab="Concordant state", 
                                     display_ylim=c(0, 1))

names <- paste0("idio_state[",idxs,"]")
temp <- lapply(tree_samples_1s[names], function(m) m[chains, , drop = FALSE])
util$plot_disc_pushforward_quantiles(temp, names, 
                                     xlab="Year", ylab="Idiosyncratic state", 
                                     display_ylim=c(0, 1))





