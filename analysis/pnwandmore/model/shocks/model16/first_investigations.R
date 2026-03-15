rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('functions/custom_functions.R', local=util)
setwd(wd)

data <- readRDS(file.path(wd, 'output/model', 'data_11mar2026_oregonTSME_standclimate_18962024.rds'))
fit <- readRDS(file.path(wd, 'output/model/model16', "fit_11mar2026_oregonTSME_18962024_model16_indGP.rds")) 

samples <- fit$draws(c('rho_ind', 'gamma_ind', 
                       'rho_sp', 'gamma_sp',
                       'rho_sh'))


names <- dimnames(samples)$variable
samples <- lapply(1:dim(samples)[3], 
                  function(k){t(matrix(samples[1:dim(samples)[1],1:dim(samples)[2],k], 
                                       nrow = dim(samples)[1], ncol = dim(samples)[2]))})
names(samples) <- names
util$check_all_expectand_diagnostics(samples)

util$plot_expectand_pushforward(samples[['rho_sp[1]']], 30)
util$plot_expectand_pushforward(samples[['rho_ind']], 30)
util$plot_expectand_pushforward(samples[['gamma_sp[1]']], 30)
util$plot_expectand_pushforward(samples[['gamma_ind']], 30)
util$plot_expectand_pushforward(samples[['rho_sh']], 30)


samples <- fit$draws(
  c("phi_sck0", "tau_phi_sck", "alpha_phi_sck","phi_sck", 
    "omega_conc_sck0", "tau_omega_conc_sck", "alpha_omega_conc_sck", "omega_conc_sck",
    "mu_logdelta_omega_nonconc_sck", "tau_logdelta_omega_nonconc_sck", "logdelta_omega_nonconc_sck", "omega_nonconc_sck"))

names <- dimnames(samples)$variable
samples <- lapply(1:dim(samples)[3], 
                  function(k){t(matrix(samples[1:dim(samples)[1],1:dim(samples)[2],k], 
                                       nrow = dim(samples)[1], ncol = dim(samples)[2]))})
names(samples) <- names
util$check_all_expectand_diagnostics(samples)


names <- paste0('phi_sck[', 1:data$N_stands, ']')
util$plot_disc_pushforward_quantiles(samples, names)
names <- paste0('omega_nonconc_sck[', 1:data$N_stand_species, ']')
util$plot_disc_pushforward_quantiles(samples, names)
names <- paste0('omega_conc_sck[', 1:data$N_stand_species, ']')
util$plot_disc_pushforward_quantiles(samples, names)


# Generate predictions 
mod_gq <- cmdstan_model(file.path(wd, 'model/stan/model16', 'model16_indGP_withGQ.stan'))
data_gq <- data
data_gq$grainsize <-  ceiling(data_gq$N_stands/8)
data_gq$uniq_tree_ids <- NULL
data_gq$uniq_species_ids <- NULL
data_gq$uniq_stand_ids <- NULL
data_gq$N_clades <- 1
data_gq$uniq_stand_lat <- NULL
data_gq$uniq_stand_lon <- NULL
data_gq$tree_pred_idxs <- 1:data$N_trees # all trees
data_gq$N_pred <- sum(sapply(data_gq$tree_pred_idxs , function(t) length(data$tree_start_idxs[t]:data$tree_end_idxs[t])))
data_gq$N_trees_pred <- length(data_gq$tree_pred_idxs )
fit_gq <- mod_gq$generate_quantities(fit, data = data_gq, seed = 5838293, parallel_chains = 4)
gc()
log_rw_pred <- fit_gq$draws('log_rw_pred')
names <- dimnames(log_rw_pred)$variable
log_rw_pred <- lapply(1:dim(log_rw_pred)[3], 
                  function(k){t(matrix(log_rw_pred[1:dim(log_rw_pred)[1],1:dim(log_rw_pred)[2],k], 
                                       nrow = dim(log_rw_pred)[1], ncol = dim(log_rw_pred)[2]))})
names(log_rw_pred) <- names
gc()
delta_sck <- fit_gq$draws('delta_sck')
names <- dimnames(delta_sck)$variable
delta_sck <- lapply(1:dim(delta_sck)[3], 
                      function(k){t(matrix(delta_sck[1:dim(delta_sck)[1],1:dim(delta_sck)[2],k], 
                                           nrow = dim(delta_sck)[1], ncol = dim(delta_sck)[2]))})
names(delta_sck) <- names
gc()
gps <- fit_gq$draws(c('f', 'f_ind'))
names <- dimnames(gps)$variable
gps <- lapply(1:dim(gps)[3], 
                    function(k){t(matrix(gps[1:dim(gps)[1],1:dim(gps)[2],k], 
                                         nrow = dim(gps)[1], ncol = dim(gps)[2]))})
names(gps) <- names


pdf(file.path(wd, 'model/shocks/model16/figures', 'orTSME.pdf'), height = 10, width = 8)
par(mfrow = c(8,2), cex.lab = 0.8, cex.axis = 0.8, mar = c(1,4,1,1))
for(t in 1:data$N_trees){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- paste0("log_rw_pred[",idxs,"]")
  
  util$plot_conn_pushforward_quantiles(log_rw_pred, names, data$years[idxs],
                                       xlab="Year", ylab="Log ring width (per mm)", 
                                       display_ylim=c(-8, 2), display_xlim = range(data$all_years))
  points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=1, col="white")
  points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=0.5, col="black")
  text(x = .5, y = -6.9, labels = data$uniq_tree_ids[t], cex = 1, adj = 0)
  text(x = .5, y = -7.7, labels = paste0('Species:', data$uniq_species_ids[data$species_idxs[t]]), cex = 1, adj = 0)
  
  names <- paste0("delta_sck[",idxs,"]")
  util$plot_conn_pushforward_quantiles(delta_sck, names, data$years[idxs],
                                       xlab="Year", ylab="Shock amplitude", 
                                       display_ylim=c(-5, 1), display_xlim = range(data$all_years))
}
dev.off()


