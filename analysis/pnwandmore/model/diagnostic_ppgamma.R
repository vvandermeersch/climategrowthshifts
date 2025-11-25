rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('functions/custom_functions.R', local=util)
setwd(wd)

data <- readRDS(file.path(wd, 'output/model', 'data_30oct2025_az_2species.rds'))
fit <- readRDS(file.path(wd, "model/shocks/output", "fit_30oct2025_az_2species_ppgamma_ncp.rds"));gc()

samples <- util$extract_expectand_vals(fit);gc()
diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)
base_samples <- util$filter_expectands(samples,
                                       c('alpha', 
                                         'beta_gdd', 'beta_sm', 'beta_vpd',
                                         'rho_sp', 'gamma_sp',
                                         'f_tilde_sh',
                                         'mu_kappa', 'tau_kappa',
                                         'rho_sh', 'kappa_sh',
                                         'delta_tilde_sck', 
                                         'inner_tau_sck', 'outer_tau_sck', 'gamma_sck',
                                         'mu_kappa_sck', 'tau_kappa_sck', 'kappa_sck',
                                         'sigma'),
                                       check_arrays=TRUE)
util$check_all_expectand_diagnostics(base_samples)

s <- 3
trees_s <- which(data$stand_idxs==s)
for(t in trees_s){
  trans <- list(1)
  names(trans) <-  paste0('tau_gamma_sck[',s,']')
  util$plot_div_pairs(paste0('alpha_tilde_gamma_sck[',t,']'), paste0('tau_gamma_sck[',s,']'),
                      samples, diagnostics, transforms = trans)
}

util$plot_disc_pushforward_quantiles(samples, paste0('mu_gamma_sck[',1:data$N_stands,']'))

par(mfrow = c(1,3), mar = c(4,4,2,2), cex.lab = 1.2)
layout(matrix(c(1, 2, 3), ncol = 3, byrow = TRUE), widths = c(2, 1, 2))
for(s in 1:3){
  trees_s <- which(data$stand_idxs==s)
  util$plot_disc_pushforward_quantiles(samples, paste0('gamma_sck[',trees_s,']'),
                                       main = paste0('Stand ', s), display_ylim = c(0.75,1))
}


s <- 2
trees_s <- which(data$stand_idxs==s)
par(mfrow = c(4,2), mar = c(4,4,2,2), cex.lab = 1.2)
for(t in sample(trees_s,8)){
  idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- sapply(idxs_t,
                  function(x) paste0('log_rw_pred[', x, ']'))
  util$plot_conn_pushforward_quantiles(
    samples, names, plot_xs = data$years[idxs_t],
    main = paste0(data$uniq_species_ids[data$species_idxs[t]], ' - ', t), ylab = 'Log ring width (mm)', display_ylim = c(-12,2))
  points(data$years[idxs_t], data$log_rw_obs[idxs_t], pch=16, cex=1.0, col="white")
  points(data$years[idxs_t], data$log_rw_obs[idxs_t], pch=16, cex=0.8, col="black")
}         

par(mfrow = c(1,2), mar = c(4,4,2,2), cex.lab = 1.2)
for(t in c(54, 41)){
  idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- sapply(idxs_t,
                  function(x) paste0('log_rw_pred[', x, ']'))
  util$plot_conn_pushforward_quantiles(
    samples, names, plot_xs = data$years[idxs_t], baseline_values = data$log_rw_obs[idxs_t],
    main = paste0('Stand ', s, ' - ', data$uniq_species_ids[data$species_idxs[t]], ' - tree ', t), ylab = 'Log ring width (mm)', display_ylim = c(-12,2))
}  

par(mfrow = c(1,2), mar = c(4,4,2,2), cex.lab = 1.2)
for(t in c(41, 54)){
  idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- sapply(idxs_t,
                  function(x) paste0('log_rw_pred[', x, ']'))
  util$plot_conn_pushforward_quantiles(
    samples, names, plot_xs = data$years[idxs_t], baseline_values = data$log_rw_obs[idxs_t],
    main = paste0('Stand ', s, ' - ', data$uniq_species_ids[data$species_idxs[t]], ' - tree ', t), ylab = 'Log ring width (mm)', 
    display_ylim = c(-4,2.5))
  text(x = 1984, y =2.5 , label = round(util$ensemble_mcmc_quantile_est(samples[[ paste0('gamma_sck[',t,']')]], c(0.5)),2))
}  


par(mfrow = c(3,2), mar = c(4,4,2,2), cex.lab = 1.2)
for(t in c(41, 54)){
  stand <- 1
  names <- sapply(data$all_years_idxs[idxs_t],
                  function(y) paste0('f_sh[', stand, ',', y, ']'))
  util$plot_realizations(samples, names, plot_xs = data$years[idxs_t], N_plots = 50,
                         main = "", ylab = 'f sh', display_ylim = c(-7,2))
  text(x = 1990, y = -4, label = paste0('median(gamma)=', round(util$ensemble_mcmc_quantile_est(samples[[ paste0('gamma_sck[',t,']')]], c(0.5)),2)))
} 
for(t in c(41, 54)){
  idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- sapply(data$all_years_idxs[idxs_t],
                  function(y) paste0('delta_sck[', t, ',', y, ']'))
  util$plot_realizations(samples, names, plot_xs = data$years[idxs_t], N_plots = 50,
                         main = "", ylab = 'delta shock', display_ylim = c(-7,2))
  text(x = 1990, y = -4, label = paste0('median(gamma)=', round(util$ensemble_mcmc_quantile_est(samples[[ paste0('gamma_sck[',t,']')]], c(0.5)),2)))
} 
for(t in c(41, 54)){
  idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- sapply(idxs_t,
                  function(x) paste0('log_rw_pred[', x, ']'))
  util$plot_conn_pushforward_quantiles(
    samples, names, plot_xs = data$years[idxs_t], 
    main = paste0('Stand ', s, ' - ', data$uniq_species_ids[data$species_idxs[t]], ' - tree ', t), ylab = 'Log ring width (mm)', 
    display_ylim = c(-7,2))
  points(data$years[idxs_t], data$log_rw_obs[idxs_t], pch=16, cex=1.0, col="white")
  points(data$years[idxs_t], data$log_rw_obs[idxs_t], pch=16, cex=0.8, col="black")
}  

