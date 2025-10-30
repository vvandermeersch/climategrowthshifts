rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('functions/custom_functions.R', local=util)
setwd(wd)

data <- readRDS(file.path(wd, 'output/model', 'data_30oct2025_az_2species.rds'))
fit <- readRDS(file.path(wd, "model/shocks/output", "fit_30oct2025_az_2species_standgamma.rds"));gc()

samples <- util$extract_expectand_vals(fit);gc()
diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)
base_samples <- util$filter_expectands(samples,
                                       c('alpha', 
                                         'mu_gdd', 'mu_sm', 'mu_vpd',
                                         'tau_gdd', 'tau_sm', 'tau_vpd',
                                         'beta_gdd', 'beta_sm', 'beta_vpd',
                                         'mu_rho', 'mu_gamma', 'tau_rho', 'tau_gamma', 
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

s <- 2
for(y in 1:data$N_all_years){
  util$plot_div_pairs(paste0('delta_sck[',s,',',y,']'), 'inner_tau_sck', 
                      samples, diagnostics, transforms = list('inner_tau_sck' = 1))
  title(data$all_years[y])
}
for(y in 1:data$N_all_years){
  util$plot_div_pairs(paste0('delta_sck[',s,',',y,']'), 'outer_tau_sck', 
                      samples, diagnostics, transforms = list('outer_tau_sck' = 1))
  title(data$all_years[y])
}
for(y in 1:data$N_all_years){
  util$plot_div_pairs(paste0('f_tilde_sh[',s,',',y,']'), 'rho_sh', 
                      samples, diagnostics)
}
for(s in 1:data$N_species){
  util$plot_div_pairs(paste0('beta_gdd[',s,']'), 'tau_gdd[1]', 
                      samples, diagnostics, transforms = list('tau_gdd[1]' = 1))
  util$plot_div_pairs(paste0('beta_vpd[',s,']'), 'tau_vpd[1]', 
                      samples, diagnostics, transforms = list('tau_vpd[1]' = 1))
  util$plot_div_pairs(paste0('beta_sm[',s,']'), 'tau_sm[1]', 
                      samples, diagnostics, transforms = list('tau_sm[1]' = 1))
  util$plot_div_pairs(paste0('rho_sp[',s,']'), 'tau_rho[1]', 
                      samples, diagnostics, transforms = list('tau_rho[1]' = 1))
  util$plot_div_pairs(paste0('gamma_sp[',s,']'), 'tau_gamma[1]', 
                      samples, diagnostics, transforms = list('tau_gamma[1]' = 1))
  util$plot_div_pairs(paste0('kappa_sh[',s,']'), 'tau_kappa[1]', 
                      samples, diagnostics, transforms = list('tau_kappa[1]' = 1))
  util$plot_div_pairs(paste0('kappa_sck[',s,']'), 'tau_kappa_sck[1]', 
                      samples, diagnostics, transforms = list('tau_kappa_sck[1]' = 1))
}

par(mfrow = c(1,1), mar = c(4,4,2,2), cex.main = 1)
util$plot_expectand_pushforward(samples[['rho_sh']], 15,
                                flim = c(0,15),
                                display_name="rho_sh")
x <- seq(0,15, 0.01)
d <- dlnorm(x, 1.7, 0.26)
lines(d ~ x, col = util$c_light_teal, lwd = 2)

util$plot_div_pairs('rho_sh', 'rho_sp[1]', 
                    samples, diagnostics, transforms = list('tau_kappa[1]' = 1))


par(mfrow = c(1,2), mar = c(4,4,2,2), cex.main = 1)
util$plot_expectand_pushforward(samples[['inner_tau_sck']], 30,
                                flim = c(0,1),
                                display_name="inner_tau_sck")
util$plot_expectand_pushforward(samples[['outer_tau_sck']], 30,
                                flim = c(0,5),
                                display_name="outer_tau_sck")

par(mfrow = c(1,1), mar = c(4,4,2,2), cex.main = 1)
util$plot_disc_pushforward_quantiles(samples, paste0('gamma_sck[',1:3,']'))

par(mfrow = c(1,3), mar = c(4,4,2,2), cex.main = 1)
util$plot_disc_pushforward_quantiles(samples, paste0('beta_sm[',1:2,']'))
util$plot_disc_pushforward_quantiles(samples, paste0('beta_gdd[',1:2,']'))
util$plot_disc_pushforward_quantiles(samples, paste0('beta_vpd[',1:2,']'))

par(mfrow = c(2,3), mar = c(4,4,2,2), cex.lab = 1.2)
for(s in 1:3){
  trees_s <- which(data$stand_idxs==s)
  minyear_s <- min(data$all_years_idxs[data$tree_start_idxs[trees_s]])
  maxyear_s <- max(data$all_years_idxs[data$tree_end_idxs[trees_s]])
  names <- sapply(minyear_s:maxyear_s,
                  function(sp) paste0('f_sh[', s, ',', sp, ']'))
  util$plot_realizations(samples, names, plot_xs = data$all_years[minyear_s:maxyear_s], N_plots = 50,
                         main = paste0('Stand ', s), ylab = 'f shortGP',
                         display_xlim = data$all_years[c(minyear_s,maxyear_s)], display_ylim = c(-5,2))
  abline(v = 2002, lty = 2, col = 'grey50')
}
for(s in 1:3){
  trees_s <- which(data$stand_idxs==s)
  minyear_s <- min(data$all_years_idxs[data$tree_start_idxs[trees_s]])
  maxyear_s <- max(data$all_years_idxs[data$tree_end_idxs[trees_s]])
  names <- sapply(minyear_s:maxyear_s,
                  function(sp) paste0('delta_sck[', s, ',', sp, ']'))
  util$plot_realizations(samples, names, plot_xs = data$all_years[minyear_s:maxyear_s], N_plots = 50,
                         main = "", ylab = 'delta shock',
                         display_xlim = data$all_years[c(minyear_s,maxyear_s)], , display_ylim = c(-5,2))
  abline(v = 2002, lty = 2, col = 'grey50')
}

s <- 3
trees_s <- which(data$stand_idxs==s)
par(mfrow = c(4,2), mar = c(4,4,2,2), cex.lab = 1.2)
for(t in sample(trees_s,8)){
  idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- sapply(idxs_t,
                  function(x) paste0('log_rw_pred[', x, ']'))
  util$plot_conn_pushforward_quantiles(
    samples, names, plot_xs = data$years[idxs_t], baseline_values = data$log_rw_obs[idxs_t],
    main = paste0(data$uniq_species_ids[data$species_idxs[t]], ' - ', t), ylab = 'Log ring width (mm)', display_ylim = c(-12,2))
}         

par(mfrow = c(1,2), mar = c(4,4,2,2), cex.lab = 1.2)
for(t in c(116, 121)){
  idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- sapply(idxs_t,
                  function(x) paste0('log_rw_pred[', x, ']'))
  util$plot_conn_pushforward_quantiles(
    samples, names, plot_xs = data$years[idxs_t], baseline_values = data$log_rw_obs[idxs_t],
    main = paste0('Stand ', s, ' - ', data$uniq_species_ids[data$species_idxs[t]], ' - tree ', t), ylab = 'Log ring width (mm)', display_ylim = c(-12,2))
}  
