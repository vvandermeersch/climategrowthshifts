rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('functions/custom_functions.R', local=util)
setwd(wd)

data <- readRDS(file.path(wd, 'output/model', 'data_28oct2025_gymno_azca.rds'))
datasets <- readRDS(file.path(wd, 'output/model', 'datasets_28oct2025_gymno_azca.rds'))
fit <- readRDS(file.path(wd, "model/shocks/output", "fit_28oct2025_gymno_azca.rds")) 

# Computational diagnostics
diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit)
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
                                         'delta_tilde_sck', 'inner_tau_sck', 'outer_tau_sck', 'gamma_sck',
                                         'sigma'),
                                       check_arrays=TRUE)
util$check_all_expectand_diagnostics(base_samples)

util$plot_div_pairs(paste0('delta_tilde_sck[67,44]'), 'inner_tau_sck', 
                    samples, diagnostics, transforms = list('inner_tau_sck' = 1))
util$plot_div_pairs(paste0('delta_tilde_sck[71,44]'), 'inner_tau_sck', 
                    samples, diagnostics, transforms = list('inner_tau_sck' = 1))
util$plot_div_pairs(paste0('delta_tilde_sck[75,44]'), 'inner_tau_sck', 
                    samples, diagnostics, transforms = list('inner_tau_sck' = 1))
util$plot_div_pairs(paste0('delta_tilde_sck[77,44]'), 'inner_tau_sck', 
                    samples, diagnostics, transforms = list('inner_tau_sck' = 1))

par(mfrow=c(1,3), mar = c(4.5,1.5,0,1))
hist(rbeta(1024, 1000, 10), breaks=seq(0.8, 1, len = 50),
     main="", col=util$c_light, border=util$c_light_highlight,
     xlab="gamma_sck", yaxt='n', ylab="", ylim=c(0, 5e2), xlim = c(0.8,1),
     mar = c(1,1,1,1))
hist(samples[['gamma_sck']], breaks=seq(0.8, 1, len = 50),
     main="", col=util$c_dark, border=util$c_dark_highlight, add=T)
hist(abs(rnorm(1024, 0, log(1.5)/ 2.57)), breaks=seq(0, .75, len = 50),
     main="", col=util$c_light, border=util$c_light_highlight,
     xlab="inner_tau_sck", yaxt='n', ylab="", ylim=c(0, 5e2), xlim = c(0,.75),
     mar = c(1,1,1,1))
hist(samples[['inner_tau_sck']], breaks=seq(0, .75, len = 50),
     main="", col=util$c_dark, border=util$c_dark_highlight, add=T)
hist(abs(rnorm(1024, 0, log(5)/ 2.57)), breaks=seq(0, 5, len = 50),
     main="", col=util$c_light, border=util$c_light_highlight,
     xlab="outer_tau_sck", yaxt='n', ylab="", ylim=c(0, 5e2), xlim = c(0,5),
     mar = c(1,1,1,1))
hist(samples[['outer_tau_sck']],  breaks=seq(0, 5, len = 50),
     main="", col=util$c_dark, border=util$c_dark_highlight, add=T)


par(mfrow=c(1,2), mar = c(4.5,0,0,0))
xs <- seq(-5,5, 0.01)
inner <- dnorm(xs, 0, log(1.5) / 2.57)
outer <- dnorm(xs, 0, log(5) / 2.57)
plot(xs, inner, lwd =2 , col = util$c_light_highlight, type = 'l', yaxt='n', frame = FALSE, xlab = '')
lines(xs, outer, lwd =2 , col = util$c_mid_highlight)

xs <- seq(-5,5, 0.01)
inner <- dnorm(xs, 0, log(1.1) / 2.57)
outer <- dnorm(xs, 0, log(20) / 2.57)
plot(xs, inner, lwd =2 , col = util$c_light_highlight, type = 'l', yaxt='n', frame = FALSE, xlab = '')
lines(xs, outer, lwd =2 , col = util$c_mid_highlight)

par(mfrow = c(1,1), mar = c(4,4,2,2))
util$plot_expectand_pushforward(samples[['rho_sh']], 200,
                                flim = c(0,5),
                                display_name="rho_sh")

# Look at all stands
par(mfrow = c(3,3), mar = c(3,5.5,2,1), cex.lab = 1.2)
for(s in 1:data$N_stands){
  trees_s <- which(data$stand_idxs==s)
  minyear_s <- min(data$all_years_idxs[data$tree_start_idxs[trees_s]])
  maxyear_s <- max(data$all_years_idxs[data$tree_end_idxs[trees_s]])
  names <- sapply(minyear_s:maxyear_s,
                  function(sp) paste0('f_sh[', s, ',', sp, ']'))
  util$plot_realizations(samples, names, plot_xs = data$all_years[minyear_s:maxyear_s], N_plots = 50,
                         main = paste0("Stand-level GP, ", data$uniq_stand_ids[s]), ylab = 'With shocks',
                         display_xlim = data$all_years[c(minyear_s,maxyear_s)], display_ylim = c(-8,3))
  names <- sapply(minyear_s:maxyear_s,
                  function(sp) paste0('f_sh[', s, ',', sp, ']'))
  util$plot_realizations(samples, names, plot_xs = data$all_years[minyear_s:maxyear_s], N_plots = 50,
                         main = "Stand-level GP (zoom)", ylab = '',
                         display_xlim = data$all_years[c(minyear_s,maxyear_s)], display_ylim = c(-3,3))
  names <- sapply(minyear_s:maxyear_s,
                  function(sp) paste0('delta_sck[', s, ',', sp, ']'))
  util$plot_realizations(samples, names, plot_xs = data$all_years[minyear_s:maxyear_s], N_plots = 50,
                         main = "Stand-level shocks", ylab = '',
                         display_xlim = data$all_years[c(minyear_s,maxyear_s)], , display_ylim = c(-10,3))
  
}

# Look only at stands with several species
par(mfrow = c(3,3), mar = c(3,5.5,2,1), cex.lab = 1.2)
stands_mutlisp <-  which(data$uniq_stand_ids %in% names(which(table(datasets$grouped_stand) > 1)))
for(s in stands_mutlisp){
  trees_s <- which(data$stand_idxs==s)
  minyear_s <- min(data$all_years_idxs[data$tree_start_idxs[trees_s]])
  maxyear_s <- max(data$all_years_idxs[data$tree_end_idxs[trees_s]])
  names <- sapply(minyear_s:maxyear_s,
                  function(sp) paste0('f_sh[', s, ',', sp, ']'))
  util$plot_realizations(samples, names, plot_xs = data$all_years[minyear_s:maxyear_s], N_plots = 50,
                         main = paste0("Stand-level GP, ", data$uniq_stand_ids[s]), ylab = 'With shocks',
                         display_xlim = data$all_years[c(minyear_s,maxyear_s)], display_ylim = c(-8,3))
  names <- sapply(minyear_s:maxyear_s,
                  function(sp) paste0('f_sh[', s, ',', sp, ']'))
  util$plot_realizations(samples, names, plot_xs = data$all_years[minyear_s:maxyear_s], N_plots = 50,
                         main = "Stand-level GP (zoom)", ylab = '',
                         display_xlim = data$all_years[c(minyear_s,maxyear_s)], display_ylim = c(-3,3))
  names <- sapply(minyear_s:maxyear_s,
                  function(sp) paste0('delta_sck[', s, ',', sp, ']'))
  util$plot_realizations(samples, names, plot_xs = data$all_years[minyear_s:maxyear_s], N_plots = 50,
                         main = "Stand-level shocks", ylab = '',
                         display_xlim = data$all_years[c(minyear_s,maxyear_s)], , display_ylim = c(-10,3))
  
}

# S59 is super interesting
# 2 species
par(mfrow = c(1,1), mar = c(2,2,2,2))
util$plot_disc_pushforward_quantiles(samples, paste0('kappa_sh[',1:data$N_species,']'), xticklabs = data$uniq_species_ids)
par(mfrow = c(1,2), mar = c(4,4,2,2))
util$plot_expectand_pushforward(samples[[paste0('kappa_sh[',which(data$uniq_species_ids == 'PIPO'),']')]], 200,
                                flim = c(0,1),
                                display_name="kappa", main = 'PIPO')
util$plot_expectand_pushforward(samples[[paste0('kappa_sh[',which(data$uniq_species_ids == 'PSME'),']')]], 200,
                                flim = c(0,1),
                                display_name="kappa", main = 'PSME')

par(mfrow = c(1,2), mar = c(3,2.5,2,1), cex.lab = 1.2)
s <- which(data$uniq_stand_ids == 'S59')
trees_s <- which(data$stand_idxs==s)
minyear_s <- min(data$all_years_idxs[data$tree_start_idxs[trees_s]])
maxyear_s <- max(data$all_years_idxs[data$tree_end_idxs[trees_s]])
names <- sapply(minyear_s:maxyear_s,
                function(sp) paste0('f_sh[', s, ',', sp, ']'))
util$plot_realizations(samples, names, plot_xs = data$all_years[minyear_s:maxyear_s], N_plots = 50,
                       main = paste0("Stand-level GP"), ylab = '',
                       display_xlim = data$all_years[c(minyear_s,maxyear_s)], display_ylim = c(-10,4))
names <- sapply(minyear_s:maxyear_s,
                function(sp) paste0('delta_sck[', s, ',', sp, ']'))
util$plot_realizations(samples, names, plot_xs = data$all_years[minyear_s:maxyear_s], N_plots = 50,
                       main = "Stand-level shocks", ylab = '',
                       display_xlim = data$all_years[c(minyear_s,maxyear_s)], , display_ylim = c(-10,4))

trees_pipo_s <- sample(which(data$stand_idxs==s & data$species_idxs == which(data$uniq_species_ids == 'PIPO')),10)
trees_psme_s <- sample(which(data$stand_idxs==s & data$species_idxs == which(data$uniq_species_ids == 'PSME')),10)
par(mfrow = c(4,5), mar = c(3,5.5,2,1), cex.lab = 1.2)
for(t in trees_pipo_s[1:5]){
  idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- sapply(idxs_t,
                  function(x) paste0('log_rw_pred[', x, ']'))
  util$plot_conn_pushforward_quantiles(
    samples, names, plot_xs = data$years[idxs_t], baseline_values = data$log_rw_obs[idxs_t],
    main = paste0("PIPO tree ",t,", ", data$uniq_stand_ids[s]), ylab = 'Log ring width (mm)', display_ylim = c(-12,2))
}
for(t in trees_pipo_s[1:5]){
  idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- sapply(idxs_t,
                  function(x) paste0('log_rw_pred[', x, ']'))
  util$plot_conn_pushforward_quantiles(
    samples, names, plot_xs = data$years[idxs_t], baseline_values = data$log_rw_obs[idxs_t],
    main = '', ylab = 'Residuals', residual = TRUE, display_ylim = c(-6,3))
}
for(t in trees_pipo_s[6:10]){
  idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- sapply(idxs_t,
                  function(x) paste0('log_rw_pred[', x, ']'))
  util$plot_conn_pushforward_quantiles(
    samples, names, plot_xs = data$years[idxs_t], baseline_values = data$log_rw_obs[idxs_t],
    main = paste0("PIPO tree ",t,", ", data$uniq_stand_ids[s]), ylab = 'Log ring width (mm)', display_ylim = c(-12,2))
}
for(t in trees_pipo_s[6:10]){
  idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- sapply(idxs_t,
                  function(x) paste0('log_rw_pred[', x, ']'))
  util$plot_conn_pushforward_quantiles(
    samples, names, plot_xs = data$years[idxs_t], baseline_values = data$log_rw_obs[idxs_t],
    main = '', ylab = 'Residuals', residual = TRUE, display_ylim = c(-6,3))
}


for(t in trees_psme_s[1:5]){
  idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- sapply(idxs_t,
                  function(x) paste0('log_rw_pred[', x, ']'))
  util$plot_conn_pushforward_quantiles(
    samples, names, plot_xs = data$years[idxs_t], baseline_values = data$log_rw_obs[idxs_t],
    main = paste0("PSME tree ",t,", ", data$uniq_stand_ids[s]), ylab = 'Log ring width (mm)', display_ylim = c(-12,2))
}
for(t in trees_psme_s[1:5]){
  idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- sapply(idxs_t,
                  function(x) paste0('log_rw_pred[', x, ']'))
  util$plot_conn_pushforward_quantiles(
    samples, names, plot_xs = data$years[idxs_t], baseline_values = data$log_rw_obs[idxs_t],
    main = '', ylab = 'Residuals', residual = TRUE, display_ylim = c(-6,3))
}
for(t in trees_psme_s[6:10]){
  idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- sapply(idxs_t,
                  function(x) paste0('log_rw_pred[', x, ']'))
  util$plot_conn_pushforward_quantiles(
    samples, names, plot_xs = data$years[idxs_t], baseline_values = data$log_rw_obs[idxs_t],
    main = paste0("PSME tree ",t,", ", data$uniq_stand_ids[s]), ylab = 'Log ring width (mm)', display_ylim = c(-12,2))
}
for(t in trees_psme_s[6:10]){
  idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- sapply(idxs_t,
                  function(x) paste0('log_rw_pred[', x, ']'))
  util$plot_conn_pushforward_quantiles(
    samples, names, plot_xs = data$years[idxs_t], baseline_values = data$log_rw_obs[idxs_t],
    main = '', ylab = 'Residuals', residual = TRUE, display_ylim = c(-6,3))
}

