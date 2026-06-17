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
fit <- readRDS(file.path(wd, "model/shocks/output", "fit_28oct2025_gymno_azca_4chains.rds")) 
fitkappa <- readRDS(file.path(wd, "model/shocks/output", "fit_28oct2025_gymno_azca_kappa.rds")) 

samples <- util$extract_expectand_vals(fit)
sampleskappa <- util$extract_expectand_vals(fitkappa)
diagnosticskappa <- util$extract_hmc_diagnostics(fitkappa)
gc()
base_samples <- util$filter_expectands(sampleskappa,
                                       c('alpha', 
                                         'mu_gdd', 'mu_sm', 'mu_vpd',
                                         'tau_gdd', 'tau_sm', 'tau_vpd',
                                         'beta_gdd', 'beta_sm', 'beta_vpd',
                                         'mu_rho', 'mu_gamma', 'tau_rho', 'tau_gamma', 
                                         'rho_sp', 'gamma_sp',
                                         # 'f_tilde_sh',
                                         'mu_kappa', 'tau_kappa',
                                         'rho_sh', 'kappa_sh',
                                         # 'delta_tilde_sck', 
                                         'inner_tau_sck', 'outer_tau_sck', 'gamma_sck',
                                         'mu_kappa_sck', 'tau_kappa_sck', 'kappa_sck',
                                         'sigma'),
                                       check_arrays=TRUE)
util$check_all_expectand_diagnostics(base_samples)

for(s in 1:data$N_species){
  util$plot_pairs_by_chain(sampleskappa[[paste0('kappa_sck[',s,']')]], paste0('kappa shock ',data$uniq_species_ids[s],''), 
                           sampleskappa[[paste0('kappa_sh[',s,']')]], paste0('kappa short GP ',data$uniq_species_ids[s],''))
}

# No change in kappa_sh
par(mfrow = c(2,2), mar = c(2,4,2,2), cex.main = 1)
util$plot_disc_pushforward_quantiles(samples, paste0('kappa_sh[',1:data$N_species,']'), xticklabs = data$uniq_species_ids,
                                     main = 'Without species scaling of shocks', ylab = 'kappa shortGP')
util$plot_disc_pushforward_quantiles(sampleskappa, paste0('kappa_sh[',1:data$N_species,']'), xticklabs = data$uniq_species_ids,
                                     main = 'With species scaling of shocks', ylab = 'kappa shortGP')
plot.new()
util$plot_disc_pushforward_quantiles(sampleskappa, paste0('kappa_sck[',1:data$N_species,']'), xticklabs = data$uniq_species_ids,
                                     ylab = 'kappa shock')

# No change in rho_sh
par(mfrow = c(1,2), mar = c(4,4,2,2), cex.main = 1)
util$plot_expectand_pushforward(samples[['rho_sh']], 200,
                                flim = c(0,5),
                                display_name="rho_sh", main = 'Without species scaling of shocks')
util$plot_expectand_pushforward(sampleskappa[['rho_sh']], 200,
                                flim = c(0,5),
                                display_name="rho_sh", main = 'With species scaling of shocks')


# Focus on one stand
par(mfrow = c(2,2), mar = c(4,4,2,2), cex.lab = 1.2)
s <- which(data$uniq_stand_ids == 'S59') # az591 and az592
trees_s <- which(data$stand_idxs==s)
minyear_s <- min(data$all_years_idxs[data$tree_start_idxs[trees_s]])
maxyear_s <- max(data$all_years_idxs[data$tree_end_idxs[trees_s]])
names <- sapply(minyear_s:maxyear_s,
                function(sp) paste0('f_sh[', s, ',', sp, ']'))
util$plot_realizations(samples, names, plot_xs = data$all_years[minyear_s:maxyear_s], N_plots = 50,
                       main = 'Without species scaling of shocks', ylab = 'f shortGP',
                       display_xlim = data$all_years[c(minyear_s,maxyear_s)], display_ylim = c(-10,4))
abline(v = 2002, lty = 2, col = 'grey50')
util$plot_realizations(sampleskappa, names, plot_xs = data$all_years[minyear_s:maxyear_s], N_plots = 50,
                       main = 'With species scaling of shocks', ylab = '',
                       display_xlim = data$all_years[c(minyear_s,maxyear_s)], display_ylim = c(-10,4))
abline(v = 2002, lty = 2, col = 'grey50')
names <- sapply(minyear_s:maxyear_s,
                function(sp) paste0('delta_sck[', s, ',', sp, ']'))
util$plot_realizations(samples, names, plot_xs = data$all_years[minyear_s:maxyear_s], N_plots = 50,
                       main = "", ylab = 'delta shock',
                       display_xlim = data$all_years[c(minyear_s,maxyear_s)], , display_ylim = c(-10,4))
abline(v = 2002, lty = 2, col = 'grey50')
util$plot_realizations(sampleskappa, names, plot_xs = data$all_years[minyear_s:maxyear_s], N_plots = 50,
                       main = "", ylab = '',
                       display_xlim = data$all_years[c(minyear_s,maxyear_s)], , display_ylim = c(-10,4))
abline(v = 2002, lty = 2, col = 'grey50')
