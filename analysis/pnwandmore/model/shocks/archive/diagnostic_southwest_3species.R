rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('functions/custom_functions.R', local=util)
setwd(wd)

data <- readRDS(file.path(wd, 'output/model', 'data_04nov2025_southwest_3species.rds'))
datasets <- readRDS(file.path(wd, 'output/model', 'datasets_04nov2025_southwest_3species.rds'))
samples <- readRDS(file.path(wd, "model/shocks/output", "samples_04nov2025_southwest_3species.rds"))

base_samples <- util$filter_expectands(samples,
                                       c('alpha', 
                                         'mu_gdd', 'mu_sm', 'mu_vpd',
                                         'tau_gdd', 'tau_sm', 'tau_vpd',
                                         'beta_gdd', 'beta_sm', 'beta_vpd',
                                         'mu_rho', 'mu_gamma', 'tau_rho', 'tau_gamma', 
                                         'rho_sp', 'gamma_sp',
                                         'mu_kappa', 'tau_kappa',
                                         'rho_sh', 'kappa_sh',
                                         'inner_tau_sck', 'outer_tau_sck', 'gamma_sck',
                                         'mu_kappa_sck', 'tau_kappa_sck', 'kappa_sck',
                                         'sigma'),
                                       check_arrays=TRUE)
util$check_all_expectand_diagnostics(base_samples)

par(mfrow = c(1,4), mar = c(4,1.5,1,1), cex.lab = 1.2, cex.main = 1)
hist(abs(rnorm(4e3, 0, log(1.1)/ 2.57)), breaks=seq(0, 1, 0.025),
     main="", col=util$c_light, border=util$c_light_highlight,
     xlab=expression(sigma), yaxt='n', ylab="", ylim=c(0, 3e3), xlim = c(0,1),
     mar = c(1,1,1,1))
hist(samples[['sigma']], breaks=seq(0, 1, 0.025),
     main="", col=util$c_dark, border=util$c_dark_highlight, add=T)

hist(abs(rlnorm(4e3, 1.7, 0.26)), breaks=seq(0, 20, 0.5),
     main="", col=util$c_light, border=util$c_light_highlight,
     xlab=expression(rho[short]), yaxt='n', ylab="", ylim=c(0, 1e3), xlim = c(0,15),
     mar = c(1,1,1,1))
hist(samples[['rho_sh']], breaks=seq(0, 20, 0.5),
     main="", col=util$c_dark, border=util$c_dark_highlight, add=T)

hist(abs(rnorm(4e3, 0, log(1.5)/ 2.57)), breaks=seq(0, 1, 0.025),
     main="", col=util$c_light, border=util$c_light_highlight,
     xlab=expression(tau[inner]), yaxt='n', ylab="", ylim=c(0, 1e3), xlim = c(0,1),
     mar = c(1,1,1,1))
hist(samples[['inner_tau_sck']], breaks=seq(0, 1, 0.025),
     main="", col=util$c_dark, border=util$c_dark_highlight, add=T)

hist(abs(rnorm(4e3, 0, log(5)/ 2.57)), breaks=seq(0, 5, 0.15),
     main="", col=util$c_light, border=util$c_light_highlight,
     xlab=expression(tau[outer]), yaxt='n', ylab="", ylim=c(0, 1e3), xlim = c(0,5),
     mar = c(1,1,1,1))
hist(samples[['outer_tau_sck']], breaks=seq(0, 5, 0.15),
     main="", col=util$c_dark, border=util$c_dark_highlight, add=T)

util$plot_pairs_by_chain(samples[['rho_sh']], expression(rho[sh]),
                         samples[['inner_tau_sck']], expression(tau[inner]))

par(mfrow = c(2,3), mar = c(4,5.5,2,1), cex.lab = 1.2, cex.main = 1)
util$plot_disc_pushforward_quantiles(samples, paste0('rho_sp[',1:data$N_species,']'),
                                     xticklabs = data$uniq_species_ids, ylab = expression(rho[species]))
util$plot_disc_pushforward_quantiles(samples, paste0('gamma_sp[',1:data$N_species,']'),
                                     xticklabs = data$uniq_species_ids, ylab = expression(gamma[species]))
util$plot_disc_pushforward_quantiles(samples, paste0('gamma_sck[',1:data$N_species,']'),
                                     xticklabs = data$uniq_species_ids, ylab = expression(phi[shock]))
util$plot_disc_pushforward_quantiles(samples, paste0('beta_gdd[',1:data$N_species,']'),
                                     xticklabs = data$uniq_species_ids, ylab = expression(beta[GDD]))
util$plot_disc_pushforward_quantiles(samples, paste0('beta_vpd[',1:data$N_species,']'),
                                     xticklabs = data$uniq_species_ids, ylab = expression(beta[VPD]))
util$plot_disc_pushforward_quantiles(samples, paste0('beta_sm[',1:data$N_species,']'),
                                     xticklabs = data$uniq_species_ids, ylab = expression(beta[SM]))


pdf(file.path(wd, 'model/shocks/figures', '18sampledstands_southwest_3species.pdf'), height = 12, width = 15)
par(mfrow = c(9,4), mar = c(3,5.5,2,1), cex.lab = 1.2, cex.main = 1)
for(s in sample(1:data$N_stands, 18)){
  
  trees_idx <- which(data$stand_idxs==s)
  minyear_s <- min(data$all_years_idxs[data$tree_start_idxs[trees_idx]])
  maxyear_s <- max(data$all_years_idxs[data$tree_end_idxs[trees_idx]])
  names <- sapply(minyear_s:maxyear_s,
                  function(y) paste0('f_sh[', s, ',', y, ']'))
  util$plot_realizations(samples, names, plot_xs = data$all_years[minyear_s:maxyear_s], N_plots = 50,
                         main = paste0("Stand-level GP - ", data$uniq_stand_ids[s], ' (', 
                                       length(unique( data$species_idxs[trees_idx])),' species)'), 
                         ylab = '',
                         display_xlim = c(1980,2020), display_ylim = c(-8,3))
  text(x = 1980, y = -7.2, 
       label = paste0('Q5-95 [', 
                      paste0(round(util$ensemble_mcmc_quantile_est(samples[['rho_sh']], 
                                                                   c(0.05,0.95)),3), collapse = '-'), ']'),
       cex = 0.7, adj = 0)
  
  stand_species_idxs <- unique(data.frame(
    species = data$species_idxs[trees_idx],
    substands = data$substand_idxs[trees_idx]
  ))
  
  for(i in 1:nrow(stand_species_idxs)){
    sp <- stand_species_idxs[i, 'species']
    sbs <- stand_species_idxs[i, 'substands'] 
    
    minyear_sbs <- min(data$all_years_idxs[data$tree_start_idxs[trees_idx[data$substand_idxs[trees_idx] == sbs]]])
    maxyear_sbs <- max(data$all_years_idxs[data$tree_end_idxs[trees_idx[data$substand_idxs[trees_idx] == sbs]]])
    
    names <- sapply(minyear_sbs:maxyear_sbs,
                    function(y) paste0('delta_sck[', sbs, ',', y, ']'))
    util$plot_realizations(samples, names, plot_xs = data$all_years[minyear_sbs:maxyear_sbs], N_plots = 50,
                           main = paste0("Species-specific stand-level shocks (",data$uniq_species_ids[sp],
                                         ", ",length(trees_idx[data$species_idxs[trees_idx] == sp])," trees)"), ylab = '',
                           display_xlim = c(1980,2020), display_ylim = c(-8,3))
    text(x = 1980, y = -7.2, 
         label = paste0('Q5-95 [', 
                        paste0(round(util$ensemble_mcmc_quantile_est(samples[[paste0('gamma_sck[',sp,']')]], 
                                                                     c(0.05,0.95)),3), collapse = '-'), ']'),
         cex = 0.7, adj = 0)
  }
  
  if(nrow(stand_species_idxs) != 3){
    for(i in 1:(3-nrow(stand_species_idxs))){
      plot.new()
    }
  }
  
  
}
dev.off()

pdf(file.path(wd, 'model/shocks/figures', '18sampledtrees_southwest_3species.pdf'), height = 15, width = 13)
par(mfrow = c(9,2), mar = c(3,5.5,2,1), cex.lab = 1.2, cex.main = 1)
for(t in sample(1:data$N_trees, 18)){
  
  tree_idx <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- sapply(tree_idx, function(i) paste0('log_rw_pred[', i,']'))
  util$plot_conn_pushforward_quantiles(samples, names, plot_xs = data$years[tree_idx], 
                                       display_xlim = c(1980,2020), display_ylim = c(-9.5,3), ylab = 'Marginal quantiles')
  abline(h = log(1e-4), lty = 2, col = util$c_light_teal)
  points(data$years[tree_idx], data$log_rw_obs[tree_idx], pch=16, cex=1.2, col="white")
  points(data$years[tree_idx], data$log_rw_obs[tree_idx], pch=16, cex=0.8, col="black")
  util$plot_conn_pushforward_quantiles(samples, names, plot_xs = data$years[tree_idx], baseline_values = data$log_rw_obs[tree_idx],
                                       display_xlim = c(1980,2020), residual = TRUE, ylab = 'Residuals')
  
}
dev.off()

