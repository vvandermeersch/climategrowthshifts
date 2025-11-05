rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('functions/custom_functions.R', local=util)
setwd(wd)

data <- readRDS(file.path(wd, 'output/model', 'data_03nov2025_aznmca_2speciesplots.rds'))

# Stand-level shocks with species scaling
fit <- readRDS(file.path(wd, "model/shocks/output", "fit_03nov2025_aznmca_2speciesplots.rds"))
samples <- util$extract_expectand_vals(fit);gc()

par(mfrow = c(9,3), mar = c(3,5.5,2,1), cex.lab = 1.2, cex.main = 1)
for(s in 1:data$N_stands){
  
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
       label = paste0('Short GP rho: Q5-95 [', 
                      paste0(round(util$ensemble_mcmc_quantile_est(samples[['rho_sh']], 
                                                                   c(0.05,0.95)),3), collapse = '-'), ']'),
       cex = 0.9, adj = 0)
  
  names <- sapply(minyear_s:maxyear_s,
                  function(y) paste0('delta_sck[', s, ',', y, ']'))
  util$plot_realizations(samples, names, plot_xs = data$all_years[minyear_s:maxyear_s], N_plots = 50,
                         main = "Stand-level shocks", ylab = '',
                         display_xlim = c(1980,2020) , display_ylim = c(-8,3))
  text(x = 1980, y = -7.2, 
       label = paste0('Stand-level gamma: Q5-95 [', 
                      paste0(round(util$ensemble_mcmc_quantile_est(samples[[paste0('gamma_sck[',s,']')]], 
                                                                   c(0.05,0.95)),3), collapse = '-'), ']'),
       cex = 0.9, adj = 0)
  
  plot.new()
}

# Species-specific stand-level shocks
fit <- readRDS(file.path(wd, "model/shocks/output", "fit_03nov2025_aznmca_2speciesplots_standsxspecies.rds"))
samples <- util$extract_expectand_vals(fit);gc()

par(mfrow = c(9,3), mar = c(3,5.5,2,1), cex.lab = 1.2, cex.main = 1)
for(s in 1:data$N_stands){
  
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
      
    names <- sapply(minyear_s:maxyear_s,
                    function(y) paste0('delta_sck[', sbs, ',', y, ']'))
    util$plot_realizations(samples, names, plot_xs = data$all_years[minyear_s:maxyear_s], N_plots = 50,
                           main = paste0("Species-specific stand-level shocks (",data$uniq_species_ids[sp],")"), ylab = '',
                           display_xlim = c(1980,2020) , display_ylim = c(-8,3))
    text(x = 1980, y = -7.2, 
         label = paste0('Q5-95 [', 
                        paste0(round(util$ensemble_mcmc_quantile_est(samples[[paste0('gamma_sck')]], 
                                                                     c(0.05,0.95)),3), collapse = '-'), ']'),
         cex = 0.7, adj = 0)
  }
}

# Species-specific stand-level shocks with species-specific gamma
fit <- readRDS(file.path(wd, "model/shocks/output", "fit_03nov2025_aznmca_2speciesplots_speciesshocks_speciesgamma.rds"))
samples <- util$extract_expectand_vals(fit);gc()

par(mfrow = c(9,3), mar = c(3,5.5,2,1), cex.lab = 1.2, cex.main = 1)
for(s in 1:data$N_stands){
  
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
    
    names <- sapply(minyear_s:maxyear_s,
                    function(y) paste0('delta_sck[', sbs, ',', y, ']'))
    util$plot_realizations(samples, names, plot_xs = data$all_years[minyear_s:maxyear_s], N_plots = 50,
                           main = paste0("Species-specific stand-level shocks (",data$uniq_species_ids[sp],")"), ylab = '',
                           display_xlim = c(1980,2020), display_ylim = c(-8,3))
    text(x = 1980, y = -7.2, 
         label = paste0('Q5-95 [', 
                        paste0(round(util$ensemble_mcmc_quantile_est(samples[[paste0('gamma_sck[',sp,']')]], 
                                                                     c(0.05,0.95)),3), collapse = '-'), ']'),
         cex = 0.7, adj = 0)
  }
}
par(mfrow = c(1,1), mar = c(2,4.5,2,2))
util$plot_disc_pushforward_quantiles(samples, paste0('gamma_sck[',1:data$N_species,']'), xticklabs = data$uniq_species_ids, ylab = expression(gamma))

par(mfrow = c(5,2), mar = c(2,4.5,2,2))
for(t in sample(1:data$N_trees, 5)){
  
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
