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

par(mfrow = c(1,3), mar = c(4,5.5,2,1), cex.lab = 1.2, cex.main = 1)
util$plot_expectand_pushforward(samples[['rho_sh']], 20)
util$plot_expectand_pushforward(samples[['inner_tau_sck']], 20)
util$plot_expectand_pushforward(samples[['outer_tau_sck']], 20)

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



par(mfrow = c(9,4), mar = c(3,5.5,2,1), cex.lab = 1.2, cex.main = 1)
for(s in sample(1:data$N_stands, 9)){
  
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
