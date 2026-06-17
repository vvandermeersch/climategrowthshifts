rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('functions/custom_functions.R', local=util)
setwd(wd)

data <- readRDS(file.path(wd, 'output/model', 'data_05nov2025_southwest_min3speciesplots.rds'))
datasets <- readRDS(file.path(wd, 'output/model', 'datasets_05nov2025_southwest_min3speciesplots.rds'))
samples <- readRDS(file.path(wd, "model/shocks/output", "samples_05nov2025_southwest_min3speciesplots.rds"))
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


par(mfrow = c(1,4), mar = c(4,5.5,2,1), cex.lab = 1.2, cex.main = 1)
util$plot_expectand_pushforward(samples[['sigma']], 20, display_name = expression(sigma))
util$plot_expectand_pushforward(samples[['rho_sh']], 20, display_name = expression(rho[sh]))
util$plot_expectand_pushforward(samples[['inner_tau_sck']], 20, display_name = expression(tau[inner]))
util$plot_expectand_pushforward(samples[['outer_tau_sck']], 20, display_name = expression(tau[outer]))


util$plot_pairs_by_chain(samples[['sigma']], expression(sigma),
                         samples[['rho_sh']], expression(rho[sh]))
util$plot_pairs_by_chain(samples[['sigma']], expression(sigma),
                         samples[['inner_tau_sck']], expression(tau[inner]))



par(mfrow = c(6,6), mar = c(3,5.5,2,1), cex.lab = 1.2, cex.main = 1)
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
    
    minyear_sbs <- min(data$all_years_idxs[data$tree_start_idxs[trees_idx[data$substand_idxs[trees_idx] == sbs]]])
    maxyear_sbs <- max(data$all_years_idxs[data$tree_end_idxs[trees_idx[data$substand_idxs[trees_idx] == sbs]]])
    
    names <- sapply(minyear_sbs:maxyear_sbs,
                    function(y) paste0('delta_sck[', sbs, ',', y, ']'))
    util$plot_realizations(samples, names, plot_xs = data$all_years[minyear_sbs:maxyear_sbs], N_plots = 50,
                           main = paste0("Species-specific shocks (",data$uniq_species_ids[sp],
                                         ", ",length(trees_idx[data$species_idxs[trees_idx] == sp])," trees)"), ylab = '',
                           display_xlim = c(1980,2020), display_ylim = c(-8,3))
    text(x = 1980, y = -7.2, 
         label = paste0('Q5-95 [', 
                        paste0(round(util$ensemble_mcmc_quantile_est(samples[[paste0('gamma_sck[',sp,']')]], 
                                                                     c(0.05,0.95)),3), collapse = '-'), ']'),
         cex = 0.7, adj = 0)
  }
  
  if(nrow(stand_species_idxs) != 5){
    for(i in 1:(5-nrow(stand_species_idxs))){
      plot.new()
    }
  }
}

