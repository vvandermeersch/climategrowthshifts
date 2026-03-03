rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(cmdstanr)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)


data <-readRDS(file.path(wd, 'output/model', 'data_30jan2025_gymnosperms_standclimate_19502024_7species.rds'))
fit <- readRDS(file.path(wd, 'output/model/model13', 'fit_16feb2026_gymnosperms_standclimate_19502024_7species_model13_updatedpriors_onlyconcordance.rds'))
fit_nonc_old <- readRDS(file.path(wd, 'output/model/model13', 'fit_30jan2025_gymnosperms_standclimate_19502024_7species_model13_updatedpriors.rds'))

# Parameter samples
param_samples <- fit$draws(variables = c("mu_alpha","mu_gdd", "mu_vpd", "mu_pre", 
                                             "tau_alpha","tau_gdd", "tau_vpd", "tau_pre",
                                             "alpha","beta_gdd", "beta_vpd", "beta_pre",
                                             
                                             "mu_tau_sck", "tau_tau_sck", "tau_sck",
                                             "phi_sck", "omega_conc_sck",
                                             
                                             "mu_rho", "tau_rho", "rho_sp",
                                             "mu_gamma", "tau_gamma", "gamma_sp",
                                             
                                             "rho_sh",
                                             "mu_kappa", "tau_kappa", "kappa_sh",
                                             
                                             "phi_sck0", "tau_phi_sck", "alpha_phi_sck",
                                             "omega_conc_sck0", "mu_omega_conc_sck", "tau_omega_conc_sck", "alpha_omega_conc_sck",
                                             # "mu_logdelta_omega_nonconc_sck", "tau_logdelta_omega_nonconc_sck", "logdelta_omega_nonconc_sck", "omega_nonconc_sck",
                                             
                                             "sigma"))
names <- dimnames(param_samples)$variable
param_samples <- lapply(1:dim(param_samples)[3], function(k){t(matrix(param_samples[1:dim(param_samples)[1],1:dim(param_samples)[2],k], 
                                                                nrow = dim(param_samples)[1], ncol = dim(param_samples)[2]))})
names(param_samples) <- names
util$check_all_expectand_diagnostics(param_samples)
gc()


standgp_samples <- fit_nonc$draws(variables = c("f_sh"))
names <- dimnames(standgp_samples)$variable
standgp_samples <- lapply(1:dim(standgp_samples)[3], function(k){t(matrix(standgp_samples[1:dim(standgp_samples)[1],1:dim(standgp_samples)[2],k], 
                                                                      nrow = dim(standgp_samples)[1], ncol = dim(standgp_samples)[2]))})
names(standgp_samples) <- names

# Generate quantities
mod_gq <- cmdstan_model(file.path(wd, 'model/stan', 'model13_updatedpriors_wGQ.stan'))
data_gq <- data
data_gq$grainsize <-  ceiling(data_gq$N_stands/8)
data_gq$uniq_tree_ids <- NULL
data_gq$uniq_species_ids <- NULL
data_gq$uniq_stand_ids <- NULL
data_gq$N_clades <- 1
data_gq$uniq_stand_lat <- NULL
data_gq$uniq_stand_lon <- NULL

data_gq$tree_pred_idxs <- which(grepl('wy067', data$uniq_tree_ids))
data_gq$tree_pred_idxs <- which(data$stand_idxs == 6)
data_gq$N_pred <- sum(sapply(data_gq$tree_pred_idxs , function(t) length(data$tree_start_idxs[t]:data$tree_end_idxs[t])))
data_gq$N_trees_pred <- length(data_gq$tree_pred_idxs )

data_gq$tree_pred_idxs <- 1:data$N_trees
data_gq$N_pred <- sum(sapply(data_gq$tree_pred_idxs , function(t) length(data$tree_start_idxs[t]:data$tree_end_idxs[t])))
data_gq$N_trees_pred <- length(data_gq$tree_pred_idxs )
  
fit_gq <- mod_gq$generate_quantities(fit_nonc, data = data_gq, seed = 5838293, parallel_chains = 1)
gc()
gq_samples <- fit_gq$draws(variables = c('sck_nonconc_state','sck_conc_state'))
gc()
names <- dimnames(gq_samples)$variable
gq_samples <- lapply(1:dim(gq_samples)[3], function(k){t(matrix(gq_samples[1:dim(gq_samples)[1],1:dim(gq_samples)[2],k], 
                                                                nrow = dim(gq_samples)[1], ncol = dim(gq_samples)[2]))})
names(gq_samples) <- names

# Look at the stand or094
tree_idxs <- which(grepl('or094', data$uniq_tree_ids))
stand <- data$stand_idxs[tree_idxs[1]]
substand <-  data$stand_species_idxs[tree_idxs[1]]
util$ensemble_mcmc_quantile_est(param_samples[[paste0('phi_sck[',stand,']')]], c(0.5))
util$ensemble_mcmc_quantile_est(param_samples[[paste0('omega_conc_sck[',substand,']')]], c(0.5))
par(mfrow = c(7,2), mar = c(2,3,0,1), mgp = c(2,0.5,0),
    cex.lab = 0.8, cex.axis = 0.7)
for(t in tree_idxs){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  stand <- data$stand_idxs[t]
  # names <- paste0("log_rw_pred[",idxs,"]")
  # 
  # util$plot_conn_pushforward_quantiles(gq_samples, names, data$years[idxs]-1979,
  #                                      xlab="Year", ylab="Log ring width (per mm)", 
  #                                      display_ylim=c(-8, 2), display_xlim = c(1980, 2024)-1979)
  # points(data$years[idxs]-1979, log(data$rw_obs[idxs]), pch=16, cex=1, col="white")
  # points(data$years[idxs]-1979, log(data$rw_obs[idxs]), pch=16, cex=0.5, col="black")
  # text(x = .5, y = -6.9, labels = data$uniq_tree_ids[t], cex = 1, adj = 0)
  # text(x = .5, y = -7.7, labels = paste0('Species:', data$uniq_species_ids[data$species_idxs[t]]), cex = 1, adj = 0)
  
  # start <- 1+(stand-1)*data$N_all_years
  # end <- stand*data$N_all_years
  # lines(1:data$N_all_years, data$pre_obs[start:end]/10, col = 'darkblue')
  # lines(1:data$N_all_years, data$vpd_obs[start:end]/20, col = 'darkred')
  
  names <- paste0("sck_state[",idxs,"]")
  util$plot_disc_pushforward_quantiles(gq_samples, names, data$years[idxs],
                                       xlab="Year", ylab="Shock state",
                                       display_ylim=c(-0.05, 1.05))
  
  # names <- paste0("delta_sck[",idxs,"]")
  # util$plot_conn_pushforward_quantiles(gq_samples, names, data$years[idxs],
  #                                      xlab="Year", ylab="Shock amplitude", 
  #                                      display_ylim=c(-5, 1), display_xlim = c(1950, 2024))
}


for(y in 1:max(data$N_years[tree_idxs])){
  gq_samples[[paste0('shockhere[', y, ']')]] <- 0
}

ny <- rep(0, max(data$N_years[tree_idxs]))
for(t in tree_idxs){
  print(t)
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  for(i in idxs){
    y <- data$years[i]-data$all_years[1]+1
    gq_samples[[paste0('shockhere[', y, ']')]]  <- 
      gq_samples[[paste0('shockhere[', y, ']')]] +
      gq_samples[[paste0("sck_state[",i,"]")]] 
    ny[y] <- ny[y]+1
  }
}

gq_samples[[paste0('shockhere_all')]] <- 0
for(y in 1:max(data$N_years[tree_idxs])){
  gq_samples[[paste0('shockhere[', y, ']')]] <-  gq_samples[[paste0('shockhere[', y, ']')]]/ny[y]*100
  gq_samples[[paste0('shockhere_all')]] <- gq_samples[[paste0('shockhere_all')]] + 
    gq_samples[[paste0('shockhere[', y, ']')]]
}
gq_samples[[paste0('shockhere_all')]] <- gq_samples[[paste0('shockhere_all')]]/max(data$N_years[tree_idxs])

util$ensemble_mcmc_quantile_est(gq_samples[[paste0('shockhere[', y, ']')]], c(0.5))

par(mfrow = c(1,1), mar = c(1.5,3,1,1))
names <- paste0("shockhere[",1:max(data$N_years[tree_idxs]),"]")
util$plot_disc_pushforward_quantiles(gq_samples, names, data$years[idxs],
                                     xlab="Year", ylab="% trees in shock",
                                     display_ylim=c(-0.05, 100.05))

util$ensemble_mcmc_quantile_est(gq_samples[['shockhere_all']], c(0.05, 0.5, 0.95))


shock_prob <- param_samples[[paste0('phi_sck[',stand,']')]]*
  param_samples[[paste0('omega_conc_sck[',substand,']')]] +
  (1-param_samples[[paste0('phi_sck[',stand,']')]])*
  param_samples[[paste0('omega_nonconc_sck[',substand,']')]]
util$ensemble_mcmc_quantile_est(shock_prob, c(0.05, 0.5, 0.95))


# Look at concordant vs. non-concordant shocks



