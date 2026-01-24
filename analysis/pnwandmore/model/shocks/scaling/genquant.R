rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(terra)
library(rnaturalearth)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)


data <- readRDS(file.path(wd, 'output/model', 'data_15jan2025_PIPO_standclimate_19802024.rds'))
fit_19802024 <- readRDS(file.path(wd, 'output/model/scaling', 'fit_15jan2025_PIPO_standclimate_19802024_customgrainsize_improved2_wgenquant.rds'))
param_samples <- fit_19802024$draws(variables = c("beta_gdd", "beta_vpd", "beta_pre", "phi_sck", "omega_conc_sck", "rho_sp", "rho_sh"))

# Stand 12 seems interesting!
s <- 12
util$ensemble_mcmc_quantile_est(t(matrix(param_samples[1:1000,1:4,paste0('phi_sck[',s,']')], nrow = 1000, ncol = 4)), c(0.5))
util$ensemble_mcmc_quantile_est(t(matrix(param_samples[1:1000,1:4,paste0('omega_conc_sck[',s,']')], nrow = 1000, ncol = 4)), c(0.5))
trees_s <- which(data$stand_idxs == 12)
data$stand_species_idxs[trees_s]
data$uniq_tree_ids[trees_s]

par(mfrow = c(2,2))
for(i in 255:256){
  t <- trees_s[i]
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- paste0("log_rw_pred[",idxs,"]")
  samples <- fit_19802024$draws(variables = names)
  samples <- lapply(1:dim(samples)[3], function(k){t(matrix(samples[1:dim(samples)[1],1:dim(samples)[2],k], nrow = dim(samples)[1], ncol = dim(samples)[2]))})
  names(samples) <- names
  
  util$plot_conn_pushforward_quantiles(samples, names, data$years[idxs],
                                       xlab="Year", ylab="Log Ring Width Per mm", 
                                       display_ylim=c(-8, 2), display_xlim = c(1980, 2024))
  points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=1.2, col="white")
  points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=0.8, col="black")
  text(x = 1983.5, y = -7, labels = data$uniq_tree_ids[t], cex = 0.8)
  
  
  names <- paste0("sck_state[",idxs,"]")
  samples <- fit_19802024$draws(variables = names)
  samples <- lapply(1:dim(samples)[3], function(k){t(matrix(samples[1:dim(samples)[1],1:dim(samples)[2],k], nrow = dim(samples)[1], ncol = dim(samples)[2]))})
  names(samples) <- names
  util$plot_conn_pushforward_quantiles(samples, names, data$years[idxs],
                                       xlab="Year", ylab="Shock state", 
                                       display_ylim=c(-0.05, 1.05), display_xlim = c(1980, 2024))
}

par(mfrow = c(1,1))
util$plot_expectand_pushforward(t(matrix(param_samples[1:1000,1:4,paste0('omega_conc_sck[',12,']')], nrow = 1000, ncol = 4)), 50,
                                flim = c(0,1), ylim = c(0,50),
                                display_name=expression(omega[shock]^concordant),
                                col = '#8f2727')


util$plot_expectand_pushforward(t(matrix(param_samples[1:1000,1:4,paste0('omega_conc_sck[',1,']')], nrow = 1000, ncol = 4)), 50,
                                flim = c(0,1), ylim = c(0,50),
                                display_name=expression(omega[shock]^concordant),
                                col = '#8f2727')
for(s in 1:data$N_stand_species){
  if(s == 12 | s == 1){next}
  util$plot_expectand_pushforward(t(matrix(param_samples[1:1000,1:4,paste0('omega_conc_sck[',s,']')], nrow = 1000, ncol = 4)), 50,
                                  flim = c(0,1), ylim = c(0,50),
                                  display_name=expression(omega[shock]^concordant),
                                  col = '#8f2727', add = TRUE)
}

