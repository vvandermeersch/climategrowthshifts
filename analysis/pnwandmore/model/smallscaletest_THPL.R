rm(list = ls())
wd <- '/home/victor/projects/climategrowthshifts/analysis/pnw'

# Small-scale test model that account fro different core lengths

# Load model fit and data
fit <- readRDS(file.path(wd, 'output', 'model', 'fit_25june2025_testTHPL.rds')) # this was run on GenOuest
data <- readRDS(file.path(wd, 'output', 'model', 'data_25june2025_testTHPL.rds'))

# Load datasets (useful if I want to plot stands)
# Load treering data
datasets <- readRDS(file = file.path(wd, 'data/itrdb', paste0('itrdb_info.rds')))
datasets <- datasets[datasets$last_year >= 1999,] # at least 20 years of observations
  
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)  

set.seed(123456789)

# Posterior quantificiation - diagnostics
diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit)
base_samples <- util$filter_expectands(samples,
                                       c('alpha', 'beta_gdd', 'beta_sm', 
                                         'rho', 'gamma',
                                         'f_tilde_sh',
                                         'rho_sh', 'gamma_sh',
                                         'kappa_sh_free',
                                         'sigma'),
                                       check_arrays=TRUE)
util$check_all_expectand_diagnostics(base_samples)


# Retrodictive check
par(mfrow=c(3, 2))
for (t in  sample(1:data$N_trees, 6)) {
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  rw_names <- sapply(idxs,
                     function(n) paste0('log_rw_pred[', n, ']'))
  util$plot_conn_pushforward_quantiles(samples, rw_names, data$years[idxs],
                                       xlab="Year", ylab="Log Ring Width Per mm", 
                                       display_ylim=c(-2, 2))
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=1.0, col="white")
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=0.8, col="black")
  #abline(v=1974, lwd=2, lty=2, col="#DDDDDD")
}

par(mfrow=c(2, 1))
pred_names <- sapply(1:data$N, function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conditional_median_quantiles(samples, pred_names, data$gdd_obs,
                                       10, 35, 1, data$log_rw_obs, 
                                       residual=FALSE,
                                       ylab = 'Marginal quantiles',
                                       xlab="GDD (x100degC)")
util$plot_conditional_median_quantiles(samples, pred_names, data$gdd_obs,
                                       10, 35, 1, data$log_rw_obs, 
                                       residual=TRUE,
                                       ylab = 'Marginal quantiles (resid.)',
                                       xlab="GDD (x100degC)")

par(mfrow=c(2, 1))
pred_names <- sapply(1:data$N, function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conditional_median_quantiles(samples, pred_names, data$sm_obs,
                                       17, 37, 1, data$log_rw_obs, 
                                       residual=FALSE,
                                       ylab = 'Marginal quantiles',
                                       xlab="Soil moisture (%)")
util$plot_conditional_median_quantiles(samples, pred_names, data$sm_obs,
                                       17, 37, 1, data$log_rw_obs, 
                                       residual=TRUE,
                                       ylab = 'Marginal quantiles (resid.)',
                                       xlab="Soil moisture (%)")

par(mfrow=c(2, 1))
pred_names <- sapply(1:data$N, function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conditional_median_quantiles(samples, pred_names, data$vpd_obs,
                                       1, 13, 1, data$log_rw_obs, 
                                       residual=FALSE,
                                       ylab = 'Marginal quantiles',
                                       xlab="VPD (hPa)")
util$plot_conditional_median_quantiles(samples, pred_names, data$vpd_obs,
                                       1, 13, 1, data$log_rw_obs, 
                                       residual=TRUE,
                                       ylab = 'Marginal quantiles (resid.)',
                                       xlab="VPD (hPa)")



# Inference
par(mfrow=c(1, 3))
util$plot_expectand_pushforward(samples[['beta_gdd']], 250,
                                flim = c(-0.1,0.1),
                                display_name="beta_gdd")
xs <- seq(-1, 1, 0.01)
ys <- dnorm(xs, 0, log(1.8) / 2.57)*50
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['beta_sm']], 250,
                                flim = c(-0.1,0.1),
                                display_name="beta_sm")
xs <- seq(-1, 1, 0.01)
ys <- dnorm(xs, 0, log(1.8) / 2.57) *10
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['beta_vpd']], 250,
                                flim = c(-0.1,0.1),
                                display_name="beta_sm")
xs <- seq(-1, 1, 0.01)
ys <- dnorm(xs, 0, log(1.8) / 2.57) *10
lines(xs, ys, lwd=2, col=util$c_light)


par(mfrow=c(1, 2))
util$plot_expectand_pushforward(samples[['rho']], 250,
                                display_name="rho")

util$plot_expectand_pushforward(samples[['rho_sh']], 250,
                                display_name="rho_sh")
