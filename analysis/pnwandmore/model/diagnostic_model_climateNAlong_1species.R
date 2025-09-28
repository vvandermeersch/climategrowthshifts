rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('functions/custom_functions.R', local=util)
setwd(wd)

# Model data
data <- readRDS(file.path(wd, 'output/model/climatena', 'data_11sept2025_long_PIPO.rds'))
# fit <- readRDS(file.path(wd, "output/model/climatena/fit_11sept2025_differentunits_partialpooling_2clades_centered_extended_onlyPIPO.rds"))
samples <- readRDS(file.path(wd, "output/model/climatena/samples_11sept2025_differentunits_partialpooling_2clades_centered_extended_onlyPIPO_constrainedrhos_updatedpriors.rds"))

diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit)
base_samples <- util$filter_expectands(samples,
                                       c('alpha', 
                                         'beta_gdd', 'beta_winterprec', 'beta_ffp',
                                         'rho_sp', 'gamma_sp',
                                         # 'f_tilde_sh',
                                         'rho_sh', 'kappa_sh',
                                         'sigma'),
                                       check_arrays=TRUE)
util$check_all_expectand_diagnostics(base_samples)

par(mfrow=c(1, 3), mar = c(4,4,1,1))
xs <- seq(-0.1, 0.3, 0.01)
ys <- dnorm(xs, 0, log(1.8) / 2.57)
util$plot_expectand_pushforward(samples[['beta_gdd']], 75,
                                flim = c(-0.1,0.1),
                                display_name=expression(beta[GDD]),
                                col = '#27278f')
util$plot_expectand_pushforward(samples[['beta_winterprec']], 75,
                                flim = c(-0.1,0.1),
                                display_name=expression(beta[winter_prec]),
                                col = '#27278f')
util$plot_expectand_pushforward(samples[['beta_ffp']], 75,
                                flim = c(-0.1,0.1),
                                display_name=expression(beta[FFP]),
                                col = '#27278f')


util$plot_expectand_pushforward(samples[['sigma']], 25,
                                flim = c(0,1),
                                display_name=expression(beta[GDD]),
                                col = '#27278f')

par(mfrow=c(1, 2), mar = c(4,4,1,1))
util$plot_expectand_pushforward(samples[['rho_sp']], 75,
                                flim = c(0,20),
                                display_name=expression(rho[sp]),
                                col = '#27278f')
util$plot_expectand_pushforward(samples[['rho_sh']], 75,
                                flim = c(0,10),
                                display_name=expression(rho[sp]),
                                col = '#27278f')



util$plot_pairs_by_chain(samples[['beta_gdd']], 'beta_gdd', samples[['beta_ffp']], 'beta_ffp')

util$plot_pairs_by_chain(samples[['beta_gdd']], 'beta_gdd', samples[['beta_ffp']], 'beta_ffp')
util$plot_pairs_by_chain(samples[['beta_gdd']], 'beta_gdd', samples[['rho_sp']], 'rho_sp')
util$plot_pairs_by_chain(samples[['rho_sh']], 'rho_sh', samples[['rho_sp']], 'rho_sp')
util$plot_pairs_by_chain(samples[['rho_sh']], 'rho_sh', samples[['alpha']], 'alpha')


par(mfrow=c(1, 1), mar = c(4,4,1,1))
util$plot_expectand_pushforward(samples[['sigma']], 25,
                                flim = c(0,1),
                                display_name=expression(beta[GDD]))
xs <- seq(0, 1, 0.01)
ys <- dnorm(xs, 0, 0.095 / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)


par(mfrow=c(3, 2), mar = c(4,4,1,1))
for (t in  sample(1:data$N_trees, 6)) {
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  rw_names <- sapply(idxs,
                     function(n) paste0('log_rw_pred[', n, ']'))
  util$plot_conn_pushforward_quantiles(samples, rw_names, data$years[idxs],
                                       xlab="Year", ylab="Log Ring Width Per mm", main = data$uniq_tree_ids[t],
                                       display_ylim=c(-3, 3), display_xlim = c(1900, 2024))
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=1.0, col="white")
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=0.8, col="black")
}
