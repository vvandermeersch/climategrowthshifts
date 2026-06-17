rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(cmdstanr)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)


folder <- '/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model'

data <- readRDS(file.path(folder,'data_18may2026_PIPO_30stands_19502022.rds'))

base_samples <- readRDS(file.path(wd, 'output/model/model21', 'basesamples_model21_PIPO_30stands_19502022_fcombined.rds'))
util$check_all_expectand_diagnostics(base_samples)

# params <- names(base_samples)
# delta_clim_smp <- util$filter_expectands(base_samples,  params[grepl('delta', params)])
# util$check_all_expectand_diagnostics(delta_clim_smp)
# 
# other_smp <- util$filter_expectands(base_samples,  params[!grepl('delta', params)])
# util$check_all_expectand_diagnostics(other_smp, min_ess_hat_per_chain= 50)

fit <- readRDS(file.path(wd, 'output/model/model21/threading', 'fit_model21_PIPO_30stands_19502022.rds'))
fit$diagnostic_summary()
fit$time()



par(mfrow = c(2,2), mar = c(4,4,1,1))
util$plot_expectand_pushforward(base_samples[['alpha']], 50,
                                display_name = 'alpha',
                                flim = c(-log(10),log(10)))
prior <- rnorm(1e6, 0, log(10)/2.32)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[['beta_gdd']], 50,
                                display_name = 'beta_gdd',
                                flim = c(-log(1.8),log(1.8)))
prior <- rnorm(1e6, 0, log(1.8)/2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[['beta_pre']], 50,
                                display_name = 'beta_pre',
                                flim = c(-log(1.8),log(1.8)))
prior <- rnorm(1e6, 0, log(1.8)/2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[['beta_vpd']], 50,
                                display_name = 'beta_vpd',
                                flim = c(-log(1.8),log(1.8)))
prior <- rnorm(1e6, 0, log(1.8)/2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)


par(mfrow = c(3,2), mar = c(4,4,1,1))
util$plot_expectand_pushforward(base_samples[['rho_sp']], 50,
                                display_name = 'rho_sp',
                                flim = c(10,80))
abline(v = data$N_all_years, lty = 2, lwd = 2)
prior <- rlnorm(1e6, 3.7, 0.35)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[['gamma_sp']], 50,
                                display_name = 'gamma_sp',
                                flim = c(0,log(10)))
prior <- rnorm(1e6, 0, log(10)/2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[['rho_sh']], 50,
                                display_name = 'rho_sh',
                                flim = c(10,80))
abline(v = data$N_all_years, lty = 2, lwd = 2)
prior <- rlnorm(1e6, 3.0, 0.42)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[['gamma_sh']], 30,
                                display_name = 'gamma_sh',
                                flim = c(0,log(3)))
prior <- rnorm(1e6, 0, log(3)/2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[['rho_ind']], 30,
                                display_name = 'rho_ind',
                                flim = c(1,10))
abline(v = data$N_all_years, lty = 2, lwd = 2)
prior <- rlnorm(1e6, 1.3, 0.38)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[['gamma_ind']], 30,
                                display_name = 'gamma_ind',
                                flim = c(0,log(3)))
prior <- rnorm(1e6, 0, log(3)/2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)

par(mfrow = c(1,2))
util$plot_expectand_pushforward(base_samples[['tau_clim']], 50,
                                display_name = 'tau_clim',
                                flim = c(0,log(2)))
prior <- rnorm(1e6, 0, log(2)/2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[['sigma']], 50,
                                display_name = 'sigma',
                                flim = c(0,log(1.3)))
prior <- rnorm(1e6, 0, log(1.3)/2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)

cols <- colorRampPalette(c(util$c_light, util$c_mid, util$c_dark))(data$N_stands)
par(mfrow = c(3,2))
for(i in 1:data$N_stands){
  util$plot_expectand_pushforward(base_samples[[paste0('phi_sck[',i,']')]], 30,
                                  display_name = 'p(stand concordance)',
                                  flim = c(0,1), add = (i != 1), col = cols[i])
}
prior <- rbeta(1e6, 2.3, 6.07)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
for(i in 1:data$N_stands){
  util$plot_expectand_pushforward(base_samples[[paste0('omega_conc_sck[',i,']')]], 30,
                                  display_name = 'p(tree shock | stand concordance)',
                                  flim = c(0,1), add = (i != 1), col = cols[i])
}
prior <- rbeta(1e6, 3.48, 1.86)
# prior <- rbeta(1e6, 4.24, 2.80)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
for(i in 1:data$N_stands){
  util$plot_expectand_pushforward(base_samples[[paste0('omega_shutdown[',i,']')]], 30,
                                  display_name = 'p(growth shutdown | tree shock)',
                                  flim = c(0,1), add = (i != 1), col = cols[i])
}
# prior <- rbeta(1e6, 1.26, 5.97)
prior <- rbeta(1e6, 1.08, 1.64)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
plot.new()
util$plot_expectand_pushforward(base_samples[[paste0('tau_conc')]], 30,
                                display_name = 'tau_conc',
                                flim = c(0,3))
prior <- rnorm(1e6, 0, 3 / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
# for(i in 1:data$N_stands){
#   util$plot_expectand_pushforward(base_samples[[paste0('mu_conc[',i,']')]], 30,
#                                   display_name = 'mu_conc',
#                                   flim = c(0,3), add = (i != 1), col = cols[i])
# }
# prior <- rlnorm(1e6, -0.65, 0.89)
# prior <- rnorm(1e6, 1.5, 1.5/2.57)
# lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
# for(i in 1:data$N_stands){
#   util$plot_expectand_pushforward(base_samples[[paste0('tau_conc[',i,']')]], 30,
#                                   display_name = 'tau_conc',
#                                   flim = c(0,3), add = (i != 1), col = cols[i])
# }
# util$plot_expectand_pushforward(base_samples[[paste0('tau_conc')]], 30,
#                                   display_name = 'tau_conc',
#                                   flim = c(0,1.5))


par(mfrow = c(2,2), mar = c(4,4,1,1))
util$plot_disc_pushforward_quantiles(base_samples, paste0('phi_sck[',1:data$N_stands,']'), display_ylim = c(0,1),
                                     ylab = "p(stand concordance)")
prior <- rbeta(1e6, 2.3, 6.07)
lines(density(prior)$y, density(prior)$x, col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_disc_pushforward_quantiles(base_samples, paste0('omega_conc_sck[',1:data$N_stands,']'), display_ylim = c(0,1),
                                     ylab = "p(tree shock | stand concordance)")
prior <- rbeta(1e6, 4.24,  2.80)
lines(density(prior)$y, density(prior)$x, col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_disc_pushforward_quantiles(base_samples, paste0('omega_shutdown[',1:data$N_stands,']'), display_ylim = c(0,1),
                                     ylab = "p(growth shutdown | tree shock)")
prior <- rbeta(1e6, 1.66, 6.86)
lines(density(prior)$y, density(prior)$x, col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[[paste0('tau_conc')]], 60,
                                display_name = 'tau_conc',
                                flim = c(0,3))
prior <- rnorm(1e6, 0, 3 / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)


norm100 <- function(x) {
  200 * (x - min(x)) /(max(x) - min(x))
}

par(mfrow = c(2,2))
thetas_baseline <- c(100, 5, 0.5)
omega_thetas <- 4
alphas = thetas_baseline / omega_thetas + c(1, 1, 1)
thetas_idio <- gtools::rdirichlet(1e6, alphas)
util$plot_disc_pushforward_quantiles(base_samples, paste0('thetas_idio[',1:data$N_trees,',1]'), display_ylim = c(0,1),
                                     ylab = "p(no idio. shock)")
lines(norm100(density(thetas_idio[,1])$y), density(thetas_idio[,1])$x, col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_disc_pushforward_quantiles(base_samples, paste0('thetas_idio[',1:data$N_trees,',2]'), display_ylim = c(0,1),
                                     ylab = "p(idio. perturb.)")
lines(norm100(density(thetas_idio[,2])$y), density(thetas_idio[,2])$x, col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_disc_pushforward_quantiles(base_samples, paste0('thetas_idio[',1:data$N_trees,',3]'), display_ylim = c(0,1),
                                     ylab = "p(idio. shutdown)")
lines(norm100(density(thetas_idio[,3])$y), density(thetas_idio[,3])$x, col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[['tau_idio']], 60,
                                display_name = 'tau_idio',
                                flim = c(0,1.5))
prior <- rnorm(1e6, 0, log(5) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)



util$plot_pairs_by_chain(base_samples[['rho_ind']], 'rho_ind',
                         base_samples[['gamma_ind']], 'gamma_ind')

util$plot_pairs_by_chain(base_samples[['rho_ind']], 'rho_ind',
                         base_samples[['rho_sp']], 'rho_sp')

util$plot_pairs_by_chain(base_samples[['rho_ind']], 'rho_ind',
                         base_samples[['sigma']], 'sigma')

util$plot_pairs_by_chain(base_samples[['rho_ind']], 'rho_ind',
                         base_samples[['rho_sh']], 'rho_sh')

util$plot_pairs_by_chain(base_samples[['rho_ind']], 'rho_ind',
                         base_samples[['gamma_sp']], 'gamma_sp')

par(mfrow = c(1,1), mar = c(4,4,1,1))
util$plot_disc_pushforward_quantiles(base_samples, paste0('alpha_stand[',1:data$N_stands,']'), display_ylim = c(-2,2), ylab = bquote(alpha[stand]), xlab = "Stand")



util$plot_pairs_by_chain(base_samples[['gamma_ind']], 'gamma_ind',
                         base_samples[['tau_clim']], 'tau_clim')

util$plot_pairs_by_chain(base_samples[['gamma_sp']], 'gamma_sp',
                         base_samples[['tau_clim']], 'tau_clim')

util$plot_pairs_by_chain(base_samples[['gamma_sh']], 'gamma_sh',
                         base_samples[['tau_clim']], 'tau_clim')

