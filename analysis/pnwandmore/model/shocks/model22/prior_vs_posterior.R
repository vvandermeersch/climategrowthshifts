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

base_samples <- readRDS(file.path(folder, 'model22', 'basesamples_model22_HGSP_PIPO_30stands_19502022_nofsh_nofind_alphastand.rds'))
old_base_samples <- readRDS(file.path(folder, 'model21', 'basesamples_model21_HGSP_PIPO_30stands_19502022_nomu_sametau_newrhos.rds'))
util$check_all_expectand_diagnostics(base_samples)

params <- names(base_samples)
base_params <- util$filter_expectands(base_samples,  params[!grepl('delta|thetas', params)])
util$check_all_expectand_diagnostics(base_params)

fit <- readRDS(file.path(folder, 'model22', 'fit_model22_HGSP_PIPO_30stands_19502022_nofsh_nofind_alphastand.rds'))
fit$time()
fit$diagnostic_summary()

util$plot_pairs_by_chain(base_samples[['sigma']], 'sigma',
                         base_samples[['tau_stand']], 'tau_stand')

util$plot_pairs_by_chain(base_samples[['tau_clim']], 'tau_clim',
                         base_samples[['tau_stand']], 'tau_stand')

cols <- colorRampPalette(c(util$c_light, util$c_mid, util$c_dark))(data$N_stands)
par(mfrow = c(2,2))
for(i in 1:data$N_stands){
  util$plot_expectand_pushforward(base_samples[[paste0('phi_sck[',i,']')]], 30,
                                  display_name = 'p(stand concordance)',
                                  flim = c(0,1), add = (i != 1), col = cols[i],
                                  border = paste0(cols[i], 10))
}
prior <- rbeta(1e6, 2.3, 6.07)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
for(i in 1:data$N_stands){
  util$plot_expectand_pushforward(base_samples[[paste0('omega_conc_sck[',i,']')]], 30,
                                  display_name = 'p(tree shock | stand concordance)',
                                  flim = c(0,1), add = (i != 1), col = cols[i],
                                  border = paste0(cols[i], 10))
}
prior <- rbeta(1e6, 3.48, 1.86)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
for(i in 1:data$N_stands){
  util$plot_expectand_pushforward(base_samples[[paste0('omega_shutdown[',i,']')]], 30,
                                  display_name = 'p(growth shutdown | tree shock)',
                                  flim = c(0,1), add = (i != 1), col = cols[i],
                                  border = paste0(cols[i], 10))
}
prior <- rbeta(1e6, 1.08, 1.64)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[[paste0('tau_conc')]], 30,
                                display_name = 'tau_conc',
                                flim = c(0,3))
prior <- rnorm(1e6, 0, 3 / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)





norm100 <- function(x) {
  100 * (x - min(x)) /(max(x) - min(x))
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




par(mfrow = c(1,2), mar = c(4,4,1,1))
util$plot_expectand_pushforward(base_samples[['rho_sp']], 50,
                                display_name = 'rho_sp',
                                flim = c(10,80))
util$plot_expectand_pushforward(old_base_samples[['rho_sp']], 50,
                                display_name = 'rho_sp',
                                flim = c(10,80), 
                                col = util$c_light, add = T)
abline(v = data$N_all_years, lty = 2, lwd = 2)
prior <- rlnorm(1e6, 3.7, 0.35)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(base_samples[['gamma_sp']], 50,
                                display_name = 'gamma_sp',
                                flim = c(0,log(10)))
util$plot_expectand_pushforward(old_base_samples[['gamma_sp']], 50,
                                display_name = 'gamma_sp',
                                flim = c(0,log(10)),
                                col = util$c_light, add = T)
prior <- rnorm(1e6, 0, log(10)/2.57)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
