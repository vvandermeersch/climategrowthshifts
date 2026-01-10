
wd <- "/home/victor/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)

library(cmdstanr)

data <- readRDS(file.path(wd, 'output/model', 'data_28dec2025_8standsPIPO.rds'))
data$uniq_tree_ids <- NULL
data$uniq_species_ids <- NULL
data$uniq_stand_ids <- NULL
data$N_clades <- 1

model <- cmdstan_model(file.path(wd, 'model/stan', 'model10_speciesshocks_improved_nopooling.stan'),
                       cpp_options = list(stan_threads = TRUE))

fit_cmd <- model$sample(data = data, threads_per_chain = 4,
                    chains = 4, parallel_chains = 4, seed = 5838294,
                    iter_warmup = 1000, iter_sampling = 1000, refresh = 10)

fit$profiles()

samples <- util$extract_expectand_vals(fit)
base_samples <- util$filter_expectands(samples,
                                       c('alpha', 
                                         'beta_gdd', 
                                         'beta_vpd', 
                                         'beta_pre',
                                         'rho_sp', 'gamma_sp',
                                         'sigma'),
                                       check_arrays=TRUE)
util$check_all_expectand_diagnostics(base_samples)



par(mfrow=c(1, 3))
util$plot_expectand_pushforward(samples[['beta_gdd[1]']], 25,
                                flim = c(-0.05,0.05),
                                display_name="beta_gdd")
xs <- seq(-1, 1, 0.01)
ys <- dnorm(xs, 0, log(1.8) / 2.57)*50
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['beta_vpd[1]']], 25,
                                flim = c(-0.05,0.05),
                                display_name="beta_vpd")
xs <- seq(-1, 1, 0.01)
ys <- dnorm(xs, 0, log(1.8) / 2.57) *10
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['beta_pre[1]']], 25,
                                flim = c(-0.05,0.05),
                                display_name="beta_pre")
xs <- seq(-1, 1, 0.01)
ys <- dnorm(xs, 0, log(1.8) / 2.57) *10
lines(xs, ys, lwd=2, col=util$c_light)
