rm(list = ls())
wd <- "/home/victor/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)

library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
options(cmdstanr_output_dir = file.path(wd, 'output/model/tmp'))

data <- readRDS(file.path(wd, 'output/model', 'data_18may2026_PIPO_10stands_19202021.rds'))
data$uniq_tree_ids <- NULL
data$uniq_species_ids <- NULL
data$uniq_stand_ids <- NULL
data$N_clades <- 1
data$uniq_stand_lat <- NULL
data$uniq_stand_lon <- NULL

data$grainsize <-  1

model <- cmdstan_model(file.path(wd, 'model/stan/model21bis', 'model21_1species_fcombined_nofind.stan'),
                       cpp_options = list(stan_threads = TRUE))

fit <- model$sample(data = data, threads_per_chain = 1,
                    chains = 4, parallel_chains = 4, seed = 487995,
                    iter_warmup = 1500, iter_sampling = 500,
                    refresh = 10)

fit$time()
fit$diagnostic_summary()
rds_file <- file.path(wd, 'output/model/model21bis', 'fit_model21bis_PIPO_10stands_19202021.rds')
fit$save_object(file = rds_file)

samples <- fit$draws(c('alpha',
                       'beta_gdd', 'beta_vpd', 'beta_pre',
                       
                       'delta_clim', 
                       'tau_clim',
                       
                       'rho_merged', 'gamma_merged',
                       
                       'rho_sh', 'gamma_sh', 
                       'f_tilde_sh','f_sh',
                       
                       'tau_conc',
                       'phi_sck',
                       'omega_conc_sck',
                       'omega_shutdown',
                       
                       'tau_idio',
                       'thetas_idio',
                       
                       'sigma'
))
names <- dimnames(samples)$variable
samples <- lapply(1:dim(samples)[3],
                  function(k){t(matrix(samples[1:dim(samples)[1],1:dim(samples)[2],k],
                                       nrow = dim(samples)[1], ncol = dim(samples)[2]))})
names(samples) <- names
util$check_all_expectand_diagnostics(samples)
rds_file <- file.path(wd, 'output/model/model21bis', 'basesamples_model21bis_PIPO_10stands_19202021.rds')
saveRDS(samples, file = rds_file)
