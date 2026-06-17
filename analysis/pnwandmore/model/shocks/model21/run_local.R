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

model <- cmdstan_model(file.path(wd, 'model/stan/model21', 'model21_1species_fcombined.stan'),
                       cpp_options = list(stan_threads = TRUE))

# Ease warmup
# fit_mle <- model$optimize(data = data, threads = 1, seed = 487995, iter = 2e5)
# mle <- fit_mle$mle() 
# params_to_init <- c("alpha", "beta_gdd", "beta_pre", "beta_vpd",
#                     "rho_sp", "rho_ind", "rho_sh")
# 
# pos_params <- c("rho_sp", "rho_ind", "rho_sh")
# jitter_init <- function() {
#   out <- lapply(params_to_init, function(nm) {
#     x <- unname(mle[nm])
#     if (nm %in% pos_params) x * exp(rnorm(1, 0, 0.10))
#     else x + rnorm(1, 0, 0.10) 
#   })
#   setNames(out, params_to_init)
# }
# set.seed(487995)
# init_list <- lapply(1:4, function(i) jitter_init())
# 

fit <- model$sample(data = data, threads_per_chain = 2,
                    chains = 4, parallel_chains = 4, seed = 487995,
                    iter_warmup = 1500, iter_sampling = 500,
                    refresh = 10)
