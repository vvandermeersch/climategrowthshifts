rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(cmdstanr)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)


library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
options(cmdstanr_output_dir = file.path(wd, 'output/model/tmp'))

data <- readRDS(file.path(wd, 'output/model', 'data_18may2026_PSME_30stands_19802024.rds'))
data$uniq_tree_ids <- NULL
data$uniq_species_ids <- NULL
data$uniq_stand_ids <- NULL
data$N_clades <- 1
data$uniq_stand_lat <- NULL
data$uniq_stand_lon <- NULL

data$grainsize <-  ceiling(data$N_stands/2)

model <- cmdstan_model(file.path(wd, 'model/stan/model22bis/hgsp', 'model22_HSGP_1species_nofsh_nofind_alphastand.stan'),
                       cpp_options = list(stan_threads = TRUE))

fit <- model$sample(data = data, threads_per_chain = 1,
                    chains = 4, parallel_chains = 4, seed = 5838293,
                    refresh = 10, iter_warmup = 100, iter_sampling = 100)

# 6823.3 with 4 within chain cores



samples <- fit$draws(c('alpha',
                       'beta_gdd', 'beta_vpd', 'beta_pre',
                       
                       'delta_clim', 'tau_clim',
                       
                       'rho_sp', 'gamma_sp',
                       
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

util$plot_pairs_by_chain(samples[['sigma']], 'sigma',
                         samples[['rho_sp']], 'rho_sp')



fit_mle <- model$optimize(
  data     = data,
  threads = 1,
  seed     = 123,
  iter     = 20000
)


mle <- fit_mle$mle() 

params_to_init <- c("alpha", "beta_gdd", "beta_pre", "beta_vpd",
             "tau_stand", "tau_clim", "tau_idio", "sigma",
             "rho_sp", "gamma_sp")

pos_params <- c("tau_stand", "tau_clim", "tau_idio", "sigma", "rho_sp", "gamma_sp")

jitter_init <- function() {
  out <- lapply(params_to_init, function(nm) {
    x <- unname(mle[nm])
    if (nm %in% pos_params) x * exp(rnorm(1, 0, 0.10))
    else x + rnorm(1, 0, 0.10) 
  })
  setNames(out, params_to_init)
}

set.seed(123)
init_list <- lapply(1:4, function(i) jitter_init())

fit <- model$sample(data = data, threads_per_chain = 1,
                    chains = 4, parallel_chains = 4, seed = 5838293,
                    init = init_list,
                    refresh = 10, iter_warmup = 1000, iter_sampling = 1000)
