library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)


setwd('~/projects/climategrowthshifts/analysis/pnwandmore')

model <- cmdstan_model('model/stan/profiling/model9_marginalshocks_zerotauinner_censoring_onlyconcordant_heterogeneous.stan')

data <- readRDS(file = file.path('model/improve_computation/output/model', 'data_05dec2025_PIPO_5stands.rds'))
data$rw_obs <- ifelse(data$log_rw_obs == log(1e-4), 0, exp(data$log_rw_obs))

data$uniq_tree_ids <- NULL
data$uniq_stand_ids <- NULL
data$uniq_species_ids <- NULL

fit <- model$sample(data = data, chains = 1, seed = 1234567)

fit$profiles()


modpar <- cmdstan_model('model/stan/profiling/model9_marginalshocks_zerotauinner_censoring_onlyconcordant_heterogeneous_reducesum.stan',
                     cpp_options = list(stan_threads = TRUE))

# data <- readRDS(file = file.path('output', 'model', 'data_15nov2025_azPIPO.rds'))

fitpar <- modpar$sample(
  data = data,
  threads_per_chain = 5,
  chains = 1,
  parallel_chains = 1,
  seed = 1234567
)

fitpar$profiles()
