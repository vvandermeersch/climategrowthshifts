rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"

library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)

data <- readRDS(file = file.path(wd, 'output/model', 'data_14dec2025_PIPO_5stands.rds'))
data$uniq_tree_ids <- NULL
data$uniq_species_ids <- NULL
data$uniq_stand_ids <- NULL

model <- cmdstan_model(file.path(wd, 'model/stan', 
                                 'model10_speciesshocks_improved.stan'),
                       cpp_options = list(stan_threads = TRUE))

fit <- model$sample(data = data, threads_per_chain = 2,
                    chains = 1, parallel_chains = 1, seed = 1234567)
fit$profiles()

stanfit <- rstan::read_stan_csv(fit$output_files())

# 1418.6 thred 1
# 971.9 thread 2
# 1101 thread 3
# 757.4, 1081.6, 1010.9 thread 4
# 1054.3 thread 5


model <- cmdstan_model(file.path(wd, 'model/stan/profiling', 
                                 'model9_marginalshocks_zerotauinner_censoring_onlyconcordant_heterogeneous_reducesum_improved_v2.stan'),
                       cpp_options = list(stan_threads = TRUE))

fit2 <- model$sample(data = data, threads_per_chain = 2,
                    chains = 1, parallel_chains = 1, seed = 1234567)
fit2$profiles()
