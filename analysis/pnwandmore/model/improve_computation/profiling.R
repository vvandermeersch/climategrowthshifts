library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)


setwd('~/projects/climategrowthshifts/analysis/pnwandmore')

model <- cmdstan_model('model/stan/profiling/model9_marginalshocks_zerotauinner_censoring_onlyconcordant_heterogeneous.stan')

data <- readRDS(file = file.path('output', 'model', 'data_15nov2025_az592_az628.rds'))
data$rw_obs <- ifelse(data$log_rw_obs == log(1e-4), 0, exp(data$log_rw_obs))

data$uniq_tree_ids <- NULL
data$uniq_stand_ids <- NULL
data$uniq_species_ids <- NULL

fit <- model$sample(data = data, chains = 1)
