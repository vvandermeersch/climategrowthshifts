rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"

library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)


data <- readRDS(file = file.path(wd, 'output/model', 'data_26dec2025_allspecies_longcoverage_improved.rds'))
data$uniq_tree_ids <- NULL
data$uniq_species_ids <- NULL
data$uniq_stand_ids <- NULL

model <- cmdstan_model(file.path(wd, 'model/stan', 
                                 'model11_speciesshocks_standpred.stan'),
                       cpp_options = list(stan_threads = TRUE))

fit <- model$sample(data = data, threads_per_chain = 1,
                    chains = 1, parallel_chains = 1, seed = 1234567,
                    iter_warmup = 1, iter_sampling = 1)



data <- readRDS(file = file.path(wd, 'output/model', 'data_26dec2025_test.rds'))
data$uniq_tree_ids <- NULL
data$uniq_species_ids <- NULL
data$uniq_stand_ids <- NULL
model <- cmdstan_model(file.path(wd, 'model/stan', 
                                 'model10_speciesshocks_improved_pooling.stan'),
                       cpp_options = list(stan_threads = TRUE))

fit2 <- model$sample(data = data, threads_per_chain = 10,
                    chains = 1, parallel_chains = 1, seed = 1234567,
                    iter_warmup = 50, iter_sampling = 50)

fit$profiles()
fit2$profiles()
