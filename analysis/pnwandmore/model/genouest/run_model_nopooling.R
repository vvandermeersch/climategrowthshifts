
rm(list = ls())
wd <- "/scratch/vvandermeersch/climategrowthshifts"
setwd(wd)

library(rstan)
options(mc.cores = 10)


datam <- readRDS(file.path(wd, 'input', 'data_11july2025.rds'))
datam$uniq_tree_ids <- NULL
datam$uniq_species_ids <- NULL
datam$uniq_stand_ids <- NULL

# file.path(wd, 'stan', 'model6_with3predictors_pnw_difflength_partialpooling_2clades_centered.stan')
fit <- stan(file=file.path(wd, 'stan', 'model6_with3predictors_pnw_difflength_nopooling.stan'),
            data=datam, seed=5838299, 
            chains = 4,
            warmup=1000, iter=2024, refresh=10)
saveRDS(fit, file = file.path(wd, 'output', 'fit_11july2025_nopooling.rds'))