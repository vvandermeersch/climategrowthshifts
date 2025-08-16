rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)

# Model data
# data <- readRDS(file.path(wd, 'output/model', 'data_11july2025.rds'))

# Posterior quantification
fit <- readRDS(file.path(wd, 'output/model', 'fit_11july2025_nopooling.rds')) # run on Margot
gc()

samples <- util$extract_expectand_vals(fit)
rm(fit);gc()
base_samples <- util$filter_expectands(samples,
                                       c('beta_gdd', 'beta_sm', 'beta_vpd'),
                                       check_arrays=TRUE)
rm(samples);gc()
saveRDS(base_samples, file.path(wd, 'output/model/shrinkage', 'slopes_11july2025_nopooling.rds'))
