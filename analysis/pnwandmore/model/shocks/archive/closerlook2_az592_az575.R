# rm(list = ls())

wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('functions/custom_functions.R', local=util)
setwd(wd)
library(rstan)
options(mc.cores = parallel::detectCores())

data <- readRDS(file.path(wd, 'model/shocks/output/tuning_nc_az592_az575', 'data_with_wsck.rds'))

# First fit
fit <- stan(file=file.path(wd, 'model', 'stan', 'model8_with3predictors_nopooling_standphi_nonconcordantshock.stan'),
             data=data, seed=5838293,
             chains = 4,
             warmup=1000, iter=2024, refresh=10)
samples <- util$extract_expectand_vals(fit)

# We prepare a new configuration
name_shocks <- paste("delta_tilde_sck[", 1:data$N, "]", sep='')
median_shocks <- sapply(name_shocks, function(x) util$ensemble_mcmc_quantile_est(samples[[x]], 0.5))
hist(median_shocks, breaks = 30, ylim = c(0,50))
abline(v = c(-4, 4))
data$w_sck <- unname(ifelse(abs(median_shocks) > 4, 1, 0))

# New fit
fit2 <- stan(file=file.path(wd, 'model', 'stan', 'model8_with3predictors_nopooling_standphi_nonconcordantshock.stan'),
             data=data, seed=5838293,
             chains = 4,
             warmup=1000, iter=2024, refresh=10)
samples2 <- util$extract_expectand_vals(fit2)
diagnostics2 <- util$extract_hmc_diagnostics(fit2)
util$check_all_hmc_diagnostics(diagnostics2)

# We keep previous configuration, and complement it
median_shocks2 <- sapply(name_shocks, function(x) util$ensemble_mcmc_quantile_est(samples2[[x]], 0.5))
median_shocks[510]
median_shocks2[510]
data$w_sck[abs(median_shocks2) > 4] <- 1

# New fit
fit3 <- stan(file=file.path(wd, 'model', 'stan', 'model8_with3predictors_nopooling_standphi_nonconcordantshock.stan'),
             data=data, seed=5838293,
             chains = 4,
             warmup=1000, iter=2024, refresh=10)
samples3 <- util$extract_expectand_vals(fit3)
diagnostics3 <- util$extract_hmc_diagnostics(fit3)
util$check_all_hmc_diagnostics(diagnostics3)

# What if we temper the non-centering?
data$w_sck[data$w_sck == 0] <- 0.2
fit4 <- stan(file=file.path(wd, 'model', 'stan', 'model8_with3predictors_nopooling_standphi_nonconcordantshock.stan'),
             data=data, seed=5838293,
             chains = 4,
             warmup=1000, iter=2024, refresh=10)
samples4 <- util$extract_expectand_vals(fit4)
diagnostics4 <- util$extract_hmc_diagnostics(fit4)
util$check_all_hmc_diagnostics(diagnostics4)

