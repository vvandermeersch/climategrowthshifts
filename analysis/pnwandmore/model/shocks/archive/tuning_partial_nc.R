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

data <- readRDS(file = file.path(wd, 'output', 'model', 'data_15nov2025_azca_pipo2022.rds'))

# Starting with a full non-centering configuration
data$w_sck <- rep(0, data$N)
fit <- stan(file=file.path(wd, 'model', 'stan', 'model8_with3predictors_nopooling_standphi_nonconcordantshock.stan'),
            data=data, seed=5838293,
            chains = 4,
            warmup=1000, iter=2024, refresh=10)
# Runtime: 2842 seconds
saveRDS(fit, file.path(wd, 'model/shocks/output/tuning_nc', 'fit.rds'))
samples <- util$extract_expectand_vals(fit)
diagnostics <- util$extract_hmc_diagnostics(fit)

util$plot_pairs_by_chain(samples[['sigma']], 'sigma', samples[['inner_tau_sck']], 'inner_tau_sck')
util$plot_div_pairs('sigma', 'inner_tau_sck', samples, diagnostics)

util$plot_div_pairs(paste("delta_tilde_sck[", 238, "]", sep=''), 'sigma', samples, diagnostics, transforms = list('sigma' = 1, 'inner_tau_sck' = 1))
util$plot_div_pairs(paste("delta_tilde_sck[", 237, "]", sep=''), 'inner_tau_sck', samples, diagnostics, transforms = list('sigma' = 1, 'inner_tau_sck' = 1))

samples[['sigma2+inner_tau_sck2']] <- samples[['sigma']]^2 + samples[['inner_tau_sck']]^2
util$plot_div_pairs(paste("delta_tilde_sck[", 238, "]", sep=''), 'sigma2+inner_tau_sck2', samples, diagnostics, transforms = list('sigma2+inner_tau_sck2' = 1))
util$plot_div_pairs(paste("delta_tilde_sck[", 237, "]", sep=''), 'sigma2+inner_tau_sck2', samples, diagnostics, transforms = list('sigma2+inner_tau_sck2' = 1))

names <-  paste("delta_sck[", 1:data$N, "]", sep='')
median_shocks <- sapply(names, function(x) util$ensemble_mcmc_quantile_est(samples[[x]], 0.5))
inner_tau_q99.9 <- round(util$ensemble_mcmc_quantile_est(samples[['inner_tau_sck']], 0.999),2)
hist(median_shocks, breaks = 30, main = "", col = util$c_light, border = util$c_light_highlight)
abline(v = 4*inner_tau_q99.9, col = util$c_dark_teal)
abline(v = -4*inner_tau_q99.9, col = util$c_dark_teal)
big_shocks <- ifelse(median_shocks > 4*inner_tau_q99.9 | median_shocks < -4*inner_tau_q99.9, TRUE, FALSE)
sum(big_shocks)/length(median_shocks)

# Centering what we consider big shocks!
data$w_sck <- ifelse(big_shocks, 1, 0)
fit2 <- stan(file=file.path(wd, 'model', 'stan', 'model8_with3predictors_nopooling_standphi_nonconcordantshock.stan'),
             data=data, seed=5838293,
             chains = 4,
             warmup=1000, iter=2024, refresh=10)
saveRDS(fit2, file.path(wd, 'model/shocks/output/tuning_nc', 'fit2.rds'))
samples2 <- util$extract_expectand_vals(fit2)
diagnostics2 <- util$extract_hmc_diagnostics(fit2)

util$plot_div_pairs('sigma', 'inner_tau_sck', samples2, diagnostics2)
abline(a = 0, b = 1, lty = 2)

samples2[['sigma2+inner_tau_sck2']] <- samples2[['sigma']]^2 + samples2[['inner_tau_sck']]^2
util$plot_div_pairs(paste("delta_tilde_sck[", 238, "]", sep=''), 'sigma2+inner_tau_sck2', samples2, diagnostics2, transforms = list('sigma2+inner_tau_sck2' = 1))
util$plot_div_pairs(paste("delta_tilde_sck[", 238, "]", sep=''), 'inner_tau_sck', samples2, diagnostics2, transforms = list('inner_tau_sck' = 1))

names <-  paste("delta_sck[", 1:data$N, "]", sep='')
median_shocks <- sapply(names, function(x) util$ensemble_mcmc_quantile_est(samples2[[x]], 0.5))
inner_tau_q99.9 <- round(util$ensemble_mcmc_quantile_est(samples2[['inner_tau_sck']], 0.999),2)
hist(median_shocks, breaks = 30, main = "", col = util$c_light, border = util$c_light_highlight)
abline(v = 4*inner_tau_q99.9, col = util$c_dark_teal)
abline(v = -4*inner_tau_q99.9, col = util$c_dark_teal)


big_shocks <- ifelse(median_shocks > 4*inner_tau_q99.9 | median_shocks < -4*inner_tau_q99.9, TRUE, FALSE)
medium_shocks <- ifelse(median_shocks > 3*inner_tau_q99.9 | median_shocks < -3*inner_tau_q99.9, TRUE, FALSE)

util$plot_div_pairs(paste("delta_tilde_sck[", which(medium_shocks & !big_shocks)[1:9], "]", sep=''), 'inner_tau_sck', samples2, diagnostics, transforms = list('sigma' = 1, 'inner_tau_sck' = 1))
sapply(paste("delta_sck[", which(medium_shocks & !big_shocks)[1:9], "]", sep=''), function(x) util$ensemble_mcmc_quantile_est(samples2[[x]], 0.5))

# 3. What if we partially centered the medium shocks?
data$w_sck <- ifelse(big_shocks, 1, ifelse(!big_shocks & medium_shocks, 0.2, 0))
fit3 <- stan(file=file.path(wd, 'model', 'stan', 'model8_with3predictors_nopooling_standphi_nonconcordantshock.stan'),
             data=data, seed=5838293,
             chains = 4,
             warmup=1000, iter=2024, refresh=10)
# Runtime: 3738 seconds
saveRDS(fit3, file.path(wd, 'model/shocks/output/tuning_nc', 'fit3.rds'))
samples3 <- util$extract_expectand_vals(fit3)
diagnostics3 <- util$extract_hmc_diagnostics(fit3)

util$plot_div_pairs('sigma', 'inner_tau_sck', samples3, diagnostics3)
abline(a = 0, b = 1, lty = 2)

util$plot_div_pairs(paste("delta_tilde_sck[", 2322, "]", sep=''), 'sigma', samples3, diagnostics3, transforms = list('sigma' = 1))
util$ensemble_mcmc_quantile_est(samples3[[paste("delta_sck[", 2322, "]", sep='')]], 0.5)

# 4. What if we totally centered the medium shocks?
median_shocks <- sapply(names, function(x) util$ensemble_mcmc_quantile_est(samples3[[x]], 0.5))
inner_tau_q99.9 <- round(util$ensemble_mcmc_quantile_est(samples3[['inner_tau_sck']], 0.999),2)
big_shocks <- ifelse(median_shocks > 4*inner_tau_q99.9 | median_shocks < -4*inner_tau_q99.9, TRUE, FALSE)
medium_shocks <- ifelse(median_shocks > 3*inner_tau_q99.9 | median_shocks < -3*inner_tau_q99.9, TRUE, FALSE)

data$w_sck <- ifelse(big_shocks, 1, ifelse(!big_shocks & medium_shocks, 1, 0))
data$w_sck[2322]
fit4 <- stan(file=file.path(wd, 'model', 'stan', 'model8_with3predictors_nopooling_standphi_nonconcordantshock.stan'),
             data=data, seed=5838293,
             chains = 4,
             warmup=1000, iter=2024, refresh=10)
# Runtime: 2338.39 seconds
saveRDS(fit4, file.path(wd, 'model/shocks/output/tuning_nc', 'fit4.rds'))
samples4 <- util$extract_expectand_vals(fit4)
diagnostics4 <- util$extract_hmc_diagnostics(fit4)

util$plot_div_pairs('sigma', 'inner_tau_sck', samples4, diagnostics4)

util$plot_div_pairs(paste("delta_tilde_sck[", 2322, "]", sep=''), 'inner_tau_sck', samples3, diagnostics3, transforms = list('inner_tau_sck' = 1))
util$plot_div_pairs(paste("delta_tilde_sck[", 2322, "]", sep=''), 'inner_tau_sck', samples4, diagnostics4, transforms = list('inner_tau_sck' = 1))

util$plot_div_pairs(paste("delta_tilde_sck[", sample(1:data$N,9), "]", sep=''), 'inner_tau_sck', samples4, diagnostics4, transforms = list('inner_tau_sck' = 1))
median_shocks <- sapply(names, function(x) util$ensemble_mcmc_quantile_est(samples4[[x]], 0.5))

util$ensemble_mcmc_quantile_est(samples3[[paste("delta_sck[", 1826, "]", sep='')]], 0.5)
util$ensemble_mcmc_quantile_est(samples4[[paste("delta_sck[", 1826, "]", sep='')]], 0.5)

# 5. To avoid extreme tails let's decrease the non-centering of non-shocks
median_shocks <- sapply(names, function(x) util$ensemble_mcmc_quantile_est(samples4[[x]], 0.5))
inner_tau_q99.9 <- round(util$ensemble_mcmc_quantile_est(samples4[['inner_tau_sck']], 0.999),2)
big_shocks <- ifelse(median_shocks > 4*inner_tau_q99.9 | median_shocks < -4*inner_tau_q99.9, TRUE, FALSE)
medium_shocks <- ifelse(median_shocks > 3*inner_tau_q99.9 | median_shocks < -3*inner_tau_q99.9, TRUE, FALSE)

data$w_sck <- ifelse(big_shocks, 1, ifelse(!big_shocks & medium_shocks, 1, 0.2))
fit5 <- stan(file=file.path(wd, 'model', 'stan', 'model8_with3predictors_nopooling_standphi_nonconcordantshock.stan'),
             data=data, seed=5838293,
             chains = 4,
             warmup=1000, iter=2024, refresh=10)
# 1851 seconds!
saveRDS(fit5, file.path(wd, 'model/shocks/output/tuning_nc', 'fit5.rds'))
samples5 <- util$extract_expectand_vals(fit5)
diagnostics5 <- util$extract_hmc_diagnostics(fit5)

util$plot_div_pairs('sigma', 'inner_tau_sck', samples5, diagnostics5)

util$plot_div_pairs(paste("delta_tilde_sck[", 1826, "]", sep=''), 'inner_tau_sck', samples4, diagnostics4, transforms = list('inner_tau_sck' = 1))
util$plot_div_pairs(paste("delta_tilde_sck[", 1826, "]", sep=''), 'inner_tau_sck', samples5, diagnostics5, transforms = list('inner_tau_sck' = 1))

util$plot_div_pairs(paste("delta_tilde_sck[", sample(1:data$N,9), "]", sep=''), 'inner_tau_sck', samples5, 
                    diagnostics5, transforms = list('inner_tau_sck' = 1))

util$plot_div_pairs(paste("delta_tilde_sck[", 1757, "]", sep=''), 'inner_tau_sck', samples5, 
                    diagnostics5, transforms = list('inner_tau_sck' = 1))


samples5[['sigma2+inner_tau_sck2']] <- samples5[['sigma']]^2 + samples5[['inner_tau_sck']]^2
util$plot_div_pairs(paste("delta_tilde_sck[", 1757, "]", sep=''), 'inner_tau_sck', samples5, 
                    diagnostics5, transforms = list('inner_tau_sck' = 1))
util$plot_div_pairs(paste("delta_tilde_sck[", 1757, "]", sep=''), 'sigma', samples5, 
                    diagnostics5, transforms = list('sigma' = 1))
util$plot_div_pairs(paste("delta_tilde_sck[", 1757, "]", sep=''), 'sigma2+inner_tau_sck2', samples5, 
                    diagnostics5, transforms = list('sigma2+inner_tau_sck2' = 1))

util$plot_div_pairs(paste("delta_tilde_sck[", sample(1:data$N,9), "]", sep=''), 'sigma2+inner_tau_sck2', samples5, 
                    diagnostics5, transforms = list('sigma2+inner_tau_sck2' = 1))

util$ensemble_mcmc_quantile_est(samples4[[paste("delta_sck[", 1757, "]", sep='')]], 0.5)
util$ensemble_mcmc_quantile_est(samples5[[paste("delta_sck[", 1757, "]", sep='')]], 0.5)

util$plot_div_pairs(paste("delta_tilde_sck[", 2123, "]", sep=''), 'inner_tau_sck', samples5, 
                    diagnostics5, transforms = list('inner_tau_sck' = 1))

# Look at divergences
# pdf(file.path(wd, 'model/shocks/figures', 'divergences_stanphi_nonconc_fit5.pdf'), height = 4, width = 5.5)
# util$plot_div_pairs(paste("delta_tilde_sck[", 1:data$N, "]", sep=''), 'sigma2+inner_tau_sck2', samples5, 
#                       diagnostics5, transforms = list('sigma2+inner_tau_sck2' = 1))
# dev.off()

# 6. Consider a larger definition of shocks? Like 2.57(median(tau))
names <-  paste("delta_sck[", 1:data$N, "]", sep='')
median_shocks <- sapply(names, function(x) util$ensemble_mcmc_quantile_est(samples5[[x]], 0.5))
inner_tau_q50 <- round(util$ensemble_mcmc_quantile_est(samples5[['inner_tau_sck']], 0.5),2)
shocks <- ifelse(median_shocks > 2.57*inner_tau_q50 | median_shocks < -2.57*inner_tau_q50, TRUE, FALSE)

base_samples <- util$filter_expectands(samples5,
                                       c('delta_tilde_sck'),
                                       check_arrays=TRUE)
util$check_all_expectand_diagnostics(base_samples)
# Could we use Right/left tail hat{xi} warnings? To define intermediate shocks?

data$w_sck <- ifelse(shocks, 1, 0.2)
fit6 <- stan(file=file.path(wd, 'model', 'stan', 'model8_with3predictors_nopooling_standphi_nonconcordantshock.stan'),
             data=data, seed=5838293,
             chains = 4,
             warmup=1000, iter=2024, refresh=10)
# 2523 seconds!
saveRDS(fit6, file.path(wd, 'model/shocks/output/tuning_nc', 'fit6.rds'))
samples6 <- util$extract_expectand_vals(fit6)
diagnostics6 <- util$extract_hmc_diagnostics(fit6)

util$plot_div_pairs(paste("delta_tilde_sck[", 2123, "]", sep=''), 'inner_tau_sck', samples5, 
                    diagnostics5, transforms = list('inner_tau_sck' = 1))
util$plot_div_pairs(paste("delta_tilde_sck[", 2123, "]", sep=''), 'inner_tau_sck', samples6, 
                    diagnostics6, transforms = list('inner_tau_sck' = 1))

util$plot_div_pairs('sigma', 'inner_tau_sck', samples6, diagnostics6)

util$plot_div_pairs(paste("delta_tilde_sck[", 1000, "]", sep=''), 'inner_tau_sck', samples6, 
                    diagnostics6, transforms = list('inner_tau_sck' = 1))

util$plot_div_pairs(paste("delta_tilde_sck[", 1388, "]", sep=''), 'inner_tau_sck', samples6, 
                    diagnostics6, transforms = list('inner_tau_sck' = 1), xlim = c(-5,5))

util$plot_div_pairs(paste("delta_tilde_sck[", 219, "]", sep=''), 'inner_tau_sck', samples6, 
                    diagnostics6, transforms = list('inner_tau_sck' = 1), xlim = c(-5,5))

# Did not work at all, I got more divergence
# Crude test:
data$w_sck <- ifelse(shocks, 1, 1)
fit7 <- stan(file=file.path(wd, 'model', 'stan', 'model8_with3predictors_nopooling_standphi_nonconcordantshock.stan'),
             data=data, seed=5838293,
             chains = 4,
             warmup=1000, iter=2024, refresh=10)
# 1897 seconds!
saveRDS(fit7, file.path(wd, 'model/shocks/output/tuning_nc', 'fit7.rds'))
samples7 <- util$extract_expectand_vals(fit7)
diagnostics7 <- util$extract_hmc_diagnostics(fit7)
util$check_all_hmc_diagnostics(diagnostics7)

util$plot_div_pairs('sigma', 'inner_tau_sck', samples7, diagnostics7)
abline(a = 0, b = 1)

util$plot_div_pairs(paste("delta_tilde_sck[", 1000, "]", sep=''), 'inner_tau_sck', samples7, 
                    diagnostics7, transforms = list('inner_tau_sck' = 1))
util$plot_div_pairs(paste("delta_tilde_sck[", 1757, "]", sep=''), 'sigma', samples7, 
                    diagnostics7, transforms = list('sigma' = 1), xlim = c(-2,2))

util$plot_pairs_by_chain(samples7[[paste("delta_tilde_sck[", 1757, "]", sep='')]], paste("delta_tilde_sck[", 1757, "]", sep=''), 
                         log(samples7[["inner_tau_sck"]]), 'log(inner_tau_sck)')
