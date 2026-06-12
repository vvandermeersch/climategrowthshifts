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

data <- readRDS(file = file.path(wd, 'output/model', 'data_15nov2025_azPIPO.rds'))

par(mfrow = c(1,2))
hist(data$log_rw_obs, col = util$c_mid, border = 'white', main ='')
hist(data$log_rw_obs, col = util$c_mid, border = 'white', main ='',
     ylim = c(0,500))

data$w_sck <- ifelse(data$log_rw_obs < -4, 1, 0)

# First fit
fit <- stan(file=file.path(wd, 'model', 'stan', 'model8_with3predictors_nopooling_standphi_nonconcordantshock.stan'),
            data=data, seed=5838293,
            chains = 4,
            warmup=1000, iter=2024, refresh=10)
saveRDS(fit, file.path(wd, 'model/shocks/output/azPIPO', 'fit.rds'))
samples <- util$extract_expectand_vals(fit)

# Quickly look at posteriors
nf <- layout(
  matrix(c(1,1,1,2,2,2,
           3,3,4,4,5,5), 
         ncol=6, byrow=TRUE)
)
util$plot_expectand_pushforward(samples[['inner_tau_sck']], 50,
                                flim = c(0,1.5), 
                                display_name=expression(tau[shock]^inner))
xs <- seq(0, 2, 0.01)
ys <- dnorm(xs, 0, log(1.2) / 2.57)
lines(xs, ys, lwd=2, col=util$c_mid_teal)

util$plot_expectand_pushforward(samples[['outer_tau_sck']], 50,
                                flim = c(0,10), 
                                display_name=expression(tau[shock]^outer))
xs <- seq(0, 10, 0.01)
ys <- dnorm(xs, 0, 10 / 2.57)
lines(xs, ys, lwd=2, col=util$c_mid_teal)
util$plot_expectand_pushforward(samples[['phi_sck[1]']], 50,
                                flim = c(0,1), 
                                display_name=expression(phi[schock]))
for(s in 2:data$N_stands){
  util$plot_expectand_pushforward(samples[[paste0('phi_sck[',s,']')]], 50,
                                  flim = c(0,1), 
                                  display_name=expression(phi[schock]), add = T)
}
xs <- seq(0, 1, 0.01)
ys <- dbeta(xs, 2, 10)
lines(xs, ys, lwd=2, col=util$c_mid_teal)

util$plot_expectand_pushforward(samples[['omega_conc_sck']], 50,
                                flim = c(0,1), 
                                display_name=expression(omega[shock]^concordant))
xs <- seq(0, 1, 0.01)
ys <- dbeta(xs, 230, 14)
lines(xs, ys, lwd=2, col=util$c_mid_teal)
util$plot_expectand_pushforward(samples[['omega_nonconc_sck']], 50,
                                flim = c(0,1), 
                                display_name=expression(omega[shock]^nonconcordant))
xs <- seq(0, 1, 0.01)
ys <- dbeta(xs, 1, 20)
lines(xs, ys, lwd=2, col=util$c_mid_teal)


util$plot_pairs_by_chain(samples[['sigma']], 'sigma', 
                         samples[['inner_tau_sck']], 'inner_tau_sck')


name_shocks <- paste("delta_tilde_sck[", 1:data$N, "]", sep='')
median_shocks <- sapply(name_shocks, function(x) util$ensemble_mcmc_quantile_est(samples[[x]], 0.5))
hist(median_shocks, breaks = 30)

par(mfrow = c(1,2))
hist(median_shocks, breaks = 50, col = util$c_mid, border = 'white', main ='')
hist(median_shocks, breaks = 50, col = util$c_mid, border = 'white', main ='',
     ylim = c(0,1000))

name_intermediate_shocks <- sub("\\..*$", "", names(which(round(median_shocks,1) == -3)))
par(mfrow = c(2,5))
for(i in name_intermediate_shocks){
  util$plot_expectand_pushforward(samples[[i]], 50,
                                  flim = c(-15,15), 
                                  display_name=i)
}

data$w_sck <- ifelse(abs(median_shocks) > 3, 1, 0)

fit2 <- stan(file=file.path(wd, 'model', 'stan', 'model8_with3predictors_nopooling_standphi_nonconcordantshock.stan'),
            data=data, seed=5838293,
            chains = 4,
            warmup=1000, iter=2024, refresh=10)
saveRDS(fit2, file.path(wd, 'model/shocks/output/azPIPO', 'fit2.rds'))
# 16400 seconds
diagnostics2 <- util$extract_hmc_diagnostics(fit2)
util$check_all_hmc_diagnostics(diagnostics2)
samples2 <- util$extract_expectand_vals(fit2)

util$plot_pairs_by_chain(log(samples2[['sigma']]), 'log(sigma)', 
                         log(samples2[['inner_tau_sck']]), 'log(inner_tau_sck)')
util$plot_div_pairs('sigma', 'inner_tau_sck', samples2, diagnostics2, transforms = list('sigma' = 1, 'inner_tau_sck' = 1))



# Quickly look at posteriors
nf <- layout(
  matrix(c(1,1,1,2,2,2,
           3,3,4,4,5,5), 
         ncol=6, byrow=TRUE)
)
util$plot_expectand_pushforward(samples2[['inner_tau_sck']], 50,
                                flim = c(0,1.5), 
                                display_name=expression(tau[shock]^inner))
xs <- seq(0, 2, 0.01)
ys <- dnorm(xs, 0, log(1.2) / 2.57)
lines(xs, ys, lwd=2, col=util$c_mid_teal)

util$plot_expectand_pushforward(samples2[['outer_tau_sck']], 50,
                                flim = c(0,10), 
                                display_name=expression(tau[shock]^outer))
xs <- seq(0, 10, 0.01)
ys <- dnorm(xs, 0, 10 / 2.57)
lines(xs, ys, lwd=2, col=util$c_mid_teal)
util$plot_expectand_pushforward(samples2[['phi_sck[1]']], 50,
                                flim = c(0,1), 
                                display_name=expression(phi[schock]))
for(s in 2:data$N_stands){
  util$plot_expectand_pushforward(samples2[[paste0('phi_sck[',s,']')]], 50,
                                  flim = c(0,1), 
                                  display_name=expression(phi[schock]), add = T)
}
xs <- seq(0, 1, 0.01)
ys <- dbeta(xs, 2, 10)
lines(xs, ys, lwd=2, col=util$c_mid_teal)

util$plot_expectand_pushforward(samples2[['omega_conc_sck']], 50,
                                flim = c(0,1), 
                                display_name=expression(omega[shock]^concordant))
xs <- seq(0, 1, 0.01)
ys <- dbeta(xs, 230, 14)
lines(xs, ys, lwd=2, col=util$c_mid_teal)
util$plot_expectand_pushforward(samples2[['omega_nonconc_sck']], 50,
                                flim = c(0,1), 
                                display_name=expression(omega[shock]^nonconcordant))
xs <- seq(0, 1, 0.01)
ys <- dbeta(xs, 1, 20)
lines(xs, ys, lwd=2, col=util$c_mid_teal)
