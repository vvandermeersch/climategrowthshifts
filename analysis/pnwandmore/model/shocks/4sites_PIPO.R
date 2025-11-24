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

data <- readRDS(file = file.path(wd, 'output/model', 'data_15nov2025_az592_az575_az585_az598.rds'))
data$w_sck <- rep(0, data$N)

# First fit
fit <- stan(file=file.path(wd, 'model', 'stan', 'model8_with3predictors_nopooling_standphi_nonconcordantshock.stan'),
            data=data, seed=5838293,
            chains = 4,
            warmup=1000, iter=2024, refresh=10)
saveRDS(fit, file.path(wd, 'model/shocks/output/4sites_PIPO', 'fit.rds'))
samples <- util$extract_expectand_vals(fit)
diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)
util$plot_div_pairs('sigma', 'inner_tau_sck', samples, diagnostics, transforms = list('sigma' = 1, 'inner_tau_sck' = 1))

# Look at shocks!
name_shocks <- paste("delta_tilde_sck[", 1:data$N, "]", sep='')
median_shocks <- sapply(name_shocks, function(x) util$ensemble_mcmc_quantile_est(samples[[x]], 0.5))
hist(median_shocks, breaks = 30, ylim = c(0,50))

# Ok, let's take +/-10 as a threshold to configure the parametrisation for our shocks!
data$w_sck <- unname(ifelse(abs(median_shocks) > 10, 1, 0))
sum(data$w_sck)

# Second fit
fit2 <- stan(file=file.path(wd, 'model', 'stan', 'model8_with3predictors_nopooling_standphi_nonconcordantshock.stan'),
            data=data, seed=5838293,
            chains = 4,
            warmup=1000, iter=2024, refresh=10)
saveRDS(fit2, file.path(wd, 'model/shocks/output/4sites_PIPO', 'fit2.rds'))
samples2 <- util$extract_expectand_vals(fit2)
diagnostics2 <- util$extract_hmc_diagnostics(fit2)
util$check_all_hmc_diagnostics(diagnostics2)
util$plot_div_pairs('sigma', 'inner_tau_sck', samples2, diagnostics2, transforms = list('sigma' = 1, 'inner_tau_sck' = 1))

# Look at shocks!
name_shocks <- paste("delta_tilde_sck[", 1:data$N, "]", sep='')
median_shocks2 <- sapply(name_shocks, function(x) util$ensemble_mcmc_quantile_est(samples2[[x]], 0.5))
hist(median_shocks2, breaks = 30, ylim = c(0,50))

# Ok, let's take +/-4 as a threshold to configure the parametrisation for our shocks!
data$w_sck <- unname(ifelse(abs(median_shocks2) > 4, 1, 0))
sum(data$w_sck)

# Third fit
fit3 <- stan(file=file.path(wd, 'model', 'stan', 'model8_with3predictors_nopooling_standphi_nonconcordantshock.stan'),
             data=data, seed=5838293,
             chains = 4,
             warmup=1000, iter=2024, refresh=10)
saveRDS(fit3, file.path(wd, 'model/shocks/output/4sites_PIPO', 'fit3.rds'))
samples3 <- util$extract_expectand_vals(fit3)
diagnostics3 <- util$extract_hmc_diagnostics(fit3)
util$check_all_hmc_diagnostics(diagnostics3)
util$plot_div_pairs('sigma', 'inner_tau_sck', samples3, diagnostics3, transforms = list('sigma' = 1, 'inner_tau_sck' = 1))

# Look at shocks!
name_shocks <- paste("delta_tilde_sck[", 1:data$N, "]", sep='')
median_shocks3 <- sapply(name_shocks, function(x) util$ensemble_mcmc_quantile_est(samples3[[x]], 0.5))
hist(median_shocks3, breaks = 30, ylim = c(0,50))

# Look at trees
par(mfrow = c(3,3))
for(t in 1:data$N_trees){
  idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- sapply(idxs_t,
                  function(x) paste0('delta_sck[', x, ']'))
  util$plot_realizations(samples3, names, plot_xs = data$years[idxs_t], N_plots = 100,
                         ylab = expression(delta[shocks]),
                         display_xlim = data$all_years[c(1,length(data$all_years))], display_ylim = c(-10,10),
                         main = paste0('Tree ', t))
}

# Look at stands
par(mfrow = c(2,2))
for(s in 1:data$N_stands){
  names <- sapply(1:data$N_all_years,
                  function(y) paste0('f_sh[', s, ',', y, ']'))
  util$plot_realizations(samples3, names, plot_xs = data$all_years[1:data$N_all_years], N_plots = 50,
                         main = '', ylab = 'f_sh',
                         display_xlim = data$all_years[c(1,data$N_all_years)], display_ylim = c(-5,2))
}
abline(v = 2002, lty = 2)

# What is happening in 2002 on stand 4?
name_y <- paste0('f_sh[', 4, ',',  which(data$all_years == 2002), ']')
t <- 81 
idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
name_x <-  paste0('delta_sck[', idxs_t[which(data$years[idxs_t] ==2002)], ']') 
util$plot_pairs_by_chain(samples3[[name_x]], name_x, 
                         samples3[[name_y]], name_y)
data$w_sck[2414] # Was not considered as a shock

# Quickly look at posteriors
nf <- layout(
  matrix(c(1,1,1,2,2,2,
           3,3,4,4,5,5), 
         ncol=6, byrow=TRUE)
)
util$plot_expectand_pushforward(samples3[['inner_tau_sck']], 50,
                                flim = c(0,1.5), 
                                display_name=expression(tau[shock]^inner))
xs <- seq(0, 2, 0.01)
ys <- dnorm(xs, 0, 2 / 2.57)
lines(xs, ys, lwd=2, col=util$c_mid_teal)

util$plot_expectand_pushforward(samples3[['outer_tau_sck']], 50,
                                flim = c(0,10), 
                                display_name=expression(tau[shock]^outer))
xs <- seq(0, 10, 0.01)
ys <- dnorm(xs, 0, 10 / 2.57)
lines(xs, ys, lwd=2, col=util$c_mid_teal)
util$plot_expectand_pushforward(samples3[['phi_sck[1]']], 50,
                                flim = c(0,1), 
                                display_name=expression(phi[schock]))
util$plot_expectand_pushforward(samples3[['phi_sck[2]']], 50,
                                flim = c(0,1), 
                                display_name=expression(phi[schock]), add = T)
util$plot_expectand_pushforward(samples3[['phi_sck[3]']], 50,
                                flim = c(0,1), 
                                display_name=expression(phi[schock]), add = T)
util$plot_expectand_pushforward(samples3[['phi_sck[4]']], 50,
                                flim = c(0,1), 
                                display_name=expression(phi[schock]), add = T)
xs <- seq(0, 1, 0.01)
ys <- dbeta(xs, 2, 10)
lines(xs, ys, lwd=2, col=util$c_mid_teal)

util$plot_expectand_pushforward(samples3[['omega_conc_sck']], 50,
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
