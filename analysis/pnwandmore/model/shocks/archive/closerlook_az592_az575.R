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

data <- readRDS(file = file.path(wd, 'output', 'model', 'data_15nov2025_az592_az575.rds'))

#--------------------------------------------------#
# Starting with a full non-centering configuration # 
#--------------------------------------------------#

data$w_sck <- rep(0, data$N)
if(FALSE){
  fit <- stan(file=file.path(wd, 'model', 'stan', 'model8_with3predictors_nopooling_standphi_nonconcordantshock.stan'),
              data=data, seed=5838293,
              chains = 4,
              warmup=1000, iter=2024, refresh=10)
  saveRDS(fit, file.path(wd, 'model/shocks/output/tuning_nc_az592_az575', 'fit.rds'))
}else{
  fit <- readRDS(file.path(wd, 'model/shocks/output/tuning_nc_az592_az575', 'fit.rds'))
}

samples <- util$extract_expectand_vals(fit)
diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

util$plot_pairs_by_chain(samples[['sigma']], 'sigma', samples[['inner_tau_sck']], 'inner_tau_sck')

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
ys <- dnorm(xs, 0, log(2) / 2.57)
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
util$plot_expectand_pushforward(samples[['phi_sck[2]']], 50,
                                flim = c(0,1), 
                                display_name=expression(phi[schock]), add = T)
xs <- seq(0, 1, 0.01)
ys <- dbeta(xs, 2, 10)
lines(xs, ys, lwd=2, col=util$c_mid_teal)

util$plot_expectand_pushforward(samples[['omega_conc_sck']], 50,
                                flim = c(0,1), 
                                display_name=expression(omega[shock]^concordant))
xs <- seq(0, 5, 0.01)
ys <- dbeta(xs, 10, 2)
lines(xs, ys, lwd=2, col=util$c_mid_teal)
util$plot_expectand_pushforward(samples[['omega_nonconc_sck']], 50,
                                flim = c(0,1), 
                                display_name=expression(omega[shock]^nonconcordant))
xs <- seq(0, 5, 0.01)
ys <- dbeta(xs, 1, 20)
lines(xs, ys, lwd=2, col=util$c_mid_teal)

par(mfrow = c(2,1))
for(s in 1:data$N_stands){
  names <- sapply(1:data$N_all_years,
                  function(sp) paste0('f_sh[', s, ',', sp, ']'))
  util$plot_realizations(samples, names, plot_xs = data$all_years[1:data$N_all_years], N_plots = 50,
                         main = '', ylab = 'f_sh',
                         display_xlim = data$all_years[c(1,data$N_all_years)], display_ylim = c(-5,2))
}

par(mfrow = c(1,1))
util$plot_expectand_pushforward(samples[['rho_sh']], 50,
                                flim = c(1,5), 
                                display_name=expression(rho[short]))

par(mfrow = c(3,3))
for(t in 1:data$N_trees){
  idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- sapply(idxs_t,
                  function(x) paste0('delta_tilde_sck[', x, ']'))
  util$plot_realizations(samples, names, plot_xs = data$years[idxs_t], N_plots = 100,
                         ylab = expression(delta[shocks]),
                         display_xlim = data$all_years[c(1,length(data$all_years))], display_ylim = c(-20,10),
                         main = paste0('Tree ', t))
}

par(mfrow = c(1,1))
t <- 6 
idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
names <- sapply(idxs_t,
                function(x) paste0('delta_sck[', x, ']'))
util$plot_realizations(samples, names, plot_xs = data$years[idxs_t], N_plots = 50,
                       ylab = expression(delta[shocks]),
                       display_xlim = data$all_years[c(1,length(data$all_years))], display_ylim = c(-10,2),
                       main = paste0('Tree ', t))
abline(v = c(2002))
name_shock <- paste("delta_tilde_sck[", idxs_t[which(data$all_years == 2002)], "]", sep='')
util$plot_div_pairs(name_shock, 'inner_tau_sck', samples, diagnostics, transforms = list('sigma' = 1, 'inner_tau_sck' = 1))

#----------------------#
# Refine configuration # 
#----------------------#
name_shocks <- paste("delta_tilde_sck[", 1:data$N, "]", sep='')
median_shocks <- sapply(name_shocks, function(x) util$ensemble_mcmc_quantile_est(samples[[x]], 0.5))
hist(median_shocks, breaks = 30, ylim = c(0,50))
data$w_sck <- unname(ifelse(median_shocks < -20, 1, 0))
if(FALSE){
  fit2 <- stan(file=file.path(wd, 'model', 'stan', 'model8_with3predictors_nopooling_standphi_nonconcordantshock.stan'),
              data=data, seed=5838293,
              chains = 4,
              warmup=1000, iter=2024, refresh=10)
  saveRDS(fit2, file.path(wd, 'model/shocks/output/tuning_nc_az592_az575', 'fit2.rds'))
}else{
  fit2 <- readRDS(file.path(wd, 'model/shocks/output/tuning_nc_az592_az575', 'fit2.rds'))
}
samples2 <- util$extract_expectand_vals(fit2)
diagnostics2 <- util$extract_hmc_diagnostics(fit2)
util$check_all_hmc_diagnostics(diagnostics2)

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
ys <- dnorm(xs, 0, log(2) / 2.57)
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
util$plot_expectand_pushforward(samples[['phi_sck[2]']], 50,
                                flim = c(0,1), 
                                display_name=expression(phi[schock]), add = T)
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

par(mfrow = c(2,1))
for(s in 1:data$N_stands){
  names <- sapply(1:data$N_all_years,
                  function(sp) paste0('f_sh[', s, ',', sp, ']'))
  util$plot_realizations(samples2, names, plot_xs = data$all_years[1:data$N_all_years], N_plots = 50,
                         main = paste0('Stand ',s), ylab = 'f_sh',
                         display_xlim = data$all_years[c(1,data$N_all_years)], display_ylim = c(-5,2))
}
util$plot_expectand_pushforward(samples2[['rho_sh']], 50,
                                flim = c(1,10), 
                                display_name=expression(rho[sh]))
util$plot_div_pairs('rho_sh', 'inner_tau_sck', samples2, diagnostics2)
util$plot_div_pairs('rho_sh', 'sigma', samples2, diagnostics2)


#-----------------------------------#
# I changed the prior on tau_innner # 
#-----------------------------------#
name_shocks <- paste("delta_sck[", 1:data$N, "]", sep='')
median_shocks <- sapply(name_shocks, function(x) util$ensemble_mcmc_quantile_est(samples2[[x]], 0.5))
hist(median_shocks, breaks = 30, ylim = c(0,50))
data$w_sck <- unname(ifelse(median_shocks < -5, 1, 0))
if(FALSE){
  fit3 <- stan(file=file.path(wd, 'model', 'stan', 'model8_with3predictors_nopooling_standphi_nonconcordantshock.stan'),
               data=data, seed=5838293,
               chains = 4,
               warmup=1000, iter=2024, refresh=10)
  saveRDS(fit3, file.path(wd, 'model/shocks/output/tuning_nc_az592_az575', 'fit3.rds'))
}else{
  fit3 <- readRDS(file.path(wd, 'model/shocks/output/tuning_nc_az592_az575', 'fit3.rds'))
}
samples3 <- util$extract_expectand_vals(fit3)
diagnostics3 <- util$extract_hmc_diagnostics(fit3)

par(mfrow = c(2,1))
for(s in 1:data$N_stands){
  names <- sapply(1:data$N_all_years,
                  function(sp) paste0('f_sh[', s, ',', sp, ']'))
  util$plot_realizations(samples3, names, plot_xs = data$all_years[1:data$N_all_years], N_plots = 50,
                         main = '', ylab = 'f_sh',
                         display_xlim = data$all_years[c(1,data$N_all_years)], display_ylim = c(-5,2))
}
par(mfrow = c(1,1))
util$plot_expectand_pushforward(samples3[['rho_sh']], 50,
                                flim = c(1,10), 
                                display_name=expression(rho[sh]))
util$plot_div_pairs('rho_sh', 'inner_tau_sck', samples3, diagnostics3)
