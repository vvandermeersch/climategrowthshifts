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

data <- readRDS(file = file.path(wd, 'output', 'model', 'data_15nov2025_az592.rds'))

data$N_stand_trees <- array(data$N_stand_trees, dim = 1)
data$N_stand_years <- array(data$N_stand_years, dim = 1)
data$stand_start_years_idxs <- array(data$stand_start_years_idxs, dim = 1)
data$stand_trees_end_idxs <- array(data$stand_trees_end_idxs, dim = 1)
data$stand_trees_start_idxs <- array(data$stand_trees_start_idxs, dim = 1)


#--------------------------------------------------#
# Starting with a full non-centering configuration # 
#--------------------------------------------------#

data$w_sck <- rep(0, data$N)
# fit <- stan(file=file.path(wd, 'model', 'stan', 'model8_with3predictors_nopooling_standphi_nonconcordantshock.stan'),
#             data=data, seed=5838293,
#             chains = 4,
#             warmup=1000, iter=2024, refresh=10)
# saveRDS(fit, file.path(wd, 'model/shocks/output/tuning_nc_az592', 'fit.rds'))
fit <- readRDS(file.path(wd, 'model/shocks/output/tuning_nc_az592', 'fit.rds'))
samples <- util$extract_expectand_vals(fit)
diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

util$plot_pairs_by_chain(samples[['sigma']], 'sigma', samples[['inner_tau_sck']], 'inner_tau_sck')
util$plot_div_pairs('sigma', 'inner_tau_sck', samples, diagnostics)

util$plot_div_pairs(paste("delta_tilde_sck[", 54, "]", sep=''), 'inner_tau_sck', samples, diagnostics, transforms = list('sigma' = 1, 'inner_tau_sck' = 1))
util$plot_div_pairs(paste("delta_tilde_sck[", 53, "]", sep=''), 'inner_tau_sck', samples, diagnostics, transforms = list('sigma' = 1, 'inner_tau_sck' = 1))

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
                                flim = c(0,5), 
                                display_name=expression(tau[shock]^outer))
xs <- seq(0, 5, 0.01)
ys <- dnorm(xs, 0, log(10) / 2.57)
lines(xs, ys, lwd=2, col=util$c_mid_teal)
util$plot_expectand_pushforward(samples[['phi_sck[1]']], 50,
                                flim = c(0,1), 
                                display_name=expression(phi[schock]))
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

# Look at one tree
par(mfrow = c(3,1))
t <- 1
idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
names <- sapply(idxs_t,
                function(x) paste0('delta_sck[', x, ']'))
util$plot_realizations(samples, names, plot_xs = data$years[idxs_t], N_plots = 50,
                       main = '', ylab = expression(delta[shocks]),
                       display_xlim = data$all_years[c(1,length(data$all_years))], display_ylim = c(-7,2))
name_shock <- paste("delta_tilde_sck[", idxs_t[which(data$years[idxs_t]==2002)], "]", sep='')
util$plot_div_pairs(name_shock, 'inner_tau_sck', samples, diagnostics, transforms = list('sigma' = 1, 'inner_tau_sck' = 1))
util$plot_pairs_by_chain(samples[[name_shock]], name_shock, log(samples[['inner_tau_sck']]), 'log(inner_tau_sck)')
util$plot_pairs_by_chain(samples[['outer_tau_sck']], 'outer_tau_sck', samples[['inner_tau_sck']], 'inner_tau_sck')

par(mfrow = c(3,1))
t <- 2
idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
names <- sapply(idxs_t,
                function(x) paste0('delta_sck[', x, ']'))
util$plot_realizations(samples, names, plot_xs = data$years[idxs_t], N_plots = 50,
                       main = '', ylab = expression(delta[shocks]),
                       display_xlim = data$all_years[c(1,length(data$all_years))], display_ylim = c(-10,2))
name_shock <- paste("delta_tilde_sck[", idxs_t[which(data$years[idxs_t]==2002)], "]", sep='')
util$plot_div_pairs(name_shock, 'inner_tau_sck', samples, diagnostics, transforms = list('sigma' = 1, 'inner_tau_sck' = 1))
util$plot_pairs_by_chain(samples[[name_shock]], name_shock, log(samples[['inner_tau_sck']]), 'log(inner_tau_sck)')
util$plot_pairs_by_chain(samples[['outer_tau_sck']], 'outer_tau_sck', samples[['inner_tau_sck']], 'inner_tau_sck')

par(mfrow = c(3,1))
t <- 3
idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
names <- sapply(idxs_t,
                function(x) paste0('delta_sck[', x, ']'))
util$plot_realizations(samples, names, plot_xs = data$years[idxs_t], N_plots = 50,
                       main = '', ylab = expression(delta[shocks]),
                       display_xlim = data$all_years[c(1,length(data$all_years))], display_ylim = c(-10,2))
name_shock <- paste("delta_sck[", idxs_t[which(data$years[idxs_t]==2002)], "]", sep='')
util$plot_div_pairs(name_shock, 'inner_tau_sck', samples, diagnostics, transforms = list('sigma' = 1, 'inner_tau_sck' = 1))
util$plot_pairs_by_chain(samples[[name_shock]], name_shock, log(samples[['inner_tau_sck']]), 'log(inner_tau_sck)')
util$plot_pairs_by_chain(samples[[name_shock]], name_shock, log(samples[['outer_tau_sck']]), 'log(outer_tau_sck)')
util$plot_pairs_by_chain(samples[[name_shock]], name_shock, samples[['outer_tau_sck']], 'outer_tau_sck')

par(mfrow = c(1,2))
s <- 1
names <- sapply(1:data$N_all_years,
                function(sp) paste0('f_sh[', s, ',', sp, ']'))
util$plot_realizations(samples, names, plot_xs = data$all_years[1:data$N_all_years], N_plots = 50,
                       main = '', ylab = 'f_sh',
                       display_xlim = data$all_years[c(1,data$N_all_years)], display_ylim = c(-5,2))
util$plot_pairs_by_chain(samples[[names[data$all_years == 2002]]], names[data$all_years == 2002], samples[['outer_tau_sck']], 'outer_tau_sck')

par(mfrow = c(7,3))
for(t in 1:data$N_trees){
  idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- sapply(idxs_t,
                  function(x) paste0('delta_sck[', x, ']'))
  util$plot_realizations(samples, names, plot_xs = data$years[idxs_t], N_plots = 50,
                         ylab = expression(delta[shocks]),
                         display_xlim = data$all_years[c(1,length(data$all_years))], display_ylim = c(-10,2),
                         main = paste0('Tree ', t))
}


#------------------------------------------------------#
# Refine configuration and increase prior of tau_outer # 
#------------------------------------------------------#
for(t in 1:data$N_trees){
  print(t)
  idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  if(t %in% c(1,14,15,20)){
    data$w_sck[idxs_t] <- 0 # no shock for these trees
    next
  }else{
    print("2002")
    data$w_sck[idxs_t[which(data$all_years == 2002)]] <- 1 # 2002 is shared across most trees
  }
  
  if(t %in% c(12,19)){
    print("2006")
    data$w_sck[idxs_t[which(data$all_years == 2006)]] <- 1 # another shock in 2006
  }
}

# fit2 <- stan(file=file.path(wd, 'model', 'stan', 'model8_with3predictors_nopooling_standphi_nonconcordantshock.stan'),
#             data=data, seed=5838293,
#             chains = 4,
#             warmup=1000, iter=2024, refresh=10)
# saveRDS(fit2, file.path(wd, 'model/shocks/output/tuning_nc_az592', 'fit2.rds'))
fit2 <- readRDS(file.path(wd, 'model/shocks/output/tuning_nc_az592', 'fit2.rds'))
samples2 <- util$extract_expectand_vals(fit2)
diagnostics2 <- util$extract_hmc_diagnostics(fit2)
util$check_all_hmc_diagnostics(diagnostics2)

util$plot_pairs_by_chain(samples2[['sigma']], 'sigma', samples2[['inner_tau_sck']], 'inner_tau_sck')
util$plot_div_pairs('sigma', 'inner_tau_sck', samples2, diagnostics2)
abline(a = 0, b = 1)

util$plot_div_pairs(paste("delta_tilde_sck[", 54, "]", sep=''), 'inner_tau_sck', samples2, diagnostics2, transforms = list('sigma' = 1, 'inner_tau_sck' = 1))
util$plot_div_pairs(paste("delta_tilde_sck[", 53, "]", sep=''), 'inner_tau_sck', samples2, diagnostics2, transforms = list('sigma' = 1, 'inner_tau_sck' = 1))

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
ys <- dnorm(xs, 0, log(150) / 2.57)
lines(xs, ys, lwd=2, col=util$c_mid_teal)
util$plot_expectand_pushforward(samples2[['phi_sck[1]']], 50,
                                flim = c(0,1), 
                                display_name=expression(phi[schock]))
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

par(mfrow = c(1,1))
util$plot_pairs_by_chain(samples2[['outer_tau_sck']], 'outer_tau_sck', samples2[['omega_conc_sck']], 'omega_conc_sck')
util$plot_pairs_by_chain(samples2[['inner_tau_sck']], 'inner_tau_sck', samples2[['omega_conc_sck']], 'omega_conc_sck')

par(mfrow = c(1,2))
s <- 1
names <- sapply(1:data$N_all_years,
                function(sp) paste0('f_sh[', s, ',', sp, ']'))
util$plot_realizations(samples2, names, plot_xs = data$all_years[1:data$N_all_years], N_plots = 50,
                       main = '', ylab = 'f_sh',
                       display_xlim = data$all_years[c(1,data$N_all_years)], display_ylim = c(-5,2))
util$plot_pairs_by_chain(samples2[[names[data$all_years == 2002]]], names[data$all_years == 2002], samples2[['outer_tau_sck']], 'outer_tau_sck')

name_shocks <- paste("delta_sck[", which(data$w_sck == 1), "]", sep='')
util$plot_div_pairs(name_shocks, 'inner_tau_sck', samples2, diagnostics2, transforms = list('sigma' = 1, 'inner_tau_sck' = 1))
util$plot_div_pairs(name_shocks, 'sigma', samples2, diagnostics2, transforms = list('sigma' = 1, 'inner_tau_sck' = 1))

par(mfrow = c(7,3))
for(t in 1:data$N_trees){
  idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- sapply(idxs_t,
                  function(x) paste0('delta_sck[', x, ']'))
  util$plot_realizations(samples2, names, plot_xs = data$years[idxs_t], N_plots = 50,
                         ylab = expression(delta[shocks]),
                         display_xlim = data$all_years[c(1,length(data$all_years))], display_ylim = c(-10,2),
                         main = paste0('Tree ', t))
  abline(v = c(2002,2006))
}

name_nonshocks <- paste("delta_sck[", sample(which(data$w_sck != 1),9), "]", sep='')
util$plot_div_pairs(name_nonshocks, 'sigma', samples2, diagnostics2, transforms = list('sigma' = 1, 'inner_tau_sck' = 1))

par(mfrow = c(2,1))
t <- 6 
idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
names <- sapply(idxs_t,
                function(x) paste0('delta_sck[', x, ']'))
util$plot_realizations(samples2, names, plot_xs = data$years[idxs_t], N_plots = 50,
                       ylab = expression(delta[shocks]),
                       display_xlim = data$all_years[c(1,length(data$all_years))], display_ylim = c(-10,2),
                       main = paste0('Tree ', t))
abline(v = c(1996,2002))
name_shock <- paste("delta_sck[", idxs_t[which(data$all_years == 2006)], "]", sep='')
util$plot_div_pairs(name_shock, 'inner_tau_sck', samples2, diagnostics2, transforms = list('sigma' = 1, 'inner_tau_sck' = 1))




#------------------------------------#
# Increase prior of tau_outer, again # 
#------------------------------------#

# fit3 <- stan(file=file.path(wd, 'model', 'stan', 'model8_with3predictors_nopooling_standphi_nonconcordantshock.stan'),
#              data=data, seed=5838293,
#              chains = 4,
#              warmup=1000, iter=2024, refresh=10)
# saveRDS(fit3, file.path(wd, 'model/shocks/output/tuning_nc_az592', 'fit3.rds'))
fit3 <- readRDS(file.path(wd, 'model/shocks/output/tuning_nc_az592', 'fit3.rds'))
samples3 <- util$extract_expectand_vals(fit3)
diagnostics3 <- util$extract_hmc_diagnostics(fit3)
util$check_all_hmc_diagnostics(diagnostics3)

util$plot_pairs_by_chain(samples3[['sigma']], 'sigma', samples3[['inner_tau_sck']], 'inner_tau_sck')
util$plot_div_pairs('sigma', 'inner_tau_sck', samples3, diagnostics3, transforms =  list('sigma' = 1, 'inner_tau_sck' = 1))

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
ys <- dnorm(xs, 0, log(2) / 2.57)
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
xs <- seq(0, 1, 0.01)
ys <- dbeta(xs, 2, 10)
lines(xs, ys, lwd=2, col=util$c_mid_teal)

util$plot_expectand_pushforward(samples3[['omega_conc_sck']], 50,
                                flim = c(0,1), 
                                display_name=expression(omega[shock]^concordant))
xs <- seq(0, 5, 0.01)
ys <- dbeta(xs, 10, 2)
lines(xs, ys, lwd=2, col=util$c_mid_teal)
util$plot_expectand_pushforward(samples3[['omega_nonconc_sck']], 50,
                                flim = c(0,1), 
                                display_name=expression(omega[shock]^nonconcordant))
xs <- seq(0, 5, 0.01)
ys <- dbeta(xs, 1, 20)
lines(xs, ys, lwd=2, col=util$c_mid_teal)

par(mfrow = c(1,1))
util$plot_pairs_by_chain(samples3[['outer_tau_sck']], 'outer_tau_sck', samples3[['omega_conc_sck']], 'omega_conc_sck')
util$plot_pairs_by_chain(samples3[['inner_tau_sck']], 'inner_tau_sck', samples3[['omega_conc_sck']], 'omega_conc_sck')

par(mfrow = c(1,1))
s <- 1
names <- sapply(1:data$N_all_years,
                function(sp) paste0('f_sh[', s, ',', sp, ']'))
util$plot_realizations(samples3, names, plot_xs = data$all_years[1:data$N_all_years], N_plots = 50,
                       main = '', ylab = 'f_sh',
                       display_xlim = data$all_years[c(1,data$N_all_years)], display_ylim = c(-5,2))
util$plot_pairs_by_chain(samples3[[names[data$all_years == 2002]]], names[data$all_years == 2002], samples3[['outer_tau_sck']], 'outer_tau_sck')
util$plot_pairs_by_chain(samples3[[names[data$all_years == 2002]]], names[data$all_years == 2002], samples3[['rho_sh']], 'rho_sh')
util$plot_pairs_by_chain(samples3[['outer_tau_sck']], 'outer_tau_sck', samples3[['rho_sh']], 'rho_sh')

name_shocks <- paste("delta_sck[", which(data$w_sck == 1), "]", sep='')
util$plot_div_pairs(name_shocks, 'inner_tau_sck', samples3, diagnostics3, transforms = list('sigma' = 1, 'inner_tau_sck' = 1))

par(mfrow = c(7,3))
for(t in 1:data$N_trees){
  idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- sapply(idxs_t,
                  function(x) paste0('delta_sck[', x, ']'))
  util$plot_realizations(samples3, names, plot_xs = data$years[idxs_t], N_plots = 50,
                         ylab = expression(delta[shocks]),
                         display_xlim = data$all_years[c(1,length(data$all_years))], display_ylim = c(-10,2),
                         main = paste0('Tree ', t))
  abline(v = c(2002,2006))
}

#-------------------------------------------------------#
# Domain expertise is inconsistent with smaller values  # 
#-------------------------------------------------------#

fit4 <- stan(file=file.path(wd, 'model', 'stan', 'model8_with3predictors_nopooling_standphi_nonconcordantshock.stan'),
             data=data, seed=5838293,
             chains = 4,
             warmup=1000, iter=2024, refresh=10)
saveRDS(fit4, file.path(wd, 'model/shocks/output/tuning_nc_az592', 'fit4.rds'))
diagnostics4 <- util$extract_hmc_diagnostics(fit4)
util$check_all_hmc_diagnostics(diagnostics4)
samples4 <- util$extract_expectand_vals(fit4)
samples4[["sigma2+tau2"]] <- samples4[['sigma']]^2 + samples4[['inner_tau_sck']]^2

util$plot_pairs_by_chain(samples4[['sigma']], 'sigma', samples4[['inner_tau_sck']], 'inner_tau_sck')
util$plot_div_pairs('sigma', 'inner_tau_sck', samples4, diagnostics4, transforms =  list('sigma' = 1, 'inner_tau_sck' = 1))

par(mfrow = c(7,3))
for(t in 1:data$N_trees){
  idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- sapply(idxs_t,
                  function(x) paste0('delta_sck[', x, ']'))
  util$plot_realizations(samples4, names, plot_xs = data$years[idxs_t], N_plots = 50,
                         ylab = expression(delta[shocks]),
                         display_xlim = data$all_years[c(1,length(data$all_years))], display_ylim = c(-10,2),
                         main = paste0('Tree ', t))
  abline(v = c(2002,2006))
}

par(mfrow = c(1,1))
t <- 11
idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
names <- sapply(idxs_t,
                function(x) paste0('delta_sck[', x, ']'))
util$plot_realizations(samples4, names, plot_xs = data$years[idxs_t], N_plots = 50,
                       ylab = expression(delta[shocks]),
                       display_xlim = data$all_years[c(1,length(data$all_years))], display_ylim = c(-10,2),
                       main = paste0('Tree ', t))
abline(v = c(2007))
name_shock <- paste("delta_sck[", idxs_t[which(data$all_years == 2007)], "]", sep='')
util$plot_div_pairs(name_shock, 'inner_tau_sck', samples4, diagnostics4, transforms = list('sigma' = 1, 'inner_tau_sck' = 1))
util$plot_div_pairs(name_shock, 'sigma2+tau2', samples4, diagnostics4, transforms = list('sigma2+tau2' = 1))
