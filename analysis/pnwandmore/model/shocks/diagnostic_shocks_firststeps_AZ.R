rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('functions/custom_functions.R', local=util)
setwd(wd)

data <- readRDS(file.path(wd, 'output/model/climatena', 'data_15oct2025_long_recentPIPOinAZ.rds'))
fit <- readRDS(file.path(wd, "model/shocks/output", "fit_15oct2025_long_recentPIPOinAZ_shocks_v1.rds")) 
fit_noshock <- readRDS(file.path(wd, "model/shocks/output", "fit_15oct2025_long_recentPIPOinAZ.rds"))
samples_noshock <- util$extract_expectand_vals(fit_noshock)

# Computational diagnostics
diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit)
base_samples <- util$filter_expectands(samples,
                                       c('alpha', 
                                         'beta_gdd', 'beta_winterprec', 
                                         'rho_sp', 'gamma_sp',
                                         'f_tilde_sh',
                                         'rho_sh', 'gamma_sh',
                                         'delta_tilde_sck', 'inner_tau_sck', 'outer_tau_sck', 'gamma_sck',
                                         'sigma'),
                                       check_arrays=TRUE)
util$check_all_expectand_diagnostics(base_samples)

# Diagnostics restricted to parameters for observed years
f_tilde_sh_obs <- unlist(sapply(1:data$N_stands, function(s){
  trees_s <- which(data$stand_idxs==s)
  minyear_s <- min(data$all_years_idxs[data$tree_start_idxs[trees_s]])
  maxyear_s <- max(data$all_years_idxs[data$tree_end_idxs[trees_s]])
  return(paste0('f_tilde_sh[', s, ',', minyear_s:maxyear_s, ']'))
}))
delta_tilde_sck_obs <- unlist(sapply(1:data$N_stands, function(s){
  trees_s <- which(data$stand_idxs==s)
  minyear_s <- min(data$all_years_idxs[data$tree_start_idxs[trees_s]])
  maxyear_s <- max(data$all_years_idxs[data$tree_end_idxs[trees_s]])
  return(paste0('delta_tilde_sck[', s, ',', minyear_s:maxyear_s, ']'))
}))
base_samples <- util$filter_expectands(samples,
                                       c('alpha', 
                                         'beta_gdd', 'beta_winterprec', 
                                         'rho_sp', 'gamma_sp',
                                         f_tilde_sh_obs,
                                         'rho_sh', 'gamma_sh',
                                         delta_tilde_sck_obs, 'inner_tau_sck', 'outer_tau_sck', 'gamma_sck',
                                         'sigma'))
util$check_all_expectand_diagnostics(base_samples)

# Look at some misbehaving parameters
for(y in misbehaving_years){
  util$plot_div_pairs(paste0('delta_sck[1,',y,']'), 'inner_tau_sck', 
                      samples, diagnostics, transforms = list('inner_tau_sck' = 1))
  text(x = util$ensemble_mcmc_quantile_est(samples[[paste0('delta_sck[1,',y,']')]], c(0.5)), 
       y = -1.2, label = data$all_years[y], xpd=NA)
}

y <- 1
util$plot_div_pairs(paste0('delta_sck[1,',y,']'), 'inner_tau_sck', 
                    samples, diagnostics, transforms = list('inner_tau_sck' = 1))

# Look at parameters
par(mfrow = c(2,4), mar = c(4.5,5,2,1), cex.lab = 1.2)
util$plot_expectand_pushforward(samples[['rho_sp']], 50,
                                flim = c(0,30),
                                display_name="rho_sh")
util$plot_expectand_pushforward(samples[['rho_sh']], 50,
                                flim = c(0,10),
                                display_name="rho_sh")
util$plot_expectand_pushforward(samples[['gamma_sh']], 50,
                                flim = c(0,1),
                                display_name="gamma_sh")
util$plot_expectand_pushforward(samples[['sigma']], 50,
                                flim = c(0,1),
                                display_name="sigma")
util$plot_expectand_pushforward(samples[['inner_tau_sck']], 50,
                                flim = c(0,1),
                                display_name="inner_tau_sck")
util$plot_expectand_pushforward(samples[['outer_tau_sck']], 50,
                                flim = c(0,4),
                                display_name="outer_tau_sck")
util$plot_expectand_pushforward(samples[['gamma_sck']], 50,
                                flim = c(0,1),
                                display_name="gamma_sck")

util$plot_pairs_by_chain(samples[['inner_tau_sck']], 'inner_tau_sck',
                         samples[['outer_tau_sck']], 'outer_tau_sck')

par(mfrow = c(1,1), mar = c(4,5.5,2,1), cex.lab = 1.2)
util$plot_expectand_pushforward(samples[['inner_tau_sck']], 50,
                                flim = c(0,1),
                                display_name="inner_tau_sck")
xs <- seq(0,1, 0.01)
ys <- dnorm(xs, 0, log(1.5) / 2.57)
lines(xs, ys, lwd=2, col=util$c_light_teal)

par(mfrow = c(1,1), mar = c(4,5.5,2,1), cex.lab = 1.2)
util$plot_expectand_pushforward(samples[['sigma']], 100,
                                flim = c(0,1),
                                display_name="sigma")
xs <- seq(0,1, 0.01)
ys <- dnorm(xs, 0, log(1.1) / 2.57)
lines(xs, ys, lwd=2, col=util$c_light_teal)

unpermuted_samples <- extract(fit, permute=FALSE)

par(mfrow=c(2, 2))
for (c in 1:4){
  plot(1:3048, unpermuted_samples[, c, 'inner_tau_sck'], type="l", lwd=1, col=util$c_dark,
       main=paste("Chain ", c, sep=""),
       xlab="Iteration",  xlim=c(1, 3048),
       ylab="inner_tau_sck", ylim=c(0, 0.1))
  
}

stand1 <- 4
stand2 <- 5
year <- 114
util$plot_pairs_by_chain(samples[[paste0('delta_sck[', stand1, ',', year, ']')]], paste0('delta_sck[', stand1, ',', year, ']'),
                         log(samples[[paste0('inner_tau_sck')]]), 'log(inner_tau_sck)')

par(mfrow=c(2, 2))
for (c in 1:4){
  plot(1:3048, unpermuted_samples[, c, paste0('delta_sck[', stand1, ',', year, ']')], type="l", lwd=1, col=util$c_dark,
       main=paste("Chain ", c, sep=""),
       xlab="Iteration",  xlim=c(1, 3048),
       ylab=paste0('delta_sck[', stand1, ',', year, ']'), ylim=c(-0.25, 0.75))
  
}


par(mfrow = c(2,3), mar = c(3,5.5,2,1), cex.lab = 1.2)
stand <- 4
names <- sapply(1:data$N_all_years,
                function(sp) paste0('f_sh[', stand, ',', sp, ']'))
util$plot_realizations(samples, names, plot_xs = data$all_years, N_plots = 50,
                       main = "Stand-level GP", ylab = 'With shocks',
                       display_xlim = c(1900,2020), display_ylim = c(-8,2))
names <- sapply(1:data$N_all_years,
                function(sp) paste0('f_sh[', stand, ',', sp, ']'))
util$plot_realizations(samples, names, plot_xs = data$all_years, N_plots = 50,
                       main = "Stand-level GP (zoom)", ylab = '',
                       display_xlim = c(1900,2020), display_ylim = c(-1.5,1.5))
names <- sapply(1:data$N_all_years,
                function(sp) paste0('delta_sck[', stand, ',', sp, ']'))
util$plot_realizations(samples, names, plot_xs = data$all_years, N_plots = 50,
                       main = "Stand-level shocks", ylab = '',
                       display_xlim = c(1900,2020), , display_ylim = c(-8,2))
names <- sapply(1:data$N_all_years,
                function(sp) paste0('f_sh[', stand, ',', sp, ']'))
util$plot_realizations(samples_noshock, names, plot_xs = data$all_years, N_plots = 50,
                       main = "Stand-level GP", ylab = 'Without shocks',
                       display_xlim = c(1900,2020), display_ylim = c(-8,2))
names <- sapply(1:data$N_all_years,
                function(sp) paste0('f_sh[', stand, ',', sp, ']'))
util$plot_realizations(samples_noshock, names, plot_xs = data$all_years, N_plots = 50,
                       main = "Stand-level GP (zoom)", ylab = '',
                       display_xlim = c(1900,2020), display_ylim = c(-1.5,1.5))

par(mfrow = c(2,8), mar = c(3,3.5,2,1), cex.lab = 1.2)
for(stand in 1:data$N_stands){
  names <- sapply(1:data$N_all_years,
                  function(sp) paste0('f_sh[', stand, ',', sp, ']'))
  util$plot_realizations(samples, names, plot_xs = data$all_years, N_plots = 50,
                         main = "Stand-level GP", ylab = '',
                         display_xlim = c(1900,2020), display_ylim = c(-1.5,1.5))
}
for(stand in 1:data$N_stands){
  names <- sapply(1:data$N_all_years,
                  function(sp) paste0('delta_sck[', stand, ',', sp, ']'))
  util$plot_realizations(samples, names, plot_xs = data$all_years, N_plots = 50,
                         main = "Stand-level shocks", ylab = '',
                         display_xlim = c(1900,2020), , display_ylim = c(-8,2))
}



par(mfrow=c(3,1), mar = c(4.5,1.5,0,1))
prior <- rbeta(12e3, 200, 5)
hist(prior, breaks=seq(0, 1, 0.0025),
     main="", col=util$c_light, border='white',
     xlab="", yaxt='n', ylab="", ylim=c(0, 3e3), xlim = c(0.9,1),
     mar = c(1,1,1,1))
abline(v =  quantile(prior, 0.005), lty = 2, col = util$c_mid_highlight, lwd = 1.2)
# legend(x = quantile(prior, 0.005), y = 1.35e3, legend = round(quantile(prior, 0.005),2), adj = 0, 
#        x.intersp = -0.5, y.intersp = 0, box.col = util$c_mid_highlight, xjust = 0.5)
abline(v = quantile(prior, 0.995), lty = 2, col = util$c_mid_highlight, lwd = 1.2)
# legend(x = quantile(prior, 0.995), y = 1.35e3, legend = round(quantile(prior, 0.995),2), adj = 0, 
#        x.intersp = -0.5, y.intersp = 0, box.col = util$c_mid_highlight, xjust = 0.5)

prior <- rbeta(12e3, 1000, 10)
hist(prior, breaks=seq(0, 1, 0.0025),
     main="", col=util$c_light, border='white',
     xlab="", yaxt='n', ylab="", ylim=c(0, 3e3), xlim = c(0.9,1),
     mar = c(1,1,1,1))
abline(v =  quantile(prior, 0.005), lty = 2, col = util$c_mid_highlight, lwd = 1.2)
# legend(x = quantile(prior, 0.005), y = 1.35e3, legend = round(quantile(prior, 0.005),2), adj = 0, 
#        x.intersp = -0.5, y.intersp = 0, box.col = util$c_mid_highlight, xjust = 0.5)
abline(v = quantile(prior, 0.995), lty = 2, col = util$c_mid_highlight, lwd = 1.2)
# legend(x = quantile(prior, 0.995), y = 1.35e3, legend = round(quantile(prior, 0.995),2), adj = 0, 
#        x.intersp = -0.5, y.intersp = 0, box.col = util$c_mid_highlight, xjust = 0.5)

prior <- rbeta(12e3, 100, 1)
hist(prior, breaks=seq(0, 1, 0.0025),
     main="", col=util$c_light, border='white',
     xlab="gamma", yaxt='n', ylab="", ylim=c(0, 3e3), xlim = c(0.9,1),
     mar = c(1,1,1,1))
abline(v =  quantile(prior, 0.005), lty = 2, col = util$c_mid_highlight, lwd = 1.2)
# legend(x = quantile(prior, 0.005), y = 1.35e3, legend = round(quantile(prior, 0.005),2), adj = 0, 
#        x.intersp = -0.5, y.intersp = 0, box.col = util$c_mid_highlight, xjust = 0.5)
abline(v = quantile(prior, 0.995), lty = 2, col = util$c_mid_highlight, lwd = 1.2)
# legend(x = quantile(prior, 0.995), y = 1.35e3, legend = round(quantile(prior, 0.995),2), adj = 0, 
#        x.intersp = -0.5, y.intersp = 0, box.col = util$c_mid_highlight, xjust = 0.5)


hist(abs(rnorm(12e3, 0, log(1.5)/ 2.57)), breaks=seq(0, 1, 0.007),
     main="", col=util$c_light, border=util$c_light_highlight,
     xlab="inner_tau_sck", yaxt='n', ylab="", ylim=c(0, 3e3), xlim = c(0,.5),
     mar = c(1,1,1,1))
hist(samples[['inner_tau_sck']], breaks=seq(0, 1, 0.007),
     main="", col=util$c_dark, border=util$c_dark_highlight, add=T)

hist(abs(rnorm(12e3, 0, log(5)/ 2.57)), breaks=seq(0, 3, 0.05),
     main="", col=util$c_light, border=util$c_light_highlight,
     xlab="outer_tau_sck", yaxt='n', ylab="", ylim=c(0, 3e3), xlim = c(0,3),
     mar = c(1,1,1,1))
hist(samples[['outer_tau_sck']], breaks=seq(0, 3, 0.05),
     main="", col=util$c_dark, border=util$c_dark_highlight, add=T)


