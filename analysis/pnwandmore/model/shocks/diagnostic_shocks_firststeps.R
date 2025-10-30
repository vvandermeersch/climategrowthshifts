rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('functions/custom_functions.R', local=util)
setwd(wd)

data <- readRDS(file.path(wd, 'output/model/climatena', 'data_15oct2025_long_az618.rds'))
data$w_sck <-  lapply(1, function(x)  ifelse(data$all_years == 2002, 1, 0.2))
fit <- readRDS(file.path(wd, "model/shocks/output", "fit_15oct2025_long_az618_log15.rds")) # v3 is c(1902, 1904, 1961, 1977, 2002, 2004)
fit_noshock <- readRDS(file.path(wd, "model/shocks/output", "fit_15oct2025_long_az618.rds"))
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

misbehaving_years <- c(1, 55, 60, 62, 75, 103, 106) # v1
misbehaving_years <- c(1, 2, 3, 55, 60, 76, 106) # v2
misbehaving_years <- c(106) # v3
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

 
par(mfrow = c(2,3), mar = c(3,5.5,2,1), cex.lab = 1.2)
stand <- 1
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

par(mfrow = c(2,1), mar = c(3,3.5,2,1), cex.lab = 1.2)
names <- sapply(1:data$N_all_years,
                function(sp) paste0('f_sh[', stand, ',', sp, ']'))
util$plot_realizations(samples, names, plot_xs = data$all_years, N_plots = 50,
                       main = "Stand-level GP", ylab = 'Without shocks',
                       display_xlim = c(1900,2020), display_ylim = c(-2,2))
abline(v = c(1902, 1904, 1961, 1977, 2002, 2004), lty = 2, col = 'grey70', lwd = 2)
names <- sapply(1:data$N_all_years,
                function(sp) paste0('delta_sck[', stand, ',', sp, ']'))
util$plot_realizations(samples, names, plot_xs = data$all_years, N_plots = 50,
                       main = "Stand-level shocks", ylab = '',
                       display_xlim = c(1900,2020), , display_ylim = c(-3,1))
abline(v = c(1902, 1904, 1961, 1977, 2002, 2004), lty = 2, col = 'grey70', lwd = 2)
abline(h = c(-0.5,0.5), lty = 5, lwd = 2, col = 'grey20')


par(mfrow = c(3,1), mar = c(3,3.5,2,1), cex.lab = 1.2)
prior <- rbeta(3e3, 1000, 10)
hist(prior, breaks=seq(0, 1, 0.002),
     main="", col=util$c_light, border='white',
     xlab="", yaxt='n', ylab="", ylim=c(0, 1e3), xlim = c(0.9,1),
     mar = c(1,1,1,1))
hist(samples[['gamma_sck']], breaks=seq(0, 1, 0.002),
     main="", col=util$c_dark, border='white', add=T)

hist(abs(rnorm(3e3, 0, log(1.5)/ 2.57)), breaks=seq(0, 1, 0.01),
     main="", col=util$c_light, border='white',
     xlab="inner_tau_sck", yaxt='n', ylab="", ylim=c(0, 1e3), xlim = c(0,.5),
     mar = c(1,1,1,1))
hist(samples[['inner_tau_sck']], breaks=seq(0, 1, 0.01),
     main="", col=util$c_dark, border='white', add=T)

hist(abs(rnorm(3e3, 0, log(5)/ 2.57)), breaks=seq(0, 3, 0.05),
     main="", col=util$c_light, border='white',
     xlab="outer_tau_sck", yaxt='n', ylab="", ylim=c(0, 1e3), xlim = c(0,3),
     mar = c(1,1,1,1))
hist(samples[['outer_tau_sck']], breaks=seq(0, 4, 0.05),
     main="", col=util$c_dark, border='white', add=T)
