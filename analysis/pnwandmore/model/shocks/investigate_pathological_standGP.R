rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('functions/custom_functions.R', local=util)
setwd(wd)

data <- readRDS(file.path(wd, 'output/model/climatena', 'data_15oct2025_long_az618.rds'))
fit <- readRDS(file.path(wd, "model/shocks/output", "fit_15oct2025_long_az618.rds"))
samples <-  util$extract_expectand_vals(fit)
fit_before2002 <- readRDS(file.path(wd, "model/shocks/output", "fit_15oct2025_long_az618_before2002.rds"))
samples_before2002 <-  util$extract_expectand_vals(fit_before2002)


par(mfrow = c(2,3), mar = c(3,5.5,2,1), cex.lab = 1.2)
stand <- 1
names <- sapply(1:data$N_all_years,
                function(sp) paste0('f_tilde_sh[', stand, ',', sp, ']'))
util$plot_realizations(samples, names, plot_xs = data$all_years, N_plots = 50,
                       main = "f_tilde", ylab = 'Function outputs\n(entire time series)',
                       display_xlim = c(1900,2020), , display_ylim = c(-5,5))
names <- sapply(1:data$N_all_years,
                function(sp) paste0('f_sh[', stand, ',', sp, ']'))
util$plot_realizations(samples, names, plot_xs = data$all_years, N_plots = 50,
                       main = "f", ylab = '',
                       display_xlim = c(1900,2020), display_ylim = c(-8,2))
abline(v = 2002, lty = 2)
util$plot_expectand_pushforward(samples[['rho_sh']], 50,
                                flim = c(0,6),
                                display_name="", main = 'rho_sh')

names <- sapply(1:100,
                function(sp) paste0('f_tilde_sh[', stand, ',', sp, ']'))
util$plot_realizations(samples_before2002, names, plot_xs = data$all_years[data$all_years<2002], N_plots = 50,
                       main = "f_tilde", ylab = 'Function outputs\n(befgore 2002)',
                       display_xlim = c(1900,2020), display_ylim = c(-5,5))
names <- sapply(1:100,
                function(sp) paste0('f_sh[', stand, ',', sp, ']'))
util$plot_realizations(samples_before2002, names, plot_xs = data$all_years[data$all_years<2002], N_plots = 50,
                       main = "f", ylab = '',
                       display_xlim = c(1900,2020), display_ylim = c(-8,2))
util$plot_expectand_pushforward(samples_before2002[['rho_sh']], 50,
                                flim = c(0,6),
                                display_name="", main = 'rho_sh')

par(mfrow = c(2,4), mar = c(3,5.5,2,1), cex.lab = 1.2)
util$plot_expectand_pushforward(samples[['rho_sp']], 50,
                                flim = c(0,30),
                                display_name="rho_sh")
util$plot_expectand_pushforward(samples[['rho_sh']], 50,
                                flim = c(0,6),
                                display_name="rho_sh")
util$plot_expectand_pushforward(samples[['gamma_sh']], 50,
                                flim = c(0,3),
                                display_name="gamma_sh")
util$plot_expectand_pushforward(samples[['sigma']], 50,
                                flim = c(0,1),
                                display_name="sigma")
util$plot_expectand_pushforward(samples_before2002[['rho_sp']], 50,
                                flim = c(0,30),
                                display_name="rho_sh")
util$plot_expectand_pushforward(samples_before2002[['rho_sh']], 50,
                                flim = c(0,6),
                                display_name="rho_sh")
util$plot_expectand_pushforward(samples_before2002[['gamma_sh']], 50,
                                flim = c(0,3),
                                display_name="gamma_sh")
util$plot_expectand_pushforward(samples_before2002[['sigma']], 50,
                                flim = c(0,1),
                                display_name="sigma")

par(mfrow = c(1,1), mar = c(3,5.5,2,1), cex.lab = 1.2)
util$plot_pairs_by_chain(samples_before2002[['gamma_sh']], 'gamma_sh', samples_before2002[['rho_sh']], 'rho_sh')
