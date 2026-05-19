rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(cmdstanr)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)


folder <- '/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model'

data <- readRDS(file.path(folder,'data_22apr2026_PSME_except2_19202024.rds'))
datasets <- readRDS(file.path(folder,'datasets_22apr2026_PSME_except2_19202024.rds'))

samples <- readRDS(file.path(folder, 'model18/basesamples_model18_HSGP_PIED_utah_pi0.rds'))
util$check_all_expectand_diagnostics(samples)

par(mfrow = c(2,2))
util$plot_expectand_pushforward(samples[['phi_sck0']], 30,
                                display_name = 'Probability of concordant state',
                                flim = c(0,1))
prior <- rbeta(1e6, 2, 13)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
plot.new()
util$plot_expectand_pushforward(samples[['omega_conc_sck0']], 30,
                                display_name = 'Probability of shock | concordant',
                                flim = c(0,1))
prior <- rbeta(1e6, 3, 12)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(samples[['omega_shutdown0']], 30,
                                display_name = 'Probability of shutdown | shock',
                                flim = c(0,1))
prior <- rbeta(1e6, 1, 25)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)




