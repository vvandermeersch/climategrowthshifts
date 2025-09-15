rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('functions/custom_functions.R', local=util)
setwd(wd)

# Model data
data <- readRDS(file.path(wd, 'output/model/climatena', 'data_11sept2025_long_PIPO.rds'))
samples <- readRDS(file.path(wd, "output/model/climatena/samples_11sept2025_differentunits_partialpooling_2clades_centered_extended_onlyPIPO.rds"))

util$check_all_expectand_diagnostics(samples)

par(mfrow=c(1, 3), mar = c(4,4,1,1))
util$plot_expectand_pushforward(samples[['beta_gdd']], 25,
                                flim = c(-0.15,0.15),
                                display_name=expression(beta[GDD]),
                                col = '#27278f')

util$plot_expectand_pushforward(samples[['sigma']], 25,
                                flim = c(0,1),
                                display_name=expression(beta[GDD]),
                                col = '#27278f')



util$plot_pairs_by_chain(samples[['beta_gdd']], 'beta_gdd', samples[['beta_winterprec']], 'beta_winterprec')
util$plot_pairs_by_chain(samples[['beta_gdd']], 'beta_gdd', samples[['beta_ffp']], 'beta_ffp')
util$plot_pairs_by_chain(samples[['beta_gdd']], 'beta_gdd', samples[['rho_sp']], 'rho_sp')
