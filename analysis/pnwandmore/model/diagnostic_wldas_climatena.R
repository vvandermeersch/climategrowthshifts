rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('functions/custom_functions.R', local=util)
setwd(wd)

# Model data
data <- readRDS(file.path(wd, 'output/model/climatena', 'data_11august2025_short.rds'))

#
fullfit <- FALSE
fullsamples <- FALSE
if(fullfit){
  
}else if (fullsamples){
  
}else{
  base_samples_climna_short <- readRDS(file.path(wd, "output/model/climatena/samples_11august2025_differentunits_partialpooling_2clades_centered_extended_short.rds"))
  base_samples_wpre <- readRDS(file.path(wd, "output/model/samples_11july2025_partialpooling_2clades_centered_extended_wpre.rds"))
}

par(mfrow=c(1, 3), mar = c(4,4,1,1))
util$plot_expectand_pushforward(base_samples_climna_short[['mu_gdd[1]']], 50,
                                flim = c(-0.15,0.15),
                                display_name=expression(beta[GDD]),
                                col = '#27278f')
util$plot_expectand_pushforward(base_samples_climna_short[['mu_gdd[2]']], 50,
                                flim = c(-0.15,0.15),
                                col = '#278f27',
                                add = TRUE)
text(x = 0.085, y = 50, label = "Gymno.", col = '#27278f', cex = 1.2)
text(x = -0.06, y = 10, label = "Angio.", col = '#278f27', cex = 1.2)

util$plot_expectand_pushforward(base_samples_climna_short[['mu_winterprec[1]']], 50,
                                flim = c(-0.15,0.15),
                                display_name=expression(beta[winterprecip]),
                                col = '#27278f')
util$plot_expectand_pushforward(base_samples_climna_short[['mu_winterprec[2]']], 50,
                                flim = c(-0.15,0.15),
                                col = '#278f27',
                                add = TRUE)

util$plot_expectand_pushforward(base_samples_climna_short[['mu_ffp[1]']], 50,
                                flim = c(-0.15,0.15),
                                display_name=expression(beta[FFP]),
                                col = '#27278f')
util$plot_expectand_pushforward(base_samples_climna_short[['mu_ffp[2]']], 50,
                                flim =c(-0.15,0.15),
                                col = '#278f27',
                                add = TRUE)

par(mfrow=c(1, 1), mar = c(4,4,1,1))
util$plot_expectand_pushforward(base_samples_wpre[['mu_gdd[1]']], 50,
                                flim = c(-0.15,0.15),
                                display_name=expression(beta[GDD]),
                                col = '#9393C7', border = NA)
util$plot_expectand_pushforward(base_samples_wpre[['mu_gdd[2]']], 50,
                                flim = c(-0.15,0.15),
                                display_name=expression(beta[GDD]),
                                col = '#93c793', border = NA,
                                add = TRUE)
util$plot_expectand_pushforward(base_samples_climna_short[['mu_gdd[1]']], 50,
                                flim = c(-0.15,0.15),
                                col = '#27278f',
                                add = TRUE)
util$plot_expectand_pushforward(base_samples_climna_short[['mu_gdd[2]']], 50,
                                flim = c(-0.15,0.15),
                                col = '#278f27',
                                add = TRUE)


util$plot_expectand_pushforward(base_samples_wpre[['mu_wpre[1]']], 50,
                                flim = c(-0.15,0.15),
                                display_name=expression(beta[GDD]),
                                col = '#9393C7', border = NA)
util$plot_expectand_pushforward(base_samples_wpre[['mu_wpre[2]']], 50,
                                flim = c(-0.15,0.15),
                                display_name=expression(beta[GDD]),
                                col = '#93c793', border = NA,
                                add = TRUE)
util$plot_expectand_pushforward(base_samples_climna_short[['mu_winterprec[1]']], 50,
                                flim = c(-0.15,0.15),
                                col = '#27278f',
                                add = TRUE)
util$plot_expectand_pushforward(base_samples_climna_short[['mu_winterprec[2]']], 50,
                                flim = c(-0.15,0.15),
                                col = '#278f27',
                                add = TRUE)


