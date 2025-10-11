rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('functions/custom_functions.R', local=util)
setwd(wd)

data <- readRDS(file.path(wd, 'output/model', 'data_24sept2025.rds'))
samples <- readRDS(file.path(wd, "output/model/samples_24sept2025_partialpooling_2clades_centered_extended.rds"))

in2mm <-25.4 # scale factor to convert inches to mm
pdf(file.path(wd, 'figures', 'newyork2025', 'model_estimates', 'GP_parameters_posterior_distributions.pdf'),
    width=200/in2mm,height=300/in2mm,paper="special")
# Species-level estimates, with clade estimates
layout(mat = matrix(c(1, 2,3,4,5,6,7,8,9), ncol = 3, byrow = T),
       widths = c(3, 4*sum(data[["clade_idxs"]]==2)/length(data[["clade_idxs"]]), 0.7))  
par(mar = c(1,5,1,1), cex.lab = 2)
names <- sapply(which(data[["clade_idxs"]]==1),
                function(sp) paste0('rho_sp[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(samples, names, clades = data$clade_idxs[data[["clade_idxs"]]==1], xlab="", ylab=expression(rho[sp]), 
                                             xticklabs = NA, display_ylim = c(0,35), line0 = FALSE)
par(mar = c(1,0,1,1))
names <- sapply(which(data[["clade_idxs"]]==2),
                function(sp) paste0('rho_sp[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(samples, names, clades = data$clade_idxs[data[["clade_idxs"]]==2], xlab="", ylab="", 
                                             xticklabs = NA, display_ylim = c(0,35), yticklabs = NA)
par(mar = c(1,2,1,1))
util$plot_expectand_pushforward_reverse(samples[['mu_rho[1]']], 50,
                                        flim = c(0,35),
                                        display_name=expression(rho[sp]),
                                        col = '#27278f')
util$plot_expectand_pushforward_reverse(samples[['mu_rho[2]']], 50,
                                        flim = c(0,35),
                                        display_name=expression(rho[sp]),
                                        col = '#278f27', add = TRUE)

par(mar = c(1,5,1,1), cex.lab = 2)
names <- sapply(which(data[["clade_idxs"]]==1),
                function(sp) paste0('gamma_sp[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(samples, names, clades = data$clade_idxs[data[["clade_idxs"]]==1], xlab="", ylab=expression(gamma[sp]), 
                                             xticklabs = NA, display_ylim = c(0,1), line0 = FALSE)
text(x = sum(data[["clade_idxs"]]==1)/2, y = 0.8, label = "Gymnosperms", col = '#27278f', cex = 1.5)
par(mar = c(1,0,1,1))
names <- sapply(which(data[["clade_idxs"]]==2),
                function(sp) paste0('gamma_sp[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(samples, names, clades = data$clade_idxs[data[["clade_idxs"]]==2], xlab="", ylab="", 
                                             xticklabs = NA, display_ylim = c(0,1), yticklabs = NA)
text(x = sum(data[["clade_idxs"]]==2)/2+0.5, y = 0.8, label = "Angiosperms", col = '#278f27', cex = 1.5)
par(mar = c(1,2,1,1))
util$plot_expectand_pushforward_reverse(samples[['mu_gamma[1]']], 50,
                                        flim = c(0,1),
                                        display_name=expression(gamma[sp]),
                                        col = '#27278f')
util$plot_expectand_pushforward_reverse(samples[['mu_gamma[2]']], 50,
                                        flim = c(0,1),
                                        display_name=expression(gamma[sp]),
                                        col = '#278f27', add = TRUE)

par(mar = c(1,5,1,1), cex.lab = 2)
names <- sapply(which(data[["clade_idxs"]]==1),
                function(sp) paste0('kappa_sh[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(samples, names, clades = data$clade_idxs[data[["clade_idxs"]]==1], xlab="", ylab=expression(kappa[sh]), 
                                             xticklabs = NA, display_ylim = c(0,1.5), line0 = FALSE)
par(mar = c(1,0,1,1))
names <- sapply(which(data[["clade_idxs"]]==2),
                function(sp) paste0('kappa_sh[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(samples, names, clades = data$clade_idxs[data[["clade_idxs"]]==2], xlab="", ylab="", 
                                             xticklabs = NA, display_ylim = c(0,1.5), yticklabs = NA)
par(mar = c(1,2,1,1))
util$plot_expectand_pushforward_reverse(samples[['mu_kappa[1]']], 50,
                                        flim = c(0,1.5),
                                        display_name=expression(kappa[sh]),
                                        col = '#27278f')
util$plot_expectand_pushforward_reverse(samples[['mu_kappa[2]']], 50,
                                        flim = c(0,1.5),
                                        display_name=expression(kappa[sh]),
                                        col = '#278f27', add = TRUE)

dev.off()

pdf(file.path(wd, 'figures', 'newyork2025', 'model_estimates', 'climate_slopes_posterior_distributions.pdf'),
    width=200/in2mm,height=300/in2mm,paper="special")
# Species-level estimates, with clade estimates
layout(mat = matrix(c(1, 2,3,4,5,6,7,8,9), ncol = 3, byrow = T),
       widths = c(3, 4*sum(data[["clade_idxs"]]==2)/length(data[["clade_idxs"]]), 0.7))  
par(mar = c(1,5,1,1), cex.lab = 2)
names <- sapply(which(data[["clade_idxs"]]==1),
                function(sp) paste0('beta_gdd[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(samples, names, clades = data$clade_idxs[data[["clade_idxs"]]==1], xlab="", ylab=expression(beta[GDD]), 
                                             xticklabs = NA, display_ylim = c(-0.25,0.25), line0 = TRUE)
par(mar = c(1,0,1,1))
names <- sapply(which(data[["clade_idxs"]]==2),
                function(sp) paste0('beta_gdd[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(samples, names, clades = data$clade_idxs[data[["clade_idxs"]]==2], xlab="", ylab="", 
                                             xticklabs = NA, display_ylim = c(-0.25,0.25),line0 = TRUE, yticklabs = NA)
par(mar = c(1,2,1,1))
util$plot_expectand_pushforward_reverse(samples[['mu_gdd[1]']], 50,
                                        flim = c(-0.25,0.25),
                                        display_name=expression(beta[gdd]),
                                        col = '#27278f')
util$plot_expectand_pushforward_reverse(samples[['mu_gdd[2]']], 50,
                                        flim = c(-0.25,0.25),
                                        display_name=expression(beta[gdd]),
                                        col = '#278f27', add = TRUE)

par(mar = c(1,5,1,1), cex.lab = 2)
names <- sapply(which(data[["clade_idxs"]]==1),
                function(sp) paste0('beta_vpd[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(samples, names, clades = data$clade_idxs[data[["clade_idxs"]]==1], xlab="", ylab=expression(beta[VPD]), 
                                             xticklabs = NA, display_ylim = c(-0.25,0.25), line0 = TRUE)
par(mar = c(1,0,1,1))
names <- sapply(which(data[["clade_idxs"]]==2),
                function(sp) paste0('beta_vpd[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(samples, names, clades = data$clade_idxs[data[["clade_idxs"]]==2], xlab="", ylab="", 
                                             xticklabs = NA, display_ylim = c(-0.25,0.25), line0 = TRUE, yticklabs = NA)
par(mar = c(1,2,1,1))
util$plot_expectand_pushforward_reverse(samples[['mu_vpd[1]']], 50,
                                        flim = c(-0.25,0.25),
                                        display_name=expression(beta[gdd]),
                                        col = '#27278f')
util$plot_expectand_pushforward_reverse(samples[['mu_vpd[2]']], 50,
                                        flim = c(-0.25,0.25),
                                        display_name=expression(beta[gdd]),
                                        col = '#278f27', add = TRUE)

par(mar = c(1,5,1,1), cex.lab = 2)
names <- sapply(which(data[["clade_idxs"]]==1),
                function(sp) paste0('beta_sm[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(samples, names, clades = data$clade_idxs[data[["clade_idxs"]]==1], xlab="", ylab=expression(beta[SM]), 
                                             xticklabs = NA, display_ylim = c(-0.25,0.25), line0 = TRUE)
par(mar = c(1,0,1,1))
names <- sapply(which(data[["clade_idxs"]]==2),
                function(sp) paste0('beta_sm[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(samples, names, clades = data$clade_idxs[data[["clade_idxs"]]==2], xlab="", ylab="", 
                                             xticklabs = NA, display_ylim = c(-0.25,0.25), line0 = TRUE, yticklabs = NA)
par(mar = c(1,2,1,1))
util$plot_expectand_pushforward_reverse(samples[['mu_sm[1]']], 50,
                                        flim = c(-0.25,0.25),
                                        display_name=expression(beta[gdd]),
                                        col = '#27278f')
util$plot_expectand_pushforward_reverse(samples[['mu_sm[2]']], 50,
                                        flim = c(-0.25,0.25),
                                        display_name=expression(beta[gdd]),
                                        col = '#278f27', add = TRUE)
dev.off()

samples_int <- readRDS(file.path(wd, "output/model/samples_24sept2025_partialpooling2clades_1interaction.rds"))

pdf(file.path(wd, 'figures', 'newyork2025', 'model_estimates', 'winteraction', 'climate_slopes_posterior_distributions.pdf'),
    width=200/in2mm,height=400/in2mm,paper="special")
# Species-level estimates, with clade estimates
layout(mat = matrix(c(1, 2,3,4,5,6,7,8,9,10,11,12,13,14,15), ncol = 3, byrow = T),
       widths = c(3, 4*sum(data[["clade_idxs"]]==2)/length(data[["clade_idxs"]]), 0.7))  
par(mar = c(1,5,1,1), cex.lab = 2)
names <- sapply(which(data[["clade_idxs"]]==1),
                function(sp) paste0('beta_gdd[', sp, ']'))
baseline_median <- sapply(names, function(n)  util$ensemble_mcmc_quantile_est(samples[[n]], 0.5))
util$plot_disc_pushforward_quantiles_2clades(samples_int, names, clades = data$clade_idxs[data[["clade_idxs"]]==1], xlab="", ylab=expression(beta[GDD]), 
                                             xticklabs = NA, display_ylim = c(-0.25,0.25), line0 = TRUE, baseline_values = baseline_median)
par(mar = c(1,0,1,1))
names <- sapply(which(data[["clade_idxs"]]==2),
                function(sp) paste0('beta_gdd[', sp, ']'))
baseline_median <- sapply(names, function(n)  util$ensemble_mcmc_quantile_est(samples[[n]], 0.5))
util$plot_disc_pushforward_quantiles_2clades(samples_int, names, clades = data$clade_idxs[data[["clade_idxs"]]==2], xlab="", ylab="", 
                                             xticklabs = NA, display_ylim = c(-0.25,0.25),line0 = TRUE, yticklabs = NA, baseline_values = baseline_median)
par(mar = c(1,2,1,1))
util$plot_expectand_pushforward_reverse(samples_int[['mu_gdd[1]']], 50,
                                        flim = c(-0.25,0.25),
                                        display_name=expression(beta[gdd]),
                                        col = '#27278f')
util$plot_expectand_pushforward_reverse(samples_int[['mu_gdd[2]']], 50,
                                        flim = c(-0.25,0.25),
                                        display_name=expression(beta[gdd]),
                                        col = '#278f27', add = TRUE)

par(mar = c(1,5,1,1), cex.lab = 2)
names <- sapply(which(data[["clade_idxs"]]==1),
                function(sp) paste0('beta_vpd[', sp, ']'))
baseline_median <- sapply(names, function(n)  util$ensemble_mcmc_quantile_est(samples[[n]], 0.5))
util$plot_disc_pushforward_quantiles_2clades(samples_int, names, clades = data$clade_idxs[data[["clade_idxs"]]==1], xlab="", ylab=expression(beta[VPD]), 
                                             xticklabs = NA, display_ylim = c(-0.25,0.25), line0 = TRUE, baseline_values = baseline_median)
par(mar = c(1,0,1,1))
names <- sapply(which(data[["clade_idxs"]]==2),
                function(sp) paste0('beta_vpd[', sp, ']'))
baseline_median <- sapply(names, function(n)  util$ensemble_mcmc_quantile_est(samples[[n]], 0.5))
util$plot_disc_pushforward_quantiles_2clades(samples_int, names, clades = data$clade_idxs[data[["clade_idxs"]]==2], xlab="", ylab="", 
                                             xticklabs = NA, display_ylim = c(-0.25,0.25), line0 = TRUE, yticklabs = NA, baseline_values = baseline_median)
par(mar = c(1,2,1,1))
util$plot_expectand_pushforward_reverse(samples_int[['mu_vpd[1]']], 50,
                                        flim = c(-0.25,0.25),
                                        display_name=expression(beta[gdd]),
                                        col = '#27278f')
util$plot_expectand_pushforward_reverse(samples_int[['mu_vpd[2]']], 50,
                                        flim = c(-0.25,0.25),
                                        display_name=expression(beta[gdd]),
                                        col = '#278f27', add = TRUE)

par(mar = c(1,5,1,1), cex.lab = 2)
names <- sapply(which(data[["clade_idxs"]]==1),
                function(sp) paste0('beta_sm[', sp, ']'))
baseline_median <- sapply(names, function(n)  util$ensemble_mcmc_quantile_est(samples[[n]], 0.5))
util$plot_disc_pushforward_quantiles_2clades(samples_int, names, clades = data$clade_idxs[data[["clade_idxs"]]==1], xlab="", ylab=expression(beta[SM]), 
                                             xticklabs = NA, display_ylim = c(-0.25,0.25), line0 = TRUE, baseline_values = baseline_median)
par(mar = c(1,0,1,1))
names <- sapply(which(data[["clade_idxs"]]==2),
                function(sp) paste0('beta_sm[', sp, ']'))
baseline_median <- sapply(names, function(n)  util$ensemble_mcmc_quantile_est(samples[[n]], 0.5))
util$plot_disc_pushforward_quantiles_2clades(samples_int, names, clades = data$clade_idxs[data[["clade_idxs"]]==2], xlab="", ylab="", 
                                             xticklabs = NA, display_ylim = c(-0.25,0.25), line0 = TRUE, yticklabs = NA, baseline_values = baseline_median)
par(mar = c(1,2,1,1))
util$plot_expectand_pushforward_reverse(samples_int[['mu_sm[1]']], 50,
                                        flim = c(-0.25,0.25),
                                        display_name=expression(beta[gdd]),
                                        col = '#27278f')
util$plot_expectand_pushforward_reverse(samples_int[['mu_sm[2]']], 50,
                                        flim = c(-0.25,0.25),
                                        display_name=expression(beta[gdd]),
                                        col = '#278f27', add = TRUE)


par(mar = c(1,5,1,1), cex.lab = 2)
names <- sapply(which(data[["clade_idxs"]]==1),
                function(sp) paste0('beta_sm_vpd[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(samples_int, names, clades = data$clade_idxs[data[["clade_idxs"]]==1], xlab="", ylab=expression(beta[SM.VPD]), 
                                             xticklabs = NA, display_ylim = c(-0.25,0.25), line0 = TRUE)
par(mar = c(1,0,1,1))
names <- sapply(which(data[["clade_idxs"]]==2),
                function(sp) paste0('beta_sm_vpd[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(samples_int, names, clades = data$clade_idxs[data[["clade_idxs"]]==2], xlab="", ylab="", 
                                             xticklabs = NA, display_ylim = c(-0.25,0.25), line0 = TRUE, yticklabs = NA)
par(mar = c(1,2,1,1))
util$plot_expectand_pushforward_reverse(samples_int[['mu_sm_vpd[1]']], 50,
                                        flim = c(-0.25,0.25),
                                        display_name=expression(beta[gdd]),
                                        col = '#27278f')
util$plot_expectand_pushforward_reverse(samples_int[['mu_sm_vpd[2]']], 50,
                                        flim = c(-0.25,0.25),
                                        display_name=expression(beta[gdd]),
                                        col = '#278f27', add = TRUE)

par(mar = c(1,5,1,1), cex.lab = 2)
names <- sapply(which(data[["clade_idxs"]]==1),
                function(sp) paste0('beta_sm_vpd[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(samples_int, names, clades = data$clade_idxs[data[["clade_idxs"]]==1], xlab="", ylab=expression(beta[SM.VPD]), 
                                             xticklabs = NA, display_ylim = c(-0.05,0.05), line0 = TRUE)
text(x = sum(data[["clade_idxs"]]==2)/2+0.5, y = 0.04, label = "Zoom", col = 'black', cex = 1.5)
par(mar = c(1,0,1,1))
names <- sapply(which(data[["clade_idxs"]]==2),
                function(sp) paste0('beta_sm_vpd[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(samples_int, names, clades = data$clade_idxs[data[["clade_idxs"]]==2], xlab="", ylab="", 
                                             xticklabs = NA, display_ylim = c(-0.05,0.05), line0 = TRUE, yticklabs = NA)
par(mar = c(1,2,1,1))
util$plot_expectand_pushforward_reverse(samples_int[['mu_sm_vpd[1]']], 50,
                                        flim = c(-0.05,0.05),
                                        display_name=expression(beta[gdd]),
                                        col = '#27278f')
util$plot_expectand_pushforward_reverse(samples_int[['mu_sm_vpd[2]']], 50,
                                        flim = c(-0.05,0.05),
                                        display_name=expression(beta[gdd]),
                                        col = '#278f27', add = TRUE)



dev.off()


