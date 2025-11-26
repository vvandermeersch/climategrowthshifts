rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('functions/custom_functions.R', local=util)
setwd(wd)

# Model data
data <- readRDS(file.path(wd, 'output/model', 'data_11july2025_soilmoisture2m.rds'))

#
fullfit <- FALSE
fullsamples <- FALSE
if(fullfit){
  
}else if (fullsamples){
  
}else{
  base_samples_1m <- readRDS(file.path(wd, "output/model/samples_11july2025_partialpooling_2clades_centered_extended.rds"))
  base_samples_2m <- readRDS(file.path(wd, "output/model/samples_11july2025_partialpooling_2clades_centered_extended_soilmoisture2m.rds"))
  base_samples_wpre <- readRDS(file.path(wd, "output/model/samples_11july2025_partialpooling_2clades_centered_extended_wpre.rds"))
  
  base_samples_2m <- readRDS(file.path(wd, "output/model/samples_24sept2025_partialpooling_2clades_centered_extended.rds"))
  base_samples_2m_gddamjjas <- readRDS(file.path(wd, "output/model/samples_07oct2025_partialpooling_2clades_gdd_amjjas.rds"))
}

util$check_all_expectand_diagnostics(base_samples_2m)
util$check_all_expectand_diagnostics(base_samples_2m_gddamjjas)

# Look at one model estimates
base_samples <- base_samples_2m_gddamjjas
par(mfrow=c(1, 1), mar = c(4,4,1,1))
names <- sapply(1:data$N_species,
                function(sp) paste0('rho_sp[', sp, ']'))
util$plot_disc_pushforward_quantiles(base_samples, names,
                                     xlab="Species",
                                     ylab="Rho")

par(mfrow=c(1, 1), mar = c(4,4,1,1))
names <- sapply(1:data$N_species,
                function(sp) paste0('gamma_sp[', sp, ']'))
util$plot_disc_pushforward_quantiles(base_samples, names,
                                     xlab="Species",
                                     ylab="Gamma")

par(mfrow=c(1, 1), mar = c(4,4,1,1))
names <- sapply(1:data$N_species,
                function(sp) paste0('kappa_sh[', sp, ']'))
util$plot_disc_pushforward_quantiles(base_samples, names,
                                     xlab="Species",
                                     ylab="Kappa")

par(mfrow=c(3, 1), mar = c(4,4,1,1))
names <- sapply(1:data$N_species,
                function(sp) paste0('beta_gdd[', sp, ']'))
util$plot_disc_pushforward_quantiles(base_samples, names,
                                     xlab="Species",
                                     ylab="beta_gdd")
names <- sapply(1:data$N_species,
                function(sp) paste0('beta_sm[', sp, ']'))
util$plot_disc_pushforward_quantiles(base_samples, names,
                                     xlab="Species",
                                     ylab="beta_sm")
names <- sapply(1:data$N_species,
                function(sp) paste0('beta_vpd[', sp, ']'))
util$plot_disc_pushforward_quantiles(base_samples, names,
                                     xlab="Species",
                                     ylab="beta_vpd")

# Clade-level estimates
par(mfrow=c(1, 3), mar = c(4,4,1,1))
util$plot_expectand_pushforward(base_samples[['mu_gdd[1]']], 50,
                                flim = c(-0.15,0.15),
                                display_name=expression(beta[GDD]),
                                col = '#27278f')
util$plot_expectand_pushforward(base_samples[['mu_gdd[2]']], 50,
                                flim = c(-0.15,0.15),
                                col = '#278f27',
                                add = TRUE)
text(x = 0.115, y = 50, label = "Gymno.", col = '#27278f', cex = 1.2)
text(x = -0.05, y = 10, label = "Angio.", col = '#278f27', cex = 1.2)

util$plot_expectand_pushforward(base_samples[['mu_sm[1]']], 50,
                                flim = c(-0.15,0.15),
                                display_name=expression(beta[SM]),
                                col = '#27278f')
util$plot_expectand_pushforward(base_samples[['mu_sm[2]']], 50,
                                flim = c(-0.15,0.15),
                                col = '#278f27',
                                add = TRUE)

util$plot_expectand_pushforward(base_samples[['mu_vpd[1]']], 50,
                                flim = c(-0.15,0.15),
                                display_name=expression(beta[VPD]),
                                col = '#27278f')
util$plot_expectand_pushforward(base_samples[['mu_vpd[2]']], 50,
                                flim =c(-0.15,0.15),
                                col = '#278f27',
                                add = TRUE)


par(mfrow=c(1, 3), mar = c(4,4,1,1))
util$plot_expectand_pushforward(base_samples[['mu_rho[1]']], 50,
                                flim = c(1,23),
                                display_name=expression(rho),
                                col = '#27278f')
util$plot_expectand_pushforward(base_samples[['mu_rho[2]']], 50,
                                flim = c(1,23),
                                col = '#278f27',
                                add = TRUE)
text(x = 5, y = 0.5, label = "Gymno.", col = '#27278f', cex = 1.2)
text(x = 17, y = 0.28, label = "Angio.", col = '#278f27', cex = 1.2)

util$plot_expectand_pushforward(base_samples[['mu_gamma[1]']], 50,
                                flim = c(0.2,0.8),
                                display_name=expression(gamma),
                                col = '#27278f')
util$plot_expectand_pushforward(base_samples[['mu_gamma[2]']], 50,
                                flim = c(0.2,0.8),
                                col = '#278f27',
                                add = TRUE)

util$plot_expectand_pushforward(base_samples[['mu_kappa[1]']], 50,
                                flim = c(0,1),
                                display_name=expression(kappa),
                                col = '#27278f')
util$plot_expectand_pushforward(base_samples[['mu_kappa[2]']], 50,
                                flim = c(0,1),
                                col = '#278f27',
                                add = TRUE)

# Species-level estimates, displayed by clade
par(mfrow=c(3, 2), mar = c(4,4,1,1))
init <- FALSE
for(s in which(data[["clade_idxs"]]==1)){
  util$plot_expectand_pushforward(base_samples[[paste0('beta_gdd[',s,']')]], 100,
                                  flim = c(-0.3,0.25), ylim = c(0,200),
                                  display_name=expression(beta[GDD]),
                                  col = '#27278f', add = init)
  init <- TRUE
}
init <- FALSE
for(s in which(data[["clade_idxs"]]==2)){
  util$plot_expectand_pushforward(base_samples[[paste0('beta_gdd[',s,']')]], 100,
                                  flim = c(-0.3,0.25), ylim = c(0,200),
                                  display_name=expression(beta[GDD]),
                                  col = '#278f27', add = init)
  init <- TRUE
}
init <- FALSE
for(s in which(data[["clade_idxs"]]==1)){
  util$plot_expectand_pushforward(base_samples[[paste0('beta_sm[',s,']')]], 100,
                                  flim = c(-0.3,0.25), ylim = c(0,120),
                                  display_name=expression(beta[SM]),
                                  col = '#27278f', add = init)
  init <- TRUE
}
init <- FALSE
for(s in which(data[["clade_idxs"]]==2)){
  util$plot_expectand_pushforward(base_samples[[paste0('beta_sm[',s,']')]], 100,
                                  flim = c(-0.3,0.25), ylim = c(0,120),
                                  display_name=expression(beta[SM]),
                                  col = '#278f27', add = init)
  init <- TRUE
}
init <- FALSE
for(s in which(data[["clade_idxs"]]==1)){
  util$plot_expectand_pushforward(base_samples[[paste0('beta_vpd[',s,']')]], 100,
                                  flim = c(-0.3,0.25), ylim = c(0,100),
                                  display_name=expression(beta[VPD]),
                                  col = '#27278f', add = init)
  init <- TRUE
}
init <- FALSE
for(s in which(data[["clade_idxs"]]==2)){
  util$plot_expectand_pushforward(base_samples[[paste0('beta_vpd[',s,']')]], 100,
                                  flim = c(-0.3,0.25), ylim = c(0,100),
                                  display_name=expression(beta[VPD]),
                                  col = '#278f27', add = init)
  init <- TRUE
}

# Species-level estimates, with clade estimates
layout(mat = matrix(c(1, 2, 3), ncol = 3),
       heights = c(2,1), widths = c(3, 4*sum(data[["clade_idxs"]]==2)/length(data[["clade_idxs"]]), 0.7))  
par(mar = c(1,4,1,1))
names <- sapply(which(data[["clade_idxs"]]==1),
                function(sp) paste0('rho_sp[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(base_samples, names, clades = data$clade_idxs[data[["clade_idxs"]]==1], xlab="", ylab=expression(rho), 
                                             xticklabs = NA, display_ylim = c(0,35), line0 = FALSE)
text(x = sum(data[["clade_idxs"]]==1)/2, y = 25, label = "Gymnosperms", col = '#27278f', cex = 1.2)
par(mar = c(1,0,1,1))
names <- sapply(which(data[["clade_idxs"]]==2),
                function(sp) paste0('rho_sp[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(base_samples, names, clades = data$clade_idxs[data[["clade_idxs"]]==2], xlab="", ylab="", 
                                             xticklabs = NA, display_ylim = c(0,35), yticklabs = NA)
text(x = sum(data[["clade_idxs"]]==2)/2, y = 25, label = "Angiosperms", col = '#278f27', cex = 1.2)
par(mar = c(1,2,1,1))
util$plot_expectand_pushforward_reverse(base_samples[['mu_rho[1]']], 50,
                                        flim = c(0,35),
                                        display_name=expression(beta[GDD]),
                                        col = '#27278f')
util$plot_expectand_pushforward_reverse(base_samples[['mu_rho[2]']], 50,
                                        flim = c(0,35),
                                        display_name=expression(beta[GDD]),
                                        col = '#278f27', add = TRUE)



par(mar = c(1,0,1,1))
names <- sapply(which(data[["clade_idxs"]]==2),
                function(sp) paste0('gamma_sp[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(base_samples, names, clades = data$clade_idxs[data[["clade_idxs"]]==2], xlab="", ylab="", 
                                             xticklabs = NA, display_ylim = c(0,1), yticklabs = NA)
par(mar = c(1,4,1,1))
names <- sapply(which(data[["clade_idxs"]]==1),
                function(sp) paste0('gamma_sp[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(base_samples, names, clades = data$clade_idxs[data[["clade_idxs"]]==1], xlab="", ylab=expression(gamma), 
                                             xticklabs = rep(NA, sum(data[["clade_idxs"]]==1)), display_ylim = c(0,1) )

par(mar = c(0,0,1,1))
names <- sapply(which(data[["clade_idxs"]]==2),
                function(sp) paste0('beta_gdd[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(base_samples, names, clades = data$clade_idxs[data[["clade_idxs"]]==2], xlab="", ylab="", 
                                             xticklabs = NA, display_ylim = c(-0.07,0.2), yticklabs = NA)
par(mar = c(0,4.5,1,1))
names <- sapply(which(data[["clade_idxs"]]==1),
                function(sp) paste0('beta_gdd[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(base_samples, names, clades = data$clade_idxs[data[["clade_idxs"]]==1], xlab="", ylab=expression(beta[GDD]), 
                                             xticklabs = rep(NA, sum(data[["clade_idxs"]]==1)), display_ylim = c(-0.07,0.2) )

par(mar = c(0,0,1,1))
names <- sapply(which(data[["clade_idxs"]]==2),
                function(sp) paste0('beta_sm[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(base_samples, names, clades = data$clade_idxs[data[["clade_idxs"]]==2], xlab="", ylab="", 
                                             xticklabs = NA, display_ylim = c(-0.12,0.25), yticklabs = NA)
par(mar = c(0,4.5,1,1))
names <- sapply(which(data[["clade_idxs"]]==1),
                function(sp) paste0('beta_sm[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(base_samples, names, clades = data$clade_idxs[data[["clade_idxs"]]==1], xlab="", ylab=expression(beta[SM]), 
                                             xticklabs = rep(NA, sum(data[["clade_idxs"]]==1)), display_ylim = c(-0.15,0.25) )

layout(mat = matrix(c(1, 2, 3), ncol = 3),
       heights = c(2,1), widths = c(3, 4*sum(data[["clade_idxs"]]==2)/length(data[["clade_idxs"]]), 0.7))  
par(mar = c(1,4.5,1,1))
names <- sapply(which(data[["clade_idxs"]]==1),
                function(sp) paste0('beta_vpd[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(base_samples, names, clades = data$clade_idxs[data[["clade_idxs"]]==1], 
                                             xlab="", ylab=expression(beta[VPD]), 
                                             xticklabs = NA, display_ylim = c(-0.3,0.1))
text(x = sum(data[["clade_idxs"]]==1)/2, y = 0.1, label = "Gymnosperms", col = '#27278f', cex = 1.2)
par(mar = c(1,0,1,1))
names <- sapply(which(data[["clade_idxs"]]==2),
                function(sp) paste0('beta_vpd[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(base_samples, names, clades = data$clade_idxs[data[["clade_idxs"]]==2], xlab="", ylab="", 
                                             xticklabs = NA, display_ylim = c(-0.3,0.1), yticklabs = NA)
text(x = sum(data[["clade_idxs"]]==2)/2, y = 0.1, label = "Angiosperms", col = '#278f27', cex = 1.2)
par(mar = c(1,2,1,1))
util$plot_expectand_pushforward_reverse(base_samples[['mu_vpd[1]']], 50,
                                        flim = c(-0.3,0.1),
                                        display_name=expression(beta[GDD]),
                                        col = '#27278f')
util$plot_expectand_pushforward_reverse(base_samples[['mu_vpd[2]']], 50,
                                        flim = c(-0.3,0.1),
                                        display_name=expression(beta[GDD]),
                                        col = '#278f27', add = TRUE)


layout(mat = matrix(c(1, 2, 3), ncol = 3),
       heights = c(2,1), widths = c(3, 4*sum(data[["clade_idxs"]]==2)/length(data[["clade_idxs"]]), 0.7))  
par(mar = c(1,4.5,1,1))
names <- sapply(which(data[["clade_idxs"]]==1),
                function(sp) paste0('beta_sm[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(base_samples, names, clades = data$clade_idxs[data[["clade_idxs"]]==1], 
                                             xlab="", ylab=expression(beta[SM]), 
                                             xticklabs = NA, display_ylim = c(-0.1,0.3))
text(x = sum(data[["clade_idxs"]]==1)/2, y = 0.3, label = "Gymnosperms", col = '#27278f', cex = 1.2)
par(mar = c(1,0,1,1))
names <- sapply(which(data[["clade_idxs"]]==2),
                function(sp) paste0('beta_sm[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(base_samples, names, clades = data$clade_idxs[data[["clade_idxs"]]==2], xlab="", ylab="", 
                                             xticklabs = NA, display_ylim = c(-0.1,0.3), yticklabs = NA)
text(x = sum(data[["clade_idxs"]]==2)/2, y = 0.3, label = "Angiosperms", col = '#278f27', cex = 1.2)
par(mar = c(1,2,1,1))
util$plot_expectand_pushforward_reverse(base_samples[['mu_sm[1]']], 50,
                                        flim = c(-0.1,0.3),
                                        display_name=expression(beta[SM]),
                                        col = '#27278f')
util$plot_expectand_pushforward_reverse(base_samples[['mu_sm[2]']], 50,
                                        flim = c(-0.1,0.3),
                                        display_name=expression(beta[SM]),
                                        col = '#278f27', add = TRUE)


# Comparison of the two soil-moisture models
baseline_1m_median <- sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples_1m[[paste0('beta_sm[', sp, ']')]], 0.5))

layout(mat = matrix(c(1, 2), ncol = 2),
       heights = c(2,1), widths = c(3, 4*sum(data[["clade_idxs"]]==2)/length(data[["clade_idxs"]])))  
par(mar = c(1,4.5,1,1))
names <- sapply(which(data[["clade_idxs"]]==1),
                function(sp) paste0('beta_sm[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(base_samples_2m, names, clades = data$clade_idxs[data[["clade_idxs"]]==1], 
                                             baseline_values = baseline_1m_median[data[["clade_idxs"]]==1],
                                             baseline_col="darkred",
                                             xlab="", ylab=expression(beta[SM]), 
                                             xticklabs = NA, display_ylim = c(-0.1,0.3))
par(mar = c(1,0,1,1))
names <- sapply(which(data[["clade_idxs"]]==2),
                function(sp) paste0('beta_sm[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(base_samples_2m, names, clades = data$clade_idxs[data[["clade_idxs"]]==2], xlab="", ylab="", 
                                             baseline_values = baseline_1m_median[data[["clade_idxs"]]==2],
                                             baseline_col="darkred",
                                             xticklabs = NA, display_ylim = c(-0.1,0.3), yticklabs = NA)


# Now look at winter precipitation
layout(mat = matrix(c(1, 2, 3), ncol = 3),
       heights = c(2,1), widths = c(3, 4*sum(data[["clade_idxs"]]==2)/length(data[["clade_idxs"]]), 0.7))  
par(mar = c(1,4.5,1,1))
names <- sapply(which(data[["clade_idxs"]]==1),
                function(sp) paste0('beta_wpre[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(base_samples_wpre, names, clades = data$clade_idxs[data[["clade_idxs"]]==1], 
                                             xlab="", ylab=expression(beta[winterprecipitation]), 
                                             xticklabs = NA, display_ylim = c(-0.1,0.3))
text(x = sum(data[["clade_idxs"]]==1)/2, y = 0.3, label = "Gymnosperms", col = '#27278f', cex = 1.2)
par(mar = c(1,0,1,1))
names <- sapply(which(data[["clade_idxs"]]==2),
                function(sp) paste0('beta_wpre[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(base_samples_wpre, names, clades = data$clade_idxs[data[["clade_idxs"]]==2], xlab="", ylab="", 
                                             xticklabs = NA, display_ylim = c(-0.1,0.3), yticklabs = NA)
text(x = sum(data[["clade_idxs"]]==2)/2, y = 0.3, label = "Angiosperms", col = '#278f27', cex = 1.2)
par(mar = c(1,2,1,1))
util$plot_expectand_pushforward_reverse(base_samples_wpre[['mu_wpre[1]']], 50,
                                        flim = c(-0.1,0.3),
                                        display_name=expression(beta[SM]),
                                        col = '#27278f')
util$plot_expectand_pushforward_reverse(base_samples_wpre[['mu_wpre[2]']], 50,
                                        flim = c(-0.1,0.3),
                                        display_name=expression(beta[SM]),
                                        col = '#278f27', add = TRUE)


# Comparison of the two GDD models
# I removed two species in the new GDD model...
# prevdata <- readRDS("/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model/data_24sept2025.rds")
data <- readRDS("/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model/data_07oct2025.rds")
# prevspecies <- which(prevdata$uniq_species_ids %in% data$uniq_species_ids)
# EDIT: but now I have run again the old GDD model without the 2 species
base_samples_gdd <- readRDS(file.path(wd, "output/model/samples_07oct2025_partialpooling_2clades.rds"))

baseline_median <- sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples_gdd[[paste0('beta_gdd[', sp, ']')]], 0.5))

in2mm <-25.4 # scale factor to convert inches to mm
pdf(file = file.path(wd, 'figures', 'newyork2025', 'model_estimates', 'gdd_window_comparisons.pdf'),
    width=250/in2mm,height=150/in2mm, paper="special")
layout(mat = matrix(c(1, 2), ncol = 2),
       heights = c(2,1), widths = c(3, 4*sum(data[["clade_idxs"]]==2)/length(data[["clade_idxs"]])))  
par(mar = c(1,4.5,1,1))
names <- sapply(which(data[["clade_idxs"]]==1),
                function(sp) paste0('beta_gdd[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(base_samples_2m_gddamjjas, names, clades = data$clade_idxs[data[["clade_idxs"]]==1], 
                                             baseline_values = baseline_median[data[["clade_idxs"]]==1],
                                             baseline_col="darkred",
                                             xlab="", ylab=expression(beta[GDD]), 
                                             xticklabs = NA, display_ylim = c(-0.1,0.3))
text(x = sum(data[["clade_idxs"]]==1)/4, y = 0.3, label = "Red: estimates from full-year GDD", col = 'darkred', cex = 1)
par(mar = c(1,0,1,1))
names <- sapply(which(data[["clade_idxs"]]==2),
                function(sp) paste0('beta_gdd[', sp, ']'))
util$plot_disc_pushforward_quantiles_2clades(base_samples_2m_gddamjjas, names, clades = data$clade_idxs[data[["clade_idxs"]]==2], xlab="", ylab="", 
                                             baseline_values = baseline_median[data[["clade_idxs"]]==2],
                                             baseline_col="darkred",
                                             xticklabs = NA, display_ylim = c(-0.1,0.3), yticklabs = NA)
dev.off()

