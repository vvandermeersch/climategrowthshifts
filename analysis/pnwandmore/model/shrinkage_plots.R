rm(list = ls())
library(ggplot2)
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)

data <- readRDS(file.path(wd, 'output/model', 'data_11july2025.rds'))
slopes_partial <- readRDS(file.path(wd, 'output/model/shrinkage', 'slopes_11july2025_partialpooling.rds'))
slopes_nopooling <- readRDS(file.path(wd, 'output/model/shrinkage', 'slopes_11july2025_nopooling.rds'))

n_obs_perspecies <- sapply(1:data$N_species, function(s) sum(data$species_idxs==s))

var <- 'beta_vpd'
shrinkage_df <- rbind(
  data.frame(
    species = 1:data$N_species,
    clade = ifelse(data$clade_idxs == 1, 'gymnosperms', 'angiosperms'),
    q10 = sapply(1:data$N_species, function(s) util$ensemble_mcmc_quantile_est(slopes_partial[[paste0(var,'[',s,']')]], c(0.1))),
    q50 = sapply(1:data$N_species, function(s) util$ensemble_mcmc_quantile_est(slopes_partial[[paste0(var,'[',s,']')]], c(0.5))),
    q90 = sapply(1:data$N_species, function(s) util$ensemble_mcmc_quantile_est(slopes_partial[[paste0(var,'[',s,']')]], c(0.9))),
    pooling = 'Partial pooling'),
  data.frame(
    species = 1:data$N_species,
    clade = ifelse(data$clade_idxs == 1, 'gymnosperms', 'angiosperms'),
    q10 = sapply(1:data$N_species, function(s) util$ensemble_mcmc_quantile_est(slopes_nopooling[[paste0(var,'[',s,']')]], c(0.1))),
    q50 = sapply(1:data$N_species, function(s) util$ensemble_mcmc_quantile_est(slopes_nopooling[[paste0(var,'[',s,']')]], c(0.5))),
    q90 = sapply(1:data$N_species, function(s) util$ensemble_mcmc_quantile_est(slopes_nopooling[[paste0(var,'[',s,']')]], c(0.9))),
    pooling = 'No pooling')
)

unc_partial <- shrinkage_df[shrinkage_df$pooling == 'Partial pooling',]$q90-shrinkage_df[shrinkage_df$pooling == 'Partial pooling',]$q10
unc_nopool <- shrinkage_df[shrinkage_df$pooling == 'No pooling',]$q90-shrinkage_df[shrinkage_df$pooling == 'No pooling',]$q10 
unc_var <- (unc_nopool-unc_partial)/unc_nopool*100
shrinkage_df$unc_var <- rep(unc_var, 2)

slop_partial <- shrinkage_df[shrinkage_df$pooling == 'Partial pooling',]$q50
slop_nopool <- shrinkage_df[shrinkage_df$pooling == 'No pooling',]$q50
slop_var <- (slop_nopool-slop_partial)
shrinkage_df$slop_var <- rep(slop_var, 2)

ggplot(data = shrinkage_df[shrinkage_df$unc_var >= 50,]) + 
  facet_wrap(~ clade) + 
  geom_vline(xintercept = 1.5, linetype = 'dotted', color = 'grey50') +
  geom_ribbon(aes(x = pooling, ymin = q10, ymax = q90, group = species),
              position = position_dodge(width = 0.5), alpha = 0.1) +
  geom_pointrange(aes(x = pooling, y = q50, ymin = q10, ymax = q90, group = species, color = as.character(species)),
                  position = position_dodge(width = 0.5)) +
  theme_classic() +
  theme(legend.position = 'none') +
  coord_cartesian(xlim = c(1.25, 1.75)) +
  labs(y = var, x = '')

ggplot(data = shrinkage_df[shrinkage_df$slop_var >= 0.1,]) + 
  facet_wrap(~ clade) + 
  geom_vline(xintercept = 1.5, linetype = 'dotted', color = 'grey50') +
  geom_ribbon(aes(x = pooling, ymin = q10, ymax = q90, group = species),
              position = position_dodge(width = 0.5), alpha = 0.1) +
  geom_pointrange(aes(x = pooling, y = q50, ymin = q10, ymax = q90, group = species, color = as.character(species)),
                  position = position_dodge(width = 0.5)) +
  theme_classic() +
  theme(legend.position = 'none') +
  coord_cartesian(xlim = c(1.25, 1.75)) +
  labs(y = var, x = '')


plot(abs(slop_var) ~ n_obs_perspecies)
