rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('functions/custom_functions.R', local=util)
setwd(wd)

# Model data
data <- readRDS(file.path(wd, 'output/model', 'data_24sept2025.rds'))

base_samples <- readRDS(file.path(wd, "output/model/samples_24sept2025_partialpooling_2clades_centered_extended_gdd_amjjas.rds"))


betas_climate <- data.frame(
  species = data$uniq_species_ids,
  clade = data$clade_idxs,
  betagdd_q5 = sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('beta_gdd[', sp, ']')]], 0.05)),
  betagdd_q50 = sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('beta_gdd[', sp, ']')]], 0.5)),
  betagdd_q95= sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('beta_gdd[', sp, ']')]], 0.95)),
  betasm_q5 = sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('beta_sm[', sp, ']')]], 0.05)),
  betasm_q50 = sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('beta_sm[', sp, ']')]], 0.5)),
  betasm_q95= sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('beta_sm[', sp, ']')]], 0.95)),
  betavpd_q5 = sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('beta_vpd[', sp, ']')]], 0.05)),
  betavpd_q50 = sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('beta_vpd[', sp, ']')]], 0.5)),
  betavpd_q95= sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('beta_vpd[', sp, ']')]], 0.95))
)



ggplot(data = betas_climate) +
  geom_vline(aes(xintercept = 0), linetype = 'dotted') + 
  geom_hline(aes(yintercept = 0), linetype = 'dotted') + 
  geom_pointrange(aes(x=betagdd_q50,
                      ymin=betasm_q5, y=betasm_q50, ymax=betasm_q95,
                      color = as.character(clade)), size = 0.2) +
  geom_pointrange(aes(xmin=betagdd_q5, x=betagdd_q50, xmax=betagdd_q95,
                      y=betasm_q50,
                      color = as.character(clade)), size = 0.2) +
  scale_color_manual(values = c('#27278f', '#278f27')) +
  # coord_equal(xlim = c(-0.1,0.3), ylim = c(-0.1,0.3)) +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(x = "GDD", y = "Soil moisture")

ggplot(data = betas_climate) +
  geom_vline(aes(xintercept = 0), linetype = 'dotted') + 
  geom_hline(aes(yintercept = 0), linetype = 'dotted') + 
  geom_pointrange(aes(x=betagdd_q50,
                      ymin=betavpd_q5, y=betavpd_q50, ymax=betavpd_q95,
                      color = as.character(clade)), size = 0.2) +
  geom_pointrange(aes(xmin=betagdd_q5, x=betagdd_q50, xmax=betagdd_q95,
                      y=betavpd_q50,
                      color = as.character(clade)), size = 0.2) +
  scale_color_manual(values = c('#27278f', '#278f27')) +
  # coord_equal(xlim = c(-0.1,0.3), ylim = c(-0.1,0.3)) +
  theme_bw() +
  #coord_cartesian(xlim = c(-0.05,0.15), ylim = c(-0.1,0.35), expand = FALSE) +
  theme(legend.position = 'none') +
  labs(x = "GDD", y = "VPD")

ggplot(data = betas_climate) +
  geom_vline(aes(xintercept = 0), linetype = 'dotted') + 
  geom_hline(aes(yintercept = 0), linetype = 'dotted') + 
  geom_pointrange(aes(x=betasm_q50,
                      ymin=betavpd_q5, y=betavpd_q50, ymax=betavpd_q95,
                      color = as.character(clade)), size = 0.2) +
  geom_pointrange(aes(xmin=betasm_q5, x=betasm_q50, xmax=betasm_q95,
                      y=betavpd_q50,
                      color = as.character(clade)), size = 0.2) +
  scale_color_manual(values = c('#27278f', '#278f27')) +
  # coord_equal(xlim = c(-0.1,0.3), ylim = c(-0.1,0.3)) +
  theme_bw() +
  #coord_cartesian(xlim = c(-0.05,0.15), ylim = c(-0.1,0.35), expand = FALSE) +
  theme(legend.position = 'none') +
  labs(x = "Soil moisture", y = "VPD")



# What if we look at the climate we observed?
betas_climate <- merge(betas_climate, species_climate_obs)

ggplot(data = betas_climate) +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') + 
  geom_pointrange(aes(x=gdd_obs_q50,
                      ymin=betagdd_q5, y=betagdd_q50, ymax=betagdd_q95,
                      color = as.character(clade)), size = 0.2) +
  geom_pointrange(aes(xmin=gdd_obs_q5, x=gdd_obs_q50, xmax=gdd_obs_q95,
                      y=betagdd_q50,
                      color = as.character(clade)), size = 0.2) +
  scale_color_manual(values = c('#27278f', '#278f27')) +
  # coord_equal(xlim = c(-0.1,0.3), ylim = c(-0.1,0.3)) +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(x = "GDD (obs. only)", y = "Beta GDD")

ggplot(data = betas_climate) +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') + 
  geom_pointrange(aes(x=sm_obs_q50,
                      ymin=betasm_q5, y=betasm_q50, ymax=betasm_q95,
                      color = as.character(clade)), size = 0.2) +
  geom_pointrange(aes(xmin=sm_obs_q5, x=sm_obs_q50, xmax=sm_obs_q95,
                      y=betasm_q50,
                      color = as.character(clade)), size = 0.2) +
  scale_color_manual(values = c('#27278f', '#278f27')) +
  # coord_equal(xlim = c(-0.1,0.3), ylim = c(-0.1,0.3)) +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(x = "SM (obs. only)", y = "Beta SM")

ggplot(data = betas_climate) +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') + 
  geom_pointrange(aes(x=vpd_obs_q50,
                      ymin=betavpd_q5, y=betavpd_q50, ymax=betavpd_q95,
                      color = as.character(clade)), size = 0.2) +
  geom_pointrange(aes(xmin=vpd_obs_q5, x=vpd_obs_q50, xmax=vpd_obs_q95,
                      y=betavpd_q50,
                      color = as.character(clade)), size = 0.2) +
  scale_color_manual(values = c('#27278f', '#278f27')) +
  # coord_equal(xlim = c(-0.1,0.3), ylim = c(-0.1,0.3)) +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(x = "VPD (obs. only)", y = "Beta VPD")


ggplot(data = betas_climate) +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') + 
  geom_pointrange(aes(x=gdd_obs_q95-gdd_obs_q5,
                      ymin=betagdd_q5, y=betagdd_q50, ymax=betagdd_q95,
                      color = as.character(clade)), size = 0.2) +
  scale_color_manual(values = c('#27278f', '#278f27')) +
  # coord_equal(xlim = c(-0.1,0.3), ylim = c(-0.1,0.3)) +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(x = "GDD (obs. only)", y = "Beta GDD")




# Climatic range by number of stands
ggplot(data = betas_climate) +
  # geom_pointrange(aes(ymin=gdd_obs_q5, y=gdd_obs_q50, ymax=gdd_obs_q95,
  #                     x=no_datasets,
  #                     color = as.character(clade)), size = 0.2) +
  geom_point(aes(y=gdd_obs_q95-gdd_obs_q5, x=no_datasets,
                 color = as.character(clade)), size = 1.5) +
  scale_color_manual(values = c('#27278f', '#278f27')) +
  # coord_equal(xlim = c(-0.1,0.3), ylim = c(-0.1,0.3)) +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(x = "No. of datasets", y = "GDD (obs. only), Q95-Q5") +
  scale_x_log10() 

ggplot(data = betas_climate) +
  geom_point(aes(y=sm_obs_q95-sm_obs_q5, x=no_datasets,
                 color = as.character(clade)), size = 1.5) +
  scale_color_manual(values = c('#27278f', '#278f27')) +
  # coord_equal(xlim = c(-0.1,0.3), ylim = c(-0.1,0.3)) +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(x = "No. of datasets", y = "SM (obs. only), Q95-Q5") +
  scale_x_log10() 

ggplot(data = betas_climate) +
  geom_point(aes(y=vpd_obs_q95-vpd_obs_q5, x=no_datasets,
                 color = as.character(clade)), size = 1.5) +
  scale_color_manual(values = c('#27278f', '#278f27')) +
  # coord_equal(xlim = c(-0.1,0.3), ylim = c(-0.1,0.3)) +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(x = "No. of datasets", y = "VPD (obs. only), Q95-Q5") +
  scale_x_log10() 
