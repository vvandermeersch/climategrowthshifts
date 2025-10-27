rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('functions/custom_functions.R', local=util)
setwd(wd)
library(ggplot2)
library(patchwork)

# Model data
data <- readRDS(file.path(wd, 'output/model', 'data_24sept2025.rds'))
samples <- readRDS(file.path(wd, 'output/model', 'samples_24sept2025_partialpooling2clades_1interaction.rds'))
samples_noint <- readRDS(file.path(wd, 'output/model', 'samples_24sept2025_partialpooling_2clades_centered_extended.rds'))

# Baselines defined in the model
sm0 <- 25
vpd0 <- 8

species <- 'QUDG'
species <- data$uniq_species_ids[21]

#### First, look at PARAMETERS! ####

par(mfrow = c(1,3))
util$plot_expectand_pushforward(samples[[paste0('beta_sm[', which(data$uniq_species_ids == species),']')]], 50,
                                flim = c(-0.2,0.2),
                                display_name="beta_sm")
util$plot_expectand_pushforward(samples[[paste0('beta_vpd[', which(data$uniq_species_ids == species),']')]], 50,
                                flim = c(-0.2,0.2),
                                display_name="beta_vpd")
util$plot_expectand_pushforward(samples[[paste0('beta_sm_vpd[', which(data$uniq_species_ids == species),']')]], 50,
                                flim = c(-0.2,0.2),
                                display_name="beta_sm_vpd")

#### MODEL WITH INTERACTION ####

# Predictions across species range
focus_trees <- which(data$species_idxs == which(data$uniq_species_ids == species))
focus_idx <- unlist(sapply(focus_trees, function(t)
  return(data$tree_start_idxs[t]:data$tree_end_idxs[t])))

low_sm <- as.numeric(quantile(data$sm_obs[focus_trees], 0.05)) # 15.54%
high_sm <- as.numeric(quantile(data$sm_obs[focus_trees], 0.95)) # 27.72%
low_vpd <- as.numeric(quantile(data$vpd_obs[focus_trees], 0.05)) # 10.27 hPa
high_vpd <- as.numeric(quantile(data$vpd_obs[focus_trees], 0.95)) # 19.04 hPa

vpd <- c(seq(range(data$vpd_obs[focus_trees])[1], range(data$vpd_obs[focus_trees])[2], 0.1), low_vpd, high_vpd)
rw_pred_lowsm <- lapply(vpd, function(v){
  logrw <- as.numeric(samples[[paste0('alpha')]]) + 
    samples[[paste0('beta_sm[', which(data$uniq_species_ids == species), ']')]]*(low_sm-sm0) +
    samples[[paste0('beta_vpd[', which(data$uniq_species_ids == species), ']')]]*(v-vpd0) +
    samples[[paste0('beta_sm_vpd[', which(data$uniq_species_ids == species), ']')]]*(low_sm-sm0)*(v-vpd0)
  rw <- exp(logrw)
  rw <- c(util$ensemble_mcmc_quantile_est(rw, c(0.1,0.5,0.9)), sm = low_sm, vpd = v)
  return(rw)
})
rw_pred_lowsm <- as.data.frame(do.call(rbind, rw_pred_lowsm))

rw_pred_highsm <- lapply(vpd, function(v){
  logrw <- as.numeric(samples[[paste0('alpha')]]) + 
    samples[[paste0('beta_sm[', which(data$uniq_species_ids == species), ']')]]*(high_sm-sm0) +
    samples[[paste0('beta_vpd[', which(data$uniq_species_ids == species), ']')]]*(v-vpd0) +
    samples[[paste0('beta_sm_vpd[', which(data$uniq_species_ids == species), ']')]]*(high_sm-sm0)*(v-vpd0)
  rw <- exp(logrw)
  rw <- c(util$ensemble_mcmc_quantile_est(rw, c(0.1,0.5,0.9)), sm = high_sm, vpd = v)
  return(rw)
})
rw_pred_highsm <- as.data.frame(do.call(rbind, rw_pred_highsm))

predictions_vpd <- rbind(rw_pred_lowsm, rw_pred_highsm)
predictions_vpd$sm_class <- ifelse(predictions_vpd$sm == low_sm, paste0('Low (', round(low_sm,1), '%)'),
                                   ifelse(predictions_vpd$sm == high_sm, paste0('High (', round(high_sm,1), '%)'), NA))

pred_withint <- ggplot() +
  geom_line(aes(x = vpd, y = `10%`, group = sm_class, color = as.character(sm_class)), 
            linetype = 'dashed', data = predictions_vpd) +
  geom_line(aes(x = vpd, y = `90%`, group = sm_class, color = as.character(sm_class)), 
            linetype = 'dashed', data = predictions_vpd) +
  geom_line(aes(x = vpd, y = `50%`, group = sm_class, color = as.character(sm_class)), data = predictions_vpd) +
  scale_color_manual(values = c('#448c81', '#c35f35'), name='Soil moisture', guide = guide_legend(order = 2)) +
  theme_classic() +
  labs(x = 'VPD (hPa)', y = 'Ring width (mm)') +
  ggtitle('Model with interaction') + 
  theme(legend.position = 'inside', legend.position.inside = c(0.15, 0.2),
        legend.background = element_blank(), plot.title = element_text(size = 11)) +
  ylim(c(0, 1))

#### MODEL WITHOUT INTERACTION ####

# Predictions across species range
focus_trees <- which(data$species_idxs == which(data$uniq_species_ids == species))
focus_idx <- unlist(sapply(focus_trees, function(t)
  return(data$tree_start_idxs[t]:data$tree_end_idxs[t])))

low_sm <- as.numeric(quantile(data$sm_obs[focus_trees], 0.05)) # 15.54%
high_sm <- as.numeric(quantile(data$sm_obs[focus_trees], 0.95)) # 27.72%
low_vpd <- as.numeric(quantile(data$vpd_obs[focus_trees], 0.05)) # 10.27 hPa
high_vpd <- as.numeric(quantile(data$vpd_obs[focus_trees], 0.95)) # 19.04 hPa

vpd <- c(seq(range(data$vpd_obs[focus_trees])[1], range(data$vpd_obs[focus_trees])[2], 0.1), low_vpd, high_vpd)
rw_pred_lowsm <- lapply(vpd, function(v){
  logrw <- as.numeric(samples_noint[[paste0('alpha')]]) + 
    samples_noint[[paste0('beta_sm[', which(data$uniq_species_ids == species), ']')]]*(low_sm-sm0) +
    samples_noint[[paste0('beta_vpd[', which(data$uniq_species_ids == species), ']')]]*(v-vpd0)
  rw <- exp(logrw)
  rw <- c(util$ensemble_mcmc_quantile_est(rw, c(0.1,0.5,0.9)), sm = low_sm, vpd = v)
  return(rw)
})
rw_pred_lowsm <- as.data.frame(do.call(rbind, rw_pred_lowsm))

rw_pred_highsm <- lapply(vpd, function(v){
  logrw <- as.numeric(samples_noint[[paste0('alpha')]]) + 
    samples_noint[[paste0('beta_sm[', which(data$uniq_species_ids == species), ']')]]*(high_sm-sm0) +
    samples_noint[[paste0('beta_vpd[', which(data$uniq_species_ids == species), ']')]]*(v-vpd0) 
  rw <- exp(logrw)
  rw <- c(util$ensemble_mcmc_quantile_est(rw, c(0.1,0.5,0.9)), sm = high_sm, vpd = v)
  return(rw)
})
rw_pred_highsm <- as.data.frame(do.call(rbind, rw_pred_highsm))

predictions_vpd <- rbind(rw_pred_lowsm, rw_pred_highsm)
predictions_vpd$sm_class <- ifelse(predictions_vpd$sm == low_sm, paste0('Low (', round(low_sm,1), '%)'),
                             ifelse(predictions_vpd$sm == high_sm, paste0('High (', round(high_sm,1), '%)'), NA))

pred_withoutint <- ggplot() +
  geom_line(aes(x = vpd, y = `10%`, group = sm_class, color = as.character(sm_class)), 
            linetype = 'dashed', data = predictions_vpd) +
  geom_line(aes(x = vpd, y = `90%`, group = sm_class, color = as.character(sm_class)), 
            linetype = 'dashed', data = predictions_vpd) +
  geom_line(aes(x = vpd, y = `50%`, group = sm_class, color = as.character(sm_class)), data = predictions_vpd) +
  scale_color_manual(values = c('#448c81', '#c35f35'), name='Soil moisture', guide = guide_legend(order = 2)) +
  theme_classic() +
  labs(x = 'VPD (hPa)', y = 'Ring width (mm)') +
  ggtitle('Model without interaction') + 
  theme(legend.position = 'none' , plot.title = element_text(size = 11)) +
  ylim(c(0, 1))

pred_withint + pred_withoutint




