rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('functions/custom_functions.R', local=util)
setwd(wd)

library(ggplot2)

data <- readRDS(file.path(wd, 'output/model', 'data_24sept2025.rds'))
datasets <- readRDS(file.path(wd, 'output/model', 'datasets_24sept2025.rds'))
base_samples <- readRDS(file.path(wd, "output/model/samples_24sept2025_partialpooling_2clades_centered_extended.rds"))

# Species
species <- unique(datasets[,c('species_code', 'species_name')])
species[species$species_code == 'PISA', 'species_name'] <- "Pinus sabiniana"
temp <- stringr::str_split_fixed(species$species_name, " ", 3)
species$species_name_simplified <- paste(temp[,1], temp[,2], sep=" ")
species <- species[species$species_code %in% data$uniq_species_ids,]

# Species-level estimate
betas_climate <- data.frame(
  species_code = data$uniq_species_ids,
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
betas_climate <- merge(species, betas_climate)

# Range-wide GDD
species_gddrange <- readRDS(file.path(wd, 'output', 'climate', 'species_gddrange_experienced.rds'))
betas_climate$rangegdd_q5 <- betas_climate$rangegdd_q50 <- betas_climate$rangegdd_q95 <- NA
for(i in 1:nrow(betas_climate)){
  range_gdd <- species_gddrange[[betas_climate[i, 'species_name_simplified']]]
  betas_climate[i, paste0('rangegdd_q', c(5,50,95))] <- as.numeric(quantile(range_gdd, c(0.05,0.5,0.95)))
}
ggplot(data = na.omit(betas_climate)) +
  geom_vline(aes(xintercept = 0), linetype = 'dotted') + 
  geom_hline(aes(yintercept = 0), linetype = 'dotted') + 
  # geom_density_2d(aes(x = rangegdd_q50, y = betagdd_q50, col = as.character(clade), group = clade), alpha = 0.5) + 
  geom_pointrange(aes(x=rangegdd_q50,
                      ymin=betagdd_q5, y=betagdd_q50, ymax=betagdd_q95,
                      color = as.character(clade)), size = 0.2) +
  geom_pointrange(aes(xmin=rangegdd_q5, x=rangegdd_q50, xmax=rangegdd_q95,
                      y=betagdd_q50,
                      color = as.character(clade)), size = 0.2) +
  scale_color_manual(values = c('#27278f', '#278f27')) +
  # coord_equal(xlim = c(-0.1,0.3), ylim = c(-0.1,0.3)) +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(x = "Range-wide GDD (Apr.-Sep., °C)", y = "GDD slope")
ggplot(data = na.omit(betas_climate)) +
  geom_vline(aes(xintercept = 0), linetype = 'dotted') + 
  geom_hline(aes(yintercept = 0), linetype = 'dotted') + 
  geom_pointrange(aes(x=rangegdd_q50,
                      ymin=betasm_q5, y=betasm_q50, ymax=betasm_q95,
                      color = as.character(clade)), size = 0.2) +
  geom_pointrange(aes(xmin=rangegdd_q5, x=rangegdd_q50, xmax=rangegdd_q95,
                      y=betasm_q50,
                      color = as.character(clade)), size = 0.2) +
  scale_color_manual(values = c('#27278f', '#278f27')) +
  # coord_equal(xlim = c(-0.1,0.3), ylim = c(-0.1,0.3)) +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(x = "Range-wide GDD (Apr.-Sep., °C)", y = "Soil moisture slope")
ggplot(data = na.omit(betas_climate)) +
  geom_vline(aes(xintercept = 0), linetype = 'dotted') + 
  geom_hline(aes(yintercept = 0), linetype = 'dotted') + 
  geom_pointrange(aes(x=rangegdd_q50,
                      ymin=betavpd_q5, y=betavpd_q50, ymax=betavpd_q95,
                      color = as.character(clade)), size = 0.2) +
  geom_pointrange(aes(xmin=rangegdd_q5, x=rangegdd_q50, xmax=rangegdd_q95,
                      y=betavpd_q50,
                      color = as.character(clade)), size = 0.2) +
  scale_color_manual(values = c('#27278f', '#278f27')) +
  # coord_equal(xlim = c(-0.1,0.3), ylim = c(-0.1,0.3)) +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(x = "Range-wide GDD (Apr.-Sep., °C)", y = "VPD slope")

# Density plots
# betas_samples <- data.frame()
# for(i in 1:nrow(betas_climate)){
#   range_gdd <- species_gddrange[[betas_climate[i, 'species_name_simplified']]]
#   betas_samples_i <- 
#     data.frame(species_code = betas_climate[i, 'species_code'], clade = betas_climate[i, 'clade'],
#                median_gdd = as.numeric(quantile(range_gdd, 0.5)),
#                betagdd = as.numeric(base_samples[[paste0('beta_gdd[', which(data$uniq_species_ids == betas_climate[i, 'species_code']), ']')]]),
#                betasm = as.numeric(base_samples[[paste0('beta_sm[', which(data$uniq_species_ids == betas_climate[i, 'species_code']), ']')]]),
#                betavpd = as.numeric(base_samples[[paste0('beta_vpd[', which(data$uniq_species_ids == betas_climate[i, 'species_code']), ']')]]))
#   betas_samples <- rbind(betas_samples, betas_samples_i)
# }
# 
# ggplot(data = na.omit(betas_samples)) +
#   geom_vline(aes(xintercept = 0), linetype = 'dotted') + 
#   geom_hline(aes(yintercept = 0), linetype = 'dotted') + 
#   geom_point(aes(x = betasm, y = betavpd, col = as.character(clade)), alpha = 0.1, size = 0.3) +
#   geom_density_2d(aes(x = betasm, y = betavpd, group = species_code), col = 'white', alpha = 1, linewidth = 2) +
#   geom_density_2d(aes(x = betasm, y = betavpd, col = as.character(clade), group = species_code), alpha = 1, linewidth = 0.5) +
#   scale_color_manual(values = c('#27278f', '#278f27'))


