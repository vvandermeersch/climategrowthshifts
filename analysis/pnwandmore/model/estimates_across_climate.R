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

data <- readRDS(file.path(wd, 'output/model', 'data_07oct2025.rds'))
datasets <- readRDS(file.path(wd, 'output/model', 'datasets_07oct2025.rds'))
base_samples <- readRDS(file.path(wd, "output/model/samples_07oct2025_partialpooling_2clades_gdd_amjjas.rds"))

# Species
species <- unique(datasets[,c('species_code', 'species_name')])
species[species$species_code == 'PISA', 'species_name'] <- "Pinus sabiniana"
temp <- stringr::str_split_fixed(species$species_name, " ", 3)
species$species_name_simplified <- paste(temp[,1], temp[,2], sep=" ")
species <- species[species$species_code %in% data$uniq_species_ids,]
species$subclade <- NA
species[grepl('Pinus', species$species_name_simplified), 'subclade'] <- 'Pinus'
species[grepl('Larix|Pseudotsuga', species$species_name_simplified), 'subclade'] <- 'Laricoideae'
species[grepl('Picea', species$species_name_simplified), 'subclade'] <- 'Picea'
species[grepl('Abies|Tsuga', species$species_name_simplified), 'subclade'] <- 'Abietoideae'
species[grepl('Juniperus|Callitropsis|Calocedrus|Chamaecyparis|Thuja|Taxodium|Sequoiadendron', species$species_name_simplified), 'subclade'] <- 'Cupressaceae'
species[grepl('Quercus', species$species_name_simplified), 'subclade'] <- 'Quercus'
species[grepl('Populus', species$species_name_simplified), 'subclade'] <- 'Populus'
species[grepl('Platanus', species$species_name_simplified), 'subclade'] <- 'Platanus'


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
species_gddrange <- readRDS(file.path(wd, 'output', 'climate', 'species_gdd_amjjas_range_experienced.rds'))
betas_climate$rangegdd_q5 <- betas_climate$rangegdd_q50 <- betas_climate$rangegdd_q95 <- NA
for(i in 1:nrow(betas_climate)){
  range_gdd <- species_gddrange[[betas_climate[i, 'species_name_simplified']]]
  betas_climate[i, paste0('rangegdd_q', c(5,50,95))] <- as.numeric(quantile(range_gdd, c(0.05,0.5,0.95)))
}
gddslopes_gdd <- ggplot(data = na.omit(betas_climate)) +
  facet_grid(~subclade, scales = 'free_x') +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') + 
  # geom_density_2d(aes(x = rangegdd_q50, y = betagdd_q50, col = as.character(clade), group = clade), alpha = 0.5) + 
  geom_pointrange(aes(x=rangegdd_q50,
                      ymin=betagdd_q5, y=betagdd_q50, ymax=betagdd_q95,
                      color = subclade), size = 0.2) +
  geom_pointrange(aes(xmin=rangegdd_q5, x=rangegdd_q50, xmax=rangegdd_q95,
                      y=betagdd_q50,
                      color = subclade), size = 0.2) +
  scale_color_manual(
    breaks = c('Pinus', 'Picea', 'Laricoideae', 'Abietoideae', 'Cupressaceae',
               'Quercus', 'Populus', 'Platanus'), 
    values = c('#2A4B73', '#6886A6', '#A0C0D9', '#B97AB0', '#8b8b5f',
               '#2B7A2B', '#66A766', '#A0CFA0')
  ) +
  # coord_equal(xlim = c(-0.1,0.3), ylim = c(-0.1,0.3)) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        panel.grid = element_blank()) +
  labs(x = "Range-wide GDD (Apr.-Sep., °C)", y = "GDD slope")

smslopes_gdd <- ggplot(data = na.omit(betas_climate)) +
  facet_grid(~subclade, scales = 'free_x') +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') + 
  # geom_density_2d(aes(x = rangegdd_q50, y = betagdd_q50, col = as.character(clade), group = clade), alpha = 0.5) + 
  geom_pointrange(aes(x=rangegdd_q50,
                      ymin=betasm_q5, y=betasm_q50, ymax=betasm_q95,
                      color = subclade), size = 0.2) +
  geom_pointrange(aes(xmin=rangegdd_q5, x=rangegdd_q50, xmax=rangegdd_q95,
                      y=betasm_q50,
                      color = subclade), size = 0.2) +
  scale_color_manual(
    breaks = c('Pinus', 'Picea', 'Laricoideae', 'Abietoideae', 'Cupressaceae',
               'Quercus', 'Populus', 'Platanus'), 
    values = c('#2A4B73', '#6886A6', '#A0C0D9', '#B97AB0', '#8b8b5f',
               '#2B7A2B', '#66A766', '#A0CFA0')
  ) +
  # coord_equal(xlim = c(-0.1,0.3), ylim = c(-0.1,0.3)) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        panel.grid = element_blank()) +
  labs(x = "Range-wide GDD (Apr.-Sep., °C)", y = "SM slope")


vpdslopes_gdd <- ggplot(data = na.omit(betas_climate)) +
  facet_grid(~subclade, scales = 'free_x') +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') + 
  # geom_density_2d(aes(x = rangegdd_q50, y = betagdd_q50, col = as.character(clade), group = clade), alpha = 0.5) + 
  geom_pointrange(aes(x=rangegdd_q50,
                      ymin=betavpd_q5, y=betavpd_q50, ymax=betavpd_q95,
                      color = subclade), size = 0.2) +
  geom_pointrange(aes(xmin=rangegdd_q5, x=rangegdd_q50, xmax=rangegdd_q95,
                      y=betavpd_q50,
                      color = subclade), size = 0.2) +
  scale_color_manual(
    breaks = c('Pinus', 'Picea', 'Laricoideae', 'Abietoideae', 'Cupressaceae',
               'Quercus', 'Populus', 'Platanus'), 
    values = c('#2A4B73', '#6886A6', '#A0C0D9', '#B97AB0', '#8b8b5f',
               '#2B7A2B', '#66A766', '#A0CFA0')
  ) +
  # coord_equal(xlim = c(-0.1,0.3), ylim = c(-0.1,0.3)) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        panel.grid = element_blank()) +
  labs(x = "Range-wide GDD (Apr.-Sep., °C)", y = "VPD slope")

comb <- gddslopes_gdd + smslopes_gdd + vpdslopes_gdd + plot_layout(ncol = 1)
ggsave(comb, file = file.path(wd, 'figures', 'newyork2025', 'model_estimates', 'estimates_againtGDD_subclades.pdf'),
       height = 200, width = 400, units = 'mm')

# Range-wide SM
species_smrange <- readRDS(file.path(wd, 'output', 'climate', 'species_sm_mjj_range_experienced.rds'))
betas_climate$rangesm_q5 <- betas_climate$rangesm_q50 <- betas_climate$rangesm_q95 <- NA
for(i in 1:nrow(betas_climate)){
  range_sm <- species_smrange[[betas_climate[i, 'species_name_simplified']]]*100
  betas_climate[i, paste0('rangesm_q', c(5,50,95))] <- as.numeric(quantile(range_sm, c(0.05,0.5,0.95)))
}
gddslopes_sm <- ggplot(data = na.omit(betas_climate)) +
  facet_grid(~subclade, scales = 'free_x') +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') + 
  # geom_density_2d(aes(x = rangegdd_q50, y = betagdd_q50, col = as.character(clade), group = clade), alpha = 0.5) + 
  geom_pointrange(aes(x=rangesm_q50,
                      ymin=betagdd_q5, y=betagdd_q50, ymax=betagdd_q95,
                      color = subclade), size = 0.2) +
  geom_pointrange(aes(xmin=rangesm_q5, x=rangesm_q50, xmax=rangesm_q95,
                      y=betagdd_q50,
                      color = subclade), size = 0.2) +
  scale_color_manual(
    breaks = c('Pinus', 'Picea', 'Laricoideae', 'Abietoideae', 'Cupressaceae',
               'Quercus', 'Populus', 'Platanus'), 
    values = c('#2A4B73', '#6886A6', '#A0C0D9', '#B97AB0', '#8b8b5f',
               '#2B7A2B', '#66A766', '#A0CFA0')
  ) +
  # coord_equal(xlim = c(-0.1,0.3), ylim = c(-0.1,0.3)) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        panel.grid = element_blank()) +
  labs(x = "Range-wide SM (May-June, %)", y = "GDD slope")

smslopes_sm <- ggplot(data = na.omit(betas_climate)) +
  facet_grid(~subclade, scales = 'free_x') +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') + 
  # geom_density_2d(aes(x = rangegdd_q50, y = betagdd_q50, col = as.character(clade), group = clade), alpha = 0.5) + 
  geom_pointrange(aes(x=rangesm_q50,
                      ymin=betasm_q5, y=betasm_q50, ymax=betasm_q95,
                      color = subclade), size = 0.2) +
  geom_pointrange(aes(xmin=rangesm_q5, x=rangesm_q50, xmax=rangesm_q95,
                      y=betasm_q50,
                      color = subclade), size = 0.2) +
  scale_color_manual(
    breaks = c('Pinus', 'Picea', 'Laricoideae', 'Abietoideae', 'Cupressaceae',
               'Quercus', 'Populus', 'Platanus'), 
    values = c('#2A4B73', '#6886A6', '#A0C0D9', '#B97AB0', '#8b8b5f',
               '#2B7A2B', '#66A766', '#A0CFA0')
  ) +
  # coord_equal(xlim = c(-0.1,0.3), ylim = c(-0.1,0.3)) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        panel.grid = element_blank()) +
  labs(x = "Range-wide SM (May-June, %)", y = "SM slope")


vpdslopes_sm <- ggplot(data = na.omit(betas_climate)) +
  facet_grid(~subclade, scales = 'free_x') +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') + 
  # geom_density_2d(aes(x = rangegdd_q50, y = betagdd_q50, col = as.character(clade), group = clade), alpha = 0.5) + 
  geom_pointrange(aes(x=rangesm_q50,
                      ymin=betavpd_q5, y=betavpd_q50, ymax=betavpd_q95,
                      color = subclade), size = 0.2) +
  geom_pointrange(aes(xmin=rangesm_q5, x=rangesm_q50, xmax=rangesm_q95,
                      y=betavpd_q50,
                      color = subclade), size = 0.2) +
  scale_color_manual(
    breaks = c('Pinus', 'Picea', 'Laricoideae', 'Abietoideae', 'Cupressaceae',
               'Quercus', 'Populus', 'Platanus'), 
    values = c('#2A4B73', '#6886A6', '#A0C0D9', '#B97AB0', '#8b8b5f',
               '#2B7A2B', '#66A766', '#A0CFA0')
  ) +
  # coord_equal(xlim = c(-0.1,0.3), ylim = c(-0.1,0.3)) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        panel.grid = element_blank()) +
  labs(x = "Range-wide SM (May-June, %)", y = "VPD slope")

comb <- gddslopes_sm + smslopes_sm + vpdslopes_sm + plot_layout(ncol = 1)
ggsave(comb, file = file.path(wd, 'figures', 'newyork2025', 'model_estimates', 'estimates_againtSM_subclades.pdf'),
       height = 200, width = 400, units = 'mm')


# Range-wide altitude
species_altrange <- readRDS(file.path(wd, 'output', 'climate', 'species_altitude_range_experienced.rds'))
betas_climate$rangealt_q5 <- betas_climate$rangealt_q50 <- betas_climate$rangealt_q95 <- NA
for(i in 1:nrow(betas_climate)){
  range_alt <- species_altrange[[betas_climate[i, 'species_name_simplified']]]
  betas_climate[i, paste0('rangealt_q', c(5,50,95))] <- as.numeric(quantile(range_alt, c(0.05,0.5,0.95)))
}
gddslopes_alt <- ggplot(data = na.omit(betas_climate)) +
  facet_grid(~subclade) +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') + 
  # geom_density_2d(aes(x = rangegdd_q50, y = betagdd_q50, col = as.character(clade), group = clade), alpha = 0.5) + 
  geom_pointrange(aes(x=rangealt_q50,
                      ymin=betagdd_q5, y=betagdd_q50, ymax=betagdd_q95,
                      color = subclade), size = 0.2) +
  geom_pointrange(aes(xmin=rangealt_q5, x=rangealt_q50, xmax=rangealt_q95,
                      y=betagdd_q50,
                      color = subclade), size = 0.2) +
  scale_color_manual(
    breaks = c('Pinus', 'Picea', 'Laricoideae', 'Abietoideae', 'Cupressaceae',
               'Quercus', 'Populus', 'Platanus'), 
    values = c('#2A4B73', '#6886A6', '#A0C0D9', '#B97AB0', '#8b8b5f',
               '#2B7A2B', '#66A766', '#A0CFA0')
  ) +
  # coord_equal(xlim = c(-0.1,0.3), ylim = c(-0.1,0.3)) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        panel.grid = element_blank()) +
  labs(x = "Range-wide altitude (m)", y = "GDD slope")

smslopes_alt <- ggplot(data = na.omit(betas_climate)) +
  facet_grid(~subclade) +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') + 
  # geom_density_2d(aes(x = rangegdd_q50, y = betagdd_q50, col = as.character(clade), group = clade), alpha = 0.5) + 
  geom_pointrange(aes(x=rangealt_q50,
                      ymin=betasm_q5, y=betasm_q50, ymax=betasm_q95,
                      color = subclade), size = 0.2) +
  geom_pointrange(aes(xmin=rangealt_q5, x=rangealt_q50, xmax=rangealt_q95,
                      y=betasm_q50,
                      color = subclade), size = 0.2) +
  scale_color_manual(
    breaks = c('Pinus', 'Picea', 'Laricoideae', 'Abietoideae', 'Cupressaceae',
               'Quercus', 'Populus', 'Platanus'), 
    values = c('#2A4B73', '#6886A6', '#A0C0D9', '#B97AB0', '#8b8b5f',
               '#2B7A2B', '#66A766', '#A0CFA0')
  ) +
  # coord_equal(xlim = c(-0.1,0.3), ylim = c(-0.1,0.3)) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        panel.grid = element_blank()) +
  labs(x = "Range-wide altitude (m)", y = "SM slope")


vpdslopes_alt <- ggplot(data = na.omit(betas_climate)) +
  facet_grid(~subclade) +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') + 
  # geom_density_2d(aes(x = rangegdd_q50, y = betagdd_q50, col = as.character(clade), group = clade), alpha = 0.5) + 
  geom_pointrange(aes(x=rangealt_q50,
                      ymin=betavpd_q5, y=betavpd_q50, ymax=betavpd_q95,
                      color = subclade), size = 0.2) +
  geom_pointrange(aes(xmin=rangealt_q5, x=rangealt_q50, xmax=rangealt_q95,
                      y=betavpd_q50,
                      color = subclade), size = 0.2) +
  scale_color_manual(
    breaks = c('Pinus', 'Picea', 'Laricoideae', 'Abietoideae', 'Cupressaceae',
               'Quercus', 'Populus', 'Platanus'), 
    values = c('#2A4B73', '#6886A6', '#A0C0D9', '#B97AB0', '#8b8b5f',
               '#2B7A2B', '#66A766', '#A0CFA0')
  ) +
  # coord_equal(xlim = c(-0.1,0.3), ylim = c(-0.1,0.3)) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        panel.grid = element_blank()) +
  labs(x = "Range-wide altitude (m)", y = "VPD slope")

comb <- gddslopes_alt + smslopes_alt + vpdslopes_alt + plot_layout(ncol = 1)
ggsave(comb, file = file.path(wd, 'figures', 'newyork2025', 'model_estimates', 'estimates_againtaltitude_subclades.pdf'),
       height = 200, width = 400, units = 'mm')
