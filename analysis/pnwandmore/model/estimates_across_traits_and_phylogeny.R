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
  betavpd_q95= sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('beta_vpd[', sp, ']')]], 0.95)),
  rhosp_q5 = sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('rho_sp[', sp, ']')]], 0.05)),
  rhosp_q50 = sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('rho_sp[', sp, ']')]], 0.5)),
  rhosp_q95 = sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('rho_sp[', sp, ']')]], 0.95)),
  gammasp_q5 = sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('gamma_sp[', sp, ']')]], 0.05)),
  gammasp_q50 = sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('gamma_sp[', sp, ']')]], 0.5)),
  gammasp_q95 = sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('gamma_sp[', sp, ']')]], 0.95)),
  kappa_q5 = sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('kappa_sh[', sp, ']')]], 0.05)),
  kappa_q50 = sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('kappa_sh[', sp, ']')]], 0.5)),
  kappa_q95 = sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('kappa_sh[', sp, ']')]], 0.95))
)
betas_climate <- merge(species, betas_climate)

betas_samples <- data.frame()
for(i in 1:nrow(betas_climate)){
  betas_samples_i <-
    data.frame(species_name_simplified = betas_climate[i, 'species_name_simplified'], clade = betas_climate[i, 'clade'],
               betagdd = as.numeric(base_samples[[paste0('beta_gdd[', which(data$uniq_species_ids == betas_climate[i, 'species_code']), ']')]]),
               betasm = as.numeric(base_samples[[paste0('beta_sm[', which(data$uniq_species_ids == betas_climate[i, 'species_code']), ']')]]),
               betavpd = as.numeric(base_samples[[paste0('beta_vpd[', which(data$uniq_species_ids == betas_climate[i, 'species_code']), ']')]]),
               rho = as.numeric(base_samples[[paste0('rho_sp[', which(data$uniq_species_ids == betas_climate[i, 'species_code']), ']')]]))
  betas_samples <- rbind(betas_samples, betas_samples_i)
}


# Traits!
mao_traits <- read.csv(file.path("~/projects/mast_trait/Data/silvicsClean.csv")) # from Mao
mao_traits <- mao_traits[mao_traits$speciesName != "", c("genusName", "speciesName", "droughtTolerance", "lifeSpanMax")]
mao_traits$species_name_simplified <- paste(mao_traits$genusName, mao_traits$speciesName, sep = ' ')
betas_climate <- merge(betas_climate, mao_traits[,c("species_name_simplified", "droughtTolerance", "lifeSpanMax")], all.x = TRUE)
betas_climate$droughtTolerance <- ifelse(betas_climate$droughtTolerance == '' | is.na(betas_climate$droughtTolerance), 'Unknown', betas_climate$droughtTolerance)
betas_climate$drought_tolerance <- factor(betas_climate$droughtTolerance , levels = c("Low", "Moderate", "High", "Unknown"))
betas_climate$lifeSpanMax <- as.numeric(ifelse(betas_climate$lifeSpanMax == '' | is.na(betas_climate$lifeSpanMax), NA, betas_climate$lifeSpanMax))
betas_samples <- merge(betas_samples, mao_traits[,c("species_name_simplified", "droughtTolerance", "lifeSpanMax")], all.x = TRUE)
betas_samples$droughtTolerance <- ifelse(betas_samples$droughtTolerance == '' | is.na(betas_samples$droughtTolerance), 'Unknown', betas_samples$droughtTolerance)
betas_samples$drought_tolerance <- factor(betas_samples$droughtTolerance , levels = c("Low", "Moderate", "High", "Unknown"))
betas_samples$lifeSpanMax <- as.numeric(ifelse(betas_samples$lifeSpanMax == '' | is.na(betas_samples$lifeSpanMax), NA, betas_samples$lifeSpanMax))
droughttrait_plot <- ggplot(data = betas_samples) +
  geom_boxplot(aes(x=drought_tolerance, y = betavpd, group = paste0(drought_tolerance, clade, species_name_simplified), color = as.character(clade)), 
               position = position_dodge(preserve = "single"), outliers = FALSE) +
  scale_color_manual(values = c('#27278f', '#278f27')) +
  labs(y = 'VDP slope', x = 'Drought tolerance') +
  theme_classic() +
  theme(legend.position = 'none')
ggsave(droughttrait_plot, file = file.path(wd, 'figures', 'newyork2025', 'model_estimates', 'vpd_estimates_droughttrait.pdf'), 
       width = 300, height = 200, units = "mm")
lifespantrait_plot <- ggplot(data = betas_samples) +
  geom_boxplot(aes(x=lifeSpanMax, y = rho, group = paste0(clade, species_name_simplified), color = as.character(clade)), 
               position = position_dodge(preserve = "single"), outliers = FALSE, width = 100) +
  scale_color_manual(values = c('#27278f', '#278f27')) +
  labs(y = 'Rho', x = 'Max. life span') +
  theme_classic() +
  theme(legend.position = 'none')
ggsave(lifespantrait_plot, file = file.path(wd, 'figures', 'newyork2025', 'model_estimates', 'rho_estimates_lifespantrait.pdf'), 
       width = 300, height = 200, units = "mm")

augusto_traits <- read.csv(file.path(wd, 'data', 'traits', 'augusto2025', 'Dataset-3_species.csv')) 
rho_samples <- data.frame()
for(i in 1:nrow(betas_climate)){
  rho_samples_i <-
    data.frame(species_name_simplified = betas_climate[i, 'species_name_simplified'], clade = betas_climate[i, 'clade'],
               rho = as.numeric(base_samples[[paste0('rho_sp[', which(data$uniq_species_ids == betas_climate[i, 'species_code']), ']')]]))
  rho_samples <- rbind(rho_samples, rho_samples_i)
}

augusto_traits <- merge(augusto_traits, rho_samples, by.x = 'Plant_Name', by.y = 'species_name_simplified') 
ggplot(data = augusto_traits) +
  geom_boxplot(aes(x=Amass, y = rho, group = Plant_Name), 
               position = position_dodge(preserve = "single"), outliers = FALSE) +
  scale_color_manual(values = c('#27278f', '#278f27')) +
  labs(y = 'Rho', x = 'Height max') +
  theme_classic() +
  theme(legend.position = 'none')



# Phylogeny!
library(ggtree)
source(file.path(wd, 'getphylo.R'))
names(sppfull) <- c('species_name', 'species_code', 'phylo_name')
sppfull <- merge(sppfull, species)
phy.plants.here$tip.label <- sppfull[match(phy.plants.here$tip.label, sppfull$phylo_name),'species_name_simplified']

betas_climate <- merge(betas_climate, sppfull[,c('species_code', 'phylo_name')])
names(betas_climate)[names(betas_climate)=='species_name_simplified'] <- 'label'
tree <- ggtree::ggtree(phy.plants.here) 
datatree <- tree$data
datatree <- merge(datatree, betas_climate, all.x = TRUE)
kippenberger <- c("#8B174DFF", "#AE2565FF", "#C1447EFF", "#D06C9BFF", "#DA9FB8FF",
                  "#ADBE7CFF", "#8BA749FF", "#6E8537FF", "#4F5F28FF", "#343D1FFF")

tree <- ggtree(datatree) +
  geom_tiplab(aes(color = betagdd_q50), size = 3) +
  scale_color_gradientn(colors = kippenberger, limits = c(-0.22,0.22)) +
  theme(legend.position = 'none') +
  coord_cartesian(xlim = c(0, 420), ylim = c(0, nrow(betas_climate)+1), expand = FALSE)

gdd <- ggplot(data = datatree) +
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey80') +
    geom_pointrange(aes(x = betagdd_q50, xmin = betagdd_q5, xmax= betagdd_q95,
                        color = betagdd_q50, y = y),
                    size = 0.2) +
    scale_color_gradientn(colors = kippenberger, limits = c(-0.22,0.22)) +
    guides(
      color = guide_colorbar(order = 0,
                             frame.colour = "grey20", ticks.colour = NA,
                             frame.linewidth = 0.2)) +
    theme_classic() +
    theme(axis.line.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), axis.title.y = element_blank(),
          axis.text.x = element_text(size = 8),
          plot.margin = margin(t = 5.5, b = 5.5, r = 12, l = 10),
          legend.position = "inside", legend.text = element_text(size = 8),
          legend.position.inside = c(0.9,0.25), legend.title = element_blank(),
          legend.direction="vertical", legend.key.width  = unit(8, "pt")) +
    coord_cartesian(xlim = c(-0.05,0.2), ylim = c(0, nrow(betas_climate)+1), expand = FALSE) +
    labs(x = 'GDD slope')

gdd_phylogeny <- tree + gdd + plot_layout(ncol = 2, width = c(0.65,0.4))
ggsave(gdd_phylogeny, file = file.path(wd, 'figures', 'newyork2025', 'model_estimates', 'gdd_estimates_phylogeny.pdf'), 
       width = 300, height = 200, units = "mm")

tree <- ggtree(datatree) +
  geom_tiplab(aes(color = betavpd_q50), size = 3) +
  scale_color_gradientn(colors = kippenberger, limits = c(-0.22,0.22)) +
  theme(legend.position = 'none') +
  coord_cartesian(xlim = c(0, 420), ylim = c(0, nrow(betas_climate)+1), expand = FALSE)

vpd <- ggplot(data = datatree) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey80') +
  geom_pointrange(aes(x = betavpd_q50, xmin = betavpd_q5, xmax= betavpd_q95,
                      color = betavpd_q50, y = y),
                  size = 0.2) +
  scale_color_gradientn(colors = kippenberger, limits = c(-0.21,0.21)) +
  guides(
    color = guide_colorbar(order = 0,
                           frame.colour = "grey20", ticks.colour = NA,
                           frame.linewidth = 0.2)) +
  theme_classic() +
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(size = 8),
        plot.margin = margin(t = 5.5, b = 5.5, r = 12, l = 10),
        legend.position = "inside", legend.text = element_text(size = 8),
        legend.position.inside = c(0.9,0.25), legend.title = element_blank(),
        legend.direction="vertical", legend.key.width  = unit(8, "pt")) +
  coord_cartesian(xlim = c(-0.25,0.1), ylim = c(0, nrow(betas_climate)+1), expand = FALSE) +
  labs(x = 'VPD slope')

vpd_phylogeny <- tree + vpd + plot_layout(ncol = 2, width = c(0.65,0.4))
ggsave(vpd_phylogeny, file = file.path(wd, 'figures', 'newyork2025', 'model_estimates', 'vpd_estimates_phylogeny.pdf'), 
       width = 300, height = 200, units = "mm")

tree <- ggtree(datatree) +
  geom_tiplab(aes(color = betasm_q50), size = 3) +
  scale_color_gradientn(colors = kippenberger, limits = c(-0.22,0.22)) +
  theme(legend.position = 'none') +
  coord_cartesian(xlim = c(0, 420), ylim = c(0, nrow(betas_climate)+1), expand = FALSE)

sm <- ggplot(data = datatree) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey80') +
  geom_pointrange(aes(x = betasm_q50, xmin = betasm_q5, xmax= betasm_q95,
                      color = betasm_q50, y = y),
                  size = 0.2) +
  scale_color_gradientn(colors = kippenberger, limits = c(-0.22,0.22)) +
  guides(
    color = guide_colorbar(order = 0,
                           frame.colour = "grey20", ticks.colour = NA,
                           frame.linewidth = 0.2)) +
  theme_classic() +
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(size = 8),
        plot.margin = margin(t = 5.5, b = 5.5, r = 12, l = 10),
        legend.position = "inside", legend.text = element_text(size = 8),
        legend.position.inside = c(0.9,0.25), legend.title = element_blank(),
        legend.direction="vertical", legend.key.width  = unit(8, "pt")) +
  coord_cartesian(xlim = c(-0.05,0.35), ylim = c(0, nrow(betas_climate)+1), expand = FALSE) +
  labs(x = 'SM slope')

sm_phylogeny <- tree + sm + plot_layout(ncol = 2, width = c(0.65,0.4))
ggsave(sm_phylogeny, file = file.path(wd, 'figures', 'newyork2025', 'model_estimates', 'sm_estimates_phylogeny.pdf'), 
       width = 300, height = 200, units = "mm")

tree <- ggtree(datatree) +
  geom_tiplab(aes(color = rhosp_q50), size = 3) +
  scale_color_gradientn(colors = c('#c4e6c3','#96d2a4','#6dbc90','#4da284','#36877a','#266b6e','#1d4f60'), limits = c(1,32)) +
  theme(legend.position = 'none') +
  coord_cartesian(xlim = c(0, 420), ylim = c(0, nrow(betas_climate)+1), expand = FALSE)

rho <- ggplot(data = datatree) +
  geom_pointrange(aes(x = rhosp_q50, xmin = rhosp_q5, xmax= rhosp_q95,
                      color = rhosp_q50, y = y),
                  size = 0.2) +
  scale_color_gradientn(colors = c('#c4e6c3','#96d2a4','#6dbc90','#4da284','#36877a','#266b6e','#1d4f60'), limits = c(1,32)) +
  guides(
    color = guide_colorbar(order = 0,
                           frame.colour = "grey20", ticks.colour = NA,
                           frame.linewidth = 0.2)) +
  theme_classic() +
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(size = 8),
        plot.margin = margin(t = 5.5, b = 5.5, r = 12, l = 10),
        legend.position = "inside", legend.text = element_text(size = 8),
        legend.position.inside = c(0.9,0.35), legend.title = element_blank(),
        legend.direction="vertical", legend.key.width  = unit(8, "pt")) +
  coord_cartesian(xlim = c(0,40), ylim = c(0, nrow(betas_climate)+1), expand = FALSE) +
  labs(x = 'Species-specific rhos')
rho_phylogeny <- tree + rho + plot_layout(ncol = 2, width = c(0.65,0.4))
ggsave(rho_phylogeny, file = file.path(wd, 'figures', 'newyork2025', 'model_estimates', 'rho_estimates_phylogeny.pdf'), 
       width = 300, height = 200, units = "mm")

tree <- ggtree(datatree) +
  geom_tiplab(aes(color = gammasp_q50), size = 3) +
  scale_color_gradientn(colors = c('#c4e6c3','#96d2a4','#6dbc90','#4da284','#36877a','#266b6e','#1d4f60'), limits = c(0.28,0.7)) +
  theme(legend.position = 'none') +
  coord_cartesian(xlim = c(0, 420), ylim = c(0, nrow(betas_climate)+1), expand = FALSE)

gamma <- ggplot(data = datatree) +
  geom_pointrange(aes(x = gammasp_q50, xmin = gammasp_q5, xmax= gammasp_q95,
                      color = gammasp_q50, y = y),
                  size = 0.2) +
  scale_color_gradientn(colors = c('#c4e6c3','#96d2a4','#6dbc90','#4da284','#36877a','#266b6e','#1d4f60'), limits = c(0.28,0.7)) +
  guides(
    color = guide_colorbar(order = 0,
                           frame.colour = "grey20", ticks.colour = NA,
                           frame.linewidth = 0.2)) +
  theme_classic() +
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(size = 8),
        plot.margin = margin(t = 5.5, b = 5.5, r = 12, l = 10),
        legend.position = "inside", legend.text = element_text(size = 8),
        legend.position.inside = c(0.9,0.35), legend.title = element_blank(),
        legend.direction="vertical", legend.key.width  = unit(8, "pt")) +
  coord_cartesian(xlim = c(0.25,0.75), ylim = c(0, nrow(betas_climate)+1), expand = FALSE) +
  labs(x = 'Species-specific gammas')
gamma_phylogeny <- tree + gamma + plot_layout(ncol = 2, width = c(0.65,0.4))
ggsave(gamma_phylogeny, file = file.path(wd, 'figures', 'newyork2025', 'model_estimates', 'gamma_estimates_phylogeny.pdf'), 
       width = 300, height = 200, units = "mm")

tree <- ggtree(datatree) +
  geom_tiplab(aes(color = rhosp_q50), size = 3) +
  scale_color_gradientn(colors = c('#c4e6c3','#96d2a4','#6dbc90','#4da284','#36877a','#266b6e','#1d4f60'), limits = c(1,32)) +
  theme(legend.position = 'none') +
  coord_cartesian(xlim = c(0, 420), ylim = c(0, nrow(betas_climate)+1), expand = FALSE)
rho_AND_gamma_phylogeny <- tree + rho + gamma + plot_layout(ncol = 3, width = c(0.65, 0.4, 0.4))
ggsave(rho_AND_gamma_phylogeny, file = file.path(wd, 'figures', 'newyork2025', 'model_estimates', 'rho_gamma_estimates_phylogeny.pdf'), 
       width = 400, height = 200, units = "mm")

tree <- ggtree(datatree) +
  geom_tiplab(aes(color = kappa_q50), size = 3) +
  scale_color_gradientn(colors = c('#c4e6c3','#96d2a4','#6dbc90','#4da284','#36877a','#266b6e','#1d4f60'), limits = c(0.1,1.05)) +
  theme(legend.position = 'none') +
  coord_cartesian(xlim = c(0, 420), ylim = c(0, nrow(betas_climate)+1), expand = FALSE)

kappa <- ggplot(data = datatree) +
  geom_pointrange(aes(x = kappa_q50, xmin = kappa_q5, xmax= kappa_q95,
                      color = kappa_q50, y = y),
                  size = 0.2) +
  scale_color_gradientn(colors = c('#c4e6c3','#96d2a4','#6dbc90','#4da284','#36877a','#266b6e','#1d4f60'), limits = c(0.1,1.05)) +
  guides(
    color = guide_colorbar(order = 0,
                           frame.colour = "grey20", ticks.colour = NA,
                           frame.linewidth = 0.2)) +
  theme_classic() +
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(size = 8),
        plot.margin = margin(t = 5.5, b = 5.5, r = 12, l = 10),
        legend.position = "inside", legend.text = element_text(size = 8),
        legend.position.inside = c(0.9,0.83), legend.title = element_blank(),
        legend.direction="vertical", legend.key.width  = unit(8, "pt")) +
  coord_cartesian(xlim = c(0,1.2), ylim = c(0, nrow(betas_climate)+1), expand = FALSE) +
  labs(x = 'Species-specific gammas')
kappa_phylogeny <- tree + kappa + plot_layout(ncol = 2, width = c(0.65,0.4))
ggsave(kappa_phylogeny, file = file.path(wd, 'figures', 'newyork2025', 'model_estimates', 'kappa_estimates_phylogeny.pdf'), 
       width = 300, height = 200, units = "mm")
