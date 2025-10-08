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
gp_params <- data.frame(
  species_code = data$uniq_species_ids,
  clade = data$clade_idxs,
  rhosp_q5 = sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('rho_sp[', sp, ']')]], 0.05)),
  rhosp_q50 = sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('rho_sp[', sp, ']')]], 0.5)),
  rhosp_q95 = sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('rho_sp[', sp, ']')]], 0.95)),
  gammasp_q5 = sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('gamma_sp[', sp, ']')]], 0.05)),
  gammasp_q50 = sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('gamma_sp[', sp, ']')]], 0.5)),
  gammasp_q95 = sapply(1:data$N_species, function(sp)  util$ensemble_mcmc_quantile_est(base_samples[[paste0('gamma_sp[', sp, ']')]], 0.95))
)
gp_params <- merge(gp_params, species)

gp_params_samples <- data.frame()
for(i in 1:nrow(species)){
  gp_params_i <-
    data.frame(species_name_simplified = species[i, 'species_name_simplified'], 
               rho = as.numeric(base_samples[[paste0('rho_sp[', which(data$uniq_species_ids == species[i, 'species_code']), ']')]]),
               gamma = as.numeric(base_samples[[paste0('gamma_sp[', which(data$uniq_species_ids == species[i, 'species_code']), ']')]]))
  gp_params_samples <- rbind(gp_params_samples, gp_params_i)
}

# GP covariate function plots!
exp_quad <- function(dx, gamma, rho){
  return(gamma^2 * exp(-1/2*(dx/rho)^2))
}
gp_plot <- data.frame()
sp <- 'Pseudotsuga menziesii'
gp_params_sp <- gp_params_samples[gp_params_samples$species_name_simplified == sp,]
for(s in sample(1:nrow(gp_params_sp),100)){
  rho <- gp_params_sp[s, 'rho']
  gamma <- gp_params_sp[s, 'gamma']
  gp_plot_s <-
    data.frame(
      sp, rho, gamma, f = sapply(seq(1,40,0.1), function(i) exp_quad(i, gamma = gamma, rho = rho)),
      year = seq(1,40,0.1)
    )
  gp_plot <- rbind(gp_plot, gp_plot_s)
}
sp <- 'Pinus edulis'
gp_params_sp <- gp_params_samples[gp_params_samples$species_name_simplified == sp,]
for(s in sample(1:nrow(gp_params_sp),100)){
  rho <- gp_params_sp[s, 'rho']
  gamma <- gp_params_sp[s, 'gamma']
  gp_plot_s <-
    data.frame(
      sp, rho, gamma, f = sapply(seq(1,40,0.1), function(i) exp_quad(i, gamma = gamma, rho = rho)),
      year = seq(1,40,0.1)
    )
  gp_plot <- rbind(gp_plot, gp_plot_s)
}
sp <- 'Pinus ponderosa'
gp_params_sp <- gp_params_samples[gp_params_samples$species_name_simplified == sp,]
for(s in sample(1:nrow(gp_params_sp),100)){
  rho <- gp_params_sp[s, 'rho']
  gamma <- gp_params_sp[s, 'gamma']
  gp_plot_s <-
    data.frame(
      sp, rho, gamma, f = sapply(seq(1,40,0.1), function(i) exp_quad(i, gamma = gamma, rho = rho)),
      year = seq(1,40,0.1)
    )
  gp_plot <- rbind(gp_plot, gp_plot_s)
}

example_cov <- ggplot(data = gp_plot) +
  geom_line(aes(x = year, y = f, group = paste0(rho, gamma), color = sp), alpha = 0.1) +
  scale_color_manual(values = c('#2A4B73', "#2B7A2B", "#6A0066"), name = '',
                     breaks = c('Pseudotsuga menziesii', 'Pinus edulis', 'Pinus ponderosa'),
                     labels = sapply(c('Pseudotsuga menziesii', 'Pinus edulis', 'Pinus ponderosa'), function(s){
                       as.expression(bquote(.(s)~'('*rho == .(round(gp_params[gp_params$species_name_simplified ==s, 'rhosp_q50'],2))*','~
                                              gamma == .(round(gp_params[gp_params$species_name_simplified ==s, 'gammasp_q50'],2))*')'))})) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  coord_cartesian(expand = FALSE) +
  theme_classic() +
  theme(legend.position = "inside", legend.position.inside = c(0.7,0.7)) +
  labs(y = 'Covariance function', x = 'Lag (year)')
ggsave(example_cov, , file = file.path(wd, 'figures', 'newyork2025', 'model_showcase', 'covariance_functions.pdf'), 
       width = 200, height = 150, units = "mm")




gp_plot <- data.frame()
for(sp in species[grepl('Pinus', species$species_name_simplified), 'species_name_simplified']){
  gp_params_sp <- gp_params_samples[gp_params_samples$species_name_simplified == sp,]
  for(s in sample(1:nrow(gp_params_sp),10)){
    rho <- gp_params_sp[s, 'rho']
    gamma <- gp_params_sp[s, 'gamma']
    gp_plot_s <-
      data.frame(
        sp, rho, gamma, f = sapply(seq(1,80,0.1), function(i) exp_quad(i, gamma = gamma, rho = rho)),
        year = seq(1,80,0.1)
      )
    gp_plot <- rbind(gp_plot, gp_plot_s)
  }
}

pinus_cov <- ggplot(data = gp_plot) +
  geom_line(aes(x = year, y = f, group = paste0(rho, gamma), color = sp), alpha = 0.3) +
  guides(color = guide_legend(override.aes = list(alpha = 1), ncol = 2)) +
  coord_cartesian(expand = FALSE) +
  theme_classic() +
  labs(y = 'Covariance function', x = 'Lag (year)')
ggsave(pinus_cov, , file = file.path(wd, 'figures', 'newyork2025', 'model_showcase', 'pinus_covariance_functions.pdf'), 
       width = 300, height = 150, units = "mm")
