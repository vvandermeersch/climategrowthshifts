rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('functions/custom_functions.R', local=util)
setwd(wd)


library(rnaturalearth)
library(terra)
library(tidyterra)
library(ggplot2)
world <- vect(ne_countries(scale = "medium", returnclass = "sf"))
westam <- vect(ne_states(country = c("Canada", "Mexico", "United States of America"), returnclass = "sf"))


# Model data
data <- readRDS(file.path(wd, 'output/model', 'data_24sept2025.rds'))

fit <- readRDS(file.path(wd, "output/model/fit_24sept2025_partialpooling_2clades_centered_extended_gdd_amjjas.rds"))
# base_samples <- readRDS(file.path(wd, "output/model/samples_24sept2025_partialpooling_2clades_centered_extended_gdd_amjjas.rds"))


par(mfrow = c(1,1), mar = c(4,4,1,1))
stand <- unique(data$stand_idxs)[1]
prednames <- sapply(1:data$N_all_years, function(n) paste0('f_sh[', stand,',', n, ']'))
samples_fsh <- rstan::extract(fit, pars = prednames, permuted=FALSE)
N <- length(prednames)
samples_fsh <- lapply(1:N, function(n) t(samples_fsh[,,n]))
names(samples_fsh) <- prednames
util$plot_realizations(samples_fsh, prednames, data$all_years, N_plots=50,
                       xlab="Year", display_ylim=c(-3, 3), display_xlim = c(1980, 2010))


fsh_pred <- data.frame()
fsh_trend <- data.frame()
for(stand in unique(data$stand_idxs)){
  
  trees_stand <- which(data$stand_idxs == stand)
  idx_trees_stand <- unlist(sapply(trees_stand, function(t) return(data$tree_start_idxs[t]:data$tree_end_idxs[t])))
  idx_years_trees_stand <- data$all_years_idxs[idx_trees_stand]
  idx_years_obs_stand <- min(idx_years_trees_stand):max(idx_years_trees_stand)
  years_trees_stand <- data$years[idx_trees_stand]
  years_obs_stand <- min(years_trees_stand):max(years_trees_stand)
  
  
  prednames <- sapply(1:data$N_all_years, function(n) paste0('f_sh[', stand,',', n, ']'))
  samples_fsh <- rstan::extract(fit, pars = prednames, permuted=FALSE)
  N <- length(prednames)
  samples_fsh <- lapply(1:N, function(n) t(samples_fsh[,,n]))
  names(samples_fsh) <- prednames
  
  prednames <- sapply(idx_years_obs_stand, function(n) paste0('f_sh[', stand,',', n, ']'))
  util$plot_realizations(samples_fsh, prednames, years_obs_stand, N_plots=50,
                         xlab="Year", display_ylim=c(-3, 3), display_xlim = c(min(years_obs_stand), max(years_obs_stand)))
  abline(v=1990, lwd=1, lty=2, col="grey50")
  
  fsh_pred <- rbind(
    fsh_pred,
    data.frame(
      standnum = stand,
      stand = data$uniq_stand_ids[stand],
      fsh_q5 = sapply(prednames, function(s){util$ensemble_mcmc_quantile_est(samples_fsh[[s]],  c(0.05))}),
      fsh_q50 = sapply(prednames, function(s){util$ensemble_mcmc_quantile_est(samples_fsh[[s]],  c(0.5))}),
      fsh_q95 = sapply(prednames, function(s){util$ensemble_mcmc_quantile_est(samples_fsh[[s]],  c(0.95))})
    )
  )
  
  prednames_before1990 <- sapply(which(idx_years_obs_stand %in% which(data$all_years <= 1990)), function(n) paste0('f_sh[', stand,',', n, ']'))
  prednames_after1990 <- sapply(which(idx_years_obs_stand %in% which(data$all_years > 1990)), function(n) paste0('f_sh[', stand,',', n, ']'))
  fsh_trend <- rbind(
    fsh_trend,
    data.frame(
      standnum = stand,
      stand = data$uniq_stand_ids[stand],
      fsh_q50_before1990 = util$ensemble_mcmc_quantile_est(sapply(prednames_before1990, function(s) samples_fsh[[s]]), c(0.5)),
      fsh_q50_after1990 = util$ensemble_mcmc_quantile_est(sapply(prednames_after1990, function(s) samples_fsh[[s]]), c(0.5))
    )
  )
  
}



fsh_trend$perc_var <- (fsh_trend$fsh_q50_after1990 - fsh_trend$fsh_q50_before1990)/fsh_trend$fsh_q50_before1990
datasets <- readRDS(file.path(wd, 'output/model', 'datasets_24sept2025.rds'))
stand_coordinates <- unique(data.frame(stand = datasets$grouped_stand, lat = round(datasets$north_lat,2), lon = round(datasets$west_lon,2)))

fsh_trend_coordinates <- merge(fsh_trend, stand_coordinates)
fsh_trend_coordinates$perc_var <- ifelse(fsh_trend_coordinates$perc_var  > 2, 2, ifelse(fsh_trend_coordinates$perc_var < -2, -2, fsh_trend_coordinates$perc_var))
ggplot(data = fsh_trend_coordinates) +
  geom_spatvector(data = westam, color = 'white', fill = 'grey90') +
  geom_point(aes(x = lon, y = lat, color = perc_var)) +
  scale_color_viridis_c(breaks = seq(-2,2,1), labels = c('<-200%', '-100%', '0', '100%', '>200%'), name = '',
                        guide = guide_colourbar(direction = "horizontal", barwidth = 10, barheight = 0.5)) +
  ggplot2::theme_minimal() +
  theme(axis.title = ggplot2::element_blank(), legend.position = "bottom") +
  ggplot2::coord_sf(xlim = c(-124.925, -98.025), ylim = c(23.065,52.925), expand = FALSE) 
  

ggplot(data = fsh_trend_coordinates) +
  geom_spatvector(data = westam, color = 'white', fill = 'grey90') +
  stat_summary_hex(aes(x = lon, y = lat, z = perc_var), fun = function(x) mean(x), bins = 20) +
  scale_fill_viridis_c(breaks = seq(-2,2,1), labels = c('<-200%', '-100%', '0', '100%', '>200%'), name = '',
                        guide = guide_colourbar(direction = "horizontal", barwidth = 10, barheight = 0.5)) +
  ggplot2::theme_minimal() +
  theme(axis.title = ggplot2::element_blank(), legend.position = "bottom") +
  ggplot2::coord_sf(xlim = c(-124.925, -98.025), ylim = c(23.065,52.925), expand = FALSE) 







