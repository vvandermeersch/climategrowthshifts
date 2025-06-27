rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnw"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)

library(rstan)
library(dplyr)

# Load treering data
datasets <- readRDS(file = file.path(wd, 'data/itrdb', paste0('itrdb_info.rds')))
datasets <- datasets[datasets$last_year >= 1999,] # at least 20 years of observations
datasets <- datasets[datasets$dataset != 'wa_149',] # temporary

# Keep only PIPO in California?
datasets <- datasets[datasets$state == 'california' & datasets$species_code == 'PIPO',]

raw_data <- data.frame()
for(d in 1:nrow(datasets)){
  
  raw_data_d <- read.csv(file.path(wd, 'data/itrdb', datasets[d, 'state'], 'total_ring_width/data', datasets[d, 'dataset'], 'cleaned_tree_data_with_flags.csv'))
  
  coltokeep <-  !grepl(colnames(raw_data_d), pattern = 'outlier')
  raw_data_d <- raw_data_d[, coltokeep]
  
  raw_data_d <- raw_data_d[raw_data_d$age_CE >= 1980 & !is.na(raw_data_d$age_CE), ]
  raw_data_d <- raw_data_d[ , colSums(is.na(raw_data_d)) < length(1980:2010)]
  
  raw_data_long <- raw_data_d %>%
    tidyr::pivot_longer(cols = -age_CE, names_to = 'core_id', values_to ='rw_core')
  names(raw_data_long)[1] <- 'year'
  raw_data_long$stand <- datasets[d, 'dataset']
  raw_data_long$species <- datasets[d, 'species_code']
  
  raw_data_long[raw_data_long$rw_core %in% c(-8, -9999), 'rw_core'] <- NA # dealing with weird values in or_105 dataset
  
  raw_data_long$tree_id <- stringr::str_split_i(raw_data_long$core_id, pattern = '_', i = 1)
  raw_data_long$tree_id <- paste0(raw_data_long$stand, "_", substr(raw_data_long$tree_id,1,nchar(raw_data_long$tree_id)-1))
  
  raw_data_d <- raw_data_long %>%
    group_by(stand, species, tree_id, year) %>%
    summarise(rw_core_ave = mean(rw_core, na.rm = TRUE), .groups = 'drop')
  
  trees_to_remove <- unique(raw_data_d[is.na(raw_data_d$rw_core_ave) | raw_data_d$rw_core_ave == 0 , 'tree_id'])
  raw_data_d <- raw_data_d[!(raw_data_d$tree_id %in% trees_to_remove$tree_id),]
  
  raw_data <- rbind(raw_data, raw_data_d)
}

# Map
library(rnaturalearth)
library(terra)
library(tidyterra)
library(ggplot2)
world <- vect(ne_countries(scale = "medium", returnclass = "sf"))
quick_map <- ggplot() +
  geom_spatvector(data = world, color = 'white', fill = 'grey90') +
  geom_point(data = datasets, aes(x = east_lon, y = north_lat, color = species_code)) +
  theme_minimal() +
  theme(axis.title = element_blank()) +
  coord_sf(xlim = c(-127, -110), ylim = c(35,52))

# Gather stands if they have the same latitude and longitude (rounded to 1e-2 degree, ie ~1.1km)
group_dataset <- datasets %>%
  group_by(round(north_lat,2), round(south_lat,2), round(east_lon,2), round(west_lon,2)) %>%
  mutate(group_dataset = paste0('S',cur_group_id())) %>%
  ungroup() %>%
  dplyr::select(dataset, group_dataset)
raw_data <- merge(raw_data, group_dataset, by.x = 'stand', by.y = 'dataset')

# Load climate data
clim_pred <- readRDS(file.path(wd, "output", "climate",  "climpredictors_ext_pnw_26june2025.rds"))
clim_pred$soilmoist_mjj <- clim_pred$soilmoist_mjj*100 # in percentage
clim_pred$gdd <- clim_pred$gdd/100 # in x100 degC
clim_pred$vpd_mjj <- clim_pred$vpd_mjj # in hPa?

# Sizes
uniq_tree_ids <- unique(raw_data$tree_id)
N_trees <- length(uniq_tree_ids)
uniq_stand_ids <- unique(raw_data$group_dataset)
N_stands <- length(unique(raw_data$group_dataset))
uniq_species_ids <- unique(raw_data$species) #  Important that JUOC is first here!
N_species <- length(uniq_species_ids)
all_years <-  min(raw_data$year):max(raw_data$year)
N_all_years <- length(all_years)

# Format data into ragged arrays
log_rw_obs <- c()
gdd_obs <- c()
sm_obs <- c()
vpd_obs <- c()
years <- c()
all_years_idxs <- c()
N_years <- c()
stand_idxs <- c()
species_idxs <- c()
idx <- 1
tree_start_idxs <- c()
tree_end_idxs <- c()
for (tid in uniq_tree_ids) {
  
  raw_data_tree <- raw_data[raw_data$tree_id == tid & raw_data$year %in% all_years,]
  
  years_tree <- raw_data_tree$year
  all_years_idxs_tree <- sapply(years_tree, function(y) which(all_years == y))
  N_years_tree <- length(years_tree)
  if(N_years_tree > 45){stop()}
  
  log_rw_obs_tree <- sapply(years_tree, 
                            function(y) 
                              log(raw_data_tree$rw_core_ave[raw_data_tree$year == y][1]))
  log_rw_obs <- c(log_rw_obs, log_rw_obs_tree)
  
  gdd_obs_tree <- sapply(years_tree, 
                         function(y) 
                           as.numeric(clim_pred$gdd[clim_pred$year == y & clim_pred$dataset == raw_data_tree$stand[1]][1]))
  gdd_obs <- c(gdd_obs, gdd_obs_tree)
  
  sm_obs_tree <- sapply(years_tree, 
                        function(y) 
                          as.numeric(clim_pred$soilmoist_mjj[clim_pred$year == y & clim_pred$dataset == raw_data_tree$stand[1]][1]))
  sm_obs <- c(sm_obs, sm_obs_tree)
  
  vpd_obs_tree <- sapply(years_tree, 
                         function(y) 
                           as.numeric(clim_pred$vpd_mjj[clim_pred$year == y & clim_pred$dataset == raw_data_tree$stand[1]][1]))
  vpd_obs <- c(vpd_obs, vpd_obs_tree)
  
  years <- c(years, years_tree)
  all_years_idxs <- c(all_years_idxs, all_years_idxs_tree)
  N_years <- c(N_years, N_years_tree)
  
  stand_tree <- which(uniq_stand_ids == raw_data_tree$group_dataset[1])
  # stand_tree <- which(uniq_stand_ids == raw_data_tree$region[1])
  stand_idxs <- c(stand_idxs, stand_tree)
  
  species_tree <- which(uniq_species_ids == raw_data_tree$species[1])
  species_idxs <- c(species_idxs, species_tree)
  
  tree_start_idxs <- c(tree_start_idxs, idx)
  idx <- idx + N_years_tree
  tree_end_idxs <- c(tree_end_idxs, idx - 1)
}

# Cross check sizes
N_trees
length(log_rw_obs)
length(gdd_obs)
length(sm_obs)
length(years)
length(all_years_idxs)
length(N_years)
length(tree_start_idxs)
length(tree_end_idxs)

# Phylogenetic matrix
# Cphy <- ape::vcv.phylo(phy.plants.here,corr=TRUE)

# Collection data into list
N <- length(years)
data <- mget(c('N', 'N_all_years', 'N_trees', 
               'log_rw_obs', 'gdd_obs', 'sm_obs', 'vpd_obs',
               'all_years', 'years', 'all_years_idxs', 'N_years', 
               'stand_idxs', 'N_stands',
               'species_idxs', 'N_species',
               'tree_start_idxs', 'tree_end_idxs'))
# saveRDS(data, file = file.path(wd, 'output/model', 'data_26june2025_californiaPIPO.rds'))


# Posterior quantificiation - diagnostics
fit <- readRDS(file.path(wd, 'output/model', 'fit_26june2025_californiaPIPO_winteract.rds')) # run on Margot
diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit)
base_samples <- util$filter_expectands(samples,
                                       c('alpha', 'beta_gdd', 'beta_sm', 
                                         'rho', 'gamma',
                                         'f_tilde_sh',
                                         'rho_sh', 'gamma_sh',
                                         'sigma'),
                                       check_arrays=TRUE)
util$check_all_expectand_diagnostics(base_samples)

# Retrodictive check
par(mfrow=c(3, 2), mar = c(4,4,1,1))
for (t in  sample(1:data$N_trees, 18)) {
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  rw_names <- sapply(idxs,
                     function(n) paste0('log_rw_pred[', n, ']'))
  util$plot_conn_pushforward_quantiles(samples, rw_names, data$years[idxs],
                                       xlab="Year", ylab="Log Ring Width Per mm", 
                                       display_ylim=c(-1.5, 3), display_xlim = c(1980, 2024))
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=1.0, col="white")
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=0.8, col="black")
  abline(v=2012, lwd=1, lty=2, col="#DDDDDD")
  abline(v=2016, lwd=1, lty=2, col="#DDDDDD")
}

par(mfrow=c(2, 2), mar = c(4,4,1,1))
indxs <-which(data$sm_obs < 24)
pred_names <- sapply(indxs, function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conditional_median_quantiles(samples, pred_names, data$vpd_obs[indxs],
                                       5, 21, 1, data$log_rw_obs[indxs], 
                                       residual=FALSE,
                                       ylab = 'Marginal quantiles',
                                       xlab="VPD (hPa)", main = 'SM < 24')
util$plot_conditional_median_quantiles(samples, pred_names, data$vpd_obs[indxs],
                                       5, 21, 1, data$log_rw_obs[indxs], 
                                       residual=TRUE,
                                       ylab = 'Marginal quantiles (resid.)',
                                       xlab="VPD (hPa)")
indxs <-which(data$sm_obs > 27)
pred_names <- sapply(indxs, function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conditional_median_quantiles(samples, pred_names, data$vpd_obs[indxs],
                                       5, 21, 1, data$log_rw_obs[indxs], 
                                       residual=FALSE,
                                       ylab = 'Marginal quantiles',
                                       xlab="VPD (hPa)", main = 'SM > 27')
util$plot_conditional_median_quantiles(samples, pred_names, data$vpd_obs[indxs],
                                       5, 21, 1, data$log_rw_obs[indxs], 
                                       residual=TRUE,
                                       ylab = 'Marginal quantiles (resid.)',
                                       xlab="VPD (hPa)")





# Inference
par(mfrow=c(1, 4))
util$plot_expectand_pushforward(samples[['beta_gdd']], 50,
                                flim = c(-0.2,0.2),
                                display_name="beta_gdd")
xs <- seq(-1, 1, 0.01)
ys <- dnorm(xs, 0, log(1.8) / 2.57)*50
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['beta_sm']], 50,
                                flim = c(-0.2,0.2),
                                display_name="beta_sm")
xs <- seq(-1, 1, 0.01)
ys <- dnorm(xs, 0, log(1.8) / 2.57) *10
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['beta_vpd']], 50,
                                flim = c(-0.2,0.2),
                                display_name="beta_vpd")
xs <- seq(-1, 1, 0.01)
ys <- dnorm(xs, 0, log(1.8) / 2.57) *10
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['beta_smvpd']], 50,
                                flim = c(-0.2,0.2),
                                display_name="beta_smvpd")
xs <- seq(-1, 1, 0.01)
ys <- dnorm(xs, 0, log(1.8) / 2.57) *10
lines(xs, ys, lwd=2, col=util$c_light)


par(mfrow=c(1, 2))
util$plot_expectand_pushforward(samples[['rho']], 25,
                                display_name="rho")

util$plot_expectand_pushforward(samples[['rho_sh']], 25,
                                display_name="rho_sh")
