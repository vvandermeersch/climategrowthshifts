
rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnw"
setwd(wd)
# util <- new.env()
# source('mcmc_analysis_tools_rstan.R', local=util)
# source('mcmc_visualization_tools.R', local=util)

library(rstan)
library(dplyr)

# Load treering data
datadir <- file.path(wd, 'data/itrdb/oregon/total_ring_width/data')
datasets <- readRDS(file = file.path(wd, 'data/itrdb', paste0('itrdb_info.rds')))
datasets <- datasets[datasets$last_year >= 2010,]

raw_data <- data.frame()
for(d in 1:nrow(datasets)){
  
  raw_data_d <- read.csv(file.path(datadir, datasets[d, 'dataset'], 'cleaned_tree_data_with_flags.csv'))
  
  coltokeep <-  !grepl(colnames(raw_data_d), pattern = 'outlier')
  raw_data_d <- raw_data_d[, coltokeep]
  
  raw_data_d <- raw_data_d[raw_data_d$age_CE >= 1980 & raw_data_d$age_CE <= 2010, ]
  raw_data_d <- raw_data_d[ , colSums(is.na(raw_data_d)) < length(1980:2010)]
  
  raw_data_long <- raw_data_d %>%
    tidyr::pivot_longer(cols = -age_CE, names_to = 'core_id', values_to ='rw_core')
  names(raw_data_long)[1] <- 'year'
  raw_data_long$stand <- datasets[d, 'dataset']
  raw_data_long$species <- datasets[d, 'species_code']
  
  raw_data_long[raw_data_long$rw_core %in% c(-8), 'rw_core'] <- NA # dealing with weird values in or_105 dataset
  
  raw_data_long$tree_id <- stringr::str_split_i(raw_data_long$core_id, pattern = '_', i = 1)
  raw_data_long$tree_id <- paste0(raw_data_long$stand, "_", substr(raw_data_long$tree_id,1,nchar(raw_data_long$tree_id)-1))
  
  raw_data_d <- raw_data_long %>%
    group_by(stand, species, tree_id, year) %>%
    summarise(rw_core_ave = mean(rw_core, na.rm = TRUE), .groups = 'drop')
  
  trees_to_remove <- unique(raw_data_d[is.na(raw_data_d$rw_core_ave) | raw_data_d$rw_core_ave == 0 , 'tree_id'])
  raw_data_d <- raw_data_d[!(raw_data_d$tree_id %in% trees_to_remove$tree_id),]
  
  raw_data <- rbind(raw_data, raw_data_d)
}

# Load climate data
clim_pred <- readRDS(file.path(wd, "output", "climate",  "climpredictors_pnw_june2025.rds"))
clim_pred$soilmoist_jja <- clim_pred$soilmoist_jja*100 # in percentage
clim_pred$gdd_ings <- clim_pred$gdd_ings/10 # in 10degC

# Sizes
uniq_tree_ids <- unique(raw_data$tree_id)
N_trees <- length(uniq_tree_ids)

uniq_stand_ids <- unique(raw_data$stand)
N_stands <- length(unique(raw_data$stand))

uniq_species_ids <- unique(raw_data$species) #  Important that Abam is first here!
N_species <- length(uniq_species_ids)

# Common years
all_years <- NULL
for (tid in uniq_tree_ids) {
  years <- raw_data[raw_data$tree_id == tid,]$year
  if (is.null(all_years))
    all_years <- years
  else
    all_years <- intersect(years, all_years)
}

N_all_years <- length(all_years)

# Format data into ragged arrays
log_rw_obs <- c()
gsl_obs <- c()
gdd_obs <- c()
sm_obs <- c()
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
  if(N_years_tree > 31){stop()}
  
  log_rw_obs_tree <- sapply(years_tree, 
                            function(y) 
                              log(raw_data_tree$rw_core_ave[raw_data_tree$year == y][1]))
  log_rw_obs <- c(log_rw_obs, log_rw_obs_tree)
  
  gdd_obs_tree <- sapply(years_tree, 
                         function(y) 
                           as.numeric(clim_pred$gdd_ings[clim_pred$year == y & clim_pred$dataset == raw_data_tree$stand[1]][1]))
  gdd_obs <- c(gdd_obs, gdd_obs_tree)
  
  sm_obs_tree <- sapply(years_tree, 
                        function(y) 
                          as.numeric(clim_pred$soilmoist_jja[clim_pred$year == y & clim_pred$dataset == raw_data_tree$stand[1]][1]))
  sm_obs <- c(sm_obs, sm_obs_tree)
  
  years <- c(years, years_tree)
  all_years_idxs <- c(all_years_idxs, all_years_idxs_tree)
  N_years <- c(N_years, N_years_tree)
  
  stand_tree <- which(uniq_stand_ids == raw_data_tree$stand[1])
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

# Collection data into list
N <- length(years)

data <- mget(c('N', 'N_all_years', 'N_trees', 
               'log_rw_obs', 'gdd_obs', 'sm_obs',
               'all_years', 'years', 'all_years_idxs', 'N_years', 
               'stand_idxs', 'N_stands',
               'tree_start_idxs', 'tree_end_idxs'))

# Posterior Quantification
fit <- stan(file=file.path(wd, 'model/stan/model6_with2predictors_pnw.stan'),
            data=data, seed=5838299, cores = 4,
            warmup=1000, iter=2024, refresh=10)


