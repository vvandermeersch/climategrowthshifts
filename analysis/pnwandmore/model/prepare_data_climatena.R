rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)

library(rstan)
library(dplyr)

# Load treering data
datasets <- readRDS(file.path(wd, 'input', 'itrdb', 'datasets_summary_all.rds'))
ringwidth_series <- readRDS(file.path(wd, 'input', 'itrdb', 'ringwidth_series_usonly_from1902.rds'))
datasets <- datasets[datasets$dataset %in% unique(ringwidth_series$dataset),]

# Temporary - missing stands due to land mask
datasets <- datasets[datasets$dataset != 'wa149',] # temporary 
datasets <- datasets[datasets$dataset != 'ca673',] # temporary
datasets <- datasets[datasets$dataset != 'az563',] # temporary

# Temporary - missing species 
datasets <- datasets[datasets$dataset != 'co689',] # temporary 
datasets <- datasets[datasets$dataset != 'co690',] # temporary
datasets <- datasets[datasets$dataset != 'co691',] # temporary

# datasets <- datasets[datasets$dataset != 'can673',] # temporary

# toremove <- c("az589","ca681","ca682","ca683","ca684","ca685","ca736","mexi060","mexi063","mt130","ut576")
# datasets <- datasets[!(datasets$dataset %in% toremove),]

# Keep only datasets within WLDAS extent
datasets <- datasets[datasets$north_lat >=25.065 & datasets$north_lat <=52.925 & datasets$east_lon <= -89.025 &  datasets$east_lon >= -124.925,]

datasets <- datasets[datasets$last_year >= 1999,] # at least 20 years of observations

# Temporary, dropping Angiosperms
angiosperms <- c('PLRA', 'QUDG', 'QULO', 'QUGA', 'PPFR', 'PPTR', 'PPDE')
# datasets <- datasets[!(datasets$species_code %in% todrop),]

# Same species
datasets[datasets$species_code == 'ABBI', c('species_name', 'species_code')] <- 
  unique(datasets[datasets$species_code == 'ABLA', c('species_name', 'species_code')])

source(file.path(wd, 'getphylo.R'))
phy.plants.here$tip.label <- sppfull[match(phy.plants.here$tip.label, sppfull$phylo.name),'shortname']

raw_data <- data.frame()
for(d in 1:nrow(datasets)){
  
  raw_data_d <- ringwidth_series[ringwidth_series$dataset == datasets[d, 'dataset'], ]
  
  raw_data_d <- raw_data_d[raw_data_d$year >= 1902 & raw_data_d$year <= 2024 & !is.na(raw_data_d$year), ]
  # raw_data_d <- raw_data_d[ , colSums(is.na(raw_data_d)) < length(1980:2010)]
  
  raw_data_d$species_code <- datasets[d, 'species_code']
  
  # raw_data_long[raw_data_long$rw_core %in% c(-8, -9999), 'rw_core'] <- NA # dealing with weird values in or_105 dataset
  
  # create a unique tree id (across all datasets)
  raw_data_d$tree_id_uniq <- paste0(raw_data_d$dataset, "_", raw_data_d$tree_id)
  
  raw_data_d <- aggregate(
    rw_mm ~ dataset + species_code + tree_id_uniq + year,
    data = raw_data_d,
    FUN = function(x) mean(x, na.rm = TRUE),
    na.action = na.pass
  )
  names(raw_data_d)[names(raw_data_d) == "rw_mm"] <- "rw_avg_mm"
  
  trees_to_remove <- unique(raw_data_d[is.na(raw_data_d$rw_avg_mm) | raw_data_d$rw_avg_mm == 0 , 'tree_id_uniq'])
  raw_data_d <- raw_data_d[!(raw_data_d$tree_id_uniq %in% trees_to_remove),]
  
  raw_data <- rbind(raw_data, raw_data_d)
}

# Gather stands if they have the same latitude and longitude (rounded to 1e-2 degree, ie ~1.1km)
group_keys <- interaction(
  round(datasets$north_lat,2),
  round(datasets$south_lat,2),
  round(datasets$east_lon,2),
  round(datasets$west_lon,2),
  drop = TRUE
)
datasets$grouped_stand <- paste0("S", as.integer(factor(group_keys)))
raw_data <- merge(raw_data,  datasets[, c("dataset", "grouped_stand")], by.x = 'dataset', by.y = 'dataset')

# Load climate data
clim_pred <- readRDS(file.path(wd, "output", "climate", 'climatena',  "climpredictors_11august2025.rds"))
clim_pred$DD5 <- clim_pred$DD5/10 # in x10 degC
clim_pred$FFP <- clim_pred$FFP # in days
PPT_prev <- clim_pred[clim_pred$year %in% 1901:2023, c('dataset', 'year','PPT11', 'PPT12')] 
PPT_prev$year <- PPT_prev$year + 1
names(PPT_prev)[3:4] <- paste0(names(PPT_prev)[3:4], '_prev')
clim_pred <- merge(clim_pred, PPT_prev, all = TRUE)
clim_pred$PPT11to04 <- clim_pred$PPT11_prev + clim_pred$PPT12_prev +
  clim_pred$PPT01 + clim_pred$PPT02 +
  clim_pred$PPT03 + clim_pred$PPT04 
clim_pred$PPT11to04 <- clim_pred$PPT11to04/10 # in cm!
clim_pred <- na.omit(clim_pred)


# Check
climpredsites <- unique(na.omit(clim_pred)$dataset)
rwsites <- unique(raw_data$dataset)
rwsites[!(rwsites %in% climpredsites)]

# Sizes
uniq_tree_ids <- unique(raw_data$tree_id_uniq)
N_trees <- length(uniq_tree_ids)

uniq_stand_ids <- unique(raw_data$grouped_stand)
N_stands <- length(unique(raw_data$grouped_stand))

uniq_species_ids <- unique(raw_data$species_code)
N_species <- length(uniq_species_ids)
clade_idxs <- as.numeric(uniq_species_ids %in% angiosperms)+1 # Gymnosperms=1, Angiosperms=2

all_years <-  min(raw_data$year):max(raw_data$year)
N_all_years <- length(all_years)

# Format data into ragged arrays
log_rw_obs <- c()
gdd_obs <- c()
winterprec_obs <- c()
ffp_obs <- c()
years <- c()
all_years_idxs <- c()
N_years <- c()
stand_idxs <- c()
species_idxs <- c()

idx <- 1
tree_start_idxs <- c()
tree_end_idxs <- c()

for(tid in uniq_tree_ids) {
  print(tid)
  
  raw_data_tree <- raw_data[raw_data$tree_id_uniq == tid & raw_data$year %in% all_years,]
  
  years_tree <- raw_data_tree$year
  all_years_idxs_tree <- sapply(years_tree, function(y) which(all_years == y))
  N_years_tree <- length(years_tree)
  if(N_years_tree > N_all_years | N_years_tree < 20){stop()}
  
  gdd_obs_tree <- sapply(years_tree, 
                         function(y) 
                           clim_pred$DD5[clim_pred$year == y & clim_pred$dataset == raw_data_tree$dataset[1]][1])
  if(any(is.na(gdd_obs_tree))){
    stop(paste0('Missing predictors for stand ', raw_data_tree$dataset[1]))
  }
  gdd_obs <- c(gdd_obs, gdd_obs_tree)
  
  log_rw_obs_tree <- sapply(years_tree, 
                            function(y) 
                              log(raw_data_tree$rw_avg_mm[raw_data_tree$year == y][1]))
  log_rw_obs <- c(log_rw_obs, log_rw_obs_tree)
  
  ffp_obs_tree <- sapply(years_tree, 
                        function(y) 
                          clim_pred$FFP[clim_pred$year == y & clim_pred$dataset == raw_data_tree$dataset[1]][1])
  ffp_obs <- c(ffp_obs, ffp_obs_tree)
  
  winterprec_obs_tree <- sapply(years_tree, 
                         function(y) 
                           clim_pred$PPT11to04[clim_pred$year == y & clim_pred$dataset == raw_data_tree$dataset[1]][1]
                         
                         )
  winterprec_obs <- c(winterprec_obs, winterprec_obs_tree)
  
  years <- c(years, years_tree)
  all_years_idxs <- c(all_years_idxs, all_years_idxs_tree)
  N_years <- c(N_years, N_years_tree)
  
  stand_tree <- which(uniq_stand_ids == raw_data_tree$grouped_stand[1])
  # stand_tree <- which(uniq_stand_ids == raw_data_tree$region[1])
  stand_idxs <- c(stand_idxs, stand_tree)
  
  species_tree <- which(uniq_species_ids == raw_data_tree$species_code[1])
  species_idxs <- c(species_idxs, species_tree)
  
  tree_start_idxs <- c(tree_start_idxs, idx)
  idx <- idx + N_years_tree
  tree_end_idxs <- c(tree_end_idxs, idx - 1)
}

# Cross check sizes
N_trees
length(log_rw_obs)
length(gdd_obs)
length(winterprec_obs)
length(ffp_obs)
length(years)
length(all_years_idxs)
length(N_years)
length(tree_start_idxs)
length(tree_end_idxs)
sum(is.na(gdd_obs)) # check clim. pred


# Collection data into list
N <- length(years)
N_clades <- 2
data <- mget(c('N', 'N_all_years', 'N_trees', 
               'log_rw_obs', 'gdd_obs', 'winterprec_obs', 'ffp_obs',
               'all_years', 'years', 'all_years_idxs', 'N_years', 
               'stand_idxs', 'N_stands',
               'species_idxs', 'N_species',
               'clade_idxs', 'N_clades',
               'tree_start_idxs', 'tree_end_idxs',
               # the following at not necessary for the model
               # but useful for plots, etc.
               'uniq_tree_ids', 'uniq_stand_ids', 'uniq_species_ids'
))
data$years <- as.numeric(data$years)
saveRDS(data, file = file.path(wd, 'output/model/climatena', 'data_11august2025_longonlyus.rds'))
saveRDS(datasets, file.path(wd, 'output/model', 'datasets_11july2025.rds'))
saveRDS(phy.plants.here, file.path(wd, 'output/model', 'phylotree_11july2025.rds'))
