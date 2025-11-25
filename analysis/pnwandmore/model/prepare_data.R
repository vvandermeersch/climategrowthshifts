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

# Defining Angiosperms
angiosperms <- c('PLRA', 'QUDG', 'QULO', 'QUGA', 'PPFR', 'PPTR', 'PPDE')

# Dropping two Mexican species whose ranges are outside WLDAS extent
todrop <- c('PIGR', 'PICU')
datasets <- datasets[!(datasets$species_code %in% todrop),]

# Dropping Angiosperms (temporary)
datasets <- datasets[!(datasets$species_code %in% angiosperms),]

# Keep only some states (temporary)
datasets <- datasets[datasets$state %in% c('az', 'nm', 'ca', 'nv', 'ut', 'co'),]

# Keep only 6 datasets (three sites)
# datasets <- datasets[datasets$dataset %in% c('az591', 'az592', 'az566', 'az567', 'az689', 'az569'),]

# Keep only some species (temporary)
datasets <- datasets[datasets$species_code %in% c('PIPO', 'PSME', 'PIED', 'JUOS', 'PIAZ'),]

# Same species
datasets[datasets$species_code == 'ABBI', c('species_name', 'species_code')] <- 
  unique(datasets[datasets$species_code == 'ABLA', c('species_name', 'species_code')])

source(file.path(wd, 'getphylo.R'))
phy.plants.here$tip.label <- sppfull[match(phy.plants.here$tip.label, sppfull$phylo.name),'shortname']

# Create a grouped_stand key to gather stands if they have the same latitude and longitude (rounded to 1e-2 degree, ie ~1.1km)
group_keys <- interaction(
  round(datasets$north_lat,2),
  round(datasets$south_lat,2),
  round(datasets$east_lon,2),
  round(datasets$west_lon,2),
  drop = TRUE
)
datasets$grouped_stand <- paste0("S", as.integer(factor(group_keys)))

# Selected a subset of sites with at least 2 species!
count_species <- aggregate(species_code ~ grouped_stand, data = datasets, FUN = function(x) length(unique(x)))
selected_stands <- count_species[count_species$species_code >= 2, 'grouped_stand']
datasets <- datasets[datasets$grouped_stand %in% selected_stands,]

# Prepare tree ring data!
ringwidth_series <- readRDS(file.path(wd, 'input', 'itrdb', 'ringwidth_series_all.rds'))

raw_data <- data.frame()
for(d in 1:nrow(datasets)){
  
  raw_data_d <- ringwidth_series[ringwidth_series$dataset == datasets[d, 'dataset'], ]
  
  raw_data_d <- raw_data_d[raw_data_d$year >= 1980 & raw_data_d$year <= 2023 & !is.na(raw_data_d$year), ]
  # raw_data_d <- raw_data_d[ , colSums(is.na(raw_data_d)) < length(1980:2010)]
  
  raw_data_d$species_code <- datasets[d, 'species_code']
  
  # raw_data_long[raw_data_long$rw_core %in% c(-8, -9999), 'rw_core'] <- NA # dealing with weird values in or_105 dataset
  
  # create a unique tree id (across all datasets)
  raw_data_d$original_tree_id <- raw_data_d$tree_id
  raw_data_d$tree_id_uniq <- paste0(raw_data_d$dataset, "_", raw_data_d$tree_id)
  
  raw_data_d <- aggregate(
    rw_mm ~ dataset + species_code + original_tree_id + tree_id_uniq + year,
    data = raw_data_d,
    FUN = function(x) mean(x, na.rm = TRUE),
    na.action = na.pass
  )
  names(raw_data_d)[names(raw_data_d) == "rw_mm"] <- "rw_avg_mm"
  
  raw_data_d <- na.omit(raw_data_d)
  raw_data_d$year <- as.numeric(raw_data_d$year)
  
  # Here, we check that there is no missing year in the individual time series
  count_bytrees <- aggregate(
    rw_avg_mm ~ tree_id_uniq,
    data = raw_data_d,
    FUN = function(x) length(x)
  )
  years_bytrees <- aggregate(
    year ~ tree_id_uniq,
    data = raw_data_d,
    FUN = function(x) max(x)-min(x)+1
  )
  check_length <- merge(count_bytrees, years_bytrees)
  trees_to_remove <- check_length[check_length$rw_avg_mm != check_length$year, 'tree_id_uniq']
  raw_data_d <- raw_data_d[!(raw_data_d$tree_id_uniq %in% trees_to_remove),]
  
  # remove less than 20 years observed
  trees_to_remove <- check_length[check_length$rw_avg_mm < 20, 'tree_id_uniq']
  raw_data_d <- raw_data_d[!(raw_data_d$tree_id_uniq %in% trees_to_remove),]
  
  # remove trees with 0mm observations - NOT ANYMORE
  # trees_to_remove <- unique(raw_data_d[raw_data_d$rw_avg_mm == 0 , 'tree_id_uniq'])
  # raw_data_d <- raw_data_d[!(raw_data_d$tree_id_uniq %in% trees_to_remove),]
  
  raw_data <- rbind(raw_data, raw_data_d)
}
raw_data <- merge(raw_data,  datasets[, c("dataset", "grouped_stand")], by.x = 'dataset', by.y = 'dataset')

# Deal with potential duplicates
# If several trees from different ITRDB datasets are on the same site and have the same ID, the same number of years
# and the same mean ringwidth, we keep only one of them
potential_duplicates <- aggregate(dataset ~ original_tree_id, data = raw_data, FUN = function(x) length(unique(x)))
potential_duplicates_groupedstand <- aggregate(grouped_stand ~ original_tree_id, data = raw_data, FUN = function(x) length(unique(x)))
potential_duplicates <- merge(potential_duplicates, potential_duplicates_groupedstand)
potential_duplicates <- potential_duplicates[potential_duplicates$dataset > 1 & potential_duplicates$grouped_stand == 1,]
cat(paste0(nrow(potential_duplicates), ' trees are potential duplicates!\n'))
for(t in potential_duplicates$original_tree_id){
  raw_data_t <- raw_data[raw_data$original_tree_id == t, ] 
  count_years <- aggregate(year ~ tree_id_uniq, data = raw_data_t,  FUN = length)
  if(length(unique(count_years$year)) == 1){
    mean_rw <- aggregate(rw_avg_mm ~ tree_id_uniq, data = raw_data_t,  FUN = mean)
    if(length(unique(mean_rw$rw_avg_mm)) == 1){
      tree_to_keep <- sample(unique(raw_data_t$tree_id_uniq), 1)
      raw_data <- raw_data[-which(raw_data$original_tree_id == t & raw_data$tree_id_uniq != tree_to_keep),]
    }else{
      print(t)
      stop()
    }
  }else if(length(unique(count_years$year)) > 1){
    print(t)
    stop()
  }
}
datasets <- datasets[datasets$dataset %in% unique(raw_data$dataset),]

# Create a substand key to gather trees from the same species on the same grouped_stand
raw_data$substand <- paste0(raw_data$grouped_stand, '_', raw_data$species_code)

# Load climate data
clim_pred <- readRDS(file.path(wd, "output", "climate",  "climpredictors_10oct2025.rds"))
clim_pred$soilmoist_mjj <- clim_pred$soilmoist_2m_mjj*100 # in percentage
clim_pred$gdd <- clim_pred$gdd/100 # in x100 degC
clim_pred$gdd_amjjas <- clim_pred$gdd_amjjas/100 # in x100 degC
clim_pred$vpd_mjj <- clim_pred$vpd_mjj # in hPa?
clim_pred$pre_jja <- clim_pred$pre_jja/10 # in cm? 
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

uniq_substand_ids <- unique(raw_data$substand)
N_substands <- length(unique(raw_data$substand))

uniq_species_ids <- unique(raw_data$species_code)
N_species <- length(uniq_species_ids)
clade_idxs <- as.numeric(uniq_species_ids %in% angiosperms)+1 # Gymnosperms=1, Angiosperms=2

all_years <-  min(raw_data$year):max(raw_data$year)
N_all_years <- length(all_years)

# Format data into ragged arrays
log_rw_obs <- c()
gdd_obs <- c()
gdd_amjjas_obs  <- c()
sm_obs <- c()
vpd_obs <- c()
pre_jja_obs <- c()
years <- c()
all_years_idxs <- c()
N_years <- c()
stand_idxs <- c()
substand_idxs <- c()
species_idxs <- c()

idx <- 1
tree_start_idxs <- c()
tree_end_idxs <- c()

aconstant <- 1e-04
for(tid in uniq_tree_ids) {
  print(tid)
  
  raw_data_tree <- raw_data[raw_data$tree_id_uniq == tid & raw_data$year %in% all_years,]
  
  years_tree <- raw_data_tree$year
  all_years_idxs_tree <- sapply(years_tree, function(y) which(all_years == y))
  N_years_tree <- length(years_tree)
  if(N_years_tree > 45 | N_years_tree < 20){stop()}
  
  gdd_obs_tree <- sapply(years_tree, 
                         function(y) 
                           as.numeric(clim_pred$gdd[clim_pred$year == y & clim_pred$dataset == raw_data_tree$dataset[1]][1]))
  if(any(is.na(gdd_obs_tree))){
    stop(paste0('Missing predictors for stand ', raw_data_tree$dataset[1]))
  }
  gdd_obs <- c(gdd_obs, gdd_obs_tree)
  
  log_rw_obs_tree <- sapply(years_tree, 
                            function(y) 
                              log(aconstant + raw_data_tree$rw_avg_mm[raw_data_tree$year == y][1]))
  log_rw_obs <- c(log_rw_obs, log_rw_obs_tree)
  
  gdd_obs_amjjas_tree <- sapply(years_tree, 
                                function(y) 
                                  as.numeric(clim_pred$gdd_amjjas[clim_pred$year == y & clim_pred$dataset == raw_data_tree$dataset[1]][1]))
  gdd_amjjas_obs <- c(gdd_amjjas_obs, gdd_obs_amjjas_tree)
  
  pre_obs_jja_tree <- sapply(years_tree, 
                                function(y) 
                                  as.numeric(clim_pred$pre_jja[clim_pred$year == y & clim_pred$dataset == raw_data_tree$dataset[1]][1]))
  pre_jja_obs <- c(pre_jja_obs, pre_obs_jja_tree)
  
  sm_obs_tree <- sapply(years_tree, 
                        function(y) 
                          as.numeric(clim_pred$soilmoist_mjj[clim_pred$year == y & clim_pred$dataset == raw_data_tree$dataset[1]][1]))
  sm_obs <- c(sm_obs, sm_obs_tree)
  
  vpd_obs_tree <- sapply(years_tree, 
                         function(y) 
                           as.numeric(clim_pred$vpd_mjj[clim_pred$year == y & clim_pred$dataset == raw_data_tree$dataset[1]][1]))
  vpd_obs <- c(vpd_obs, vpd_obs_tree)
  
  years <- c(years, years_tree)
  all_years_idxs <- c(all_years_idxs, all_years_idxs_tree)
  N_years <- c(N_years, N_years_tree)
  
  stand_tree <- which(uniq_stand_ids == raw_data_tree$grouped_stand[1])
  # stand_tree <- which(uniq_stand_ids == raw_data_tree$region[1])
  stand_idxs <- c(stand_idxs, stand_tree)
  
  substand_tree <- which(uniq_substand_ids == raw_data_tree$substand[1])
  substand_idxs <- c(substand_idxs, substand_tree)
  
  species_tree <- which(uniq_species_ids == raw_data_tree$species_code[1])
  species_idxs <- c(species_idxs, species_tree)
  
  tree_start_idxs <- c(tree_start_idxs, idx)
  idx <- idx + N_years_tree
  tree_end_idxs <- c(tree_end_idxs, idx - 1)
}

substand_species_idxs <- c()
for(sid in uniq_substand_ids){
  substand_species_idxs <- c(
    substand_species_idxs, 
    which(uniq_species_ids == raw_data[raw_data$substand == sid, 'species_code'][1])
  )
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
sum(is.na(gdd_obs)) # check clim. pred


# Collection data into list
N <- length(years)
N_clades <- 1
data <- mget(c('N', 'N_all_years', 'N_trees', 
               'log_rw_obs', 'gdd_obs', 'gdd_amjjas_obs', 
               'sm_obs', 'vpd_obs', 'pre_jja_obs',
               'all_years', 'years', 'all_years_idxs', 'N_years', 
               'stand_idxs', 'N_stands', 'substand_idxs', 'N_substands', 'substand_species_idxs',
               'species_idxs', 'N_species',
               'clade_idxs', 'N_clades',
               'tree_start_idxs', 'tree_end_idxs',
               # the following at not necessary for the model
               # but useful for plots, etc.
               'uniq_tree_ids', 'uniq_stand_ids', 'uniq_substand_ids', 'uniq_species_ids'
               ))
data$years <- as.numeric(data$years)
saveRDS(data, file = file.path(wd, 'output/model', 'data_05nov2025_southwest_5species_min2speciesplots.rds'))
saveRDS(datasets, file.path(wd, 'output/model', 'datasets_05nov2025_southwest_5species_min2speciesplots.rds'))
