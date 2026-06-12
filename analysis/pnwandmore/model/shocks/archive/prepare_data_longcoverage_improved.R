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

# Keep only US (because of PRISM)
datasets <- datasets[!(datasets$state %in% c("can", "mexi")), ]

# Defining Angiosperms
angiosperms <- c('PLRA', 'QUDG', 'QULO', 'QUGA', 'PPFR', 'PPTR', 'PPDE')

# Temporary - missing species 
datasets <- datasets[datasets$dataset != 'co689',] # temporary 
datasets <- datasets[datasets$dataset != 'co690',] # temporary
datasets <- datasets[datasets$dataset != 'co691',] # temporary

# Same species
datasets[datasets$species_code == 'ABBI', c('species_name', 'species_code')] <- 
  unique(datasets[datasets$species_code == 'ABLA', c('species_name', 'species_code')])

# Create a grouped_stand key to gather stands if they have the same latitude and longitude (rounded to 1e-2 degree, ie ~1.1km)
group_keys <- interaction(
  round(datasets$north_lat,2),
  round(datasets$south_lat,2),
  round(datasets$east_lon,2),
  round(datasets$west_lon,2),
  drop = TRUE
)
datasets$grouped_stand <- paste0("S", as.integer(factor(group_keys)))

# Temporary
# datasets <- datasets[datasets$species_code == 'PIPO',]
# datasets <- datasets[datasets$dataset %in% c('co591', 'mt128', 'wa136', 'ut539',
#                                              'or085', 'wy032', 'id018', 'az599'),]
# datasets <- datasets[datasets$dataset %in% c('co591', 'mt128', 'wa136', 'ut539',
#                                              'or085', 'wy032', 'id018', 'az599',
#                                              'ca729', 'co623', 'id033', 'az608',
#                                              'az601', 'nm610', 'mt167', 'id015'),]

# Prepare tree ring data!
ringwidth_series <- readRDS(file.path(wd, 'input', 'itrdb', 'ringwidth_series_usonly_from1896.rds'))

raw_data <- data.frame()
for(d in 1:nrow(datasets)){
  
  raw_data_d <- ringwidth_series[ringwidth_series$dataset == datasets[d, 'dataset'], ]
  
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
length( unique(datasets$grouped_stand))

# Deal with potential duplicates
# If several trees from different ITRDB datasets are on the same site and have the same ID, the same number of years
# and the same mean ringwidth, we keep only one of them
potential_duplicates <- aggregate(dataset ~ original_tree_id + species_code, data = raw_data, FUN = function(x) length(unique(x)))
potential_duplicates_groupedstand <- aggregate(grouped_stand ~ original_tree_id + species_code, data = raw_data, FUN = function(x) length(unique(x)))
potential_duplicates <- merge(potential_duplicates, potential_duplicates_groupedstand)
potential_duplicates <- potential_duplicates[potential_duplicates$dataset > 1 & potential_duplicates$grouped_stand == 1,]
cat(paste0(nrow(potential_duplicates), ' trees are potential duplicates!\n'))
for(t in potential_duplicates$original_tree_id){
  sp_t <- potential_duplicates[potential_duplicates$original_tree_id == t, 'species_code']
  raw_data_t <- raw_data[raw_data$original_tree_id == t & raw_data$species_code == sp_t, ] 
  count_years <- aggregate(year ~ tree_id_uniq, data = raw_data_t,  FUN = length)
  if(length(unique(count_years$year)) == 1){
    mean_rw <- aggregate(rw_avg_mm ~ tree_id_uniq, data = raw_data_t,  FUN = mean)
    if(length(unique(mean_rw$rw_avg_mm)) == 1){
      tree_to_keep <- sample(unique(raw_data_t$tree_id_uniq), 1)
      raw_data <- raw_data[-which(raw_data$original_tree_id == t & raw_data$tree_id_uniq != tree_to_keep),]
    }else{
      if(t == 'pic382'){
        # keep both trees because they are different
      }else{
        print(t)
        stop()
      }
    }
  }else if(length(unique(count_years$year)) > 1){
    if(t == 'at2022'){
      tree_to_keep <- 'co589_at2022' # longer time series
      raw_data <- raw_data[-which(raw_data$original_tree_id == t & raw_data$tree_id_uniq != tree_to_keep),]
    }else{
      print(t)
      stop()
    }
    
  }
}
datasets <- datasets[datasets$dataset %in% unique(raw_data$dataset),]
length( unique(datasets$grouped_stand))

# Load climate data
clim_pred <- readRDS(file.path(wd, "output", "climate", "prism",  "climpredictors_20dec2025.rds"))
clim_pred$gdd_all <- clim_pred$gdd_all/100 # in x100 degC
clim_pred$vpd_mjja <- clim_pred$vpd_mjja # in hPa?
clim_pred$ppt_ndjfma <- clim_pred$ppt_ndjfma/100 # in dm? 
clim_pred <- na.omit(clim_pred)

# Check
climpredsites <- unique(na.omit(clim_pred)$dataset)
rwsites <- unique(raw_data$dataset)
rwsites[!(rwsites %in% climpredsites)]

# Sizes
uniq_tree_ids <- unique(raw_data$tree_id_uniq)
N_trees <- length(uniq_tree_ids)

uniq_stand_ids <- unique(raw_data$grouped_stand)
N_stands <- length(uniq_stand_ids)

uniq_substand_ids <- unique(paste(raw_data$grouped_stand, raw_data$species_code, sep = '_'))
N_substands <- length(uniq_substand_ids)

uniq_species_ids <- unique(raw_data$species_code)
N_species <- length(uniq_species_ids)
clade_idxs <- as.numeric(uniq_species_ids %in% angiosperms)+1 # Gymnosperms=1, Angiosperms=2

all_years <-  min(raw_data$year):max(raw_data$year)
N_all_years <- length(all_years)

# Format data into ragged arrays
rw_obs <- c()
years <- c()
all_years_idxs <- c()
N_years <- c()
stand_idxs <- c()
stand_species_idxs <- c()
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
  # if(N_years_tree > 45 | N_years_tree < 20){stop()}
  
  # gdd_obs_tree <- sapply(years_tree, 
  #                        function(y) 
  #                          as.numeric(clim_pred$gdd_all[clim_pred$year == y & clim_pred$dataset == raw_data_tree$dataset[1]][1]))
  # if(any(is.na(gdd_obs_tree))){
  #   stop(paste0('Missing predictors for stand ', raw_data_tree$dataset[1]))
  # }
  # gdd_obs <- c(gdd_obs, gdd_obs_tree)
  
  rw_obs_tree <- sapply(years_tree, 
                        function(y) 
                          raw_data_tree$rw_avg_mm[raw_data_tree$year == y][1])
  rw_obs <- c(rw_obs, rw_obs_tree)
  
  # pre_obs_tree <- sapply(years_tree, 
  #                            function(y) 
  #                              as.numeric(clim_pred$ppt_ndjfma[clim_pred$year == y & clim_pred$dataset == raw_data_tree$dataset[1]][1]))
  # pre_obs <- c(pre_obs, pre_obs_tree)
  # 
  # vpd_obs_tree <- sapply(years_tree, 
  #                        function(y) 
  #                          as.numeric(clim_pred$vpd_mjja[clim_pred$year == y & clim_pred$dataset == raw_data_tree$dataset[1]][1]))
  # vpd_obs <- c(vpd_obs, vpd_obs_tree)
  
  years <- c(years, years_tree)
  all_years_idxs <- c(all_years_idxs, all_years_idxs_tree)
  N_years <- c(N_years, N_years_tree)
  
  stand_tree <- which(uniq_stand_ids == raw_data_tree$grouped_stand[1])
  stand_idxs <- c(stand_idxs, stand_tree)
  
  # stand_species_tree <- which(uniq_stand_species_ids == paste0(raw_data_tree$grouped_stand[1], '_', raw_data_tree$species_code[1]))
  # stand_species_idxs <- c(stand_species_idxs, stand_species_tree)
  
  species_tree <- which(uniq_species_ids == raw_data_tree$species_code[1])
  species_idxs <- c(species_idxs, species_tree)
  
  tree_start_idxs <- c(tree_start_idxs, idx)
  idx <- idx + N_years_tree
  tree_end_idxs <- c(tree_end_idxs, idx - 1)
}

N_trees_stand <- c()
N_years_stand <- c()
stand_start_years_idxs <- c()
idx <- 1
idx_sbs <- 1
stand_trees_start_idxs <- c()
stand_trees_end_idxs <- c()

gdd_obs  <- c()
pre_obs <- c()
vpd_obs <- c()

N_substands_stand <- c()
stand_substand_start_idxs <- c()
stand_substand_end_idxs <- c()
substand_species_idxs <- c()

for(s in uniq_stand_ids){

  N_trees_stand_s <- sum(stand_idxs == which(uniq_stand_ids == s))
  N_trees_stand <- c(N_trees_stand, N_trees_stand_s)
  
  idx_prev <- idx
  stand_trees_start_idxs <- c(stand_trees_start_idxs, idx)
  idx <- idx + N_trees_stand_s
  stand_trees_end_idxs <- c(stand_trees_end_idxs, idx - 1)
  
  min_year <- min(all_years_idxs[tree_start_idxs[idx_prev:(idx-1)]])
  max_year <- max(all_years_idxs[tree_end_idxs[idx_prev:(idx-1)]])
  
  N_years_stand <- c(N_years_stand, length(min_year:max_year))
  stand_start_years_idxs <- c(stand_start_years_idxs, min_year)
  
  dataset <- datasets[datasets$grouped_stand ==s, 'dataset'][1]
  
  gdd_obs_stand <- sapply(all_years, 
                          function(y) 
                            as.numeric(clim_pred$gdd_all[clim_pred$year == y & clim_pred$dataset == dataset]))
  gdd_obs <- c(gdd_obs, gdd_obs_stand)
  
  pre_obs_stand <- sapply(all_years, 
                          function(y) 
                            as.numeric(clim_pred$ppt_ndjfma[clim_pred$year == y & clim_pred$dataset == dataset]))
  pre_obs <- c(pre_obs, pre_obs_stand)
  
  vpd_obs_stand <- sapply(all_years, 
                          function(y) 
                            as.numeric(clim_pred$vpd_mjja[clim_pred$year == y & clim_pred$dataset == dataset]))
  vpd_obs <- c(vpd_obs, vpd_obs_stand)
  
  N_substands_stand_s <- length(unique(species_idxs[stand_idxs == which(uniq_stand_ids == s)]))
  N_substands_stand <- c(N_substands_stand, N_substands_stand_s)
  
  stand_substand_start_idxs <- c(stand_substand_start_idxs, idx_sbs)
  idx_sbs <- idx_sbs + N_substands_stand_s
  stand_substand_end_idxs <- c(stand_substand_end_idxs, idx_sbs - 1)
  
  stand_species <- unique(species_idxs[stand_idxs == which(uniq_stand_ids == s)])
  substand_species_idxs <- c(substand_species_idxs, stand_species)
}

substand_idxs <- c()
substand_rel_idxs <- c()
for(tid in uniq_tree_ids){
  
  raw_data_tree <- raw_data[raw_data$tree_id_uniq == tid & raw_data$year %in% all_years,]
  
  substand <- which(uniq_substand_ids == unique(paste(raw_data_tree$grouped_stand, raw_data_tree$species_code, sep = '_')))
  substand_idxs <- c(substand_idxs, substand)
  
  stand <- stand_idxs[uniq_tree_ids == tid]
  stand_species <- substand_species_idxs[stand_substand_start_idxs[stand]:stand_substand_end_idxs[stand]]
  species <- species_idxs[uniq_tree_ids == tid]
  if(!(species %in% stand_species)){stop()} # safety check
  substand_rel_idxs <- c(substand_rel_idxs, which(stand_species == species))
  
}



# Cross check sizes
N_trees
length(rw_obs)
length(gdd_obs)
length(vpd_obs)
length(pre_obs)
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
               'rw_obs', 'gdd_obs', 
               'vpd_obs', 'pre_obs',
               'all_years', 'years', 'all_years_idxs', 'N_years', 
               'stand_idxs', 'N_stands', 
               'N_trees_stand', 'stand_trees_start_idxs', 'stand_trees_end_idxs',
               'N_years_stand', 'stand_start_years_idxs',
               
               'N_substands', 'N_substands_stand', 
               'substand_idxs', 'substand_species_idxs', 'substand_rel_idxs',
               'stand_substand_start_idxs', 'stand_substand_end_idxs',
               
               'species_idxs', 'N_species',
               'clade_idxs', 'N_clades',
               'tree_start_idxs', 'tree_end_idxs',
               # the following at not necessary for the model
               # but useful for plots, etc.
               'uniq_tree_ids', 'uniq_stand_ids', 'uniq_species_ids'
))
data$years <- as.numeric(data$years)
saveRDS(data, file = file.path(wd, 'output/model', 'data_28dec2025_standpred_improved.rds'))
saveRDS(datasets, file.path(wd, 'output/model', 'datasets_28dec2025_stanpred_improved.rds'))


