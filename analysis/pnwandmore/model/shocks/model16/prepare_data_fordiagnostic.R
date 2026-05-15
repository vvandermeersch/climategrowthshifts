rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)

library(rstan)
# library(dplyr)

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
  drop = TRUE,
  sep = '/'
)
datasets$grouped_stand <- paste0("S", as.integer(factor(group_keys)))
group_coordinates <- data.frame(group_keys, grouped_stand = paste0("S", as.integer(factor(group_keys))))

# Remove the angiosperms
datasets <- datasets[!(datasets$species_code %in% angiosperms),]

# Remove species that occur at only one stand
toremove <- names(table(datasets$species_code))[table(datasets$species_code) == 1]
datasets <- datasets[!(datasets$species_code %in% toremove),]

# Temporary
# datasets <- datasets[datasets$species_code == 'PIPO',]
# datasets <- datasets[datasets$dataset %in% c('co591', 'mt128', 'wa136', 'ut539',
#                                              'or085', 'wy032', 'id018', 'az599'),]
# datasets <- datasets[datasets$dataset %in% c('co591', 'mt128', 'wa136', 'ut539',
#                                              'or085', 'wy032', 'id018', 'az599',
#                                              'ca729', 'co623', 'id033', 'az608',
#                                              'az601', 'nm610', 'mt167', 'id015'),]
# set.seed(1234)
# datasets32 <- sample(datasets$dataset, 32)
# datasets <- datasets[datasets$dataset %in% datasets32,]

# Temporary
# datasets <- datasets[datasets$species_code %in% c("JUOC", "PICO", "PIFL", "THPL", "TSME", 'SEGI') |
#                        datasets$dataset %in% c('co591', 'mt128', 'wa136', 'ut539', 'or085', 'wy032', 'id018', 'az599'),]
# datasets <- datasets[datasets$species_code %in% c("TSME"),]
# 
# datasets <- datasets[datasets$state %in% c("or"),]


datasets <- datasets[datasets$species_code %in% c('ABCO', 'ABLA', 'ABMA', 'CADE', 'CANO', 'JUOS', 'JUSC', 'LALY', 
                                                  'PCEN', 'PIAL', 'PIAR', 'PIAZ', 'PIJE', 'PILA', 'PILO', 'PISF', 'SEGI'),]

datasets <- datasets[datasets$species_code %in% c('ABCO', 'ABLA', 'ABMA', 'CADE', 'CANO', 'JUOS', 'JUSC', 'LALY', 'PCEN', 'PIAL'),]
                                                  
                                                  
                                                  # Prepare tree ring data!
ringwidth_series <- readRDS(file.path(wd, 'input', 'itrdb', 'ringwidth_series_usonly_from1896.rds'))

all_years <-  1896:2024
# all_years <- min(ringwidth_series$year): max(ringwidth_series$year)

raw_data <- data.frame()
for(d in 1:nrow(datasets)){
  
  raw_data_d <- ringwidth_series[ringwidth_series$dataset %in% datasets[d, 'dataset'], ]
  
  raw_data_d <- raw_data_d[raw_data_d$year %in% all_years, ]
  
  raw_data_d$species_code <- datasets[d, 'species_code']
  
  # temporary
  # raw_data_d <- raw_data_d[raw_data_d$year >= 1950 & raw_data_d$year <= 2025,]
  
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
  
  cat(paste0(length(unique(raw_data_d$tree_id_uniq)), ' unique trees!\n'))
  
  if(length(unique(raw_data_d$tree_id_uniq)) > 10){
    randomtrees <- sample(unique(raw_data_d$tree_id_uniq), 10)
    raw_data_d <- raw_data_d[raw_data_d$tree_id_uniq %in% randomtrees, ]
  }
  
  raw_data <- rbind(raw_data, raw_data_d)
}
raw_data <- merge(raw_data,  datasets[, c("dataset", "grouped_stand")], by.x = 'dataset', by.y = 'dataset')
length( unique(datasets$grouped_stand))

all_years <- min(raw_data$year): max(raw_data$year)

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

uniq_stand_species_ids <- unique(paste(raw_data$grouped_stand, raw_data$species_code, sep = '_'))
N_stand_species <- length(uniq_stand_species_ids)

uniq_species_ids <- unique(raw_data$species_code)
N_species <- length(uniq_species_ids)
clade_idxs <- as.numeric(uniq_species_ids %in% angiosperms)+1 # Gymnosperms=1, Angiosperms=2

# all_years <-  min(raw_data$year):max(raw_data$year)
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
  
  rw_obs_tree <- sapply(years_tree, 
                        function(y) 
                          raw_data_tree$rw_avg_mm[raw_data_tree$year == y][1])
  rw_obs <- c(rw_obs, rw_obs_tree)
  
  years <- c(years, years_tree)
  all_years_idxs <- c(all_years_idxs, all_years_idxs_tree)
  N_years <- c(N_years, N_years_tree)
  
  stand_tree <- which(uniq_stand_ids == raw_data_tree$grouped_stand[1])
  stand_idxs <- c(stand_idxs, stand_tree)
  
  stand_species_tree <- which(uniq_stand_species_ids == paste0(raw_data_tree$grouped_stand[1], '_', raw_data_tree$species_code[1]))
  stand_species_idxs <- c(stand_species_idxs, stand_species_tree)
  
  species_tree <- which(uniq_species_ids == raw_data_tree$species_code[1])
  species_idxs <- c(species_idxs, species_tree)
  
  tree_start_idxs <- c(tree_start_idxs, idx)
  idx <- idx + N_years_tree
  tree_end_idxs <- c(tree_end_idxs, idx - 1)
}

N_stand_trees <- c()
N_stand_years <- c()
stand_start_years_idxs <- c()
idx <- 1
stand_trees_start_idxs <- c()
stand_trees_end_idxs <- c()
stand_tree_idxs <- c()

gdd_obs  <- c()
pre_obs <- c()
vpd_obs <- c()

# for cheap comparison
clim_pred_climatena <- readRDS(file.path(wd, "output", "climate", "climatena",  "climpredictors_11august2025.rds"))
PPT_prev <- clim_pred_climatena[clim_pred_climatena$year %in% 1901:2023, c('dataset', 'year','PPT11', 'PPT12')] 
PPT_prev$year <- PPT_prev$year + 1
names(PPT_prev)[3:4] <- paste0(names(PPT_prev)[3:4], '_prev')
clim_pred_climatena <- merge(clim_pred_climatena, PPT_prev, all = TRUE)
clim_pred_climatena$PPT11to04 <- clim_pred_climatena$PPT11_prev + clim_pred_climatena$PPT12_prev +
  clim_pred_climatena$PPT01 + clim_pred_climatena$PPT02 +
  clim_pred_climatena$PPT03 + clim_pred_climatena$PPT04 
clim_pred_climatena$PPT11to04 <- clim_pred_climatena$PPT11to04/100 # in dm!
clim_pred_climatena <- na.omit(clim_pred_climatena)
pre_obs_climatena <- c() 

for(s in uniq_stand_ids){
  # if(s == 'S30'){stop()}
  N_stand_trees_s <- sum(stand_idxs == which(uniq_stand_ids == s))
  N_stand_trees <- c(N_stand_trees, N_stand_trees_s)
  
  idx_prev <- idx
  stand_trees_start_idxs <- c(stand_trees_start_idxs, idx)
  idx <- idx + N_stand_trees_s
  stand_trees_end_idxs <- c(stand_trees_end_idxs, idx - 1)
  
  stand_tree_idxs_here <- which(stand_idxs == which(uniq_stand_ids == s))
  stand_tree_idxs <- c(stand_tree_idxs, stand_tree_idxs_here)
  
  min_year <- min(all_years_idxs[tree_start_idxs[stand_tree_idxs_here]])
  max_year <- max(all_years_idxs[tree_end_idxs[stand_tree_idxs_here]])
  print(max(idx_prev:(idx-1)))
  
  N_stand_years <- c(N_stand_years, length(min_year:max_year))
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
  
  # cheap comparison
  pre_obs_stand_climatena <- sapply(all_years, 
                          function(y) 
                            as.numeric(clim_pred_climatena$PPT11to04[clim_pred_climatena$year == y & clim_pred_climatena$dataset == dataset]))
  pre_obs_climatena <- c(pre_obs_climatena, unlist(pre_obs_stand_climatena))
}

# Cheap comparison
par(mfrow = c(1,1))
hist((pre_obs-pre_obs_climatena)/pre_obs_climatena*100)
quantile((pre_obs-pre_obs_climatena)/pre_obs_climatena*100, c(0.05,0.95))


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
N_clades <- 1
data <- mget(c('N', 'N_all_years', 'N_trees', 
               'rw_obs', 'gdd_obs', 
               'vpd_obs', 'pre_obs',
               'all_years', 'years', 'all_years_idxs', 'N_years', 
               'stand_idxs', 'N_stands', 
               'N_stand_trees', 'stand_tree_idxs', 
               'stand_trees_start_idxs', 'stand_trees_end_idxs',
               'N_stand_years', 'stand_start_years_idxs',
               'stand_species_idxs', 'N_stand_species', 
               'species_idxs', 'N_species',
               'clade_idxs', 'N_clades',
               'tree_start_idxs', 'tree_end_idxs',
               # the following at not necessary for the model
               # but useful for plots, etc.
               'uniq_tree_ids', 'uniq_stand_ids', 'uniq_species_ids'
))
data$years <- as.numeric(data$years)

group_coordinates <- unique(group_coordinates)
data$uniq_stand_lat <- stringr::str_split_i(group_coordinates[group_coordinates$grouped_stand %in% data$uniq_stand_ids, "group_keys"], '/', 1)
data$uniq_stand_lon <- stringr::str_split_i(group_coordinates[group_coordinates$grouped_stand %in% data$uniq_stand_ids, "group_keys"], '/', 3)


saveRDS(data, file = file.path(wd, 'output/model', 'data_15mar2026_10species_10trees_18962024.rds'))
saveRDS(datasets, file.path(wd, 'output/model', 'datasets_15mar2026_10species_10trees_18962024.rds'))


