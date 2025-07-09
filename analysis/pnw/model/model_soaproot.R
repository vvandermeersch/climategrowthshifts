
rm(list = ls())
wd <- '/home/victor/projects/climategrowthshifts/analysis/pnw'
library(dplR)
library(rstan)
rstan_options(auto_write = TRUE)

# Read RWL file
years <- years <- c(1983:2025)
raw_data <- dplR::read.rwl(file.path(wd, 'data', 'others', 'soaproot_saddle', 'PonderosaMaster.rwl'))
idrows <- which(rownames(raw_data) %in% as.character(years))
raw_data <- raw_data[idrows,]
raw_data$year <- rownames(raw_data)
rwdat_long <- reshape(raw_data, varying = names(raw_data)[names(raw_data) != "year"],
                      v.names = "rw_mm", timevar = "core_id",
                      times = names(raw_data)[names(raw_data) != "year"],
                      direction = "long")
rwdat_long$tree_id <- sub("^(.*\\d)[^\\d]*$", "\\1", rwdat_long$core_id)
rwdat_long <- rwdat_long[, c("year", 'tree_id', "core_id", "rw_mm")]

raw_data <- aggregate(
  rw_mm ~ tree_id + year,
  data = rwdat_long,
  FUN = function(x) mean(x, na.rm = TRUE),
  na.action = na.pass
)
names(raw_data)[names(raw_data) == "rw_mm"] <- "rw_avg_mm"
raw_data <- na.omit(raw_data)

treecounts <- aggregate(year ~ tree_id, data = raw_data, FUN = length)
treestoremove <- c(treecounts[treecounts$year < 20, 'tree_id'], 'SS04-16')
raw_data <- raw_data[!(raw_data$tree_id %in% treestoremove),]


# Load climate data
clim_pred <- readRDS(file.path(wd, 'output', 'climate', 'climpredictors_soaproot.rds'))
clim_pred$soilmoist_mjj <- clim_pred$soilmoist_1m_mjj*100 # in percentage
clim_pred$gdd <- clim_pred$gdd/100 # in x100 degC
clim_pred$vpd_mjj <- clim_pred$vpd_mjj # in hPa?
clim_pred <- na.omit(clim_pred)

# Sizes
uniq_tree_ids <- unique(raw_data$tree_id)
N_trees <- length(uniq_tree_ids)

all_years <-  min(raw_data$year):max(raw_data$year)
N_all_years <- length(all_years)

# Format data into ragged arrays
log_rw_obs <- c()
gdd_obs <- c()
sm_obs <- c()
vpd_obs <- c()
vpd_obs_n1 <- c()
vpd_obs_n2 <- c()
vpd_obs_n3 <- c()
years <- c()
all_years_idxs <- c()
N_years <- c()

idx <- 1
tree_start_idxs <- c()
tree_end_idxs <- c()

for(tid in uniq_tree_ids) {
  print(tid)
  
  raw_data_tree <- raw_data[raw_data$tree_id == tid & raw_data$year %in% all_years,]
  
  years_tree <- raw_data_tree$year
  all_years_idxs_tree <- sapply(years_tree, function(y) which(all_years == y))
  N_years_tree <- length(years_tree)
  if(N_years_tree > 45 | N_years_tree < 20){stop()}
  
  gdd_obs_tree <- sapply(years_tree, 
                         function(y) 
                           as.numeric(clim_pred$gdd[clim_pred$year == y][1]))
  if(any(is.na(gdd_obs_tree))){
    stop(paste0('Missing predictors for stand ', raw_data_tree$dataset[1]))
  }
  gdd_obs <- c(gdd_obs, gdd_obs_tree)
  
  log_rw_obs_tree <- sapply(years_tree, 
                            function(y) 
                              log(raw_data_tree$rw_avg_mm[raw_data_tree$year == y][1]))
  log_rw_obs <- c(log_rw_obs, log_rw_obs_tree)
  
  
  
  sm_obs_tree <- sapply(years_tree, 
                        function(y) 
                          as.numeric(clim_pred$soilmoist_mjj[clim_pred$year == y][1]))
  sm_obs <- c(sm_obs, sm_obs_tree)
  
  vpd_obs_tree <- sapply(years_tree, 
                         function(y) 
                           as.numeric(clim_pred$vpd_mjj[clim_pred$year == y][1]))
  vpd_obs <- c(vpd_obs, vpd_obs_tree)
  
  vpd_obs_tree_n1 <- sapply(years_tree, 
                         function(y) 
                           as.numeric(clim_pred$vpd_mjj[clim_pred$year == as.character(as.numeric(y)-1)][1]))
  vpd_obs_n1 <- c(vpd_obs_n1, vpd_obs_tree_n1)
  
  vpd_obs_tree_n2 <- sapply(years_tree, 
                         function(y) 
                           as.numeric(clim_pred$vpd_mjj[clim_pred$year == as.character(as.numeric(y)-2)][1]))
  vpd_obs_n2 <- c(vpd_obs_n2, vpd_obs_tree_n2)
  
  vpd_obs_tree_n3 <- sapply(years_tree, 
                         function(y) 
                           as.numeric(clim_pred$vpd_mjj[clim_pred$year == as.character(as.numeric(y)-3)][1]))
  vpd_obs_n3 <- c(vpd_obs_n3, vpd_obs_tree_n3)
  
  years <- c(years, years_tree)
  all_years_idxs <- c(all_years_idxs, all_years_idxs_tree)
  N_years <- c(N_years, N_years_tree)
  
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
sum(is.na(gdd_obs)) # check clim. pred

# Collection data into list
N <- length(years)
data <- mget(c('N', 'N_all_years', 'N_trees', 
               'log_rw_obs', 'gdd_obs', 'sm_obs', 
               'vpd_obs', 'vpd_obs_n1', 'vpd_obs_n2', 'vpd_obs_n3',
               'all_years', 'years', 'all_years_idxs', 'N_years', 
               'tree_start_idxs', 'tree_end_idxs'
))
data$years <- as.numeric(data$years)

# Posterior quantification
fit <- stan(file=file.path(wd, 'model/stan/model6_with3predictors_pnw_difflength_onespecies_legacies.stan'),
            data=data, seed=5838299, 
            chains = 4, cores = 4,
            warmup=1000, iter=2024, refresh=10)
