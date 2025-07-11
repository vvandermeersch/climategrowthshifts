
rm(list = ls())
wd <- '/home/victor/projects/climategrowthshifts/analysis/pnw'
library(dplR)
library(rstan)
rstan_options(auto_write = TRUE)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)

# Read RWL file
years <- years <- c(1983:2025)
raw_data1 <- dplR::read.rwl(file.path(wd, 'data', 'others', 'soaproot_saddle', 'PonderosaMaster.rwl'))
idrows <- which(rownames(raw_data1) %in% as.character(years))
raw_data1 <- raw_data1[idrows,]
raw_data1$year <- rownames(raw_data1)
rwdat_long <- reshape(raw_data1, varying = names(raw_data1)[names(raw_data1) != "year"],
                      v.names = "rw_mm", timevar = "core_id",
                      times = names(raw_data1)[names(raw_data1) != "year"],
                      direction = "long")
rwdat_long$tree_id <- sub("^(.*\\d)[^\\d]*$", "\\1", rwdat_long$core_id)
rwdat_long <- rwdat_long[, c("year", 'tree_id', "core_id", "rw_mm")]
raw_data1 <- aggregate(
  rw_mm ~ tree_id + year,
  data = rwdat_long,
  FUN = function(x) mean(x, na.rm = TRUE),
  na.action = na.pass
)
names(raw_data1)[names(raw_data1) == "rw_mm"] <- "rw_avg_mm"
raw_data1 <- na.omit(raw_data1)

raw_data2 <- dplR::read.rwl(file.path(wd, 'data', 'others', 'soaproot_saddle', 'PonderosaMaster_full.rwl'))
idrows <- which(rownames(raw_data2) %in% as.character(years))
raw_data2 <- raw_data2[idrows,]
raw_data2$year <- rownames(raw_data2)
rwdat_long <- reshape(raw_data2, varying = names(raw_data2)[names(raw_data2) != "year"],
                      v.names = "rw_mm", timevar = "core_id",
                      times = names(raw_data2)[names(raw_data2) != "year"],
                      direction = "long")
rwdat_long$tree_id <- sub("^(.*\\d)[^\\d]*$", "\\1", rwdat_long$core_id)
rwdat_long <- rwdat_long[, c("year", 'tree_id', "core_id", "rw_mm")]
raw_data2 <- aggregate(
  rw_mm ~ tree_id + year,
  data = rwdat_long,
  FUN = function(x) mean(x, na.rm = TRUE),
  na.action = na.pass
)
names(raw_data2)[names(raw_data2) == "rw_mm"] <- "rw_avg_mm"
raw_data2 <- na.omit(raw_data2)

raw_data <- rbind(raw_data1, raw_data2)

treecounts <- aggregate(year ~ tree_id, data = raw_data, FUN = length)
treestoremove <- c(treecounts[treecounts$year < 20, 'tree_id'], 'SS04-16', 'SRA13')
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

diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)
samples <- util$extract_expectand_vals(fit)
base_samples <- util$filter_expectands(samples,
                                       c('alpha', 
                                         'beta_gdd', 'beta_sm', 
                                         'beta_vpd', 'rho_beta_vpd',
                                         'rho', 'gamma',
                                         'f_tilde_sh',
                                         'rho_sh', 'kappa',
                                         'sigma'),
                                       check_arrays=TRUE)
util$check_all_expectand_diagnostics(base_samples)

util$plot_pairs_by_chain(samples[[paste0('rho')]], paste0('rho'), samples[[paste0('gamma')]], paste0('gamma'))
util$plot_pairs_by_chain(samples[[paste0('rho_sh')]], paste0('rho_sh'), samples[[paste0('kappa')]], paste0('kappa'))
util$plot_pairs_by_chain(log(samples[[paste0('rho_beta_vpd')]]), paste0('log(rho_beta_vpd)'), samples[[paste0('beta_vpd')]], paste0('beta_vpd'))
util$plot_pairs_by_chain(samples[[paste0('beta_gdd')]], paste0('beta_gdd'), samples[[paste0('beta_sm')]], paste0('beta_sm'))

par(mfrow=c(1, 1), mar = c(4,4,1,1))
names <- sapply(1:4,
                function(sp) paste0('beta_vpd_leg[', sp, ']'))
util$plot_disc_pushforward_quantiles(samples, names,
                                     xlab="Lag",
                                     xticklabs=c(0:3),
                                     ylab="VPD effect")
