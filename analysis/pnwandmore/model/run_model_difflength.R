rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnw"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)

library(rstan)
library(dplyr)

# Load treering data
datasets <- readRDS(file.path(wd, 'input', 'itrdb', 'datasets_summary_all.rds'))

# Keep only datasets within WLDAS extent
datasets <- datasets[datasets$north_lat >=25.065 & datasets$north_lat <=52.925 & datasets$east_lon <= -89.025 &  datasets$east_lon >= -124.925,]

# datasets <- datasets[datasets$last_year >= 1999,] # at least 20 years of observations
datasets <- datasets[datasets$dataset != 'wa149',] # temporary 
datasets <- datasets[datasets$dataset != 'ca673',] # temporary
datasets <- datasets[datasets$dataset != 'az563',] # temporary

# Temporary, dropping Angiosperms
todrop <- c('QUGA', 'QULO', 'PPTR', 'PPFR', 'PLRA', 'QUDG')
datasets <- datasets[!(datasets$species_code %in% todrop),]

# Same species
datasets[datasets$species_code == 'ABBI', c('species_name', 'species_code')] <- 
  unique(datasets[datasets$species_code == 'ABLA', c('species_name', 'species_code')])

source(file.path(wd, 'getphylo.R'))

ringwidth_series <- readRDS(file.path(wd, 'input', 'itrdb', 'ringwidth_series_all.rds'))
raw_data <- data.frame()
for(d in 1:nrow(datasets)){
  
  raw_data_d <- ringwidth_series[ringwidth_series$dataset == datasets[d, 'dataset'], ]
  
  raw_data_d <- raw_data_d[raw_data_d$year >= 1980 & raw_data_d$year < 2024 & !is.na(raw_data_d$year), ]
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
clim_pred <- readRDS(file.path(wd, "output", "climate",  "climpredictors_ext_pnw_02july2025.rds"))
clim_pred$soilmoist_mjj <- clim_pred$soilmoist_mjj*100 # in percentage
clim_pred$gdd <- clim_pred$gdd/100 # in x100 degC
clim_pred$vpd_mjj <- clim_pred$vpd_mjj # in hPa?

# Sizes
uniq_tree_ids <- unique(raw_data$tree_id_uniq)
N_trees <- length(uniq_tree_ids)

uniq_stand_ids <- unique(raw_data$grouped_stand)
N_stands <- length(unique(raw_data$grouped_stand))

uniq_species_ids <- unique(raw_data$species_code) #  Important that JUOC is first here!
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
                              log(raw_data_tree$rw_avg_mm[raw_data_tree$year == y][1]))
  log_rw_obs <- c(log_rw_obs, log_rw_obs_tree)
  
  
  
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
length(sm_obs)
length(years)
length(all_years_idxs)
length(N_years)
length(tree_start_idxs)
length(tree_end_idxs)
sum(is.na(gdd_obs)) # check clim. pred

# Phylogenetic matrix
Cphy <- ape::vcv.phylo(phy.plants.here,corr=TRUE)

# Collection data into list
N <- length(years)

data <- mget(c('N', 'N_all_years', 'N_trees', 
               'log_rw_obs', 'gdd_obs', 'sm_obs', 'vpd_obs',
               'all_years', 'years', 'all_years_idxs', 'N_years', 
               'stand_idxs', 'N_stands',
               'species_idxs', 'N_species',
               'tree_start_idxs', 'tree_end_idxs'))
data$years <- as.numeric(data$years)
saveRDS(data, file = file.path(wd, 'output/model', 'data_02july2025.rds'))

# Posterior Quantification
# fit <- stan(file=file.path(wd, 'model/stan/model6_with2predictors_pnw_samefsh.stan'),
#             data=data, seed=5838299, cores = 4,
#             warmup=1000, iter=2024, refresh=10)
fit <- stan(file=file.path(wd, 'model/stan/model6_with3predictors_pnw_fullphylo_difflength_test.stan'),
            data=data, seed=5838299, 
            chains = 4, cores = 4,
            warmup=10, iter=20, refresh=10)
saveRDS(fit, file = file.path(wd, 'output/model', 'fit_20june2025.rds'))

diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit)
base_samples <- util$filter_expectands(samples,
                                       c('alpha', 'beta_gdd', 'beta_sm', 
                                         'rho', 'gamma',
                                         'f_tilde_sh',
                                         'rho_sh', 'gamma_sh',
                                         'kappa_sh_free',
                                         'sigma'),
                                       check_arrays=TRUE)
util$check_all_expectand_diagnostics(base_samples)

# Retrodictive Check
par(mfrow=c(3, 2))
for (t in  50*1:6) {
  idxs <- tree_start_idxs[t]:tree_end_idxs[t]
  rw_names <- sapply(idxs,
                     function(n) paste0('log_rw_pred[', n, ']'))
  util$plot_conn_pushforward_quantiles(samples, rw_names, data$years[idxs],
                                       xlab="Year", ylab="Log Ring Width Per mm", 
                                       display_ylim=c(-4, 2),
                                       main=paste0("Stand ", uniq_stand_ids[stand_idxs[t]], 
                                                   ", Tree ", uniq_tree_ids[t]))
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=1.0, col="white")
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=0.8, col="black")
  #abline(v=1974, lwd=2, lty=2, col="#DDDDDD")
}

par(mfrow=c(3, 2))
for (t in  sample(which(data$stand_idxs == which(uniq_stand_ids == "or_093")), 6)) {
  idxs <- tree_start_idxs[t]:tree_end_idxs[t]
  rw_names <- sapply(idxs,
                     function(n) paste0('log_rw_pred[', n, ']'))
  util$plot_conn_pushforward_quantiles(samples, rw_names, data$years[idxs],
                                       xlab="Year", ylab="Log Ring Width Per mm", 
                                       display_ylim=c(-6, 2),
                                       main=paste0("Stand ", uniq_stand_ids[stand_idxs[t]], 
                                                   ", Tree ", uniq_tree_ids[t]))
  abline(v=1992, lwd=1, lty=2, col='darkgrey')
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=1.0, col="white")
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=0.8, col="black")
  
}

# Inference
par(mfrow=c(3, 2))
for (t in  10 * (1:6)) {
  idxs <- tree_start_idxs[t]:tree_end_idxs[t]
  rw_names <- sapply(idxs,
                     function(n) paste0('mu2[', n, ']'))
  util$plot_conn_pushforward_quantiles(samples, rw_names, data$years[idxs],
                                       xlab="Year", ylab="Log Ring Width Per mm", 
                                       display_ylim=c(-4, 2),
                                       main=paste0("Stand ", uniq_stand_ids[stand_idxs[t]], 
                                                   ", Tree ", uniq_tree_ids[t]))
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=1.0, col="white")
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=0.8, col="black")
  #abline(v=1974, lwd=2, lty=2, col="#DDDDDD")
}

par(mfrow=c(1, 2))
util$plot_expectand_pushforward(samples[['beta_gdd']], 250,
                                flim = c(0,0.01),
                                display_name="beta_gdd")
xs <- seq(0, 1, 0.01)
ys <- dnorm(xs, 0, log(1.8) / 2.57)*50
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['beta_sm']], 250,
                                flim = c(0,0.1),
                                display_name="beta_sm")
xs <- seq(0, 1, 0.01)
ys <- dnorm(xs, 0, log(1.8) / 2.57) *10
lines(xs, ys, lwd=2, col=util$c_light)


par(mfrow=c(1, 1), mar = c(4,4,1,1))
util$plot_expectand_pushforward(samples[[paste0('rho_sh')]], 100,
                                display_name="rho_sh",
                                flim = c(0,5))

par(mfrow=c(1, 1), mar = c(4,4,1,1))
names <- sapply(1:data$N_species,
                function(sp) paste0('kappa_sh[', sp, ']'))
util$plot_disc_pushforward_quantiles(samples, names,
                                     xlab="Species",
                                     xticklab=uniq_species_ids,
                                     ylab="Climate Response Scaling")



par(mfrow=c(2, 2), mar = c(4,4,1,1))
for (s in 1:data$N_species) {
  util$plot_expectand_pushforward(samples[[ paste0('rho[', s, ']')]], 250,
                                  flim = c(0,50),
                                  display_name= paste0("rho ", uniq_species_ids[s]))
  xs <- seq(0, 60, 0.01)
  ys <- dlnorm(xs, 3.55, 0.24)
  lines(xs, ys, lwd=2, col=util$c_light)
}

par(mfrow=c(2, 2))
for (s in 1:data$N_species) {
  util$plot_expectand_pushforward(samples[[ paste0('kappa_sh[', s, ']')]], 250,
                                  flim = c(0,50),
                                  display_name= paste0("kappa_sh_free ", uniq_species_ids[s]))
  xs <- seq(0, 60, 0.01)
  ys <- dlnorm(xs, 3.55, 0.24)
  lines(xs, ys, lwd=2, col=util$c_light)
}

par(mfrow=c(1, 1), mar = c(4,4.5,1,1))
names <- sapply(1:data$N_species,
                function(sp) paste0('kappa_sh[', sp, ']'))
util$plot_disc_pushforward_quantiles(samples, names,
                                     xlab="Species",
                                     xticklab=uniq_species_ids,
                                     ylab=expression(kappa[short]))

par(mfrow=c(1, 1), mar = c(4,4.5,1,1))
names <- sapply(1:data$N_species,
                function(sp) paste0('beta_sm[', sp, ']'))
util$plot_disc_pushforward_quantiles(samples, names,
                                     xlab="Species",
                                     xticklab=uniq_species_ids,
                                     ylab=expression(beta[SMoist]))

par(mfrow=c(1, 1), mar = c(4,4.5,1,1))
names <- sapply(1:data$N_species,
                function(sp) paste0('beta_gdd[', sp, ']'))
util$plot_disc_pushforward_quantiles(samples, names,
                                     xlab="Species",
                                     xticklab=uniq_species_ids,
                                     ylab=expression(beta[GDD]))

par(mfrow=c(1, 1), mar = c(4,4.5,1,1))
names <- sapply(1:data$N_species,
                function(sp) paste0('beta_vpd[', sp, ']'))
util$plot_disc_pushforward_quantiles(samples, names,
                                     xlab="Species",
                                     xticklab=uniq_species_ids,
                                     ylab=expression(beta[VPD]))

spind <- 4
start <- data$tree_start_idxs[min(which(data$species_idxs == spind))]
end <- data$tree_end_idxs[max(which(data$species_idxs == spind))]
pred_names <- sapply(start:end, function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conditional_median_quantiles(samples, pred_names, data$gdd_obs[start:end],
                                       0, 32, 1, data$log_rw_obs[start:end], 
                                       residual=FALSE,
                                       xlab="GDD (x100degC)", main = paste0('Species ', uniq_species_ids[spind]))

pred_names <- sapply(start:end, function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conditional_median_quantiles(samples, pred_names, data$sm_obs[start:end],
                                       15, 36, 1, data$log_rw_obs[start:end], 
                                       residual=FALSE,
                                       xlab="Soil moisture (%)", main = paste0('Species ', uniq_species_ids[spind]))

pred_names <- sapply(start:end, function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conditional_median_quantiles(samples, pred_names, data$vpd_obs[start:end],
                                       0, 20, 1, data$log_rw_obs[start:end], 
                                       residual=FALSE,
                                       xlab="VPD (hPa)", main = paste0('Species ', uniq_species_ids[spind]))


plot(c(samples[["beta_gdd[1]"]]) ~ c(samples[["beta_vpd[1]"]]),
     xlab = expression(beta[VPD]), ylab = expression(beta[GDD]))

plot(c(samples[["beta_sm[1]"]]) ~ c(samples[["beta_vpd[1]"]]),
     xlab = expression(beta[VPD]), ylab = expression(beta[SMoist]))

plot(c(samples[["beta_gdd[1]"]]) ~ c(samples[["beta_sm[1]"]]),
     xlab = expression(beta[SMoist]), ylab = expression(beta[GDD]))
