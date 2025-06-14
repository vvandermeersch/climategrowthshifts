
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

# Gather some stands (TEMPORARY, we should find a better method for this)
stands_gather <-
  data.frame(stand = datasets$dataset,
             standgathered = c("or_092", "or_093", "or_094", "or_095", rep("or_098", 7), rep('or_105', 3), rep('or_093', 3),
                               "or_107", "or_116", "or_129", "or_130", "or_131"))
raw_data <- merge(raw_data, stands_gather)


# Load climate data
clim_pred <- readRDS(file.path(wd, "output", "climate",  "climpredictors_pnw_june2025.rds"))
clim_pred$soilmoist_jja <- clim_pred$soilmoist_jja*100 # in percentage
clim_pred$gdd_ings <- clim_pred$gdd_ings/10 # in 10degC

# Sizes
uniq_tree_ids <- unique(raw_data$tree_id)
N_trees <- length(uniq_tree_ids)

uniq_stand_ids <- unique(raw_data$standgathered)
N_stands <- length(unique(raw_data$standgathered))
# altitude <- datasets[,c('dataset', 'altitude')]
# names(altitude)[1] <- 'stand'
# raw_data <- merge(raw_data, altitude)
# raw_data$region <- ifelse(raw_data$altitude < 1000, "A", "B")
# uniq_stand_ids <- unique(raw_data$region)
# N_stands <- length(unique(raw_data$region))

uniq_species_ids <- unique(raw_data$species) #  Important that JUOC is first here!
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
  
  stand_tree <- which(uniq_stand_ids == raw_data_tree$standgathered[1])
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

# Collection data into list
N <- length(years)

data <- mget(c('N', 'N_all_years', 'N_trees', 
               'log_rw_obs', 'gdd_obs', 'sm_obs',
               'all_years', 'years', 'all_years_idxs', 'N_years', 
               'stand_idxs', 'N_stands',
               'species_idxs', 'N_species',
               'tree_start_idxs', 'tree_end_idxs'))

# Posterior Quantification
# fit <- stan(file=file.path(wd, 'model/stan/model6_with2predictors_pnw_samefsh.stan'),
#             data=data, seed=5838299, cores = 4,
#             warmup=1000, iter=2024, refresh=10)
fit <- stan(file=file.path(wd, 'model/stan/model6_with2predictors_pnw_differentbetas.stan'),
            data=data, seed=5838299, cores = 4,
            warmup=1000, iter=2024, refresh=10)

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
for (t in  sample(which(data$stand_idxs == which(uniq_stand_ids == "or_106")), 6)) {
  idxs <- tree_start_idxs[t]:tree_end_idxs[t]
  rw_names <- sapply(idxs,
                     function(n) paste0('log_rw_pred[', n, ']'))
  util$plot_conn_pushforward_quantiles(samples, rw_names, data$years[idxs],
                                       xlab="Year", ylab="Log Ring Width Per mm", 
                                       display_ylim=c(-4, 2),
                                       main=paste0("Stand ", uniq_stand_ids[stand_idxs[t]], 
                                                   ", Tree ", uniq_tree_ids[t]))
  abline(v=1992, lwd=1, lty=2, col='darkgrey')
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=1.0, col="white")
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=0.8, col="black")
  
}

# Inference
par(mfrow=c(3, 2))
for (t in  1 * (1:6)) {
  idxs <- tree_start_idxs[t]:tree_end_idxs[t]
  rw_names <- sapply(idxs,
                     function(n) paste0('mu1[', n, ']'))
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

par(mfrow=c(1, 1), mar = c(4,4,1,1))
names <- sapply(1:data$N_species,
                function(sp) paste0('kappa_sh[', sp, ']'))
util$plot_disc_pushforward_quantiles(samples, names,
                                     xlab="Species",
                                     xticklab=uniq_species_ids,
                                     ylab="Climate Response Scaling")

par(mfrow=c(1, 1), mar = c(4,4,1,1))
names <- sapply(1:data$N_species,
                function(sp) paste0('beta_sm[', sp, ']'))
util$plot_disc_pushforward_quantiles(samples, names,
                                     xlab="Species",
                                     xticklab=uniq_species_ids,
                                     ylab="Climate Response Scaling")
