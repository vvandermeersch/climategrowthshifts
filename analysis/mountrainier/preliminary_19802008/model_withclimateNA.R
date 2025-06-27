### Created by Victor ###
### Edited by Mao on 25 June 2025 ###

rm(list = ls())
if(length(grep("victor", getwd()) > 0)) {
  setwd("~/projects/climategrowthshifts/analysis/mountrainier")
} else if(length(grep("Xiaomao", getwd()) > 0)) {
  setwd("C:/PhD/Project/climategrowthshifts/analysis/mountrainier")
} else if(length(grep("lizzie", getwd()) > 0)) {
  setwd("~/Documents/git/projects/grephon/climategrowthshifts/analysis/mountrainier")
}

util <- new.env()
runmodel <- FALSE
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)

library(rstan)

# Data exploration
# raw_data <- read.csv(file.path(wd, 'input/treerings/mosttrees_allspecies_allplots_19802008.csv'))
raw_data <- read.csv(file.path('input/treerings/mosttrees_allspecies_allplots_19802008.csv'))
raw_data <- raw_data[raw_data$species != 'Xano',]

clim_pred <- read.csv(file.path('data/climate/climateNAMORAStands.csv'))

# Average ring widths when both cores are present
robust_ave <- function(x, y) {
  ifelse(is.na(y), x, mean(c(x, y)))
  if (is.na(x)) {
    y
  } else if(is.na(y)) {
    x
  } else {
    mean(c(x, y))
  }
}

raw_data$rw_core_ave..mm. <- as.vector(by(raw_data, 
                                          seq_len(nrow(raw_data)), 
                                          function(r) robust_ave(r[[6]], r[[7]])))


# Reduce data size (temporary)
#raw_data <- raw_data[raw_data$tree_id %in% sample(raw_data$tree_id, 40),]

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
ffp_obs <- c()
gdd_obs <- c()
pas_obs <- c()
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
  
  log_rw_obs_tree <- sapply(years_tree, 
                            function(y) 
                              log(raw_data_tree$rw_core_ave..mm.[raw_data_tree$year == y][1]))
  log_rw_obs <- c(log_rw_obs, log_rw_obs_tree)
  
  ffp_obs_tree <- sapply(years_tree, 
                         function(y) 
                           as.numeric(clim_pred$FFP[clim_pred$year == y & clim_pred$stand_id == raw_data_tree$stand[1]][1]))
  ffp_obs <- c(ffp_obs, ffp_obs_tree)
  
  gdd_obs_tree <- sapply(years_tree, 
                         function(y) 
                           as.numeric(clim_pred$DD5[clim_pred$year == y & clim_pred$stand_id == raw_data_tree$stand[1]][1]))
  gdd_obs <- c(gdd_obs, gdd_obs_tree)
  
  pas_obs_tree <- sapply(years_tree, 
                         function(y) 
                           as.numeric(clim_pred$PAS[clim_pred$year == y & clim_pred$stand_id == raw_data_tree$stand[1]][1]))
  pas_obs <- c(pas_obs, pas_obs_tree)
  
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
length(ffp_obs)
length(gdd_obs)
length(pas_obs)
length(years)
length(all_years_idxs)
length(N_years)

length(tree_start_idxs)
length(tree_end_idxs)

# Plot all trees using new data format
# par(mfrow=c(5, 2))
# 
# for (t in 50*(1:10)) {
#   idxs <- tree_start_idxs[t]:tree_end_idxs[t]
#   
#   year <- years[idxs]
#   rw_ave <- exp(log_rw_obs[idxs])
#   
#   plot(year, rw_ave, pch=16, cex=1.0, 
#        xlab="Year", ylab="Ring Width Per mm", ylim=c(0, 10),
#        main=paste("Stand", uniq_stand_ids[stand_idxs[t]],
#                   ", Species", uniq_species_ids[species_idxs[t]],
#                   ", Tree", uniq_tree_ids[t]))
# }

# Cross check
par(mfrow=c(2, 1))

t <- 10

id <- uniq_tree_ids[t]
year <- raw_data$year[raw_data$tree_id == id]
rw_ave <- raw_data$rw_core_ave..mm.[raw_data$tree_id == id]

plot(year[year %in% all_years], rw_ave[year %in% all_years], pch=16, cex=1.0, 
     xlab="Year", ylab="Ring Width (mm)", ylim=c(0, 3.5),
     main=paste("Tree", uniq_tree_ids[t]))

idxs <- tree_start_idxs[t]:tree_end_idxs[t]
year <- years[idxs]
rw_ave <- exp(log_rw_obs[idxs])

plot(year, rw_ave, pch=16, cex=1.0, 
     xlab="Year", ylab="Ring Width Per mm", ylim=c(0, 3.5),
     main=paste("Tree", uniq_tree_ids[t]))

par(mfrow=c(2, 1))
stand_tree <- raw_data$stand[raw_data$tree_id == id][1]

ffp_ave <- as.numeric(clim_pred$FFP[clim_pred$stand_id == stand_tree & clim_pred$year %in% all_years])

plot(year[year %in% all_years], ffp_ave, pch=16, cex=1.0, 
     xlab="Year", ylab="Frost Free Period (day)", ylim=c(0, 250),
     main=paste("Tree", uniq_tree_ids[t]))

idxs <- tree_start_idxs[t]:tree_end_idxs[t]
year <- years[idxs]
ffp_ave <- ffp_obs[idxs]

plot(year, ffp_ave, pch=16, cex=1.0, 
     xlab="Year", ylab="Frost Free Period (day)", ylim=c(0, 250),
     main=paste("Tree", uniq_tree_ids[t]))


# Get the phylogenetic matrix (path to modify, I don't know how to use setwd)
print(uniq_species_ids)
sps.list <- c('Abies_amabilis', 'Tsuga_mertensiana', 'Pseudotsuga_menziesii', 'Thuja_plicata',  'Tsuga_heterophylla')
print(sps.list)
phy.plants <- read.tree(file.path('~/projects/climategrowthshifts/analysis/pnw', "input/ALLMB.tre"))
phy.plants.here <-  drop.tip(phy.plants,
                             which(!phy.plants$tip.label %in% sps.list))
Cphy <- ape::vcv.phylo(phy.plants.here,corr=TRUE)

# Collection data into list
N <- length(years)

data <- mget(c('N', 'N_all_years', 'N_trees', 
               'log_rw_obs', 'ffp_obs', 'gdd_obs', 'pas_obs',
               'all_years', 'years', 'all_years_idxs', 'N_years', 
               'stand_idxs', 'N_stands',
               'species_idxs', 'N_species',
               'tree_start_idxs', 'tree_end_idxs', 'Cphy'))

# Posterior Quantification
if(runmodel){
  fit <- stan(file=file.path('stan/model6_with3predictorsClimateNAspp.stan'),
              data=data, seed=5838299, cores = 4,
              warmup=1000, iter=2024, refresh=10)
  saveRDS(fit, "output/model6_with3predictorsClimateNAspp.rds")
}

if(!runmodel){
  fit <- readRDS("output/model6_with3predictorsClimateNA.rds")
  samples <- readRDS("output/model6_with3predictorsClimateNAsamples.rds")
}

diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit)
base_samples <- util$filter_expectands(samples,
                                       c('alpha', 
                                         'rho', 'gamma',
                                         'f_tilde_sh',
                                         'rho_sh', 'gamma_sh',
                                         'kappa_sh_free',
                                         'sigma'),
                                       check_arrays=TRUE)
util$check_all_expectand_diagnostics(base_samples)

saveRDS(samples, "output/model6_with3predictorsClimateNAsamples.rds")

# Retrodictive check

## Plot 4 random trees
par(mfrow=c(2, 2))
for (t in  sample(1:N_trees, 4)) {
  idxs <- tree_start_idxs[t]:tree_end_idxs[t]
  rw_names <- sapply(idxs,
                     function(n) paste0('log_rw_pred[', n, ']'))
  util$plot_conn_pushforward_quantiles(samples, rw_names, data$years[idxs],
                                       xlab="Year", ylab="Log Ring Width Per mm", 
                                       display_ylim=c(-4, 4),
                                       main=paste0("Stand ", uniq_stand_ids[stand_idxs[t]], 
                                                   ", species ", uniq_species_ids[species_idxs[t]]))
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=1.0, col="white")
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=0.8, col="black")
  #abline(v=1974, lwd=2, lty=2, col="#DDDDDD")
}

## Look at GDD per stands
par(mfrow=c(2, 2), mar = c(4,4.5,2.5,1))
for(s in 1:N_stands){
  start <- data$tree_start_idxs[min(which(data$stand_idxs == s))]
  end <- data$tree_end_idxs[max(which(data$stand_idxs == s))]
  pred_names <- sapply(start:end, function(n) paste0('log_rw_pred[', n, ']'))
  util$plot_conditional_median_quantiles(samples, pred_names, data$gdd_obs[start:end],
                                         0, 32, 1, data$log_rw_obs[start:end], 
                                         residual=FALSE,
                                         ylab = 'Marginal quantiles',
                                         xlab="GDD (degC)", main = paste0('Species ', uniq_species_ids[spind]))
}
## Idem,but residuals
for(s in 1:N_stands){
  start <- data$tree_start_idxs[min(which(data$stand_idxs == s))]
  end <- data$tree_end_idxs[max(which(data$stand_idxs == s))]
  pred_names <- sapply(start:end, function(n) paste0('log_rw_pred[', n, ']'))
  util$plot_conditional_median_quantiles(samples, pred_names, data$gdd_obs[start:end],
                                         0, 32, 1, data$log_rw_obs[start:end], 
                                         residual=TRUE,
                                         ylab = 'Marginal quantiles',
                                         xlab="GDD (degC)", main = paste0('Species ', uniq_species_ids[spind]))
}

par(mfrow=c(1, 3))
util$plot_expectand_pushforward(samples[['beta_ffp']], 25,
                                display_name="beta_ffp")
xs <- seq(0, 1, 0.001)
ys <- dnorm(xs, 0, log(1.8) / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['beta_gdd']], 25,
                                display_name="beta_gdd")
xs <- seq(0, 1, 0.01)
ys <- dnorm(xs, 0, log(1.8) / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['beta_pas']], 25,
                                display_name="beta_pas")
xs <- seq(0, 1, 0.01)
ys <- dnorm(xs, 0, log(1.8) / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)


par(mfrow=c(3, 1))

pred_names <- sapply(1:data$N, function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conditional_median_quantiles(samples, pred_names, data$ffp_obs,
                                       60, 200, 5, data$log_rw_obs, 
                                       xlab="FFP (days)")

util$plot_conditional_median_quantiles(samples, pred_names, data$gdd_obs,
                                       520, 2200, 50, data$log_rw_obs, 
                                       xlab="GDD (degC)")

util$plot_conditional_median_quantiles(samples, pred_names, data$pas_obs,
                                       50, 2100, 50, data$log_rw_obs, 
                                       xlab="PAS (mm)")
