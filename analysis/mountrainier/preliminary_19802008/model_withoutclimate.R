
rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/mountrainier"
setwd(wd)
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)

library(rstan)

# Data exploration
# raw_data <- read.csv(file.path(wd, 'input/treerings/mosttrees_allspecies_allplots_19802008.csv'))
raw_data <- read.csv(file.path(wd, 'input/treerings/mosttrees_allspecies_allplots_19192008.csv'))

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
length(years)
length(all_years_idxs)
length(N_years)

length(tree_start_idxs)
length(tree_end_idxs)

# Plot all trees using new data format
par(mfrow=c(5, 2))

for (t in 50*(1:10)) {
  idxs <- tree_start_idxs[t]:tree_end_idxs[t]
  
  year <- years[idxs]
  rw_ave <- exp(log_rw_obs[idxs])
  
  plot(year, rw_ave, pch=16, cex=1.0, 
       xlab="Year", ylab="Ring Width Per mm", ylim=c(0, 10),
       main=paste("Stand", uniq_stand_ids[stand_idxs[t]],
                  ", Species", uniq_species_ids[species_idxs[t]],
                  ", Tree", uniq_tree_ids[t]))
}

# Cross check
par(mfrow=c(2, 1))

t <- 12

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

# Collection data into list
N <- length(years)

data <- mget(c('N', 'N_all_years', 'N_trees', 
               'log_rw_obs', 
               'all_years', 'years', 'all_years_idxs', 'N_years', 
               'stand_idxs', 'N_stands',
               'species_idxs', 'N_species',
               'tree_start_idxs', 'tree_end_idxs'))

# Posterior Quantification
fit <- stan(file=file.path(wd, 'stan/model6.stan'),
            data=data, seed=5838299,
            warmup=1000, iter=2024, refresh=10)


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

# Retrodictive Check
par(mfrow=c(3, 2))

for (t in  10 * (1:6)) {
  idxs <- tree_start_idxs[t]:tree_end_idxs[t]
  rw_names <- sapply(idxs,
                     function(n) paste0('log_rw_pred[', n, ']'))
  util$plot_conn_pushforward_quantiles(samples, rw_names, data$years[idxs],
                                       xlab="Year", ylab="Log Ring Width Per mm", 
                                       display_ylim=c(-2, 2),
                                       main=paste("Stand", uniq_stand_ids[stand_idxs[t]], 
                                                  ", Tree", uniq_tree_ids[t]))
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=1.0, col="white")
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=0.8, col="black")
  #abline(v=1974, lwd=2, lty=2, col="#DDDDDD")
}

# Posterior inference
par(mfrow=c(3, 2))
for(i in 1:6){
  util$plot_expectand_pushforward(samples[[paste0('rho[',i,']')]], 250,
                                  display_name="rho", main = uniq_species_ids[i],
                                  flim = c(0,50))
  xs <- seq(0, 60, 0.01)
  ys <- dlnorm(xs, 3.55, 0.24)
  lines(xs, ys, lwd=2, col=util$c_light)
}

par(mfrow=c(1, 1))
util$plot_expectand_pushforward(samples[[paste0('rho_sh')]], 250,
                                  display_name="rho_sh",
                                  flim = c(0,10))
xs <- seq(0, 10, 0.01)
ys <- dlnorm(xs, 1.7, 0.26)
lines(xs, ys, lwd=2, col=util$c_light)


