
rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/modeling"
setwd(wd)
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)

library(rstan)
library(dplyr)

# Data exploration
raw_data <- read.csv(file.path(wd, 'input/ut_550/cleaned_tree_data_with_flags.csv'))
raw_data <- raw_data[raw_data$age_CE >= 1980, ]
coltokeep <-  !grepl(colnames(raw_data), pattern = 'outlier')
raw_data <- raw_data[, coltokeep]

clim_pred <- readRDS(file.path(wd, "output", "climate",  "climpredictors_ut550.rds"))
clim_pred$soilmoist_jja <- clim_pred$soilmoist_jja*100 # in percentage
clim_pred$gdd_ings <- clim_pred$gdd_ings/10 # in 10degC

# Check where to stop
for(endyear in seq(2010, 2014,1)){
  cat(paste0(endyear, '\n'))
  temp <- raw_data[raw_data$age_CE <= endyear, ]
  temp <- temp[ , colSums(is.na(temp)) == 0]
  cat(paste0(ncol(temp)-1, '\n'))
}
raw_data <- raw_data[raw_data$age_CE <= 2012, ]
raw_data <- raw_data[ , colSums(is.na(raw_data)) < length(1980:2012)]

raw_data_long <- raw_data %>%
  tidyr::pivot_longer(cols = starts_with("RCB"), names_to = 'core_id', values_to ='rw_core')
names(raw_data_long)[1] <- 'year'
raw_data_long$stand <- 1

raw_data_long$tree_id <- stringr::str_split_i(raw_data_long$core_id, pattern = '_', i = 1)
raw_data_long$tree_id <- substr(raw_data_long$tree_id,1,nchar(raw_data_long$tree_id)-1)

raw_data <- raw_data_long %>%
  group_by(stand, tree_id, year) %>%
  summarise(rw_core_ave = mean(rw_core, na.rm = TRUE))

trees_to_remove <- unique(raw_data[is.na(raw_data$rw_core_ave) | raw_data$rw_core_ave == 0 , 'tree_id'])
raw_data <- raw_data[!(raw_data$tree_id %in% trees_to_remove$tree_id),]


# Sizes
uniq_tree_ids <- unique(raw_data$tree_id)
N_trees <- length(uniq_tree_ids)

uniq_stand_ids <- unique(raw_data$stand)
N_stands <- length(unique(raw_data$stand))

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
                              log(raw_data_tree$rw_core_ave[raw_data_tree$year == y][1]))
  log_rw_obs <- c(log_rw_obs, log_rw_obs_tree)
  
  gdd_obs_tree <- sapply(years_tree, 
                         function(y) 
                           as.numeric(clim_pred$gdd_ings[clim_pred$year == y & clim_pred$ID == raw_data_tree$stand[1]][1]))
  gdd_obs <- c(gdd_obs, gdd_obs_tree)
  
  sm_obs_tree <- sapply(years_tree, 
                         function(y) 
                           as.numeric(clim_pred$soilmoist_jja[clim_pred$year == y & clim_pred$ID == raw_data_tree$stand[1]][1]))
  sm_obs <- c(sm_obs, sm_obs_tree)
  
  years <- c(years, years_tree)
  all_years_idxs <- c(all_years_idxs, all_years_idxs_tree)
  N_years <- c(N_years, N_years_tree)
  
  stand_tree <- which(uniq_stand_ids == raw_data_tree$stand[1])
  stand_idxs <- c(stand_idxs, stand_tree)
  
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

# Cross check
par(mfrow=c(2, 1))

t <- 10

id <- uniq_tree_ids[t]
year <- raw_data$year[raw_data$tree_id == id]
rw_ave <- raw_data$rw_core_ave[raw_data$tree_id == id]

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
               'log_rw_obs', 'gdd_obs', 'sm_obs',
               'all_years', 'years', 'all_years_idxs', 'N_years', 
               'stand_idxs', 'N_stands',
               'tree_start_idxs', 'tree_end_idxs'))

# Posterior Quantification
fit <- stan(file=file.path(wd, 'stan/model6_with2predictors_bristlecone.stan'),
            data=data, seed=5838299, cores = 4,
            warmup=1000, iter=2024, refresh=10)



diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit)
base_samples <- util$filter_expectands(samples,
                                       c('alpha', 
                                         'rho', 'gamma',
                                         'f_tilde_sh',
                                         'rho_sh', 'gamma_sh',
                                         'sigma'),
                                       check_arrays=TRUE)
util$check_all_expectand_diagnostics(base_samples)

# Retrodictive Check
par(mfrow=c(3, 2))

for (t in  1 * (1:6)) {
  idxs <- tree_start_idxs[t]:tree_end_idxs[t]
  rw_names <- sapply(idxs,
                     function(n) paste0('f[', n, ']'))
  util$plot_conn_pushforward_quantiles(samples, rw_names, data$years[idxs],
                                       xlab="Year", ylab="Log Ring Width Per mm", 
                                       display_ylim=c(-4, 0),
                                       main=paste0("Stand ", uniq_stand_ids[stand_idxs[t]], 
                                                  ", Tree ", uniq_tree_ids[t]))
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=1.0, col="white")
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=0.8, col="black")
  #abline(v=1974, lwd=2, lty=2, col="#DDDDDD")
}

# Posterior inference
par(mfrow=c(1, 2))

util$plot_expectand_pushforward(samples[['beta_gdd']], 25,
                                display_name="beta_gdd")
xs <- seq(0, 1, 0.01)
ys <- dnorm(xs, 0, log(1.8) / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['beta_sm']], 25,
                                display_name="beta_sm")
xs <- seq(0, 1, 0.01)
ys <- dnorm(xs, 0, log(1.8) / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)


par(mfrow=c(2, 1))
pred_names <- sapply(1:data$N, function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conditional_median_quantiles(samples, pred_names, data$gdd_obs,
                                       1500/10, 2060/10, 2.5, data$log_rw_obs, 
                                       xlab="GDD (x10 degC)")

util$plot_conditional_median_quantiles(samples, pred_names, data$sm_obs,
                                       0.17*100, 0.27*100, 0.5, data$log_rw_obs, 
                                       xlab="Soil moisture (%)")



par(mfrow=c(1, 2))
util$plot_expectand_pushforward(samples[[paste0('rho_sh')]], 100,
                                display_name="rho_sh",
                                flim = c(0,5))
xs <- seq(0, 10, 0.01)
ys <- dlnorm(xs, 1.7, 0.26)
lines(xs, ys, lwd=2, col=util$c_light)
util$plot_expectand_pushforward(samples[[paste0('rho')]], 80,
                                display_name="rho",
                                flim = c(5,35))
xs <- seq(0, 60, 0.01)
ys <- dlnorm(xs, 3.55, 0.24)
lines(xs, ys, lwd=2, col=util$c_light)

