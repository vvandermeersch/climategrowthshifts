rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(cmdstanr)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)


data <-readRDS(file.path(wd, 'output/model', 'data_30jan2025_gymnosperms_standclimate_19502024_7species.rds'))
fit <- readRDS(file.path(wd, 'output/model/model13', 'fit_16feb2026_gymnosperms_standclimate_19502024_7species_model13_updatedpriors.rds'))

# Standalone generated quantities
# mod_gq <- cmdstan_model(file.path(wd, 'model/stan', 'model13_updatedpriors_wGQ.stan'))
# data_gq <- data
# data_gq$grainsize <-  ceiling(data_gq$N_stands/8)
# data_gq$uniq_tree_ids <- NULL
# data_gq$uniq_species_ids <- NULL
# data_gq$uniq_stand_ids <- NULL
# data_gq$N_clades <- 1
# data_gq$uniq_stand_lat <- NULL
# data_gq$uniq_stand_lon <- NULL
# data_gq$tree_pred_idxs <- 1:data$N_trees # all
# data_gq$N_pred <- sum(sapply(data_gq$tree_pred_idxs , function(t) length(data$tree_start_idxs[t]:data$tree_end_idxs[t])))
# data_gq$N_trees_pred <- length(data_gq$tree_pred_idxs )
# fit_gq <- mod_gq$generate_quantities(fit, data = data_gq, seed = 5838293, parallel_chains = 1)
# gc()
# rds_file <- file.path(wd, 'output', 'fit_16feb2026_gymnosperms_standclimate_19502024_7species_model13_updatedpriors_SGQ.rds')
# fit_gq$save_object(file = rds_file) # too much memory for my computer...
# sck_conc_states <- fit_gq$draws(variables = c('sck_conc_state')) # 'sck_nonconc_state', 'delta_sck'
# saveRDS(sck_conc_states, file.path(wd, 'output/model/model13', 'sck_conc_states_16feb2026_gymnosperms_standclimate_19502024_7species_model13_updatedpriors.rds'))
# sck_nonconc_states <- fit_gq$draws(variables = c('sck_nonconc_state')) # 'sck_nonconc_state', 'delta_sck'
# saveRDS(sck_nonconc_states, file.path(wd, 'output/model/model13', 'sck_nonconc_states_16feb2026_gymnosperms_standclimate_19502024_7species_model13_updatedpriors.rds'))
# delta_sck <- fit_gq$draws(variables = c('delta_sck')) # 'sck_nonconc_state', 'delta_sck'
# saveRDS(delta_sck, file.path(wd, 'output/model/model13', 'delta_sck_16feb2026_gymnosperms_standclimate_19502024_7species_model13_updatedpriors.rds'))

delta_sck <- readRDS(file.path(wd, 'output/model/model13', 'delta_sck_16feb2026_gymnosperms_standclimate_19502024_7species_model13_updatedpriors.rds'))
names <- dimnames(delta_sck)$variable
delta_sck <- lapply(1:dim(delta_sck)[3], function(k){t(matrix(delta_sck[1:dim(delta_sck)[1],1:dim(delta_sck)[2],k], 
                                                                nrow = dim(delta_sck)[1], ncol = dim(delta_sck)[2]))})
names(delta_sck) <- names
gc()

sck_conc_states <- readRDS(file.path(wd, 'output/model/model13', 'sck_conc_states_16feb2026_gymnosperms_standclimate_19502024_7species_model13_updatedpriors.rds'))
names <- dimnames(sck_conc_states)$variable
sck_conc_states <- lapply(1:dim(sck_conc_states)[3], function(k){t(matrix(sck_conc_states[1:dim(sck_conc_states)[1],1:dim(sck_conc_states)[2],k], 
                                                              nrow = dim(sck_conc_states)[1], ncol = dim(sck_conc_states)[2]))})
names(sck_conc_states) <- names
gc()

sck_nonconc_states <- readRDS(file.path(wd, 'output/model/model13', 'sck_nonconc_states_16feb2026_gymnosperms_standclimate_19502024_7species_model13_updatedpriors.rds'))
names <- dimnames(sck_nonconc_states)$variable
sck_nonconc_states <- lapply(1:dim(sck_nonconc_states)[3], function(k){t(matrix(sck_nonconc_states[1:dim(sck_nonconc_states)[1],1:dim(sck_nonconc_states)[2],k], 
                                                                          nrow = dim(sck_nonconc_states)[1], ncol = dim(sck_nonconc_states)[2]))})
names(sck_nonconc_states) <- names
gc()


# Separate amplitude by type of shocks
delta_conc_sck <- delta_sck
names(delta_conc_sck) <- paste0("delta_conc_sck[", 1:data$N, "]")
for(i in 1:data$N){
  concname <- paste0("sck_conc_state[", i, "]")
  name <- paste0("delta_conc_sck[", i, "]")
  nsamples <- length(delta_conc_sck[[name]])
  delta_conc_sck[[name]][1:nsamples] <- ifelse(sck_conc_states[[concname]][1:nsamples] == 0, 0, delta_conc_sck[[name]][1:nsamples])
}

delta_nonconc_sck <- delta_sck
names(delta_nonconc_sck) <- paste0("delta_nonconc_sck[", 1:data$N, "]")
for(i in 1:data$N){
  nonconcname <- paste0("sck_nonconc_state[", i, "]")
  name <- paste0("delta_nonconc_sck[", i, "]")
  nsamples <- length(delta_nonconc_sck[[name]])
  delta_nonconc_sck[[name]][1:nsamples] <- ifelse(sck_nonconc_states[[nonconcname]][1:nsamples] == 0, 0, delta_nonconc_sck[[name]][1:nsamples])
}

# Average amplitude of concordant shocks in 1950-1980
newname1 <- 'n_conc_sck_19501980'
sck_conc_states[[newname1]] <- 0
newname2 <- 'avg_delta_conc_sck_19501980'
delta_conc_sck[[newname2]] <- 0
idxs <- which(data$all_years_idxs %in% which(data$all_years %in% c(1950:1980)))
for(i in idxs){
  name1 <- paste0("sck_conc_state[", i, "]")
  sck_conc_states[[newname1]] <- sck_conc_states[[newname1]] + sck_conc_states[[name1]]
  
  name2 <- paste0("delta_conc_sck[", i, "]")
  delta_conc_sck[[newname2]] <- delta_conc_sck[[newname2]] + delta_conc_sck[[name2]]
}
delta_conc_sck[[newname2]] <- delta_conc_sck[[newname2]]/sck_conc_states[[newname1]]
util$ensemble_mcmc_quantile_est(delta_conc_sck[['avg_delta_conc_sck_19501980']], c(0.05, 0.5, 0.95))

# Average amplitude of concordant shocks in 1981-today
newname1 <- 'n_conc_sck_19812020'
sck_conc_states[[newname1]] <- 0
newname2 <- 'avg_delta_conc_sck_19812020'
delta_conc_sck[[newname2]] <- 0
idxs <- which(data$all_years_idxs %in% which(data$all_years %in% c(1981:2020)))
for(i in idxs){
  name1 <- paste0("sck_conc_state[", i, "]")
  sck_conc_states[[newname1]] <- sck_conc_states[[newname1]] + sck_conc_states[[name1]]
  
  name2 <- paste0("delta_conc_sck[", i, "]")
  delta_conc_sck[[newname2]] <- delta_conc_sck[[newname2]] + delta_conc_sck[[name2]]
}
delta_conc_sck[[newname2]] <- delta_conc_sck[[newname2]]/sck_conc_states[[newname1]]
util$ensemble_mcmc_quantile_est(delta_conc_sck[['avg_delta_conc_sck_19812020']], c(0.05, 0.5, 0.95))

util$plot_expectand_pushforward(delta_conc_sck[['avg_delta_conc_sck_19501980']], B = 30,
                                flim = c(-0.3, 0), ylim = c(0,20),
                                col = util$c_mid_teal, display_name = 'Average amplitude of concordant shocks\n(log-scale)')
util$plot_expectand_pushforward(delta_conc_sck[['avg_delta_conc_sck_19812020']], B = 30,
                                flim = c(-0.3, 0), , ylim = c(0,20),
                                add = T, col = util$c_mid)
text(x = -0.11, y = 17, labels = '1981-2020', col = util$c_mid)
text(x = -0.25, y = 10, labels = '1950-1980', col = util$c_mid_teal)

# Average amplitude of non-concordant shocks in 1950-1980
newname1 <- 'n_nonconc_sck_19501980'
sck_nonconc_states[[newname1]] <- 0
newname2 <- 'avg_delta_nonconc_sck_19501980'
delta_nonconc_sck[[newname2]] <- 0
idxs <- which(data$all_years_idxs %in% which(data$all_years %in% c(1950:1980)))
for(i in idxs){
  name1 <- paste0("sck_nonconc_state[", i, "]")
  sck_nonconc_states[[newname1]] <- sck_nonconc_states[[newname1]] + sck_nonconc_states[[name1]]
  
  name2 <- paste0("delta_nonconc_sck[", i, "]")
  delta_nonconc_sck[[newname2]] <- delta_nonconc_sck[[newname2]] + delta_nonconc_sck[[name2]]
}
delta_nonconc_sck[[newname2]] <- delta_nonconc_sck[[newname2]]/sck_nonconc_states[[newname1]]
util$ensemble_mcmc_quantile_est(delta_nonconc_sck[['avg_delta_nonconc_sck_19501980']], c(0.05, 0.5, 0.95))

# Average amplitude of non-concordant shocks in 1981-today
newname1 <- 'n_nonconc_sck_19812020'
sck_nonconc_states[[newname1]] <- 0
newname2 <- 'avg_delta_nonconc_sck_19812020'
delta_nonconc_sck[[newname2]] <- 0
idxs <- which(data$all_years_idxs %in% which(data$all_years %in% c(1981:2020)))
for(i in idxs){
  name1 <- paste0("sck_nonconc_state[", i, "]")
  sck_nonconc_states[[newname1]] <- sck_nonconc_states[[newname1]] + sck_nonconc_states[[name1]]
  
  name2 <- paste0("delta_nonconc_sck[", i, "]")
  delta_nonconc_sck[[newname2]] <- delta_nonconc_sck[[newname2]] + delta_nonconc_sck[[name2]]
}
delta_nonconc_sck[[newname2]] <- delta_nonconc_sck[[newname2]]/sck_nonconc_states[[newname1]]
util$ensemble_mcmc_quantile_est(delta_nonconc_sck[['avg_delta_nonconc_sck_19812020']], c(0.05, 0.5, 0.95))

util$plot_expectand_pushforward(delta_nonconc_sck[['avg_delta_nonconc_sck_19501980']], B = 30,
                                flim = c(-0.3, 0), ylim = c(0,20),
                                col = util$c_mid_teal, display_name = 'Average amplitude of non-concordant shocks\n(log-scale)')
util$plot_expectand_pushforward(delta_nonconc_sck[['avg_delta_nonconc_sck_19812020']], B = 30,
                                flim = c(-0.3, 0), , ylim = c(0,20),
                                add = T, col = util$c_mid)
text(x = -0.105, y = 17, labels = '1981-2020', col = util$c_mid)
text(x = -0.25, y = 10, labels = '1950-1980', col = util$c_mid_teal)


idxs <- which(data$all_years_idxs %in% which(data$all_years %in% c(1950:1980)))
util$plot_expectand_pushforward(sck_conc_states[['n_conc_sck_19501980']]/length(idxs), B = 30,
                                flim = c(0, 0.1), 
                                col = util$c_mid_teal, display_name = 'Number of shocks / number of observations')
util$plot_expectand_pushforward(sck_nonconc_states[['n_nonconc_sck_19501980']]/length(idxs), B = 30,
                                flim = c(0, 0.1), ylim = c(0,0.002),
                                col = '#b5c6c6', add = T)
text(x = 0.019, y = 50, labels = 'Non-concordant', col = '#b5c6c6')
text(x = 0.085, y = 45, labels = 'Concordant', col = util$c_mid_teal)
text(x = 0, y = 100, adj = 0, labels = '1950-1980', col = 'black')

idxs <- which(data$all_years_idxs %in% which(data$all_years %in% c(1981:2020)))
util$plot_expectand_pushforward(sck_conc_states[['n_conc_sck_19812020']]/length(idxs), B = 30,
                                flim = c(0, 0.1), 
                                col = util$c_mid, display_name = 'Number of shocks / number of observations')
util$plot_expectand_pushforward(sck_nonconc_states[['n_nonconc_sck_19812020']]/length(idxs), B = 30,
                                flim = c(0, 0.1), ylim = c(0,0.002),
                                col = util$c_light, add = T)
text(x = 0.022, y = 50, labels = 'Non-concordant', col = util$c_light)
text(x = 0.09, y = 45, labels = 'Concordant', col = util$c_mid)
text(x = 0, y = 90, adj = 0, labels = '1981-2020', col = 'black')






