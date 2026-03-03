rm(list = ls())



# Generate quantities
mod_gq <- cmdstan_model(file.path(wd, 'model/stan', 'model13_updatedpriors_wGQ.stan'))
data_gq <- data
data_gq$grainsize <-  ceiling(data_gq$N_stands/8)
data_gq$uniq_tree_ids <- NULL
data_gq$uniq_species_ids <- NULL
data_gq$uniq_stand_ids <- NULL
data_gq$N_clades <- 1
data_gq$uniq_stand_lat <- NULL
data_gq$uniq_stand_lon <- NULL

data_gq$tree_pred_idxs <- which(grepl('wy067', data$uniq_tree_ids))
data_gq$tree_pred_idxs <- which(data$stand_idxs == 6)
data_gq$N_pred <- sum(sapply(data_gq$tree_pred_idxs , function(t) length(data$tree_start_idxs[t]:data$tree_end_idxs[t])))
data_gq$N_trees_pred <- length(data_gq$tree_pred_idxs )

data_gq$tree_pred_idxs <- 1:data$N_trees
data_gq$N_pred <- sum(sapply(data_gq$tree_pred_idxs , function(t) length(data$tree_start_idxs[t]:data$tree_end_idxs[t])))
data_gq$N_trees_pred <- length(data_gq$tree_pred_idxs )

fit_gq <- mod_gq$generate_quantities(fit, data = data_gq, seed = 5838293, parallel_chains = 1)
gc()
gq_samples <- fit_gq$draws(variables = c('sck_conc_state'))
gc()
names <- dimnames(gq_samples)$variable
gq_samples <- lapply(1:dim(gq_samples)[3], function(k){t(matrix(gq_samples[1:dim(gq_samples)[1],1:dim(gq_samples)[2],k], 
                                                                nrow = dim(gq_samples)[1], ncol = dim(gq_samples)[2]))})
names(gq_samples) <- names


for(y in 1:data$N_all_years){
  name <- paste0('sum_sck_conc_states[', y, ']')
  gq_samples[[name]] <- 0
  obs_y <- which(data$all_years_idxs == y)
  for(i in obs_y){
    name_i <- paste0('sck_conc_state[', i, ']')
    gq_samples[[name]] <- gq_samples[[name]] + gq_samples[[name_i]]
  }
  propname <- paste0('prop_sck_conc_states[', y, ']')
  gq_samples[[propname]] <- gq_samples[[name]]/length(obs_y)
}

par(mfrow = c(2,1), mar = c(2,6,1,1), cex.axis = 0.8, cex.lab = 0.85)
names <- paste0('prop_sck_conc_states[', 1:data$N_all_years, ']')
util$plot_disc_pushforward_quantiles(gq_samples, names,
                                     xlab="Year", 
                                     ylab="Proportion of trees\nexperiencing a concordant shock",
                                     display_ylim=c(0, 0.15), xticklabs=rep(NA, data$N_all_years))
axis(1, 
     at = which(data$all_years %in% seq(1950, 2020, 10)),
     labels = seq(1950, 2020, 10),
     tck = 0.05)

stand_years <- sapply(1:data$N_stands, function(s) 
  data$stand_start_years_idxs[s]:(data$stand_start_years_idxs[s]+ (data$N_stand_years[s]-1)))

for(y in 1:data$N_all_years){
  print(y)
  
  propname <- paste0('average_stand_prop_sck_conc_states[', y, ']')
  gq_samples[[propname]] <- 0
  stands <- as.numeric(na.omit(sapply(1:data$N_stands, function(s) ifelse(y %in% stand_years[[s]], s, NA))))
  
  for(s in stands){
    obs <- unlist(sapply(which(data$stand_idxs == s), function(i) data$tree_start_idxs[i]:data$tree_end_idxs[i]))
    obs_y <- which(data$all_years_idxs == y)
    obs_y <- intersect(obs, obs_y)
    # print(obs_y)
    
    namehere <- paste0('prop_sck_conc_states_stand')
    gq_samples[[namehere]] <- 0
    for(i in obs_y){
      name_i <- paste0('sck_conc_state[', i, ']')
      gq_samples[[namehere]] <- gq_samples[[namehere]] + gq_samples[[name_i]]
    }
    gq_samples[[namehere]] <- gq_samples[[namehere]]/length(obs_y)
    
    gq_samples[[propname]] <-  gq_samples[[propname]] + gq_samples[[namehere]]
  }
  gq_samples[[propname]]  <- gq_samples[[propname]]/length(stands)
  
}

newname <- paste0('average_stand_prop_sck_conc_states_allyears')
gq_samples[[newname]] <- 0
for(y in 1:data$N_all_years){
  propname <- paste0('average_stand_prop_sck_conc_states[', y, ']')
  gq_samples[[newname]] <-  gq_samples[[newname]] + gq_samples[[propname]]
}
gq_samples[[newname]] <- gq_samples[[newname]]/data$N_all_years

names <- paste0('average_stand_prop_sck_conc_states[', 1:data$N_all_years, ']')
util$plot_disc_pushforward_quantiles(gq_samples, names,
                                     xlab="Year", 
                                     ylab="Average proportion of trees within a stand\nexperiencing a concordant shock",
                                     display_ylim=c(0, 0.15), xticklabs=rep(NA, data$N_all_years))
axis(1, 
     at = which(data$all_years %in% seq(1950, 2020, 10)),
     labels = seq(1950, 2020, 10),
     tck = 0.05)

name <- paste0('average_stand_prop_sck_conc_states_allyears')
abline(h = util$ensemble_mcmc_quantile_est(gq_samples[[name]], c(0.5)), lty = 5, col = 'grey60')
abline(h = util$ensemble_mcmc_quantile_est(gq_samples[[name]], c(0.05, 0.95)), lty = 3, col = 'grey80')
text(x = (which(data$all_years == 1988) + which(data$all_years == 1992))/2, y = 0.11, labels = '1988-1990\n(1992?)', cex = 0.8)
segments(x0 = which(data$all_years == 1988), x1 = which(data$all_years == 1992), y0 = 0.097)
segments(x0 = which(data$all_years == 1988), y0 = 0.097, y1 = 0.094)
segments(x0 = which(data$all_years == 1992), y0 = 0.097, y1 = 0.094)
text(x = which(data$all_years == 1977), y = 0.11, labels = '1977', cex = 0.8, adj = 0.5)
text(x = which(data$all_years == 2002), y = 0.085, labels = '2002', cex = 0.8, adj = 0.5)
text(x = which(data$all_years == 2004), y = 0.077, labels = '2004', cex = 0.8, adj = 0.5)



