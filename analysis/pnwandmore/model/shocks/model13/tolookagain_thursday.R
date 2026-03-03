
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


for(y in 1:data$N_all_years){
  name <- paste0('sum_sck_nonconc_states[', y, ']')
  gq_samples[[name]] <- 0
  obs_y <- which(data$all_years_idxs == y)
  for(i in obs_y){
    name_i <- paste0('sck_nonconc_state[', i, ']')
    gq_samples[[name]] <- gq_samples[[name]] + gq_samples[[name_i]]
  }
  propname <- paste0('prop_sck_nonconc_states[', y, ']')
  gq_samples[[propname]] <- gq_samples[[name]]/length(obs_y)
}

range(data$all_years_idxs[unlist(sapply(which(data$stand_idxs == s), 
                                        function(i) data$tree_start_idxs[i]:data$tree_end_idxs[i]))])


s <- which(data$uniq_stand_ids == 'S30')
sum(data$stand_idxs == s)

data$N_stand_trees[s]




for(y in 1:data$N_all_years){
  name <- paste0('sum_sck_states[', y, ']')
  gq_samples[[name]] <- 
    gq_samples[[paste0('sum_sck_nonconc_states[', y, ']')]] +
    gq_samples[[paste0('sum_sck_conc_states[', y, ']')]] 
  obs_y <- which(data$all_years_idxs == y)
  
  propname <- paste0('prop_sck_states[', y, ']')
  gq_samples[[propname]] <- gq_samples[[name]]/length(obs_y)
}

stand_years <- sapply(1:data$N_stands, function(s) 
  data$stand_start_years_idxs[s]:(data$stand_start_years_idxs[s]+ (data$N_stand_years[s]-1)))
for(y in 1:data$N_all_years){
  
  propname <- paste0('average_stand_prop_sck_nonconc_states[', y, ']')
  gq_samples[[propname]] <- 0
  stands <- na.omit(sapply(1:data$N_stands, function(s) ifelse(y %in% stand_years[[s]], s, NA)))
  
  for(s in stands){
    obs <- unlist(sapply(which(data$stand_idxs == s), function(i) data$tree_start_idxs[i]:data$tree_end_idxs[i]))
    obs_y <- which(data$all_years_idxs == y)
    obs_y <- intersect(obs, obs_y)
    
    namestand <- paste0('prop_sck_nonconc_states_stand')
    gq_samples[[namestand]] <- 0
    for(i in obs_y){
      name_i <- paste0('sck_nonconc_state[', i, ']')
      gq_samples[[namestand]] <- gq_samples[[namestand]] + gq_samples[[name_i]]
    }
    gq_samples[[namestand]] <- gq_samples[[namestand]]/length(obs_y)
    
    gq_samples[[propname]] <-  gq_samples[[propname]] + gq_samples[[namestand]]
  }
  gq_samples[[propname]]  <- gq_samples[[propname]]/length(stands)

}

for(y in 1:data$N_all_years){
  
  propname <- paste0('average_stand_prop_sck_conc_states[', y, ']')
  gq_samples[[propname]] <- 0
  stands <- na.omit(sapply(1:data$N_stands, function(s) ifelse(y %in% stand_years[[s]], s, NA)))
  
  for(s in stands){
    obs <- unlist(sapply(which(data$stand_idxs == s), function(i) data$tree_start_idxs[i]:data$tree_end_idxs[i]))
    obs_y <- which(data$all_years_idxs == y)
    obs_y <- intersect(obs, obs_y)
    print(obs_y)
    
    namestand <- paste0('prop_sck_conc_states_stand')
    gq_samples[[namestand]] <- 0
    for(i in obs_y){
      name_i <- paste0('sck_conc_state[', i, ']')
      gq_samples[[namestand]] <- gq_samples[[namestand]] + gq_samples[[name_i]]
    }
    gq_samples[[namestand]] <- gq_samples[[namestand]]/length(obs_y)
    
    gq_samples[[propname]] <-  gq_samples[[propname]] + gq_samples[[namestand]]
  }
  gq_samples[[propname]]  <- gq_samples[[propname]]/length(stands)
  
}

par(mfrow = c(2,1), mar = c(2,6,1,1))
names <- paste0('average_stand_prop_sck_nonconc_states[', 1:data$N_all_years, ']')
util$plot_disc_pushforward_quantiles(gq_samples, names,
                                     xlab="Year", 
                                     ylab="Average proportion of trees within a stand\nexperiencing a non-concordant shock",
                                     display_ylim=c(0, 0.15), xticklabs=rep(NA, data$N_all_years))
axis(1, 
     at = which(data$all_years %in% seq(1950, 2020, 10)),
     labels = seq(1950, 2020, 10),
     tck = 0.05)


names <- paste0('average_stand_prop_sck_conc_states[', 1:data$N_all_years, ']')
util$plot_disc_pushforward_quantiles(gq_samples, names,
                                     xlab="Year", 
                                     ylab="Average proportion of trees within a stand\nexperiencing a concordant shock",
                                     display_ylim=c(0, 0.15), xticklabs=rep(NA, data$N_all_years))
axis(1, 
     at = which(data$all_years %in% seq(1950, 2020, 10)),
     labels = seq(1950, 2020, 10),
     tck = 0.05)



par(mfrow = c(2,1), mar = c(2,6,1,1))
names <- paste0('prop_sck_nonconc_states[', 1:data$N_all_years, ']')
util$plot_disc_pushforward_quantiles(gq_samples, names,
                                     xlab="Year", 
                                     ylab="Proportion of trees experiencing\na non-concordant shock",
                                     display_ylim=c(0, 0.15), xticklabs=rep(NA, data$N_all_years))
axis(1, 
     at = which(data$all_years %in% seq(1950, 2020, 10)),
     labels = seq(1950, 2020, 10),
     tck = 0.05)


names <- paste0('sum_sck_nonconc_states[', 1:data$N_all_years, ']')
util$plot_disc_pushforward_quantiles(gq_samples, names,
                                     xlab="Year", 
                                     ylab="Number of trees experiencing\na non-concordant shock",
                                     display_ylim=c(0, 150), xticklabs=rep(NA, data$N_all_years))
axis(1, 
     at = which(data$all_years %in% seq(1950, 2020, 10)),
     labels = seq(1950, 2020, 10),
     tck = 0.05)

par(mfrow = c(2,1), mar = c(2,6,1,1))
names <- paste0('prop_sck_states[', 1:data$N_all_years, ']')
util$plot_disc_pushforward_quantiles(gq_samples, names,
                                     xlab="Year", 
                                     ylab="Proportion of trees experiencing\na shock (both types)",
                                     display_ylim=c(0, 0.25), xticklabs=rep(NA, data$N_all_years))
axis(1, 
     at = which(data$all_years %in% seq(1950, 2020, 10)),
     labels = seq(1950, 2020, 10),
     tck = 0.05)

names <- paste0('sum_sck_states[', 1:data$N_all_years, ']')
util$plot_disc_pushforward_quantiles(gq_samples, names,
                                     xlab="Year", 
                                     ylab="Number of trees experiencing\na shock (both types)",
                                     display_ylim=c(0, 250), xticklabs=rep(NA, data$N_all_years))
axis(1, 
     at = which(data$all_years %in% seq(1950, 2020, 10)),
     labels = seq(1950, 2020, 10),
     tck = 0.05)


par(mfrow = c(1,1), mar = c(5,5,1,1))
plot(x = NULL, y = NULL,
     xlim = c(0, 0.15), ylim = c(0, 0.15),
     xlab = "Trees in concordant shocks (%)",
     ylab = "Trees in non-concordant shocks (%)")
abline(a = 0, b = 1, lty = 2, col = 'grey70')
for(i in 1:data$N_all_years){
  namex <- paste0('prop_sck_conc_states[', i, ']')
  x <- util$ensemble_mcmc_quantile_est(gq_samples[[namex]], c(0.05, 0.5, 0.95))
  
  namey <- paste0('prop_sck_nonconc_states[', i, ']')
  y <- util$ensemble_mcmc_quantile_est(gq_samples[[namey]], c(0.05, 0.5, 0.95))
  
  segments(x0= x[2], y0 = y[1], y1 = y[3])
  segments(x0= x[1], x1 = x[3], y0 = y[2])
  points(x= x[2], y = y[2], pch = 20)
}



