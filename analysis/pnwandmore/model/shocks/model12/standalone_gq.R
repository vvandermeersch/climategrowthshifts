


library(cmdstanr)

mod_gq <- cmdstan_model(file.path(wd, 'model/stan', 'model12.stan'))

data$grainsize <-  ceiling(data$N_stands/8)
data$uniq_tree_ids <- NULL
data$uniq_species_ids <- NULL
data$uniq_stand_ids <- NULL
data$N_clades <- 1
data$uniq_stand_lat <- NULL
data$uniq_stand_lon <- NULL
fit_gq <- mod_gq$generate_quantities(fit, data = data, seed = 5838293, parallel_chains = 4)

gq_samples <- fit_gq$draws()
names <- dimnames(gq_samples)$variable
gq_samples <- lapply(1:dim(gq_samples)[3], function(k){t(matrix(gq_samples[1:dim(gq_samples)[1],1:dim(gq_samples)[2],k], 
                                                                nrow = dim(gq_samples)[1], ncol = dim(gq_samples)[2]))})
names(gq_samples) <- names

# Stand 63 seems very shocky
c(1:data$N_stands)[sapply(1:data$N_stands, function(s){util$ensemble_mcmc_quantile_est(samples[[paste0('phi_sck[',s,']')]], c(0.5))}) > 0.3]
s <- 62
trees_stand <- which(data$stand_idxs == s)

par(mfrow = c(1,3), cex.lab = 1.2)
for(t in trees_stand){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- paste0("log_rw_pred[",idxs,"]")

  util$plot_conn_pushforward_quantiles(gq_samples, names, data$years[idxs]-1989,
                                       xlab="Year", ylab="Log ring width (per mm)", 
                                       display_ylim=c(-8, 2), display_xlim = c(1990, 2010)-1989)
  points(data$years[idxs]-1989, log(data$rw_obs[idxs]), pch=16, cex=1, col="white")
  points(data$years[idxs]-1989, log(data$rw_obs[idxs]), pch=16, cex=0.5, col="black")
  text(x = .5, y = -6.9, labels = data$uniq_tree_ids[t], cex = 1, adj = 0)
  text(x = .5, y = -7.7, labels = paste0('Species:', data$uniq_species_ids[data$species_idxs[t]]), cex = 1, adj = 0)
  
  
  names <- paste0("sck_state[",idxs,"]")
  util$plot_disc_pushforward_quantiles(gq_samples, names, data$years[idxs],
                                       xlab="Year", ylab="Shock state", 
                                       display_ylim=c(-0.05, 1.05))
  
  names <- paste0("delta_sck[",idxs,"]")
  util$plot_conn_pushforward_quantiles(gq_samples, names, data$years[idxs]-1989,
                                       xlab="Year", ylab="Shock amplitude", 
                                       display_ylim=c(-2, 1))
}


fsh_samples <- fit_upd$draws(variables = c("f_sh"))
names <- dimnames(fsh_samples)$variable
fsh_samples <- lapply(1:dim(fsh_samples)[3], function(k){t(matrix(fsh_samples[1:dim(fsh_samples)[1],1:dim(fsh_samples)[2],k], 
                                                                nrow = dim(fsh_samples)[1], ncol = dim(fsh_samples)[2]))})
names(fsh_samples) <- names

par(mfrow = c(1,2), cex.lab = 0.8)
s <- 29
idxs <- seq(1, data$N_all_years, 1)
names <- paste0("f_sh[",s,",",idxs,"]")
util$plot_conn_pushforward_quantiles(fsh_samples, names, 1:data$N_all_years,
                                     xlab="Year", ylab="Log ring width (per mm)", 
                                     display_ylim=c(-3, 3), display_xlim = c(1990, 2010)-1989)
text(x = .7, y = 2.9, labels = paste0('Stand ', s), cex = 0.7, adj = 0)
s <- 63
idxs <- seq(1, data$N_all_years, 1)
names <- paste0("f_sh[",s,",",idxs,"]")
util$plot_conn_pushforward_quantiles(fsh_samples, names, 1:data$N_all_years,
                                     xlab="Year", ylab="Log ring width (per mm)", 
                                     display_ylim=c(-3, 3), display_xlim = c(1990, 2010)-1989)
text(x = .7, y = 2.9, labels = paste0('Stand ', s), cex = 0.7, adj = 0)