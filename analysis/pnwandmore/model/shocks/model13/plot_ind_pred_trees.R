

par(mfrow = c(3,3), mar = c(4,4,1,1))
idx <- 1
for(i in 1:data_gq$N_trees_pred){
  
  t <- data_gq$tree_pred_idxs[i]
  start <- idx 
  idx <- idx + data_gq$N_years[t] - 1
  end <- idx
  
  local_idxs <- start:end
  global_idxs <- data_gq$tree_start_idxs[t]:data_gq$tree_end_idxs[t]
  
  names <- paste0("log_rw_pred[",local_idxs,"]")
  # 
  util$plot_conn_pushforward_quantiles(gq_samples, names, data_gq$all_years_idxs[global_idxs],
                                       xlab="Year", ylab="Log ring width (per mm)",
                                       display_ylim=c(-8, 2))
  # points(data$years[idxs]-1979, log(data$rw_obs[idxs]), pch=16, cex=1, col="white")
  # points(data$years[idxs]-1979, log(data$rw_obs[idxs]), pch=16, cex=0.5, col="black")
  # text(x = .5, y = -6.9, labels = data$uniq_tree_ids[t], cex = 1, adj = 0)
  # text(x = .5, y = -7.7, labels = paste0('Species:', data$uniq_species_ids[data$species_idxs[t]]), cex = 1, adj = 0)
  abline(v = c(10:12, 31:33), lty = 2)
  
  names <- paste0("sck_conc_state[",local_idxs,"]")
  util$plot_disc_pushforward_quantiles(gq_samples, names,
                                       xlab="Year", ylab="Concordant shock",
                                       display_ylim=c(-0.05, 1.05))
  abline(v = c(10:12, 31:33), lty = 2)
  
  names <- paste0("sck_nonconc_state[",local_idxs,"]")
  util$plot_disc_pushforward_quantiles(gq_samples, names,
                                       xlab="Year", ylab="Non-concordant shock",
                                       display_ylim=c(-0.05, 1.05))
  abline(v = c(10:12, 31:33), lty = 2)
  
  idx <- idx + 1
}


stand <- unique(data$stand_idxs[which(grepl('wy067', data$uniq_tree_ids))])
stand <- 6
data_gq$tree_pred_idxs <- which(data$stand_idxs == 6)
data$species_idxs[data_gq$tree_pred_idxs]

par(mfrow = c(2,1))
names <- paste0("f_sh[", stand,',',1:data$N_stand_years[stand],"]")
util$plot_conn_pushforward_quantiles(standgp_samples, names, 1:data$N_stand_years[stand],
                                     xlab="Year", ylab="Stand GP",
                                     display_ylim=c(-8, 2))
abline(v = c(10:12, 31:33), lty = 2)

for(n in names){
  newn <- paste0('kappa_', n)
  standgp_samples[[newn]] <- standgp_samples[[n]]*param_samples[[paste0('kappa_sh[',4,']')]]
}

names <- paste0("kappa_f_sh[", stand,',',1:data$N_stand_years[stand],"]")
util$plot_conn_pushforward_quantiles(standgp_samples, names, 1:data$N_stand_years[stand],
                                     xlab="Year", ylab=bquote(kappa[species] %*% "Stand GP"),
                                     display_ylim=c(-8, 2))
abline(v = c(10:12, 31:33), lty = 2)

