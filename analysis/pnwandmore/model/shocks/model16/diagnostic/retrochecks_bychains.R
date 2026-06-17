
gq_samples <- fit_gq$draws(c('log_rw_pred', 'f', 'f_ind', 'delta_conc_sck', 'delta_nonconc_sck'))
gc()
names <- dimnames(gq_samples)$variable


chains <- c(3)


gq_samples_subchains <- lapply(1:dim(gq_samples)[3], 
                               function(k){t(matrix(gq_samples[1:dim(gq_samples)[1],chains,k], 
                                                    nrow = dim(gq_samples)[1], ncol =  length(chains)))})
names(gq_samples_subchains) <- names

stand_samples <- fit$draws(c('f_sh'))
gc()
names <- dimnames(stand_samples)$variable
stand_samples_subchains <- lapply(1:dim(stand_samples)[3], 
                                  function(k){t(matrix(gq_samples[1:dim(gq_samples)[1],chains,k], 
                                                       nrow = dim(gq_samples)[1], ncol =  length(chains)))})
names(stand_samples_subchains) <- names



pdf(file.path(wd, 'model/shocks/model16/figures', 'oregonTSME_newpriors_chain3.pdf'), height = 10, width = 19.5)
par(mfrow = c(8,6), cex.lab = 0.8, cex.axis = 0.8, mar = c(1.5,4,0.5,1), mgp = c(1.5, 0.2, 0),tck = -0.03)
for(t in 1:data$N_trees){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- paste0("log_rw_pred[",idxs,"]")
  
  util$plot_conn_pushforward_quantiles(gq_samples_subchains, names, data$years[idxs],
                                       xlab="Year", ylab="Log ring width (per mm)", 
                                       display_ylim=c(-8, 2), display_xlim = range(data$all_years))
  points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=0.5, col="white")
  points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=0.25, col="black")
  text(x = .5, y = -6.9, labels = data$uniq_tree_ids[t], cex = 1, adj = 0)
  text(x = .5, y = -7.7, labels = paste0('Species:', data$uniq_species_ids[data$species_idxs[t]]), cex = 1, adj = 0)
  
  names <- paste0("delta_conc_sck[",idxs,"]")
  util$plot_conn_pushforward_quantiles(gq_samples_subchains, names, data$years[idxs],
                                       xlab="Year", ylab="Concordant\nshock amplitude", 
                                       display_ylim=c(-5, 1), display_xlim = range(data$all_years))
  
  names <- paste0("delta_nonconc_sck[",idxs,"]")
  util$plot_conn_pushforward_quantiles(gq_samples_subchains, names, data$years[idxs],
                                       xlab="Year", ylab="Non-concordant\nshock amplitude", 
                                       display_ylim=c(-5, 1), display_xlim = range(data$all_years))
  
  names <- paste0("f[",idxs,"]")
  util$plot_conn_pushforward_quantiles(gq_samples_subchains, names, data$years[idxs],
                                       xlab="Year", ylab="Long-term GP", 
                                       display_ylim=c(-5, 2), display_xlim = range(data$all_years))
  
  names <- paste0("f_ind[",idxs,"]")
  util$plot_conn_pushforward_quantiles(gq_samples_subchains, names, data$years[idxs],
                                       xlab="Year", ylab="Short-term GP", 
                                       display_ylim=c(-5, 2), display_xlim = range(data$all_years))
  
  s <- unique(data$stand_idxs[t])
  names <- paste0("f_sh[",s,',',1:data$N_all_years,"]")
  util$plot_conn_pushforward_quantiles(stand_samples_subchains, names, data$all_years,
                                       xlab="Year", ylab="Stand GP", 
                                       display_ylim=c(-8, 2), display_xlim = range(data$all_years))
  
  
}
dev.off()
