

for(s in 1:data$N_stands){
  newname <- paste0('avg_thetha_idio[2][', s,']')
  base_samples[[newname]] <- 0
  
  trees <- which(data$stand_idxs %in% s)
  for(t in which(data$stand_idxs %in% s)){
    name <- paste0('thetas_idio[',t,',2]')
    base_samples[[newname]] <- base_samples[[newname]] + base_samples[[name]]
  }
  base_samples[[newname]] <- base_samples[[newname]]/length(trees)
}

par(mfrow = c(1,1), mar = c(0,4,1,1), cex.lab = 0.9, cex.axis = 0.9)
util$plot_disc_pushforward_quantiles(base_samples, 
                                     paste0('avg_thetha_idio[2][', order(data$uniq_stand_lat),']'), 
                                     display_ylim = c(0,0.8),
                                     ylab = "Average p(idiosyncratic shock), by stand")

util$plot_conn_pushforward_quantiles(base_samples, 
                                     paste0('avg_thetha_idio[2][', order(data$uniq_stand_lat),']'), 
                                     plot_xs = as.numeric(data$uniq_stand_lat[order(data$uniq_stand_lat)]),
                                     display_ylim = c(0,0.8),
                                     ylab = "Average p(idiosyncratic shock), by stand")

par(mfrow = c(1,1), mar = c(4,4,1,1), cex.lab = 0.9, cex.axis = 0.9)
plot(x = NULL, y = NULL,
     xlim = c(0,1), ylim = c(0,0.5),
     xlab = 'p(stand concordance)', 
     ylab = "Average p(idiosyncratic shock)",
     bty = 'n')
for(s in 1:data$N_stands){
  name1 <- paste0('avg_thetha_idio[2][', s,']')
  name2 <-  paste0('phi_sck[', s,']')
  segments(x0 = util$ensemble_mcmc_quantile_est(base_samples[[name2]], c(0.5)),
           y0 = util$ensemble_mcmc_quantile_est(base_samples[[name1]], c(0.05)),
           y1 = util$ensemble_mcmc_quantile_est(base_samples[[name1]], c(0.95)),
           lwd = 0.5,  col = util$c_light)
  segments(x0 = util$ensemble_mcmc_quantile_est(base_samples[[name2]], c(0.05)),
           x1 = util$ensemble_mcmc_quantile_est(base_samples[[name2]], c(0.95)),
           y0 = util$ensemble_mcmc_quantile_est(base_samples[[name1]], c(0.5)),
           lwd = 0.5,  col = util$c_light)
  points(y = util$ensemble_mcmc_quantile_est(base_samples[[name1]], c(0.5)),
         x =  util$ensemble_mcmc_quantile_est(base_samples[[name2]], c(0.5)),
         pch = 20, cex = 0.5, col = util$c_mid)
}



for(s in 1:data$N_stand_species){
  newname <- paste0('avg_thetha_idio[2][', s,']')
  base_samples[[newname]] <- 0
  
  trees <- which(data$stand_species_idxs %in% s)
  for(t in which(data$stand_species_idxs %in% s)){
    name <- paste0('thetas_idio[',t,',2]')
    base_samples[[newname]] <- base_samples[[newname]] + base_samples[[name]]
  }
  base_samples[[newname]] <- base_samples[[newname]]/length(trees)
}


par(mfrow = c(1,1), mar = c(4,4,1,1), cex.lab = 0.9, cex.axis = 0.9)
plot(x = NULL, y = NULL,
     xlim = c(0,1), ylim = c(0,0.5),
     xlab = 'p(tree shock | stand concordance)', 
     ylab = "Average p(idiosyncratic shock)",
     bty = 'n')
for(s in 1:data$N_stand_species){
  name1 <- paste0('avg_thetha_idio[2][', s,']')
  name2 <-  paste0('omega_conc_sck[', s,']')
  segments(x0 = util$ensemble_mcmc_quantile_est(base_samples[[name2]], c(0.5)),
           y0 = util$ensemble_mcmc_quantile_est(base_samples[[name1]], c(0.05)),
           y1 = util$ensemble_mcmc_quantile_est(base_samples[[name1]], c(0.95)),
           lwd = 0.5,  col = util$c_light)
  segments(x0 = util$ensemble_mcmc_quantile_est(base_samples[[name2]], c(0.05)),
           x1 = util$ensemble_mcmc_quantile_est(base_samples[[name2]], c(0.95)),
           y0 = util$ensemble_mcmc_quantile_est(base_samples[[name1]], c(0.5)),
           lwd = 0.5,  col = util$c_light)
  points(y = util$ensemble_mcmc_quantile_est(base_samples[[name1]], c(0.5)),
         x =  util$ensemble_mcmc_quantile_est(base_samples[[name2]], c(0.5)),
         pch = 20, cex = 0.5, col = util$c_mid)
}


par(mfrow = c(1,1), mar = c(4,4,1,1), cex.lab = 0.9, cex.axis = 0.9)
plot(x = NULL, y = NULL,
     xlim = c(0,1), ylim = c(0,0.5),
     xlab = 'p(stand concordance)*p(tree shock | stand concordance)', 
     ylab = "Average p(idiosyncratic shock)",
     bty = 'n')
for(s in 1:data$N_stand_species){
  name1 <- paste0('avg_thetha_idio[2][', s,']')
  name2 <-  paste0('omega_conc_sck[', s,']')
  
  stand <- unique(data$stand_idxs[which(data$stand_species_idxs == s)])
  name3 <-  paste0('phi_sck[', stand,']')
  segments(x0 = util$ensemble_mcmc_quantile_est(base_samples[[name2]]*base_samples[[name3]], c(0.5)),
           y0 = util$ensemble_mcmc_quantile_est(base_samples[[name1]], c(0.05)),
           y1 = util$ensemble_mcmc_quantile_est(base_samples[[name1]], c(0.95)),
           lwd = 0.5,  col = util$c_light)
  segments(x0 = util$ensemble_mcmc_quantile_est(base_samples[[name2]]*base_samples[[name3]], c(0.05)),
           x1 = util$ensemble_mcmc_quantile_est(base_samples[[name2]]*base_samples[[name3]], c(0.95)),
           y0 = util$ensemble_mcmc_quantile_est(base_samples[[name1]], c(0.5)),
           lwd = 0.5,  col = util$c_light)
  points(y = util$ensemble_mcmc_quantile_est(base_samples[[name1]], c(0.5)),
         x =  util$ensemble_mcmc_quantile_est(base_samples[[name2]]*base_samples[[name3]], c(0.5)),
         pch = 20, cex = 0.5, col = util$c_mid)
}



par(mfrow = c(1,1), mar = c(4,5,1,1), cex.lab = 0.9, cex.axis = 0.9)
plot(x = NULL, y = NULL,
     xlim = c(30,50), ylim = c(0,0.5),
     xlab = 'Latitude', 
     ylab = ("p(stand concordance) × \np(tree shock | stand concordance)"),
     bty = 'n')
for(s in 1:data$N_stand_species){
  name2 <-  paste0('omega_conc_sck[', s,']')
  stand <- unique(data$stand_idxs[which(data$stand_species_idxs == s)])
  name3 <-  paste0('phi_sck[', stand,']')
  segments(x0 = as.numeric(data$uniq_stand_lat[stand]),
           y0 = util$ensemble_mcmc_quantile_est(base_samples[[name2]]*base_samples[[name3]], c(0.05)),
           y1 = util$ensemble_mcmc_quantile_est(base_samples[[name2]]*base_samples[[name3]], c(0.95)),
           lwd = 0.5,  col = util$c_light)
  points(x = as.numeric(data$uniq_stand_lat[stand]),
         y =  util$ensemble_mcmc_quantile_est(base_samples[[name2]]*base_samples[[name3]], c(0.5)),
         pch = 20, cex = 0.5, col = util$c_mid)
}

phylo <- c('ABCO', 'ABMA', 'ABLA', 'TSME', 'LALY', 'PCEN', 
                 'PILO', 'PIAR', 'PIAL', 'PIFL', 'PISF', 'PICO', 'PIJE', 'PIAZ', 
                 'SEGI', 'CADE', 'THPL', 'CANO', 'JUOC', 'JUOS', 'JUSC')
phylo_order <- match(phylo, data$uniq_species_ids)

par(mfrow = c(1,1), mar = c(4,5,1,1))
for(s in 1:data$N_species){
  newname <- paste0('avg_omega_conc[', s,']')
  base_samples[[newname]] <- 0
  
  trees <- which(data$species_idxs %in% s)
  stand_species <- unique(data$stand_species_idxs[trees])
  
  for(stsp in stand_species){
    
    name1 <-  paste0('omega_conc_sck[', stsp,']')
    
    base_samples[[newname]] <- base_samples[[newname]] + (base_samples[[name1]])
  }
  base_samples[[newname]] <- base_samples[[newname]]/length(stand_species)
}
util$plot_disc_pushforward_quantiles_sign(base_samples, paste0('avg_omega_conc[',phylo_order,']'), display_ylim = c(-0.075,1),
                                          xticklabs =  data$uniq_species_ids[phylo_order], 
                                          ylab = "Average p(tree shock | stand concordance)", 
                                          xlab = '', ignore_sign = TRUE)
for(s in data$uniq_species_ids[phylo_order]){
  text(x = which(data$uniq_species_ids[phylo_order] == s), y = -0.07, 
       labels = length(unique(data$stand_idxs[data$species_idxs == which(data$uniq_species_ids == s)])), 
       xpd = NA, cex = 0.75, col = 'grey30')
}
text(x = -0.5, y = -0.07, labels = '# stands:', 
     xpd = NA, cex = 0.75, col = 'grey30')


par(mfrow = c(1,1), mar = c(4,5,1,1))
for(s in 1:data$N_species){
  newname <- paste0('avg_concshock[', s,']')
  base_samples[[newname]] <- 0
  
  trees <- which(data$species_idxs %in% s)
  stand_species <- unique(data$stand_species_idxs[trees])
  
  for(stsp in stand_species){
    
    name1 <-  paste0('omega_conc_sck[', stsp,']')
    st <- unique(data$stand_idxs[data$stand_species_idxs == stsp])
    name2 <-  paste0('phi_sck[', st,']')
    
    base_samples[[newname]] <- base_samples[[newname]] + (base_samples[[name1]]*base_samples[[name2]])
  }
  base_samples[[newname]] <- base_samples[[newname]]/length(stand_species)
}
util$plot_disc_pushforward_quantiles_sign(base_samples, paste0('avg_concshock[',phylo_order,']'), display_ylim = c(0,0.4),
                                          xticklabs = data$uniq_species_ids[phylo_order], 
                                          ylab = "Average p(stand concordance) × \np(tree shock | stand concordance)", 
                                          xlab = '', ignore_sign = TRUE)

for(s in 1:data$N_species){
  newname <- paste0('avg_thetha_idio[2][', s,']')
  base_samples[[newname]] <- 0
  
  trees <- which(data$species_idxs %in% s)
  for(t in which(data$species_idxs %in% s)){
    name <- paste0('thetas_idio[',t,',2]')
    base_samples[[newname]] <- base_samples[[newname]] + base_samples[[name]]
  }
  base_samples[[newname]] <- base_samples[[newname]]/length(trees)
}
util$plot_disc_pushforward_quantiles_sign(base_samples, paste0('avg_thetha_idio[2][',phylo_order,']'), display_ylim = c(0,0.4),
                                          xticklabs = data$uniq_species_ids[phylo_order], ylab = "Average p(idiosyncratic shock)", 
                                          xlab = '', ignore_sign = TRUE)



par(mfrow = c(2,1), mar = c(2,4,1,1))
sp <- which(data$uniq_species_ids == 'JUOS')
st <- unique(data$stand_idxs[which(data$species_idxs == sp)])
util$plot_disc_pushforward_quantiles(base_samples, paste0('phi_sck[',st,']'), display_ylim = c(0,1),
                                     ylab = "p(stand concordance)", xlab = 'Stands x species')
text(x = 0.5, y = 0.9, labels = 'Stands with TSME', adj = 0)
stsp <- unique(data$stand_species_idxs[which(data$species_idxs == sp)])
util$plot_disc_pushforward_quantiles(base_samples, paste0('omega_conc_sck[',stsp,']'), display_ylim = c(0,1),
                                          ylab = "p(tree shock | stand concordance)", xlab = 'Stands x species')



