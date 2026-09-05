

par(mfrow = c(1,1), mar = c(4,4,1,1))
plot(y = NULL, x = NULL, xlab = 'Winter prec. (mm)', ylab = 'log(ring width)', bty = 'n',
     xlim = c(0, 4000), ylim =  c(-2,2))

species <- which(data$uniq_species_ids %in% c('JUOC', 'PIED'))
# species <- 1:data$N_species
for(s in species){
  
  probs <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
  gdd0 = 10
  pre0 = 5
  vpd0 = 23
  
  muolds <- c()
  mus <- c()
  
  stands <- unique(data$stand_idxs[data$species_idxs == s])
   
  idxs <- c(sapply(stands, function(st) seq(1+(st-1)*data$N_all_years, st*data$N_all_years, 1)))
  pres <- data$pre_obs[idxs]
  
  g <- mean(data$gdd_obs[idxs])
  v <- mean(data$vpd_obs[idxs])
  
  pres <- seq(min(pres), max(pres), length.out = 40)
  for(p in pres){
    muold <- old_samples[['mu_alpha']] +
      old_samples[[paste0('beta_gdd[',s,']')]] * (g - gdd0) +
      old_samples[[paste0('beta_pre[',s,']')]] * (p - pre0) +
      old_samples[[paste0('beta_vpd[',s,']')]] * (v - vpd0) 
    muolds <- rbind(muolds, util$ensemble_mcmc_quantile_est(muold, probs))
    
    mu <- new_samples[['mu_alpha']] +
      new_samples[[paste0('beta_gdd[',s,']')]] * (g - gdd0)/10 +
      new_samples[[paste0('beta_pre[',s,']')]] * (p - pre0)/10 +
      new_samples[[paste0('beta_pre2[',s,']')]] * ((p - pre0)/10)^2 +
      new_samples[[paste0('beta_vpd[',s,']')]] * (v - vpd0)/10 + 
      new_samples[[paste0('beta_vpd2[',s,']')]] * ((v - vpd0)/10)^2 +
      2 * new_samples[[paste0('beta_prexvpd[',s,']')]] * ((v - vpd0)/10) * ((p - pre0)/10)
    mus <- rbind(mus, util$ensemble_mcmc_quantile_est(mu, probs))
  }
  
  pres <- pres*100
  # polygon(c(pres, rev(pres)), c(muolds[,'5%'], rev(muolds[,'95%'])),
  #         col = adjustcolor(util$c_light_teal, 0.3), border = NA)
  # lines(x = pres, y = muolds[,'95%'], col = util$c_mid_teal, lty = 2)
  # lines(x = pres, y = muolds[,'5%'], col = util$c_mid_teal, lty = 2)
  # lines(x = pres, y = muolds[,'50%'], col = util$c_mid_teal)
  
  polygon(c(pres, rev(pres)), c(mus[,'5%'], rev(mus[,'95%'])),
          col = adjustcolor(util$c_light, 0.3), border = NA)
  # lines(x = pres, y = mus[,'95%'], col = util$c_mid, lty = 2)
  # lines(x = pres, y = mus[,'5%'], col = util$c_mid, lty = 2)
  lines(x = pres, y = mus[,'50%'], col = util$c_mid)
  
  text(x = mean(pres), y = max(mus[,'95%']) + 0.1, labels = data$uniq_species_ids[s],
       col = util$c_dark)
  
  
}
