

plots_altitude <- plots[order(plots$Elevation),c('Plot', 'Elevation')]

pdf(file.path(wd, 'figures', 'growthGDD_acrosselevation.pdf'), height = 26)

par(mfrow=c(9, 1))
for(i in 1:9){
  
  p <- plots_altitude[i,]
  id <- which(uniq_stand_ids == p$Plot)
  
  nmin <- data$tree_start_idxs[min(which(data$stand_idxs == id))]
  nmax <- data$tree_end_idxs[max(which(data$stand_idxs == id))]
  
  pred_names <- sapply(nmin:nmax, function(n) paste0('log_rw_pred[', n, ']'))
  
  util$plot_conditional_median_quantiles(samples, pred_names, data$gdd_obs[nmin:nmax],
                                         250/10, 1800/10, 5, data$log_rw_obs[nmin:nmax], 
                                         display_ylim = c(-0.7, 1.2),
                                         ylab= "Marginal quantiles",
                                         xlab="GDD (x10 degC)", main = paste0(p$Plot, ', altitude = ', round(p$Elevation), 'm'))
  
  
}
dev.off()

