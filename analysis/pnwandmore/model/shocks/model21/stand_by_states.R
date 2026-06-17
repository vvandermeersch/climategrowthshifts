par(mfrow = c(3,4), mar = c(2,4,1,1))

states <- unique(substr(data$uniq_tree_ids, 1,2))

for(st in states){
  
  init <- 1
  stands <- unique(data$stand_idxs[grepl(st, data$uniq_tree_ids)])
  
  
  for(s in stands){
    
    idxs <- seq(data$stand_start_years_idxs[s], length.out = data$N_stand_years[s])
    plot_xs <- data$all_years[idxs]
    names <- paste0("f_sh[",s,",",idxs,"]")
    
    
    # Construct quantiles for bin contents
    N <- length(names)
    probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
    
    calc <- function(n) {
      util$ensemble_mcmc_quantile_est(samples[[names[n]]], probs)
    }
    plot_quantiles <- sapply(1:N, calc)
    
    if(init == 1){
      plot(1, type="n", main= st,
           ylim=c(-2, 2), xlab='',
           xlim=range(data$all_years), ylab='Stand GP')
    }
    
    
    polygon(c(plot_xs, rev(plot_xs)),
            c(plot_quantiles[1,], rev(plot_quantiles[9,])),
            col = paste0(util$c_light,'30') , border = NA)
    polygon(c(plot_xs, rev(plot_xs)),
            c(plot_quantiles[2,], rev(plot_quantiles[8,])),
            col = paste0(util$c_light_highlight,'30'), border = NA)
    polygon(c(plot_xs, rev(plot_xs)),
            c(plot_quantiles[3,], rev(plot_quantiles[7,])),
            col = paste0(util$c_mid,'30'), border = NA)
    polygon(c(plot_xs, rev(plot_xs)),
            c(plot_quantiles[4,], rev(plot_quantiles[6,])),
            col = paste0(util$c_mid_highlight,'30'), border = NA)
    lines(plot_xs, plot_quantiles[5,], col=util$c_dark, lwd=2)
    
    init <- 0
  }
  
}





