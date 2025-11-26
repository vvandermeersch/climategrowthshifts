rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('functions/custom_functions.R', local=util)
setwd(wd)

data <- readRDS(file.path(wd, 'output/model/climatena', 'data_15oct2025_long_recentPIPOinAZ.rds'))
fit <- readRDS(file.path(wd, "model/shocks/output", "fit_15oct2025_long_recentPIPOinAZ.rds"))
samples <-  util$extract_expectand_vals(fit)

par(mfrow = c(8,2), mar = c(3,4.5,1,1), cex.lab = 1.2)
for(stand in 1:data$N_stands){
  names <- sapply(1:data$N_all_years,
                  function(sp) paste0('f_tilde_sh[', stand, ',', sp, ']'))
  util$plot_realizations(samples, names, plot_xs = data$all_years, N_plots = 50,
                         main = "f_tilde", ylab = 'Function outputs',
                         display_xlim = c(1900,2020), , display_ylim = c(-5,5))
  names <- sapply(1:data$N_all_years,
                  function(sp) paste0('f_sh[', stand, ',', sp, ']'))
  util$plot_realizations(samples, names, plot_xs = data$all_years, N_plots = 50,
                         main = "f", ylab = '',
                         display_xlim = c(1900,2020), display_ylim = c(-8,4))
  if(stand == 1){
    abline(v = c(1902), lty = 2, col = 'grey70', lwd = 2)
  }
  if(stand %in% c(1, 2,3,4,5,7,8)){
    abline(v = c(2002), lty = 2, col = 'grey70', lwd = 2)
  }
  # if(stand == 2){
  #   abline(v = c(1920), lty = 2, col = 'grey70', lwd = 2)
  # }
  if(stand %in% c(5,6)){
    abline(v = c(1977), lty = 2, col = 'grey70', lwd = 2)
  }
  if(stand == 5){
    abline(v = c(1904, 1959, 2004, 2012), lty = 2, col = 'grey70', lwd = 2)
  }
  if(stand == 4){
    abline(v = c(1977), lty = 2, col = 'grey70', lwd = 2)
  }
  
}
