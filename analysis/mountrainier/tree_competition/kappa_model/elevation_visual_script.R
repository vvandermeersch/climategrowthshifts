
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('confronting_with_real_data.R')

k_samples <- lapply(samplefit, function(fit){
  samples <- util$extract_expectand_vals(fit)
  samples <- util$filter_expectands(samples, c('k'))
  return(samples[['k']])
})

k_samples3 <- lapply(samplefit3, function(fit3){
  samples <- util$extract_expectand_vals(fit3)
  samples <- util$filter_expectands(samples, c('k'))
  return(samples[['k']])
})
names(k_samples) <- paste0('k', 1:length(samplefit))
names(k_samples3) <- paste0('k', 1:length(samplefit3))

util$plot_disc_pushforward_quantiles(k_samples, names(k_samples))

folder <- 'D:/ubc_study/udergrad_research/temporal ecology lab/climategrowthshifts/analysis/mountrainier/data/treerings_ailene'
elevation <- read.csv(file.path(folder, 'tree_plot_climate_temp.csv'))
elevation <- elevation[,c('Plot', 'Elevation')]
stand_names <- c('AG05','AX15', 'AV06', 'TB13','TO04')
stand_elevations <- elevation[which(elevation$Plot %in% stand_names), 'Elevation']
order(stand_elevations)

par(mfrow=c(1,2))
par(cex.lab = 0.8, cex.axis = 0.8, mar = c(4,4,1,1))
util$plot_disc_pushforward_quantiles(k_samples, paste0('k', order(stand_elevations)),
                                     xticklabs = stand_names[order(stand_elevations)],ylab="kappa(tshe)",
                                     main="2008")
mtext(text = paste0(round(stand_elevations[order(stand_elevations)],0),'m'), side = 1,
      at = 1:length(stand_names), line = 1.7, cex = 0.7, col = 'grey50')

util$plot_disc_pushforward_quantiles(k_samples3, paste0('k', order(stand_elevations)),
                                     xticklabs = stand_names[order(stand_elevations)],ylab="kappa(tshe)",
                                     main="avg(2006,2007,2008)")
mtext(text = paste0(round(stand_elevations[order(stand_elevations)],0),'m'), side = 1,
      at = 1:length(stand_names), line = 1.7, cex = 0.7, col = 'grey50')


beta_samples <- lapply(samplefit, function(fit){
  samples <- util$extract_expectand_vals(fit)
  samples <- util$filter_expectands(samples, c('beta[1]'))
  return(samples[['beta[1]']])
})
names(beta_samples) <- paste0('beta', 1:length(samplefit))

beta_samples3 <- lapply(samplefit3, function(fit3){
  samples <- util$extract_expectand_vals(fit3)
  samples <- util$filter_expectands(samples, c('beta[1]'))
  return(samples[['beta[1]']])
})
names(beta_samples3) <- paste0('beta', 1:length(samplefit3))

par(cex.lab = 0.8, cex.axis = 0.8, cex.main = 0.9)
util$plot_disc_pushforward_quantiles(beta_samples, paste0('beta', order(stand_elevations)),
                                     xticklabs = stand_names[order(stand_elevations)],
                                     main = '2008', ylab = 'beta(tshe)')
mtext(text = paste0(round(stand_elevations[order(stand_elevations)],0),'m'), side = 1,
      at = 1:length(stand_names), line = 1.7, cex = 0.7, col = 'grey50')
util$plot_disc_pushforward_quantiles(beta_samples3, paste0('beta', order(stand_elevations)),
                                     xticklabs = stand_names[order(stand_elevations)],
                                     main = 'avg(2006,2007,2008)', ylab = 'beta(tshe)')
mtext(text = paste0(round(stand_elevations[order(stand_elevations)],0),'m'), side = 1,
      at = 1:length(stand_names), line = 1.7, cex = 0.7, col = 'grey50')

gamma_samples <- lapply(samplefit, function(fit){
  samples <- util$extract_expectand_vals(fit)
  samples <- util$filter_expectands(samples, c('r[1]'))
  return(samples[['r[1]']])
})
names(gamma_samples) <- paste0('gamma', 1:length(samplefit))

gamma_samples3 <- lapply(samplefit3, function(fit3){
  samples <- util$extract_expectand_vals(fit3)
  samples <- util$filter_expectands(samples, c('r[1]'))
  return(samples[['r[1]']])
})
names(gamma_samples3) <- paste0('gamma', 1:length(samplefit3))

par(cex.lab = 0.8, cex.axis = 0.8, cex.main = 0.9)
util$plot_disc_pushforward_quantiles(gamma_samples, paste0('gamma', order(stand_elevations)),
                                     xticklabs = stand_names[order(stand_elevations)],
                                     main = '2008', ylab = 'gamma(tshe)')
mtext(text = paste0(round(stand_elevations[order(stand_elevations)],0),'m'), side = 1,
      at = 1:length(stand_names), line = 1.7, cex = 0.7, col = 'grey50')
util$plot_disc_pushforward_quantiles(gamma_samples3, paste0('gamma', order(stand_elevations)),
                                     xticklabs = stand_names[order(stand_elevations)],
                                     main = 'avg(2006,2007,2008)', ylab = 'gamma(tshe)')
mtext(text = paste0(round(stand_elevations[order(stand_elevations)],0),'m'), side = 1,
      at = 1:length(stand_names), line = 1.7, cex = 0.7, col = 'grey50')
