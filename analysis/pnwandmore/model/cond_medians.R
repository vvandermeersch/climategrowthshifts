rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnw"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)

# Model data
data <- readRDS(file.path(wd, 'output/model', 'data_28june2025.rds'))

# Posterior quantification - diagnostics
# fit <- readRDS(file.path(wd, 'output/model', 'fit_28june2025_nopooling.rds')) # run on Margot
# diagnostics <- util$extract_hmc_diagnostics(fit)
# util$check_all_hmc_diagnostics(diagnostics)
# samples <- util$extract_expectand_vals(fit)
# base_samples <- util$filter_expectands(samples,
#                                        c('log_rw_pred'),
#                                        check_arrays=TRUE)
# rm(fit, samples);gc(verbose=T)
# saveRDS(base_samples, file.path(wd, 'output/model', 'logrwpred_28june2025_nopooling.rds'))
base_samples <- readRDS(file.path(wd, 'output/model', 'logrwpred_28june2025_nopooling.rds')) # need much less memory, as GC is not very efficient in R


pdf(file = file.path(wd, 'figures', 'model28june_nopooling', 'conditional_medians.pdf'), width = 8.27, height = 11.69)
par(mfrow = c(3,1), mar = c(4,4,1,1))
pred_names <- sapply(1:data$N, function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conditional_median_quantiles(base_samples, pred_names, data$gdd_obs,
                                       2, 51, 1, data$log_rw_obs, 
                                       residual=FALSE,
                                       xlab="GDD (x100degC)")

util$plot_conditional_median_quantiles(base_samples, pred_names, data$sm_obs,
                                       9, 38, 1, data$log_rw_obs, 
                                       residual=FALSE,
                                       xlab="Soil moisture (%)")

util$plot_conditional_median_quantiles(base_samples, pred_names, data$vpd_obs,
                                       1, 30, 1, data$log_rw_obs, 
                                       residual=FALSE,
                                       xlab="VPD (hPa)")
dev.off()

pdf(file = file.path(wd, 'figures', 'model28june_nopooling', 'conditional_medians_residuals.pdf'), width = 8.27, height = 11.69)
par(mfrow = c(3,1), mar = c(4,4,1,1))
pred_names <- sapply(1:data$N, function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conditional_median_quantiles(base_samples, pred_names, data$gdd_obs,
                                       2, 51, 1, data$log_rw_obs, 
                                       residual=TRUE,
                                       xlab="GDD (x100degC)")

util$plot_conditional_median_quantiles(base_samples, pred_names, data$sm_obs,
                                       9, 38, 1, data$log_rw_obs, 
                                       residual=TRUE,
                                       xlab="Soil moisture (%)")

util$plot_conditional_median_quantiles(base_samples, pred_names, data$vpd_obs,
                                       1, 30, 1, data$log_rw_obs, 
                                       residual=TRUE,
                                       xlab="VPD (hPa)")
dev.off()
