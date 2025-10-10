rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('functions/custom_functions.R', local=util)
setwd(wd)

data <- readRDS(file.path(wd, 'output/model', 'data_24sept2025.rds'))
# fit <- readRDS(file.path(wd, "output/model/fit_24sept2025_partialpooling_2clades_centered_extended.rds"))
# 
pred_names <- sapply(1:data$N, function(n) paste0('log_rw_pred[', n, ']'))
# predictions <- rstan::extract(fit, pars = pred_names, permuted=FALSE)
# saveRDS(predictions, file.path(wd, "output/model/logrwpred_24sept2025_partialpooling_2clades_centered_extended.rds"))
predictions <- readRDS(file.path(wd, "output/model/logrwpred_24sept2025_partialpooling_2clades_centered_extended.rds"))


# # Below, from Mike's function
N <- length(pred_names)
predictions <- lapply(1:N, function(n) t(predictions[,,n]))
names(predictions) <- pred_names

in2mm <-25.4 # scale factor to convert inches to mm
pdf(file.path(wd, 'figures', 'newyork2025', 'posterior_predictive_behaviors', 'allspecies_against_gddobs_14bins.pdf'),
    width=200/in2mm,height=300/in2mm,paper="special")
par(mfrow=c(2, 1))
pred_names <- sapply(1:data$N, function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conditional_median_quantiles(predictions, pred_names, data$gdd_obs,
                                       0, 72, 5, data$log_rw_obs, 
                                       xlab="GDD (x100 degC)")

pred_names <- sapply(1:data$N, function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conditional_median_quantiles(predictions, pred_names, data$gdd_obs,
                                       0, 72, 2, data$log_rw_obs, 
                                       residual=TRUE, xlab="GDD (x100 degC)")
dev.off()

pdf(file.path(wd, 'figures', 'newyork2025', 'posterior_predictive_behaviors', 'allspecies_against_gddobs_36bins.pdf'),
    width=200/in2mm,height=300/in2mm,paper="special")
par(mfrow=c(2, 1))
pred_names <- sapply(1:data$N, function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conditional_median_quantiles(predictions, pred_names, data$gdd_obs,
                                       0, 72, 2, data$log_rw_obs, 
                                       xlab="GDD (x100 degC)")

pred_names <- sapply(1:data$N, function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conditional_median_quantiles(predictions, pred_names, data$gdd_obs,
                                       0, 72, 2, data$log_rw_obs, 
                                       residual=TRUE, xlab="GDD (x100 degC)")
dev.off()


species <- "PSME"
focustrees <- which(data$species_idxs == which(data$uniq_species_ids == species))
species_obs_idx <- unlist(sapply(focustrees, function(t) return(data$tree_start_idxs[t]:data$tree_end_idxs[t])))
pdf(file.path(wd, 'figures', 'newyork2025', 'posterior_predictive_behaviors', 'PSME_against_gddobs_36bins.pdf'),
    width=200/in2mm,height=300/in2mm,paper="special")
par(mfrow=c(2, 1))
pred_names <- sapply(species_obs_idx, function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conditional_median_quantiles(predictions, pred_names, data$gdd_obs[species_obs_idx],
                                       0, 72, 2, data$log_rw_obs[species_obs_idx], 
                                       xlab="VPD (hPa)")

pred_names <- sapply(species_obs_idx, function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conditional_median_quantiles(predictions, pred_names, data$gdd_obs[species_obs_idx],
                                       0, 72, 2, data$log_rw_obs[species_obs_idx], 
                                       residual=TRUE, xlab="VPD (hPa)")
dev.off()
