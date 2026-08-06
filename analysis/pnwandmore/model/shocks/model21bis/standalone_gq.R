library(cmdstanr)

wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('mcmc_custom_functions.R', local = util)
setwd(wd)

mod_gq <- cmdstan_model(file.path(wd, 'model/stan/model21bis/hgsp/old', 'model21bis_HSGP_multispecies_nofind_pooled_climateshocks_withGQ.stan'))

folder <- '/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model'
data <- readRDS(file.path(folder,'data_10july2026_24species_365stands_19502022.rds'))
data$grainsize <-  1
data$uniq_tree_ids <- NULL
data$uniq_species_ids <- NULL
data$uniq_stand_ids <- NULL
data$N_clades <- 1
data$uniq_stand_lat <- NULL
data$uniq_stand_lon <- NULL

fit <- readRDS(file.path(wd, 'output/model/model21bis', 'fit_model21bis_HGSP_24species_365stands_withinit_7threads_climateshocks.rds'))
gc()
fit_gq <- mod_gq$generate_quantities(fit, data = data, seed = 1234567, parallel_chains = 4,
                                     output_dir = '/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model/model21bis/tmp')

gq_samples <- fit_gq$draws(variables = c('shock_change_noshutdown', 'shock_growth_change', 'shutdown_perc'))
names <- dimnames(gq_samples)$variable
gq_samples <- lapply(1:dim(gq_samples)[3],
                       function(k){t(matrix(gq_samples[1:dim(gq_samples)[1],1:dim(gq_samples)[2],k],
                                            nrow = dim(gq_samples)[1], ncol = dim(gq_samples)[2]))})
names(gq_samples) <- names
gc()


st <- 1
par(mfrow = c(3,1), mar = c(4,5.5,1,1), cex.axis = 0.9, cex.lab = 1.2)
util$plot_conn_pushforward_quantiles_missingrings(gq_samples, paste0('shock_change_noshutdown[', st,',',1:data$N_stand_years[st],']'), 
                                                  plot_xs = data$all_years[1:data$N_stand_years[st]],
                                                  display_ylim = c(-100,20),
                                                  ylab = "Average growth variation\ndue to concordant shocks (%)")
text(x = 1955, y = 15, labels = 'Without shutdowns', xpd = TRUE, srt = 0, font = 3, cex = 1.3)
abline(v = 1974, lty = 2)
abline(v = 1996, lty = 2)
abline(v = 2002, lty = 2)
util$plot_conn_pushforward_quantiles_missingrings(gq_samples, paste0('shock_growth_change[', st,',',1:data$N_stand_years[st],']'), 
                                                  plot_xs = data$all_years[1:data$N_stand_years[st]],
                                                  display_ylim = c(-101,20),
                                                  ylab = "Average growth variation\ndue to concordant shocks (%)")
text(x = 1955, y = 15, labels = 'With shutdowns', xpd = TRUE, srt = 0, font = 3, cex = 1.3)
abline(v = 1974, lty = 2)
abline(v = 1996, lty = 2)
abline(v = 2002, lty = 2)
util$plot_disc_pushforward_quantiles(gq_samples, paste0('shutdown_perc[',st,',',1:data$N_stand_years[st],']'), 
                                     display_ylim = c(0,100),
                                     ylab = "Concordant shutdown trees (%)")
abline(v = which(data$all_years == 2002), lty = 2)






trees <- which(data$stand_idxs == 2)
idxs <- c(sapply(trees, function(t) data$tree_start_idxs[t]:data$tree_end_idxs[t]))


tree_gq_samples <- fit_gq$draws(variables = c(paste0('log_rw_pred[', idxs, ']'), paste0('delta_sck[', idxs, ']')))
names <- dimnames(tree_gq_samples)$variable
tree_gq_samples <- lapply(1:dim(tree_gq_samples)[3],
                     function(k){t(matrix(tree_gq_samples[1:dim(tree_gq_samples)[1],1:dim(tree_gq_samples)[2],k],
                                          nrow = dim(tree_gq_samples)[1], ncol = dim(tree_gq_samples)[2]))})
names(tree_gq_samples) <- names

par(mfrow = c(2,1), mar = c(2,5,1,1))
t <- trees[11]
idxs_tree <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
names <- paste0('log_rw_pred[', idxs_tree, ']')
util$plot_conn_pushforward_quantiles_missingrings(tree_gq_samples, names, data$years[idxs_tree], ylab = 'log_rw_pred')
abline(v = data$all_years[22], lty = 2)
abline(v = data$all_years[47], lty = 2)

names <- paste0('delta_sck[', idxs_tree, ']')
util$plot_conn_pushforward_quantiles_missingrings(tree_gq_samples, names, data$years[idxs_tree], ylab = 'delta_shock\n(conc. + idio)')
abline(v = data$all_years[22], lty = 2)
abline(v = data$all_years[47], lty = 2)






