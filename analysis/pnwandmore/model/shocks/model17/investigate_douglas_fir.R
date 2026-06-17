rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(cmdstanr)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)


out_dir <- '/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model'


data <- readRDS(file.path(out_dir,'data_22apr2026_PSME_except2_19202024.rds'))
datasets <- readRDS(file.path(out_dir,'datasets_22apr2026_PSME_except2_19202024.rds'))

log_rw_pred_samples <- readRDS(file.path(out_dir, 'model17', 'treesamples_logrwpred_model17_HSGP_PSME_except2_pi0.rds'))
f_samples <- readRDS(file.path(out_dir, 'model17', 'treesamples_f_model17_HSGP_PSME_except2_pi0.rds'))
find_samples <- readRDS(file.path(out_dir, 'model17', 'treesamples_find_model17_HSGP_PSME_except2_pi0.rds'))
delta_sck_samples <- readRDS(file.path(out_dir, 'model17', 'treesamples_deltasck_model17_HSGP_PSME_except2_pi0.rds'))
fsh_samples <- readRDS(file.path(out_dir, 'model17', 'standsamples_fsh_model17_HSGP_PSME_except2_pi0.rds'))


# Plot one tree! And separate chains
par(mfcol = c(5,2), cex.lab = 0.8, cex.axis = 0.8, mar = c(1.5,4,2.5,1), mgp = c(1.5, 0.2, 0),tck = -0.03)

t <- sample(1:data$N_trees, 1)

#t=902, 903
#t = 2629
# t <- 1182, 3,4... missing shocks (stand GP?)
# t <- 1149
# t <- which(data$uniq_tree_ids == sample(data$uniq_tree_ids[which(data$stand_idxs == 87)], 1))

chains <- 1
idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
names <- paste0("log_rw_pred[",idxs,"]")
temp <- lapply(log_rw_pred_samples[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Log ring width (per mm)", 
                                     display_ylim=c(-5, 2), display_xlim = range(data$all_years),
                                     main = paste0('Chain ', chains, collapse = ' '))
points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=1, col="white")
points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=0.5, col="black")
text(x = 1920, y = 1.5, labels = data$uniq_tree_ids[t], cex = 0.8, adj = 0, col = "#474747")
abline(v = data$years[idxs[which(data$rw_obs[idxs] == 0)]], col = 'grey80', lty = 2)

names <- paste0("f[",idxs,"]")
temp <- lapply(f_samples[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Short-term GP", 
                                     display_ylim=c(-5, 2), display_xlim = range(data$all_years))

names <- paste0("f_ind[",idxs,"]")
temp <- lapply(find_samples[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Short-term GP", 
                                     display_ylim=c(-5, 2), display_xlim = range(data$all_years))

names <- paste0("delta_sck[",idxs,"]")
temp <- lapply(delta_sck_samples[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Shock amplitude", 
                                     display_ylim=c(-5, 2), display_xlim = range(data$all_years))
abline(v = data$years[idxs[which(data$rw_obs[idxs] == 0)]], col = 'grey80', lty = 2)

stand <- data$stand_idxs[t]
stand_idxs <- data$all_years_idxs[idxs]
names <- paste0("f_sh[",stand, ',', stand_idxs,"]")
temp <- lapply(fsh_samples[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Stand-level GP", 
                                     display_ylim=c(-5, 2), display_xlim = range(data$all_years))
abline(v = data$years[idxs[which(data$rw_obs[idxs] == 0)]], col = 'grey80', lty = 2)


chains <- 2:4
idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
names <- paste0("log_rw_pred[",idxs,"]")
temp <- lapply(log_rw_pred_samples[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Log ring width (per mm)", 
                                     display_ylim=c(-5, 2), display_xlim = range(data$all_years),
                                     main = paste0('Chain ', chains, collapse = ' '))
points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=1, col="white")
points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=0.5, col="black")
text(x = 1920, y = 1.5, labels = data$uniq_tree_ids[t], cex = 0.8, adj = 0, col = "#474747")
abline(v = data$years[idxs[which(data$rw_obs[idxs] == 0)]], col = 'grey80', lty = 2)

names <- paste0("f[",idxs,"]")
temp <- lapply(f_samples[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Short-term GP", 
                                     display_ylim=c(-5, 2), display_xlim = range(data$all_years))

names <- paste0("f_ind[",idxs,"]")
temp <- lapply(find_samples[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Short-term GP", 
                                     display_ylim=c(-5, 2), display_xlim = range(data$all_years))

names <- paste0("delta_sck[",idxs,"]")
temp <- lapply(delta_sck_samples[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Shock amplitude", 
                                     display_ylim=c(-5, 2), display_xlim = range(data$all_years))
abline(v = data$years[idxs[which(data$rw_obs[idxs] == 0)]], col = 'grey80', lty = 2)

stand <- data$stand_idxs[t]
stand_idxs <- data$all_years_idxs[idxs]
names <- paste0("f_sh[",stand, ',', stand_idxs,"]")
temp <- lapply(fsh_samples[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Stand-level GP", 
                                     display_ylim=c(-5, 2), display_xlim = range(data$all_years))
abline(v = data$years[idxs[which(data$rw_obs[idxs] == 0)]], col = 'grey80', lty = 2)
