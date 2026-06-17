rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(cmdstanr)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)


folder <- '/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model'

data <- readRDS(file.path(folder,'data_18may2026_PIPO_1stand_19502008.rds'))

fit <- readRDS(file.path(folder, 'model19', 'fit_model19_HSGP_PIPO_10stands_19502013_fullidio.rds'))
samples <- fit$draws(c('f_sh'))
names <- dimnames(samples)$variable
samples <- lapply(1:dim(samples)[3],
                  function(k){t(matrix(samples[1:dim(samples)[1],1:dim(samples)[2],k],
                                       nrow = dim(samples)[1], ncol = dim(samples)[2]))})
names(samples) <- names

treesamples <- readRDS(file.path(folder, 'model19', 'treesamples_model19_HSGP_PIPO_10stands_19502013_fullidio.rds'))

par(mfrow=c(2,3))
t <- 20
idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
names <- paste0("log_rw_pred[",idxs,"]")
util$plot_conn_pushforward_quantiles(treesamples, names, data$years[idxs],
                                     xlab="Year", ylab="Log ring width (per mm)", 
                                     display_ylim=c(-5, 2), display_xlim = range(data$all_years))
points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=1, col="white")
points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=0.5, col="black")
text(x = 1920, y = 1.5, labels = data$uniq_tree_ids[t], cex = 0.8, adj = 0, col = "#474747")

names <- paste0("delta_sck[",idxs,"]")
util$plot_conn_pushforward_quantiles(treesamples, names, data$years[idxs],
                                     xlab="Year", ylab="Shock amplitude", 
                                     display_ylim=c(-5, 2), display_xlim = range(data$all_years))

names <- paste0("shutdown[",idxs,"]")
util$plot_disc_pushforward_quantiles(treesamples, names, 
                                     xlab="Year", ylab="Shutdown", 
                                     display_ylim=c(0, 1))

names <- paste0("conc_state[",idxs,"]")
util$plot_disc_pushforward_quantiles(treesamples, names, 
                                     xlab="Year", ylab="Concordant state", 
                                     display_ylim=c(0, 1))

names <- paste0("idio_state[",idxs,"]")
util$plot_disc_pushforward_quantiles(treesamples, names, 
                                     xlab="Year", ylab="Idiosyncratic state", 
                                     display_ylim=c(0, 1))

names <- paste0("f_sh[1,",1:data$N_all_years,"]")
util$plot_conn_pushforward_quantiles(samples, names, data$years[idxs],
                                     xlab="Year", ylab="Stand-level GP", 
                                     display_ylim=c(-5, 1))
abline(v=c(2002,2007))

