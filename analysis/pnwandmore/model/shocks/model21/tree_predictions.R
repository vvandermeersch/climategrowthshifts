rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(cmdstanr)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)


folder <- '/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model'

data <- readRDS(file.path(folder,'data_18may2026_PIPO_1stand_ut539_19502005.rds'))

fit <- readRDS(file.path(folder, 'model21', 'fit_model21_HGSP_PIPO_1stand_ut539_19502005_nomu_sametau.rds'))
samples <- fit$draws(c('f_sh', 'delta_clim'))
names <- dimnames(samples)$variable
samples <- lapply(1:dim(samples)[3],
                  function(k){t(matrix(samples[1:dim(samples)[1],1:dim(samples)[2],k],
                                       nrow = dim(samples)[1], ncol = dim(samples)[2]))})
names(samples) <- names

treesamples <- readRDS(file.path(folder, 'model21', 'treesamples_model21_HGSP_PIPO_1stand_ut539_19502005_nomu_sametau.rds'))

## Focus on shocks
chains <- c(1)
par(mfrow=c(2,4))
t <- 56
idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
names <- paste0("log_rw_pred[",idxs,"]")

util$plot_conn_pushforward_quantiles_missingrings(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Log ring width (per mm)", 
                                     display_ylim=c(-3, 2), display_xlim = range(data$all_years),
                                     main = c('Chain(s) ', paste0(chains, collapse = ", ")))
points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=1, col="white")
points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=0.5, col="black")
text(x = 1950, y = 2, labels = data$uniq_tree_ids[t], cex = 0.8, adj = 0, col = "#474747")

names <- paste0("f_sh[",data$stand_idxs[t],",",which(data$all_years %in% data$years[idxs]),"]")
temp <- lapply(samples[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Stand GP", 
                                     display_ylim=c(-3, 2), display_xlim = range(data$all_years))

names <- paste0("regr_clim[",data$stand_idxs[t],",",which(data$all_years %in% data$years[idxs]),"]")
temp <- lapply(treesamples[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Explained climate", 
                                     display_ylim=c(-3, 2), display_xlim = range(data$all_years))

names <- paste0("delta_clim[",data$stand_idxs[t],",",which(data$all_years %in% data$years[idxs]),"]")
temp <- lapply(samples[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Missing climate", 
                                     display_ylim=c(-3, 2), display_xlim = range(data$all_years))



names <- paste0("tree_conc_state[",idxs,"]")
temp <- lapply(treesamples[names], function(m) m[chains, , drop = FALSE])
util$plot_disc_pushforward_quantiles(temp, names, 
                                     xlab="Year", ylab="Concordant state", 
                                     display_ylim=c(0, 1))

names <- paste0("tree_idio_state[",idxs,"]")
temp <- lapply(treesamples[names], function(m) m[chains, , drop = FALSE])
util$plot_disc_pushforward_quantiles(temp, names, 
                                     xlab="Year", ylab="Idiosyncratic state", 
                                     display_ylim=c(0, 1))

names <- paste0("delta_sck[",idxs,"]")
temp <- lapply(treesamples[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Shock amplitude", 
                                     display_ylim=c(-3, 2), display_xlim = range(data$all_years))

names <- paste0("shutdown[",idxs,"]")
temp <- lapply(treesamples[names], function(m) m[chains, , drop = FALSE])
util$plot_disc_pushforward_quantiles(temp, names, 
                                     xlab="Year", ylab="Shutdown state", 
                                     display_ylim=c(0, 1))


## Compare climatic components
par(mfrow=c(2,1), mar = c(2,4,1,1))
s <- 3
names <- paste0("regr_clim[",s,",",data$stand_start_years_idxs[s]:data$N_stand_years[s],"]")
temp <- lapply(treesamples[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Explained climate", 
                                     display_ylim=c(-1, 2), display_xlim = c(1980,1995))
abline(v = c(1985, 1989), lty = 2)

names <- paste0("delta_clim[",s,",",data$stand_start_years_idxs[s]:data$N_stand_years[s],"]")
temp <- lapply(samples[names], function(m) m[chains, , drop = FALSE])
util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                     xlab="Year", ylab="Missing climate", 
                                     display_ylim=c(-1, 2), display_xlim = c(1980,1995))
abline(v = c(1985, 1989), lty = 2)

## Full model

pdf(file.path(wd, 'model/shocks/model21', 'trees_for_Lizzie.pdf'), width = 11, height = 8.5)

chains <- c(1:4)
par(mfrow=c(3,4), mar = c(2,4,1,1), oma = c(1,1,1,1))


for(t in c(10,12,14,16)){
  
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- paste0("log_rw_pred[",idxs,"]")
  temp <- lapply(treesamples[names], function(m) m[chains, , drop = FALSE])
  util$plot_conn_pushforward_quantiles_missingrings(temp, names, data$years[idxs],
                                                    xlab="Year", ylab="Log ring width (per mm)", 
                                                    display_ylim=c(-3, 2), display_xlim = range(data$all_years),
                                                    # main = c('Chain(s) ', paste0(chains, collapse = ", "))
  )
  points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=1, col="white")
  points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=0.5, col="black")
  text(x = 1950, y = 2, labels = data$uniq_tree_ids[t], cex = 0.8, adj = 0, col = "#474747")
  
  plot.new()
  
  names <- paste0("f[",idxs,"]")
  temp <- lapply(treesamples[names], function(m) m[chains, , drop = FALSE])
  util$plot_conn_pushforward_quantiles_missingrings(temp, names, data$years[idxs],
                                                    xlab="Year", ylab="Long-term GP (allometry)", 
                                                    display_ylim=c(-3, 2), display_xlim = range(data$all_years))
  
  names <- paste0("f_ind[",idxs,"]")
  temp <- lapply(treesamples[names], function(m) m[chains, , drop = FALSE])
  util$plot_conn_pushforward_quantiles_missingrings(temp, names, data$years[idxs],
                                                    xlab="Year", ylab="Short-term GP (microenv., canopy)", 
                                                    display_ylim=c(-3, 2), display_xlim = range(data$all_years))
  
  
  plot.new()
  
  names <- paste0("f_sh[",data$stand_idxs[t],",",which(data$all_years %in% data$years[idxs]),"]")
  temp <- lapply(samples[names], function(m) m[chains, , drop = FALSE])
  util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                       xlab="Year", ylab="Stand GP (macroenv.)", 
                                       display_ylim=c(-3, 2), display_xlim = range(data$all_years))
  
  names <- paste0("regr_clim[",data$stand_idxs[t],",",which(data$all_years %in% data$years[idxs]),"]")
  temp <- lapply(treesamples[names], function(m) m[chains, , drop = FALSE])
  util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                       xlab="Year", ylab="Explained climate", 
                                       display_ylim=c(-3, 2), display_xlim = range(data$all_years))
  
  names <- paste0("delta_clim[",data$stand_idxs[t],",",which(data$all_years %in% data$years[idxs]),"]")
  temp <- lapply(samples[names], function(m) m[chains, , drop = FALSE])
  util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                       xlab="Year", ylab="Missing macroenvironment", 
                                       display_ylim=c(-3, 2), display_xlim = range(data$all_years))
  
  names <- paste0("tree_conc_state[",idxs,"]")
  temp <- lapply(treesamples[names], function(m) m[chains, , drop = FALSE])
  util$plot_disc_pushforward_quantiles(temp, names, 
                                       xlab="Year", ylab="Concordant state", 
                                       display_ylim=c(0, 1))
  
  names <- paste0("tree_idio_state[",idxs,"]")
  temp <- lapply(treesamples[names], function(m) m[chains, , drop = FALSE])
  util$plot_disc_pushforward_quantiles(temp, names, 
                                       xlab="Year", ylab="Idiosyncratic state", 
                                       display_ylim=c(0, 1))
  
  names <- paste0("delta_sck[",idxs,"]")
  temp <- lapply(treesamples[names], function(m) m[chains, , drop = FALSE])
  util$plot_conn_pushforward_quantiles(temp, names, data$years[idxs],
                                       xlab="Year", ylab="Shock amplitude", 
                                       display_ylim=c(-3, 2), display_xlim = range(data$all_years))
  
  names <- paste0("shutdown[",idxs,"]")
  temp <- lapply(treesamples[names], function(m) m[chains, , drop = FALSE])
  util$plot_disc_pushforward_quantiles(temp, names, 
                                       xlab="Year", ylab="Shutdown state", 
                                       display_ylim=c(0, 1))
  
}

dev.off()