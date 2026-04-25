rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(cmdstanr)
library(rnaturalearth)
library(terra)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)


folder <- '/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model/'

data <- readRDS(file.path(folder,'data_26mar2026_PIPOexcept1_19202024.rds'))
fit <- readRDS(file.path(folder, 'model17/fit_model17_HSGP_PIPOexcept1.rds'))

# ----------------- #
## 1. Make a map!
# ----------------- #

world <- vect(ne_countries(scale = "medium", returnclass = "sf"))
westam <- vect(ne_states(country = c("Canada", "Mexico", "United States of America"), returnclass = "sf"))
westam <- vect(ne_states(country = c("United States of America"), returnclass = "sf"))
westam <- crop(westam, ext(c(-125, -104, 31.8, 49)))
mtn <- crop(vect(file.path(wd, 'data/usgs/physio_shp', 'physio.shp')),
            ext(c(-125, -104, 31.8, 49)))

pdf('/home/victor/projects/climategrowthshifts/docs/misc/egu26/figures/testmap.pdf',
    height = 35 / 2.54, width = 35 / 2.54)
par(mfrow = c(1,1), cex.main = 0.9, mar = c(1,1,1,1))
plot(westam, xlim = c(-125, -104), ylim = c(31.8, 49), 
     xaxs = "i", yaxs = "i", axes = FALSE, 
     box = FALSE, clip = FALSE, col = '#eeeeee', border = 'white',
     lwd = 1.5)

plot(aggregate(buffer(mtn[mtn$PROVINCE %in% c('COLORADO PLATEAUS', 'SOUTHERN ROCKY MOUNTAINS')], width = 1)), col = '#cf644820', border = 'white', add = TRUE, lwd = 2.5)
plot(aggregate(buffer(mtn[mtn$PROVINCE %in% c('NORTHERN ROCKY MOUNTAINS', 'MIDDLE ROCKY MOUNTAINS') |
                            mtn$SECTION %in% c('NORTHERN CASCADE MOUNTAINS', 'MIDDLE CASCADE MOUNTAINS','SOUTHERN CASCADE MOUNTAINS')], width = 1)), 
     col = '#cf644820', border = 'white', add = TRUE, lwd = 2.5)
plot(aggregate(buffer(mtn[mtn$SECTION == 'SIERRA NEVADA'], width = 1)), col = '#cf644820', border = 'white', add = TRUE, lwd = 2.5)

plot(aggregate(buffer(mtn[mtn$SECTION %in% c('SONORAN DESERT', 'MEXICAN HIGHLAND','SACRAMENTO')], width = 1)), 
     col = '#cf644820', border = 'white', add = TRUE, lwd = 2.5)

lines(aggregate(westam), col = 'white', lwd = 10)
points(x = data$uniq_stand_lon, y = data$uniq_stand_lat,
       pch = 19, col = 'white', cex = 2.5)
points(x = data$uniq_stand_lon, y = data$uniq_stand_lat,
       pch = 19, col = '#cf6448', cex = 1.7)
dev.off()


# --------------------------------------------------- #
## 2. Show one (non-shocky) tree to explain the model
# --------------------------------------------------- #


par(mfrow = c(1,1), cex.main = 0.9, mar = c(2,4,1,1))
t <- 70
idxs <- c(data$tree_start_idxs[t]+24):data$tree_end_idxs[t]
names <- paste0("delta_sck[",idxs,"]")
gq_samples <- fit_gq$draws(names)
names <- dimnames(gq_samples)$variable
gq_samples <- lapply(1:dim(gq_samples)[3],
                  function(k){t(matrix(gq_samples[1:dim(gq_samples)[1],1:dim(gq_samples)[2],k],
                                       nrow = dim(gq_samples)[1], ncol = dim(gq_samples)[2]))})
names(gq_samples) <- names
gc()
util$plot_conn_pushforward_quantiles(gq_samples, names, data$years[idxs],
                                     xlab="Year", ylab="Log ring width (per mm)", 
                                     display_ylim=c(-10, 2), display_xlim = range(data$years[idxs]))

t <- 1
idxs <- c(data$tree_start_idxs[t]+24):data$tree_end_idxs[t]
pdf('/home/victor/projects/climategrowthshifts/docs/misc/egu26/figures/f.pdf',
    height = 10 / 2.54, width = 26 / 2.54)
par(mfrow = c(1,1), cex.main = 0.9, mar = c(2.2,5,0.8,0.8), cex.axis = 1.8, cex.lab = 2)
names <- paste0("alpha_f[",idxs,"]")
gq_samples <- fit_gq$draws(names)
names <- dimnames(gq_samples)$variable
gq_samples <- lapply(1:dim(gq_samples)[3],
                     function(k){t(matrix(gq_samples[1:dim(gq_samples)[1],1:dim(gq_samples)[2],k],
                                          nrow = dim(gq_samples)[1], ncol = dim(gq_samples)[2]))})
names(gq_samples) <- names
gc()
util$plot_conn_pushforward_quantiles(gq_samples, names, data$years[idxs],
                                     xlab="Year", ylab="log(ring width)", 
                                     display_ylim=c(-1, 2), display_xlim = range(data$years[idxs]))
dev.off()



pdf('/home/victor/projects/climategrowthshifts/docs/misc/egu26/figures/find.pdf',
    height = 10 / 2.54, width = 26 / 2.54)
par(mfrow = c(1,1), cex.main = 0.9, mar = c(2.2,5,0.8,0.8), cex.axis = 1.8, cex.lab = 2)
names <- paste0("f_ind[",idxs,"]")
gq_samples <- fit_gq$draws(names)
names <- dimnames(gq_samples)$variable
gq_samples <- lapply(1:dim(gq_samples)[3],
                     function(k){t(matrix(gq_samples[1:dim(gq_samples)[1],1:dim(gq_samples)[2],k],
                                          nrow = dim(gq_samples)[1], ncol = dim(gq_samples)[2]))})
names(gq_samples) <- names
gc()
util$plot_conn_pushforward_quantiles(gq_samples, names, data$years[idxs],
                                     xlab="Year", ylab="log(ring width)", 
                                     display_ylim=c(-1, 2), display_xlim = range(data$years[idxs]))
dev.off()



pdf('/home/victor/projects/climategrowthshifts/docs/misc/egu26/figures/clim.pdf',
    height = 10 / 2.54, width = 26 / 2.54)
par(mfrow = c(1,1), cex.main = 0.9, mar = c(2.2,5,0.8,0.8), cex.axis = 1.8, cex.lab = 2)
names <- paste0("alpha_clim[",idxs,"]")
gq_samples <- fit_gq$draws(names)
names <- dimnames(gq_samples)$variable
gq_samples <- lapply(1:dim(gq_samples)[3],
                     function(k){t(matrix(gq_samples[1:dim(gq_samples)[1],1:dim(gq_samples)[2],k],
                                          nrow = dim(gq_samples)[1], ncol = dim(gq_samples)[2]))})
names(gq_samples) <- names
gc()
util$plot_conn_pushforward_quantiles(gq_samples, names, data$years[idxs],
                                     xlab="Year", ylab="log(ring width)", 
                                     display_ylim=c(-1, 2), display_xlim = range(data$years[idxs]))
dev.off()

pdf('/home/victor/projects/climategrowthshifts/docs/misc/egu26/figures/pred.pdf',
    height = 10 / 2.54, width = 26 / 2.54)
par(mfrow = c(1,1), cex.main = 0.9, mar = c(2.2,5,0.8,0.8), cex.axis = 1.8, cex.lab = 2)
names <- paste0("log_rw_pred[",idxs,"]")
gq_samples <- fit_gq$draws(names)
names <- dimnames(gq_samples)$variable
gq_samples <- lapply(1:dim(gq_samples)[3],
                     function(k){t(matrix(gq_samples[1:dim(gq_samples)[1],1:dim(gq_samples)[2],k],
                                          nrow = dim(gq_samples)[1], ncol = dim(gq_samples)[2]))})
names(gq_samples) <- names
gc()
util$plot_conn_pushforward_quantiles(gq_samples, names, data$years[idxs],
                                     xlab="Year", ylab="log(ring width)", 
                                     display_ylim=c(-1, 2), display_xlim = range(data$years[idxs]))
points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=1.3, col="white")
points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=0.8, col="black")
dev.off()



# Relative reduction of growth due to climate predictors
gq_samples <- readRDS(file.path(folder, 'model17/gq_samples_model17_HSGP_PIPOexcept1.rds'))
data <- readRDS(file.path(folder,'data_26mar2026_PIPOexcept1_19202024.rds'))

for(i in 1:data$N_all_years){
  gq_samples[[paste0('avg_exp_m1_clim[',i,']')]]<- 0
}

for(s in 1:data$N_stands){
  idxs <- seq(1 + data$N_all_years * (s-1), data$N_all_years * s, 1)
  idxs <- idxs[data$stand_start_years_idxs[s]:data$N_stand_years[s]]
  for(i in data$stand_start_years_idxs[s]:data$N_stand_years[s]){
    gq_samples[[paste0('avg_exp_m1_clim[',i,']')]] <- gq_samples[[paste0('avg_exp_m1_clim[',i,']')]] +
      gq_samples[[paste0('avg_exp_m1_clim_stdlvl[',idxs[i],']')]]/data$N_stands
  }
}

for(i in 1:data$N_all_years){
  gq_samples[[paste0('avg_exp_m1_clim[',i,']')]] <- gq_samples[[paste0('avg_exp_m1_clim[',i,']')]]*100
}

pdf('/home/victor/projects/climategrowthshifts/docs/misc/egu26/figures/climvar.pdf',
    height = 14 / 2.54, width = 35 / 2.54)
par(mfrow = c(1,1), cex.main = 0.9, mar = c(2.2,6.5,0.8,0.8), cex.axis = 1.8, cex.lab = 2)
names <- paste0('avg_exp_m1_clim[',1:data$N_all_years,']')
util$plot_disc_pushforward_quantiles(gq_samples, names, 
                                     ylab = 'Average relative change in growth\ndue to climate predictors (%)', 
                                     display_ylim = c(-20, 20),
                                     xticklabs = ifelse(data$all_years[1:data$N_all_years]%in% seq(1920,2000,20), data$all_years[1:data$N_all_years], ''))
abline(h = 0, lty = 2, lwd = 1.5, col = 'grey30')
dev.off()

years  <- 1:data$N_all_years
periods <- ifelse(years %in% which(data$all_years %in% 1930:1940), '30s',
                  ifelse(years %in% which(data$all_years %in% 1980:1989), '80s',
                         ifelse(years %in% which(data$all_years %in% 1990:1999), '90s', 
                                ifelse(years %in% which(data$all_years >= 2000), 'post2000', NA))))

for(p in na.omit(unique(periods))){
  print(p)
  gq_samples[[paste0('avg_exp_m1_clim_period[',p,']')]] <- 0
  print(gq_samples[[paste0('avg_exp_m1_clim_period[',p,']')]] )
}

for(i in 1:data$N_all_years){
  if(is.na(periods[i])){next}
  gq_samples[[paste0('avg_exp_m1_clim_period[',periods[i],']')]] <- gq_samples[[paste0('avg_exp_m1_clim_period[',periods[i],']')]] +
    gq_samples[[paste0('avg_exp_m1_clim[',i,']')]]/sum(na.omit(periods == periods[20]))
}

util$plot_expectand_pushforward(gq_samples[[paste0('avg_exp_m1_clim_period[','30s',']')]], 200,
                                flim = c(-5,10))
util$plot_expectand_pushforward(gq_samples[[paste0('avg_exp_m1_clim_period[','80s',']')]], 200,
                                flim = c(-5,10), add = T, col = "#6CD883")
util$plot_expectand_pushforward(gq_samples[[paste0('avg_exp_m1_clim_period[','90s',']')]], 200,
                                flim = c(-5,10), add = T)
util$plot_expectand_pushforward(gq_samples[[paste0('avg_exp_m1_clim_period[','post2000',']')]], 200,
                                flim = c(-5,10), add = T)






# Relative reduction of growth due to climate predictors

for(i in 1:data$N_all_years){
  gq_samples[[paste0('avg_exp_m1_delta_sck[',i,']')]]<- 0
}

for(s in 1:data$N_stands){
  idxs <- seq(1 + data$N_all_years * (s-1), data$N_all_years * s, 1)
  idxs <- idxs[data$stand_start_years_idxs[s]:data$N_stand_years[s]]
  for(i in data$stand_start_years_idxs[s]:data$N_stand_years[s]){
    gq_samples[[paste0('avg_exp_m1_delta_sck[',i,']')]] <- gq_samples[[paste0('avg_exp_m1_delta_sck[',i,']')]] +
      gq_samples[[paste0('avg_exp_m1_delta_sck_stdlvl[',idxs[i],']')]]/data$N_stands
  }
}

for(i in 1:data$N_all_years){
  gq_samples[[paste0('avg_exp_m1_delta_sck[',i,']')]] <- gq_samples[[paste0('avg_exp_m1_delta_sck[',i,']')]]*100
}


pdf('/home/victor/projects/climategrowthshifts/docs/misc/egu26/figures/shockvar.pdf',
    height = 11 / 2.54, width = 35 / 2.54)
par(mfrow = c(1,1), cex.main = 0.9, mar = c(2.2,6.5,0.8,0.8), cex.axis = 1.8, cex.lab = 2)
names <- paste0('avg_exp_m1_delta_sck[',1:data$N_all_years,']')
util$plot_disc_pushforward_quantiles(gq_samples, names, 
                                     ylab = 'Average relative change in growth\ndue extreme years (%)', 
                                     display_ylim = c(-20, 0.05),
                                     xticklabs = ifelse(data$all_years[1:data$N_all_years]%in% seq(1920,2000,20), data$all_years[1:data$N_all_years], ''))
abline(h = 0, lty = 2, lwd = 1.5, col = 'grey30')
dev.off()


# What the model's missing

pdf('/home/victor/projects/climategrowthshifts/docs/misc/egu26/figures/pred2_noshock.pdf',
    height = 10 / 2.54, width = 26 / 2.54)
t <- 70
idxs <- c(data$tree_start_idxs[t]+24):data$tree_end_idxs[t]
par(mfrow = c(1,1), cex.main = 0.9, mar = c(2.2,5,0.8,0.8), cex.axis = 1.8, cex.lab = 2)
names <- c(paste0("log_rw_pred_noshock[",idxs,"]"), paste0("log_rw_pred[",idxs,"]"))
gq_samples <- fit_gq$draws(names)
names <- dimnames(gq_samples)$variable
gq_samples <- lapply(1:dim(gq_samples)[3],
                     function(k){t(matrix(gq_samples[1:dim(gq_samples)[1],1:dim(gq_samples)[2],k],
                                          nrow = dim(gq_samples)[1], ncol = dim(gq_samples)[2]))})
names(gq_samples) <- names
gc()
names <- paste0("log_rw_pred_noshock[",idxs,"]")
util$plot_conn_pushforward_quantiles(gq_samples, names, data$years[idxs],
                                     xlab="Year", ylab="log(ring width)", 
                                     display_ylim=c(-4, 2), display_xlim = range(data$years[idxs]))
points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=1.3, col="white")
points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=0.8, col="black")
# names <- paste0("log_rw_pred[",idxs,"]")
# util$plot_conn_pushforward_quantiles(gq_samples, names, data$years[idxs],
#                                      xlab="Year", ylab="log(ring width)", 
#                                      display_ylim=c(-4, 2), display_xlim = range(data$years[idxs]))
# points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=1.3, col="white")
# points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=0.8, col="black")
dev.off()


pdf('/home/victor/projects/climategrowthshifts/docs/misc/egu26/figures/pred2_shock.pdf',
    height = 10 / 2.54, width = 26 / 2.54)
t <- 70
idxs <- c(data$tree_start_idxs[t]+24):data$tree_end_idxs[t]
par(mfrow = c(1,1), cex.main = 0.9, mar = c(2.2,5,0.8,0.8), cex.axis = 1.8, cex.lab = 2)
names <- c(paste0("log_rw_pred_noshock[",idxs,"]"), paste0("log_rw_pred[",idxs,"]"))
gq_samples <- fit_gq$draws(names)
names <- dimnames(gq_samples)$variable
gq_samples <- lapply(1:dim(gq_samples)[3],
                     function(k){t(matrix(gq_samples[1:dim(gq_samples)[1],1:dim(gq_samples)[2],k],
                                          nrow = dim(gq_samples)[1], ncol = dim(gq_samples)[2]))})
names(gq_samples) <- names
gc()
names <- paste0("log_rw_pred[",idxs,"]")
util$plot_conn_pushforward_quantiles(gq_samples, names, data$years[idxs],
                                     xlab="Year", ylab="log(ring width)", 
                                     display_ylim=c(-4, 2), display_xlim = range(data$years[idxs]))
points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=1.3, col="white")
points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=0.8, col="black")
# names <- paste0("log_rw_pred[",idxs,"]")
# util$plot_conn_pushforward_quantiles(gq_samples, names, data$years[idxs],
#                                      xlab="Year", ylab="log(ring width)", 
#                                      display_ylim=c(-4, 2), display_xlim = range(data$years[idxs]))
# points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=1.3, col="white")
# points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=0.8, col="black")
dev.off()






pdf('/home/victor/projects/climategrowthshifts/docs/misc/egu26/figures/pred2_t69.pdf',
    height = 4 / 2.54, width = 14 / 2.54)
t <- 69
idxs <- c(data$tree_start_idxs[t]+24):data$tree_end_idxs[t]
par(mfrow = c(1,1), cex.main = 0.9, mar = c(0.2,0.2,0.2,0.2), cex.axis = 1.8, cex.lab = 2)
names <- paste0("log_rw_pred[",idxs,"]")
gq_samples <- fit_gq$draws(names)
names <- dimnames(gq_samples)$variable
gq_samples <- lapply(1:dim(gq_samples)[3],
                     function(k){t(matrix(gq_samples[1:dim(gq_samples)[1],1:dim(gq_samples)[2],k],
                                          nrow = dim(gq_samples)[1], ncol = dim(gq_samples)[2]))})
names(gq_samples) <- names
gc()
names <- paste0("log_rw_pred[",idxs,"]")
util$plot_conn_pushforward_quantiles(gq_samples, names, data$years[idxs],
                                     xlab="Year", ylab="log(ring width)", 
                                     display_ylim=c(-8, 2), display_xlim = range(data$years[idxs]))
points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=0.6, col="white")
points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=0.4, col="black")
# names <- paste0("log_rw_pred[",idxs,"]")
# util$plot_conn_pushforward_quantiles(gq_samples, names, data$years[idxs],
#                                      xlab="Year", ylab="log(ring width)", 
#                                      display_ylim=c(-4, 2), display_xlim = range(data$years[idxs]))
# points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=1.3, col="white")
# points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=0.8, col="black")
dev.off()

