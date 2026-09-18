rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(cmdstanr)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('mcmc_custom_functions.R', local = util)
setwd(wd)

datasets <- readRDS(file.path(wd, 'output/model', 'datasets_25august2026_25species_376stands_18962024.rds'))
data <- readRDS(file.path(wd, 'output/model', 'data_25august2026_25species_376stands_18962024.rds'))

csv_files <- list.files(file.path(wd, 'output/model/model21bis', 'tmp/since1896'), full.names = T) 
fit <- cmdstanr::read_cmdstan_csv(csv_files)
gc()
params <- c("logit_phi_sck", "beta_phi_vpd", "beta_phi_pre")
base_samples <- fit$post_warmup_draws[, , grepl(paste(params, collapse = "|"), dimnames(fit$post_warmup_draws)[[3]]), drop = FALSE]
names <- dimnames(base_samples)$variable
base_samples <- lapply(1:dim(base_samples)[3],
                       function(k){t(matrix(base_samples[1:dim(base_samples)[1],1:dim(base_samples)[2],k],
                                            nrow = dim(base_samples)[1], ncol = dim(base_samples)[2]))})
names(base_samples) <- names
rm(fit);gc()


pre0 <- 5
vpd0 <- 23

all_years_incfut <- c(data$all_years, 2025:2100)

# Observations (PRISM)
for(i in 1:data$N_all_years){
  base_samples[[paste0('mean_phi[',i,']')]] <- 0
}
for(s in 1:data$N_stands){
  
  years <- seq(1+data$N_all_years*(s-1), data$N_all_years*s, 1)
  for(i in 1:length(years)){
    
    logit_phi <- base_samples[[paste0('logit_phi_sck[', s, ']')]]+
      (data$vpd_obs[years[i]]-vpd0)*base_samples[['beta_phi_vpd']]+
      (data$pre_obs[years[i]]-pre0)*base_samples[['beta_phi_pre']]
    phi <- boot::inv.logit(logit_phi)
    
    base_samples[[paste0('mean_phi[',i,']')]] <- base_samples[[paste0('mean_phi[',i,']')]] + phi/data$N_stands
    
  }
}

# par(mfrow = c(1,1), mar = c(2,4,1,1))
# names <- paste0('mean_phi[',which(data$all_years %in% hist_years),']')
# util$plot_conn_pushforward_quantiles(base_samples, names, hist_years,
#                                      ylab = 'Avg. p(concordant event)',
#                                      display_ylim = c(0.05, 0.22))
# 
# before2000 <- mean(sapply(which(data$all_years %in% 1951:1999), 
#                           function(y) util$ensemble_mcmc_quantile_est(base_samples[[paste0('mean_phi[',y,']')]], c(0.5))))
# after2000 <- mean(sapply(which(data$all_years %in% 2000:2014), 
#                          function(y) util$ensemble_mcmc_quantile_est(base_samples[[paste0('mean_phi[',y,']')]], c(0.5))))
# segments(x0 = 1951, x1 = 1999, y0 = before2000, lty = 2, col = util$c_dark)
# segments(x0 = 2000, x1 = 2013, y0 = after2000, lty = 2, col = util$c_dark)
# text(x = 1950, y = 0.21, labels = 'PRISM data', adj = 0)


# NEX projections: compare historical period
nex_df <- readRDS(file.path(wd, 'output', 'climate', 'nex', 'climprojections_17sept2026.rds'))
gcms <- unique(nex_df$gcm)

nex_samples <- list()
hist_years <- which(data$all_years %in% c(1951:2014))
for(i in hist_years){
  nex_samples[[paste0('mean_phi[',i,']')]] <- NULL
}

for(m in gcms){
  
  gcm_samples <- list()
  for(i in hist_years){
    gcm_samples[[paste0('mean_phi[',i,']')]] <- 0
  }
  
  for(s in 1:data$N_stands){
    
    for(i in hist_years){
      clim_y <- nex_df[nex_df$grouped_stand == data$uniq_stand_ids[s] & 
                         nex_df$year == data$all_years[i] &
                         nex_df$gcm == m, ]
      
      logit_phi <- base_samples[[paste0('logit_phi_sck[', s, ']')]]+
        (clim_y$vpd_mjja-vpd0)*base_samples[['beta_phi_vpd']]+
        (clim_y$pre_ndjfma/100-pre0)*base_samples[['beta_phi_pre']]
      
      phi <- boot::inv.logit(logit_phi)
      gcm_samples[[paste0('mean_phi[',i,']')]] <- gcm_samples[[paste0('mean_phi[',i,']')]] + phi/data$N_stands
    }
  }
  
  for(i in hist_years){
    nex_samples[[paste0('mean_phi[',i,']')]] <- cbind(nex_samples[[paste0('mean_phi[',i,']')]], 
                                                      gcm_samples[[paste0('mean_phi[',i,']')]])
  }
}
nex_samples_hist <- nex_samples



nex_samples <- list()
fut_years <- which(all_years_incfut %in% c(2015:2100))
for(i in fut_years){
  nex_samples[[paste0('mean_phi[',i,']')]] <- NULL
}

for(m in gcms){
  
  gcm_samples <- list()
  for(i in fut_years){
    gcm_samples[[paste0('mean_phi[',i,']')]] <- 0
  }
  
  for(s in 1:data$N_stands){
    
    for(i in fut_years){
      clim_y <- nex_df[nex_df$grouped_stand == data$uniq_stand_ids[s] & 
                         nex_df$year == all_years_incfut[i] &
                         nex_df$gcm == m, ]
      
      logit_phi <- base_samples[[paste0('logit_phi_sck[', s, ']')]]+
        (clim_y$vpd_mjja-vpd0)*base_samples[['beta_phi_vpd']]+
        (clim_y$pre_ndjfma/100-pre0)*base_samples[['beta_phi_pre']]
      
      phi <- boot::inv.logit(logit_phi)
      gcm_samples[[paste0('mean_phi[',i,']')]] <- gcm_samples[[paste0('mean_phi[',i,']')]] + phi/data$N_stands
    }
  }
  
  for(i in fut_years){
    nex_samples[[paste0('mean_phi[',i,']')]] <- cbind(nex_samples[[paste0('mean_phi[',i,']')]], 
                                                      gcm_samples[[paste0('mean_phi[',i,']')]])
  }
}
nex_samples_fut <- nex_samples



# par(mfrow = c(1,1), mar = c(2,4,1,1))
# names <- paste0('mean_phi[',hist_years,']')
# util$plot_conn_pushforward_quantiles(nex_samples, names, data$all_years[hist_years],
#                                      ylab = 'Avg. p(concordant event)', 
#                                      display_ylim = c(0.05, 0.22))
# 
# before2000 <- mean(sapply(which(data$all_years %in% 1951:1999), 
#                           function(y) util$ensemble_mcmc_quantile_est(ipsl_samples[[paste0('mean_phi[',y,']')]], c(0.5))))
# after2000 <- mean(sapply(which(data$all_years %in% 2000:2013), 
#                          function(y) util$ensemble_mcmc_quantile_est(ipsl_samples[[paste0('mean_phi[',y,']')]], c(0.5))))
# segments(x0 = 1951, x1 = 1999, y0 = before2000, lty = 2, col = util$c_dark)
# segments(x0 = 2000, x1 = 2013, y0 = after2000, lty = 2, col = util$c_dark)
# text(x = 1950, y = 0.21, labels = 'IPSL projections', adj = 0)



names_obs <- paste0('mean_phi[',which(all_years_incfut %in% 1896:2024),']')
names_hist <- paste0('mean_phi[',which(all_years_incfut %in% 1951:2014),']')
names_fut <- paste0('mean_phi[',which(all_years_incfut %in% 2015:2100),']')
plot_conn_quantiles_projvsobs(
  base_samples, nex_samples_hist, nex_samples_fut,
  names_obs, names_hist, names_fut,
  data$all_years[which(all_years_incfut %in% 1896:2024)],
  data$all_years[which(all_years_incfut %in% 1951:2014)],
  all_years_incfut[which(all_years_incfut %in% 2015:2100)],
  ylab = 'Average p(stand-level extreme event)', 
  display_ylim = c(0, 0.3), display_xlim = c(1896, 2100)
)
legend(lwd = 10, col = c('#D9D9D990', "#AFC7E080", "#f6932070"), 
       legend = c('', '', ''),
       x = 1895, y = 0.3, 
       bg = NA, box.col = NA, cex = 0.9)
legend(lwd = 1.6, col = c('#333333', "#4C78A8", "#f69320"), 
       legend = c('PRISM', 'NEX-DCP30-CMIP6 (Historical)', 'NEX-DCP30-CMIP6 (SSP2-4.5)'),
       x = 1895, y = 0.3, box.col = NA, bg = NA, cex = 0.9)





# Compare both pre- and post-2000
# base_samples[[('mean_phi_bef2000')]] <- 0
# years <-  which(data$all_years %in% 1951:1999)
# for(i in years){
#   base_samples[[('mean_phi_bef2000')]] <-  base_samples[[('mean_phi_bef2000')]] + 
#     base_samples[[paste0('mean_phi[',i,']')]]/length(years)
# }
# base_samples[[('mean_phi_aft2000')]] <- 0
# years <- which(data$all_years %in% 2000:2013)
# for(i in years){
#   base_samples[[('mean_phi_aft2000')]] <-  base_samples[[('mean_phi_aft2000')]] + 
#     base_samples[[paste0('mean_phi[',i,']')]]/length(years)
# }
# nex_samples[[paste0('mean_phi_bef2000')]] <- 0
# years <-  which(data$all_years %in% 1951:1999)
# for(i in years){
#   nex_samples[[('mean_phi_bef2000')]] <-  nex_samples[[('mean_phi_bef2000')]] + 
#     nex_samples[[paste0('mean_phi[',i,']')]]/length(years)
# }
# nex_samples[[('mean_phi_aft2000')]] <- 0
# years <- which(data$all_years %in% 2000:2013)
# for(i in years){
#   nex_samples[[('mean_phi_aft2000')]] <-  nex_samples[[('mean_phi_aft2000')]] + 
#     nex_samples[[paste0('mean_phi[',i,']')]]/length(years)
# }
# 
# par(mfrow = c(1,2), mar = c(4,1,1,1), cex.lab = 0.95, cex.axis = 0.9)
# util$plot_expectand_pushforward(base_samples[[('mean_phi_bef2000')]], 50, flim = c(0.05,0.17), col = "#333333", 
#                                 display_name = 'Avg. p(concordant event), before 2000', ylim = c(0,170))
# util$plot_expectand_pushforward(nex_samples[[('mean_phi_bef2000')]], 50, flim = c(0.05,0.17), add = T, col =  "#4C78A8")
# util$plot_expectand_pushforward(base_samples[[('mean_phi_aft2000')]], 50, flim = c(0.05,0.17), col = "#333333", 
#                                 display_name = 'Avg. p(concordant event), after 2000', ylim = c(0,170))
# util$plot_expectand_pushforward(nex_samples[[('mean_phi_aft2000')]], 50, flim = c(0.05,0.17), add = T, col =  "#4C78A8")
# 
# # Compare climate directly
# prism_df <- data.frame()
# for(s in 1:data$N_stands){
#   years <- seq(1+data$N_all_years*(s-1), data$N_all_years*s, 1)
#   
#   prism_df <- rbind(prism_df,
#                     data.frame(stand = s, year = hist_years, 
#                                vpd_obs = data$vpd_obs[years][which(data$all_years %in% hist_years)], 
#                                pre_obs = data$pre_obs[years][which(data$all_years %in% hist_years)]*100))
# }
# 
# plot(x = NULL, y = NULL, xlim = c(1951,2014), ylim = c(0,52), bty = 'n',
#      ylab = 'VPD (hPa)', cex.axis = 0.85, cex.lab = 0.9)
# vpd_q <- aggregate(vpd_mjja ~ year, data = nex_df, FUN = quantile, c(0.01, 0.5, 0.99))
# vpd_q <- vpd_q[order(vpd_q$year),]
# polygon(x = c(vpd_q$year, rev(vpd_q$year)), y = c(vpd_q$vpd_mjja[,'1%'], rev(vpd_q$vpd_mjja[,'99%'])),
#         col = "#DCBCBC80", border = NA)
# 
# 
# vpdobs_q <- aggregate(vpd_obs ~ year, data = prism_df, FUN = quantile, c(0.01, 0.5, 0.99))
# vpdobs_q <- vpdobs_q[order(vpdobs_q$year),]
# polygon(x = c(vpdobs_q$year, rev(vpdobs_q$year)), y = c(vpdobs_q$vpd_obs[,'1%'], rev(vpdobs_q$vpd_obs[,'99%'])),
#         col = "#6B8E8E80", border = NA)
# lines(x = vpd_q$year, y = vpd_q$vpd_mjja[,'50%'], col =  "#B97C7C", lwd = 2)
# lines(x = vpdobs_q$year, y = vpdobs_q$vpd_obs[,'50%'], col =  "#6B8E8E", lwd = 2)
# legend(lwd = 2, col = c('#6B8E8E', "#C79999"), legend = c('PRISM', 'IPSL-CM6A-LR'),
#        x = 1950, y = 50, box.col = NA, cex = 0.9)
# 
# 
# plot(x = NULL, y = NULL, xlim = c(1951,2014), ylim = c(0,4000), bty = 'n',
#      ylab = 'Precipitation (mm)', cex.axis = 0.85, cex.lab = 0.9)
# pre_q <- aggregate(pre_ndjfma ~ year, data = nex_df, FUN = quantile, c(0.01, 0.5, 0.99))
# pre_q <- pre_q[order(pre_q$year),]
# polygon(x = c(pre_q$year, rev(pre_q$year)), y = c(pre_q$pre_ndjfma[,'1%'], rev(pre_q$pre_ndjfma[,'99%'])),
#         col = "#DCBCBC80", border = NA)
# 
# preobs_q <- aggregate(pre_obs ~ year, data = prism_df, FUN = quantile, c(0.01, 0.5, 0.99))
# preobs_q <- preobs_q[order(preobs_q$year),]
# polygon(x = c(preobs_q$year, rev(preobs_q$year)), y = c(preobs_q$pre_obs[,'1%'], rev(preobs_q$pre_obs[,'99%'])),
#         col = "#6B8E8E80", border = NA)
# lines(x = pre_q$year, y = pre_q$pre_ndjfma[,'50%'], col =  "#B97C7C", lwd = 2)
# lines(x = preobs_q$year, y = preobs_q$pre_obs[,'50%'], col =  "#6B8E8E", lwd = 2)
# legend(lwd = 2, col = c('#6B8E8E', "#C79999"), legend = c('PRISM', 'IPSL-CM6A-LR'),
#        x = 1950, y = 4000, box.col = NA, cex = 0.9)
# 
# 
# # Compare climate only on 2000-2013, but by sites
# vpd_rmse <- data.frame()
# for(s in 1:data$N_stands){
#   years <- seq(1+data$N_all_years*(s-1), data$N_all_years*s, 1)
#   proj <- nex_df[nex_df$grouped_stand == data$uniq_stand_ids[s],]
#   
#   vpd_obs <- data$vpd_obs[years][which(data$all_years %in% c(2000:2013))]
#   vpd_proj <- proj[proj$year %in% c(2000:2013), 'vpd_mjja']
#   
#   vpd_rmse <- rbind(vpd_rmse,
#                     data.frame(stand = s, lat = data$uniq_stand_lat[s], lon = data$uniq_stand_lon[s],
#                                vpd_rmse = ModelMetrics::rmse(vpd_obs, vpd_proj),
#                                delta = mean(vpd_obs-vpd_proj)))
# }
# 
# extent <- ext(c(-126, -102, 30,49))
# us <- crop(geodata::gadm(country = "USA",level = 1, resolution = 1, path = file.path(wd, 'data')), extent)
# 
# par(mfrow = c(1,1))
# plot(us, ext = extent, clip = FALSE, lwd = 0.5, 
#      box = FALSE, buffer = FALSE, axes = FALSE,
#      mar = c(0,2,0.5,0), col = 'grey90', border = 'white')
# 
# cols <- hcl.colors(100, "YlOrRd", rev = TRUE)
# z <- pmax(1, pmin(6, vpd_rmse$vpd_rmse))
# points(vpd_rmse$lon, vpd_rmse$lat,
#        col = cols[1 + floor((z - 1) / 5 * 99)],
#        pch = 20)
# 
# legend(x = -105, y = 48,
#        legend = seq(1, 6, 1),
#        pt.bg = cols[1 + (seq(0, 5) / 5 * 99)],
#        col = 'black',
#        pch = 21,
#        title = "RMSE",
#        cex = 0.8)



