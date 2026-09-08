rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(cmdstanr)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('mcmc_custom_functions.R', local = util)
setwd(wd)

datasets <- readRDS(file.path(wd, 'output/model', 'datasets_10july2026_24species_365stands_19502022.rds'))
data <- readRDS(file.path(wd, 'output/model', 'data_10july2026_24species_365stands_19502022.rds'))
fit <- readRDS(file.path(wd, 'output/model/model21bis', 'fit_model21bis_HGSP_24species_365stands_withinit_4threads_climateshocks.rds'))

params <- c(
  "mu_alpha", "sigma_alpha", "alpha",
  "sigma_alpha_stand", "alpha_stand",
  "mu_beta_gdd", "sigma_beta_gdd", "beta_gdd",
  "mu_beta_pre", "sigma_beta_pre", "beta_pre",
  "mu_beta_vpd", "sigma_beta_vpd", "beta_vpd",
  "delta_clim",
  "tau_clim", # "kappa_clim_free", "kappa_clim",
  'log_kappa_clim', 'kappa_clim',
  "f_tilde_ind",
  "mu_log_rho", "sigma_log_rho", "rho_merged",
  "mu_log_gamma", "sigma_log_gamma", "log_gamma_merged", 'gamma_merged',
  "mu_log_tau_conc", "sigma_log_tau_conc", "tau_conc",
  "mu_phi", "sigma_phi", "logit_phi_sck", "phi_sck", "beta_phi_vpd", "beta_phi_pre",
  "mu_omega_conc", "sigma_omega_conc", "logit_omega_conc_sck", "omega_conc_sck",
  "mu_omega_shutdown", "sigma_omega_shutdown", "logit_omega_shutdown", "omega_shutdown",
  "thetas_idio", "tau_idio",
  # 'thetas_baseline', 'omega_thetas', 'thetas_idio',
  "sigma"
)

base_samples <- fit$draws(params)
names <- dimnames(base_samples)$variable
base_samples <- lapply(1:dim(base_samples)[3],
                       function(k){t(matrix(base_samples[1:dim(base_samples)[1],1:dim(base_samples)[2],k],
                                            nrow = dim(base_samples)[1], ncol = dim(base_samples)[2]))})
names(base_samples) <- names

pre0 <- 5
vpd0 <- 23
hist_years <- 1951:2013

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

par(mfrow = c(1,1), mar = c(2,4,1,1))
names <- paste0('mean_phi[',which(data$all_years %in% hist_years),']')
util$plot_conn_pushforward_quantiles(base_samples, names, hist_years,
                                     ylab = 'Avg. p(concordant event)',
                                     display_ylim = c(0.05, 0.22))

before2000 <- mean(sapply(which(data$all_years %in% 1951:1999), 
                          function(y) util$ensemble_mcmc_quantile_est(base_samples[[paste0('mean_phi[',y,']')]], c(0.5))))
after2000 <- mean(sapply(which(data$all_years %in% 2000:2013), 
                         function(y) util$ensemble_mcmc_quantile_est(base_samples[[paste0('mean_phi[',y,']')]], c(0.5))))
segments(x0 = 1951, x1 = 1999, y0 = before2000, lty = 2, col = util$c_dark)
segments(x0 = 2000, x1 = 2013, y0 = after2000, lty = 2, col = util$c_dark)
text(x = 1950, y = 0.21, labels = 'PRISM data', adj = 0)


# NEX projections: compare historical period
nex_df <- readRDS(file.path(wd, 'output', 'climate', 'nex', 'climprojections_08sept2026.rds'))

ipsl_samples <- list()
for(i in 1:data$N_all_years){
  ipsl_samples[[paste0('mean_phi[',i,']')]] <- 0
}
for(s in 1:data$N_stands){
  
  for(y in hist_years){
    
    clim_y <- nex_df[nex_df$grouped_stand == data$uniq_stand_ids[s] &
                      nex_df$year == y, ]
    
    logit_phi <- base_samples[[paste0('logit_phi_sck[', s, ']')]]+
      (clim_y$vpd_mjja-vpd0)*base_samples[['beta_phi_vpd']]+
      (clim_y$pre_ndjfma/100-pre0)*base_samples[['beta_phi_pre']]
    phi <- boot::inv.logit(logit_phi)
    
    i <- which(data$all_years == y)
    ipsl_samples[[paste0('mean_phi[',i,']')]] <- ipsl_samples[[paste0('mean_phi[',i,']')]] + phi/data$N_stands
    
  }
}

par(mfrow = c(1,1), mar = c(2,4,1,1))
names <- paste0('mean_phi[',which(data$all_years %in% hist_years),']')
util$plot_conn_pushforward_quantiles(ipsl_samples, names, hist_years,
                                     ylab = 'Avg. p(concordant event)', 
                                     display_ylim = c(0.05, 0.22))

before2000 <- mean(sapply(which(data$all_years %in% 1951:1999), 
                          function(y) util$ensemble_mcmc_quantile_est(ipsl_samples[[paste0('mean_phi[',y,']')]], c(0.5))))
after2000 <- mean(sapply(which(data$all_years %in% 2000:2013), 
                         function(y) util$ensemble_mcmc_quantile_est(ipsl_samples[[paste0('mean_phi[',y,']')]], c(0.5))))
segments(x0 = 1951, x1 = 1999, y0 = before2000, lty = 2, col = util$c_dark)
segments(x0 = 2000, x1 = 2013, y0 = after2000, lty = 2, col = util$c_dark)
text(x = 1950, y = 0.21, labels = 'IPSL projections', adj = 0)


# Compare both pre- and post-2000
base_samples[[('mean_phi_bef2000')]] <- 0
years <-  which(data$all_years %in% 1951:1999)
for(i in years){
  base_samples[[('mean_phi_bef2000')]] <-  base_samples[[('mean_phi_bef2000')]] + 
    base_samples[[paste0('mean_phi[',i,']')]]/length(years)
}
base_samples[[('mean_phi_aft2000')]] <- 0
years <- which(data$all_years %in% 2000:2013)
for(i in years){
  base_samples[[('mean_phi_aft2000')]] <-  base_samples[[('mean_phi_aft2000')]] + 
    base_samples[[paste0('mean_phi[',i,']')]]/length(years)
}
ipsl_samples[[paste0('mean_phi_bef2000')]] <- 0
years <-  which(data$all_years %in% 1951:1999)
for(i in years){
  ipsl_samples[[('mean_phi_bef2000')]] <-  ipsl_samples[[('mean_phi_bef2000')]] + 
    ipsl_samples[[paste0('mean_phi[',i,']')]]/length(years)
}
ipsl_samples[[('mean_phi_aft2000')]] <- 0
years <- which(data$all_years %in% 2000:2013)
for(i in years){
  ipsl_samples[[('mean_phi_aft2000')]] <-  ipsl_samples[[('mean_phi_aft2000')]] + 
    ipsl_samples[[paste0('mean_phi[',i,']')]]/length(years)
}

par(mfrow = c(1,2), mar = c(4,1,1,1), cex.lab = 0.95, cex.axis = 0.9)
util$plot_expectand_pushforward(base_samples[[('mean_phi_bef2000')]], 30, flim = c(0.05,0.17), col = "#6B8E8E", 
                                display_name = 'Avg. p(concordant event), before 2000', ylim = c(0,150))
util$plot_expectand_pushforward(ipsl_samples[[('mean_phi_bef2000')]], 30, flim = c(0.05,0.17), add = T, col =  "#B97C7C")
util$plot_expectand_pushforward(base_samples[[('mean_phi_aft2000')]], 30, flim = c(0.05,0.17), col = "#6B8E8E", 
                                display_name = 'Avg. p(concordant event), after 2000', ylim = c(0,150))
util$plot_expectand_pushforward(ipsl_samples[[('mean_phi_aft2000')]], 30, flim = c(0.05,0.17), add = T, col =  "#B97C7C")

# Compare climate directly
prism_df <- data.frame()
for(s in 1:data$N_stands){
  years <- seq(1+data$N_all_years*(s-1), data$N_all_years*s, 1)
  
  prism_df <- rbind(prism_df,
                    data.frame(stand = s, year = hist_years, 
                               vpd_obs = data$vpd_obs[years][which(data$all_years %in% hist_years)], 
                               pre_obs = data$pre_obs[years][which(data$all_years %in% hist_years)]*100))
}

plot(x = NULL, y = NULL, xlim = c(1951,2014), ylim = c(0,52), bty = 'n',
     ylab = 'VPD (hPa)', cex.axis = 0.85, cex.lab = 0.9)
vpd_q <- aggregate(vpd_mjja ~ year, data = nex_df, FUN = quantile, c(0.01, 0.5, 0.99))
vpd_q <- vpd_q[order(vpd_q$year),]
polygon(x = c(vpd_q$year, rev(vpd_q$year)), y = c(vpd_q$vpd_mjja[,'1%'], rev(vpd_q$vpd_mjja[,'99%'])),
        col = "#DCBCBC80", border = NA)


vpdobs_q <- aggregate(vpd_obs ~ year, data = prism_df, FUN = quantile, c(0.01, 0.5, 0.99))
vpdobs_q <- vpdobs_q[order(vpdobs_q$year),]
polygon(x = c(vpdobs_q$year, rev(vpdobs_q$year)), y = c(vpdobs_q$vpd_obs[,'1%'], rev(vpdobs_q$vpd_obs[,'99%'])),
        col = "#6B8E8E80", border = NA)
lines(x = vpd_q$year, y = vpd_q$vpd_mjja[,'50%'], col =  "#B97C7C", lwd = 2)
lines(x = vpdobs_q$year, y = vpdobs_q$vpd_obs[,'50%'], col =  "#6B8E8E", lwd = 2)
legend(lwd = 2, col = c('#6B8E8E', "#C79999"), legend = c('PRISM', 'IPSL-CM6A-LR'),
       x = 1950, y = 50, box.col = NA, cex = 0.9)


plot(x = NULL, y = NULL, xlim = c(1951,2014), ylim = c(0,4000), bty = 'n',
     ylab = 'Precipitation (mm)', cex.axis = 0.85, cex.lab = 0.9)
pre_q <- aggregate(pre_ndjfma ~ year, data = nex_df, FUN = quantile, c(0.01, 0.5, 0.99))
pre_q <- pre_q[order(pre_q$year),]
polygon(x = c(pre_q$year, rev(pre_q$year)), y = c(pre_q$pre_ndjfma[,'1%'], rev(pre_q$pre_ndjfma[,'99%'])),
        col = "#DCBCBC80", border = NA)

preobs_q <- aggregate(pre_obs ~ year, data = prism_df, FUN = quantile, c(0.01, 0.5, 0.99))
preobs_q <- preobs_q[order(preobs_q$year),]
polygon(x = c(preobs_q$year, rev(preobs_q$year)), y = c(preobs_q$pre_obs[,'1%'], rev(preobs_q$pre_obs[,'99%'])),
        col = "#6B8E8E80", border = NA)
lines(x = pre_q$year, y = pre_q$pre_ndjfma[,'50%'], col =  "#B97C7C", lwd = 2)
lines(x = preobs_q$year, y = preobs_q$pre_obs[,'50%'], col =  "#6B8E8E", lwd = 2)
legend(lwd = 2, col = c('#6B8E8E', "#C79999"), legend = c('PRISM', 'IPSL-CM6A-LR'),
       x = 1950, y = 4000, box.col = NA, cex = 0.9)


# Compare climate only on 2000-2013, but by sites

