rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(cmdstanr)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('mcmc_custom_functions.R', local = util)
setwd(wd)

folder <- '/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model'

data <- readRDS(file.path(folder,'data_10july2026_24species_365stands_19502022.rds'))

states <- sapply(1:data$N_stands, function(s){
  unique(substr(data$uniq_tree_ids[data$stand_idxs == s], 1, 2))
})

fit <- readRDS(file.path(wd, 'output/model/model21bis', 'fit_model21bis_HGSP_24species_365stands_withinit_7threads_climateshocks.rds'))
fit$time()
fit$diagnostic_summary()

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

# base_samples <- readRDS(file.path(wd, 'output/model/model21bis', 'basesamples_model21bis_HGSP_11species_99stands_withinit.rds'))

params <- c(
  "mu_alpha", "sigma_alpha", "alpha",
  "sigma_alpha_stand", "alpha_stand",
  "mu_beta_gdd", "sigma_beta_gdd", "beta_gdd",
  "mu_beta_pre", "sigma_beta_pre", "beta_pre",
  "mu_beta_vpd", "sigma_beta_vpd", "beta_vpd",
  # "delta_clim",
  "tau_clim", 
  'log_kappa_clim', 'kappa_clim',
  # "f_tilde_ind",
  "mu_log_rho", "sigma_log_rho", "rho_merged",
  "mu_log_gamma", "sigma_log_gamma", "log_gamma_merged", 'gamma_merged',
  "mu_log_tau_conc", "sigma_log_tau_conc", "tau_conc",
  "mu_phi", "sigma_phi", "logit_phi_sck", "beta_phi_vpd", "beta_phi_pre",
  "mu_omega_conc", "sigma_omega_conc", "logit_omega_conc_sck",
  "mu_omega_shutdown", "sigma_omega_shutdown", "logit_omega_shutdown",
  "thetas_idio", "tau_idio",
  # 'thetas_baseline', 'omega_thetas', 'thetas_idio',
  "sigma"
)
# params <- c('delta_clim')

for(p in params){
  message(paste0('\n\n---------\nParameter(s):',p))
  rest <- util$filter_expectands(base_samples, p, 
                                 check_arrays = TRUE)
  util$check_all_expectand_diagnostics(rest)
}

par(mfrow = c(1,2), mar = c(4,4,1,1))
util$plot_expectand_pushforward(base_samples[['beta_phi_vpd']], 50,
                                display_name = bquote(beta[VPD]^phi),
                                flim = c(-0.5,0.5))
util$plot_expectand_pushforward(base_samples[['beta_phi_pre']], 50,
                                display_name = bquote(beta[pre]^phi),
                                flim = c(-0.5,0.5))

par(mfrow = c(1,1))
util$plot_disc_pushforward_quantiles(base_samples, paste0('phi_sck[',1:data$N_stands,']'), display_ylim = c(0,1),
                                     ylab = "p(stand concordance)")


pre0 <- 5
vpd0 <- 23
pre <- 7

phi_quant <- data.frame()
for(vpd in seq(10, 25, 0.5)){
  logit_phi <- base_samples[['mu_phi']]+
    (vpd-vpd0)*base_samples[['beta_phi_vpd']]+
    (pre-pre0)*base_samples[['beta_phi_vpd']]
  phi <- boot::inv.logit(logit_phi)
  
  phi_quant <- rbind(
    phi_quant,
    data.frame(vpd, t(util$ensemble_mcmc_quantile_est(phi, c(0.05, 0.5, 0.95))))
  )
}
names(phi_quant) <- c('vpd', 'q5', 'q50', 'q95')

plot(x = NULL, y = NULL,
     xlim = c(10,25), ylim = c(0,0.2),
     xlab = 'VPD (hPa)', ylab = 'p(concordant event)',
     bty = 'n')
polygon(
  x = c(phi_quant$vpd, rev(phi_quant$vpd)),
  y = c(phi_quant$q95, rev(phi_quant$q5)),
  col = util$c_light,
  border = NA
)
lines(x = phi_quant$vpd, y = phi_quant$q50, col = util$c_mid_highlight, lwd = 2)


util$ensemble_mcmc_quantile_est(base_samples[['phi_sck[176]']], c(0.05, 0.5, 0.95))
s <- 176
range(data$vpd_obs[seq(1+data$N_all_years*(s-1), data$N_all_years*s, 1)])
range(data$pre_obs[seq(1+data$N_all_years*(s-1), data$N_all_years*s, 1)])

phi_quant <- data.frame()
for(i in seq(1+data$N_all_years*(s-1), data$N_all_years*s, 1)){
  logit_phi <- base_samples[[paste0('logit_phi_sck[', s, ']')]]+
    (data$vpd_obs[i]-vpd0)*base_samples[['beta_phi_vpd']]+
    (data$pre_obs[i]-pre0)*base_samples[['beta_phi_pre']]
  
  phi <- boot::inv.logit(logit_phi)
  
  phi_quant <- rbind(
    phi_quant,
    data.frame(data$vpd_obs[i], data$pre_obs[i], 
               t(util$ensemble_mcmc_quantile_est(phi, c(0.05, 0.5, 0.95))))
  )
}
names(phi_quant) <- c('vpd', 'pre', 'q5', 'q50', 'q95')


plot(x = NULL, y = NULL,
     xlim = c(7,17), ylim = c(8,34),
     xlab = 'VPD (hPa)', ylab = 'Winter prec. (dm)',
     bty = 'n', main = paste0('Stand ', s),
     cex.main = 1)
symbols(x = c(phi_quant$vpd, phi_quant$vpd),
        y = c(phi_quant$pre, logit_phi_sck),
        circles = c(phi_quant$q95, phi_quant$q5), 
        inches = 0.1, add = T,
        bg = c(rep(util$c_light, length(phi_quant$vpd)),
               rep(util$c_mid, length(phi_quant$vpd))),
        fg = NA)


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
names <- paste0('mean_phi[',1:data$N_all_years,']')
util$plot_conn_pushforward_quantiles(base_samples, names, data$all_years,
                                     ylab = 'Avg. p(concordant event)')

before2000 <- mean(sapply(which(data$all_years < 2000), 
                          function(y) util$ensemble_mcmc_quantile_est(base_samples[[paste0('mean_phi[',y,']')]], c(0.5))))
after2000 <- mean(sapply(which(data$all_years >= 2000), 
                          function(y) util$ensemble_mcmc_quantile_est(base_samples[[paste0('mean_phi[',y,']')]], c(0.5))))
segments(x0 = 1950, x1 = 1999, y0 = before2000, lty = 2, col = util$c_dark)
segments(x0 = 2000, x1 = 2024, y0 = after2000, lty = 2, col = util$c_dark)

base_samples[[paste0('mean_phi_bef1980')]] <- 0
years <-  which(data$all_years < 2000)
for(i in years){
  base_samples[[paste0('mean_phi_bef1980')]] <-  base_samples[[paste0('mean_phi_bef1980')]] + 
    base_samples[[paste0('mean_phi[',i,']')]]/length(years)
}

base_samples[[paste0('mean_phi_aft1980')]] <- 0
years <-  which(data$all_years >= 2000)
for(i in years){
  base_samples[[paste0('mean_phi_aft1980')]] <-  base_samples[[paste0('mean_phi_aft1980')]] + 
    base_samples[[paste0('mean_phi[',i,']')]]/length(years)
}

par(mar = c(4,4,1,1), mfrow = c(2,1))
util$plot_expectand_pushforward(base_samples[['mean_phi_bef1980']], 30,
                                flim = c(0, 0.2), display_name = 'p(concordant event), before 2000',
                                col = util$c_mid_teal)
util$plot_expectand_pushforward(base_samples[['mean_phi_aft1980']], 30,
                                flim = c(0, 0.2), display_name = 'p(concordant event), after 2000',
                                col = util$c_mid, add = F)




par(mfrow = c(1,1), mar = c(4,5,1,1))
for(s in 1:data$N_species){
  newname <- paste0('avg_omega_conc[', s,']')
  base_samples[[newname]] <- 0
  
  trees <- which(data$species_idxs %in% s)
  stand_species <- unique(data$stand_species_idxs[trees])
  
  for(stsp in stand_species){
    
    name1 <-  paste0('omega_conc_sck[', stsp,']')
    
    base_samples[[newname]] <- base_samples[[newname]] + (base_samples[[name1]])
  }
  base_samples[[newname]] <- base_samples[[newname]]/length(stand_species)
}
util$plot_disc_pushforward_quantiles_sign(base_samples, paste0('avg_omega_conc[',phylo_order,']'), display_ylim = c(-0.075,1),
                                          xticklabs =  data$uniq_species_ids[phylo_order], 
                                          ylab = "Average p(tree shock | stand concordance)", 
                                          xlab = '', ignore_sign = TRUE)

par(mfrow = c(1,1), mar = c(4,5,1,1))
for(s in 1:data$N_species){
  newname <- paste0('avg_omega_shutdown[', s,']')
  base_samples[[newname]] <- 0
  
  trees <- which(data$species_idxs %in% s)
  stand_species <- unique(data$stand_species_idxs[trees])
  
  for(stsp in stand_species){
    
    name1 <-  paste0('omega_shutdown[', stsp,']')
    
    base_samples[[newname]] <- base_samples[[newname]] + (base_samples[[name1]])
  }
  base_samples[[newname]] <- base_samples[[newname]]/length(stand_species)
}
util$plot_disc_pushforward_quantiles_sign(base_samples, paste0('avg_omega_shutdown[',phylo_order,']'), display_ylim = c(-0.075,1),
                                          xticklabs =  data$uniq_species_ids[phylo_order], 
                                          ylab = "Average p(shutdown | tree shock, stand concordance)", 
                                          xlab = '', ignore_sign = TRUE)
