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

data <- readRDS(file.path(folder,'data_07august2026_23species_153stands_2010andmore_19502024.rds'))

oldfit <- readRDS(file.path(wd, 'output/model/model21bis', 'fit_model21bis_HGSP_23species_153stands_2010andmore.rds'))
# fit <- readRDS(file.path(wd, 'output/model/model21bis', 'fit_model21bis_HGSP_vpd2_pre2_rescaled_23species_153stands_2010andmore.rds'))
# fit$time()
# fit$diagnostic_summary()

params <- c(
  "mu_alpha", 
  "sigma_alpha_stand", "alpha_stand",
  "mu_beta_gdd", "sigma_beta_gdd", "beta_gdd",
  "mu_beta_pre", "sigma_beta_pre", "beta_pre",
  "mu_beta_vpd", "sigma_beta_vpd", "beta_vpd",
  "mu_beta_vpd2", "sigma_beta_vpd2", "beta_vpd2",
  "mu_beta_pre2", "sigma_beta_pre2", "beta_pre2",
  "mu_beta_prexvpd", "sigma_beta_prexvpd", "beta_prexvpd",
  'beta_phi_vpd', 'beta_phi_pre',
  'tau_conc'
)

# new_samples <- fit$draws(params)
# names <- dimnames(new_samples)$variable
# new_samples <- lapply(1:dim(new_samples)[3],
#                        function(k){t(matrix(new_samples[1:dim(new_samples)[1],1:dim(new_samples)[2],k],
#                                             nrow = dim(new_samples)[1], ncol = dim(new_samples)[2]))})
# names(new_samples) <- names

new_samples <- readRDS(file.path(wd, 'output/model/model21bis', 'subsamples_model21bis_HGSP_2ndtaylor_rescaled_23species_153stands_2010andmore.rds'))

for(p in params){
  message(paste0('\n\n---------\nParameter(s):',p))
  rest <- util$filter_expectands(new_samples, p, 
                                 check_arrays = TRUE)
  util$check_all_expectand_diagnostics(rest, min_ess_hat_per_chain=50)
}

params <- c(
  "mu_alpha", 
  "sigma_alpha_stand", "alpha_stand",
  "mu_beta_gdd", "sigma_beta_gdd", "beta_gdd",
  "mu_beta_pre", "sigma_beta_pre", "beta_pre",
  "mu_beta_vpd", "sigma_beta_vpd", "beta_vpd",
  'beta_phi_vpd', 'beta_phi_pre',
  'tau_conc'
)

old_samples <- oldfit$draws(params)
names <- dimnames(old_samples)$variable
old_samples <- lapply(1:dim(old_samples)[3],
                       function(k){t(matrix(old_samples[1:dim(old_samples)[1],1:dim(old_samples)[2],k],
                                            nrow = dim(old_samples)[1], ncol = dim(old_samples)[2]))})
names(old_samples) <- names



util$plot_pairs_by_chain(new_samples[['beta_vpd[1]']], 'beta_vpd[1]',
                         new_samples[['beta_vpd2[1]']], 'beta_vpd2[1]')        

util$plot_pairs_by_chain(new_samples[['mu_beta_vpd2']], 'mu_beta_vpd2',
                         new_samples[['mu_beta_vpd']], 'mu_beta_vpd')  

util$plot_pairs_by_chain(old_samples[['mu_beta_vpd']], 'mu_beta_vpd (no VPD2)',
                         new_samples[['mu_beta_vpd']], 'mu_beta_vpd')  

util$plot_pairs_by_chain(old_samples[['beta_phi_vpd']], 'beta_phi_vpd (no VPD2)',
                         new_samples[['beta_phi_vpd']], 'beta_phi_vpd')  

util$plot_pairs_by_chain(old_samples[['beta_phi_pre']], 'beta_phi_pre (no VPD2)',
                         new_samples[['beta_phi_pre']], 'beta_phi_pre')  


probs <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
gdd0 = 10
pre0 = 5
vpd0 = 23
g <- mean(data$gdd_obs)
p <- mean(data$pre_obs)
muolds <- c()
mus <- c()
vpds <- seq(quantile(data$vpd_obs, 0.05), quantile(data$vpd_obs, 0.95), length.out = 40)
for(v in vpds){
  muold <- old_samples[['mu_alpha']] +
    old_samples[['mu_beta_gdd']] * (g - gdd0) +
    old_samples[['mu_beta_pre']] * (p - pre0) +
    old_samples[['mu_beta_vpd']] * (v - vpd0) 
  muolds <- rbind(muolds, util$ensemble_mcmc_quantile_est(muold, probs))
  
  mu <- new_samples[['mu_alpha']] +
    new_samples[['mu_beta_gdd']] * (g - gdd0)/10 +
    new_samples[['mu_beta_pre']] * (p - pre0)/10 +
    new_samples[['mu_beta_pre2']] * ((p - pre0)/10)^2 +
    new_samples[['mu_beta_vpd']] * (v - vpd0)/10 + 
    new_samples[['mu_beta_vpd2']] * ((v - vpd0)/10)^2 +
    2 * new_samples[['mu_beta_prexvpd']] * ((v - vpd0)/10) * ((p - pre0)/10)
  mus <- rbind(mus, util$ensemble_mcmc_quantile_est(mu, probs))
}

par(mfrow = c(1,1), mar = c(4,4,1,1))
plot(y = NULL, x = NULL, xlab = 'VPD (hPa)', ylab = 'log(ring width)', bty = 'n',
     xlim = c(1, 50), ylim =  c(-2,0.5))
polygon(c(vpds, rev(vpds)), c(muolds[,'5%'], rev(muolds[,'95%'])),
        col = adjustcolor(util$c_light_teal, 0.3), border = NA)
lines(x = vpds, y = muolds[,'95%'], col = util$c_mid_teal, lty = 2)
lines(x = vpds, y = muolds[,'5%'], col = util$c_mid_teal, lty = 2)
lines(x = vpds, y = muolds[,'50%'], col = util$c_mid_teal)

polygon(c(vpds, rev(vpds)), c(mus[,'5%'], rev(mus[,'95%'])),
        col = adjustcolor(util$c_light, 0.3), border = NA)
lines(x = vpds, y = mus[,'95%'], col = util$c_mid, lty = 2)
lines(x = vpds, y = mus[,'5%'], col = util$c_mid, lty = 2)
lines(x = vpds, y = mus[,'50%'], col = util$c_mid)

h <- hist(data$vpd_obs, breaks = 30, plot = FALSE)

counts <- -2 + h$counts / max(h$counts) 

for (i in seq_along(counts)) {
  rect(h$breaks[i], -2,
       h$breaks[i + 1], counts[i],
       col = adjustcolor("grey70", 0.3),
       border = adjustcolor("grey70", 0.5))
}





probs <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
gdd0 = 10
pre0 = 5
vpd0 = 23
g <- mean(data$gdd_obs)
v <- mean(data$vpd_obs)
muolds <- c()
mus <- c()
pres <- seq(0, 40, 1)
pres <- seq(quantile(data$pre_obs, 0.05), quantile(data$pre_obs, 0.95), length.out = 40)
for(p in pres){
  muold <- old_samples[['mu_alpha']] +
    old_samples[['mu_beta_gdd']] * (g - gdd0) +
    old_samples[['mu_beta_pre']] * (p - pre0) +
    old_samples[['mu_beta_vpd']] * (v - vpd0) 
  muolds <- rbind(muolds, util$ensemble_mcmc_quantile_est(muold, probs))
  
  mu <- new_samples[['mu_alpha']] +
    new_samples[['mu_beta_gdd']] * (g - gdd0)/10 +
    new_samples[['mu_beta_pre']] * (p - pre0)/10 +
    new_samples[['mu_beta_pre2']] * ((p - pre0)/10)^2 +
    new_samples[['mu_beta_vpd']] * (v - vpd0)/10 + 
    new_samples[['mu_beta_vpd2']] * ((v - vpd0)/10)^2 +
    2 * new_samples[['mu_beta_prexvpd']] * ((v - vpd0)/10) * ((p - pre0)/10)
  mus <- rbind(mus, util$ensemble_mcmc_quantile_est(mu, probs))
}

pres <- pres*100
par(mfrow = c(1,1), mar = c(4,4,1,1))
plot(y = NULL, x = NULL, xlab = 'Winter prec. (mm)', ylab = 'log(ring width)', bty = 'n',
     xlim = c(0, 4000), ylim =  c(-2,2))
polygon(c(pres, rev(pres)), c(muolds[,'5%'], rev(muolds[,'95%'])),
        col = adjustcolor(util$c_light_teal, 0.3), border = NA)
lines(x = pres, y = muolds[,'95%'], col = util$c_mid_teal, lty = 2)
lines(x = pres, y = muolds[,'5%'], col = util$c_mid_teal, lty = 2)
lines(x = pres, y = muolds[,'50%'], col = util$c_mid_teal)

polygon(c(pres, rev(pres)), c(mus[,'5%'], rev(mus[,'95%'])),
        col = adjustcolor(util$c_light, 0.3), border = NA)
lines(x = pres, y = mus[,'95%'], col = util$c_mid, lty = 2)
lines(x = pres, y = mus[,'5%'], col = util$c_mid, lty = 2)
lines(x = pres, y = mus[,'50%'], col = util$c_mid)

h <- hist(data$pre_obs*100, breaks = 30, plot = FALSE)

counts <- -2 + h$counts / max(h$counts) 

for (i in seq_along(counts)) {
  rect(h$breaks[i], -2,
       h$breaks[i + 1], counts[i],
       col = adjustcolor("grey70", 0.3),
       border = adjustcolor("grey70", 0.5))
}


par(mfrow = c(1,2), mar = c(4,4,1,1), cex.lab = 0.9, cex.axis = 0.85)
util$plot_expectand_pushforward(new_samples[['beta_phi_vpd']], 30, 'beta_phi_vpd', flim = c(0, 0.25))
util$plot_expectand_pushforward(old_samples[['beta_phi_vpd']], 30, 'beta_phi_vpd', flim = c(0, 0.25),
                                col = util$c_mid_teal, add = T)

util$plot_expectand_pushforward(new_samples[['beta_phi_pre']], 30, 'beta_phi_pre', flim = c(-0.25, 0))
util$plot_expectand_pushforward(old_samples[['beta_phi_pre']], 30, 'beta_phi_pre', flim = c(-0.25, 0),
                                col = util$c_mid_teal, add = T)

par(mfrow = c(2,1), mar = c(2,4.5,1,1))
util$plot_disc_pushforward_quantiles(new_samples, paste0('tau_conc[',1:data$N_species,']'),
                                     ylab = bquote(tau[shock]))
util$plot_disc_pushforward_quantiles(old_samples, paste0('tau_conc[',1:data$N_species,']'),
                                     ylab = bquote(tau[shock] ~'(no VPD²)'))

par(mfrow = c(1,2), mar = c(4,4,1,1))
util$plot_expectand_pushforward(new_samples[['mu_beta_vpd']], 30, 'mu_beta_vpd', flim = c(-0.75, 0.75))
util$plot_expectand_pushforward(old_samples[['mu_beta_vpd']]*10, 30, 'mu_beta_vpd', flim = c(-0.75, 0.75),
                                col = util$c_mid_teal, add = T)

util$plot_expectand_pushforward(new_samples[['mu_beta_pre']], 30, 'mu_beta_pre', flim = c(-0.75, 0.75))
util$plot_expectand_pushforward(old_samples[['mu_beta_pre']], 30, 'mu_beta_pre', flim = c(-0.75, 0.75),
                                col = util$c_mid_teal, add = T)


# Tension between prior and likelihood

util$plot_expectand_pushforward(new_samples[['mu_beta_vpd']], 30, 'mu_beta_vpd', flim = c(-0.75, 0.75))
prior <- rnorm(1e6, 0, log(1.3) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
prior <- rnorm(1e6, 0, log(1.75) / 2.57)
lines(density(prior), col = util$c_mid_teal, lwd = 2, lty = 5)
util$plot_expectand_pushforward(new_samples[['mu_beta_vpd2']], 30, 'mu_beta_vpd2', flim = c(-0.75, 0.75))
prior <- rnorm(1e6, 0, log(1.3) / 2.57 / 2)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
prior <- rnorm(1e6, 0, log(1.5) / 2.57)
lines(density(prior), col = util$c_mid_teal, lwd = 2, lty = 5)

util$plot_expectand_pushforward(new_samples[['sigma_beta_vpd']], 30, 'sigma_beta_vpd', flim = c(0, 0.75))
prior <- rnorm(1e6, 0, log(1.3) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
prior <- rnorm(1e6, 0, log(1.75) / 2.57)
lines(density(prior), col = util$c_mid_teal, lwd = 2, lty = 5)
util$plot_expectand_pushforward(new_samples[['sigma_beta_vpd2']], 30, 'sigma_beta_vpd2', flim = c(0, 0.75))
prior <- rnorm(1e6, 0, log(1.3) / 2.57 / 2)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
prior <- rnorm(1e6, 0, log(1.5) / 2.57)
lines(density(prior), col = util$c_mid_teal, lwd = 2, lty = 5)


util$plot_expectand_pushforward(new_samples[['mu_beta_pre']], 30, 'mu_beta_pre', flim = c(-0.75, 0.75))
prior <- rnorm(1e6, 0, log(1.3) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
prior <- rnorm(1e6, 0, log(1.75) / 2.57)
lines(density(prior), col = util$c_mid_teal, lwd = 2, lty = 5)
util$plot_expectand_pushforward(new_samples[['mu_beta_pre2']], 30, 'mu_beta_pre2', flim = c(-0.75, 0.75))
prior <- rnorm(1e6, 0, log(1.3) / 2.57 / 2)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
prior <- rnorm(1e6, 0, log(1.5) / 2.57)
lines(density(prior), col = util$c_mid_teal, lwd = 2, lty = 5)

util$plot_expectand_pushforward(new_samples[['sigma_beta_pre']], 30, 'sigma_beta_pre', flim = c(0, 0.75))
prior <- rnorm(1e6, 0, log(1.3) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
prior <- rnorm(1e6, 0, log(1.75) / 2.57)
lines(density(prior), col = util$c_mid_teal, lwd = 2, lty = 5)
util$plot_expectand_pushforward(new_samples[['sigma_beta_pre2']], 30, 'sigma_beta_pre2', flim = c(0, 0.75))
prior <- rnorm(1e6, 0, log(1.3) / 2.57 / 2)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
prior <- rnorm(1e6, 0, log(1.5) / 2.57)
lines(density(prior), col = util$c_mid_teal, lwd = 2, lty = 5)

par(mfrow = c(2,1))
util$plot_disc_pushforward_quantiles_sign(new_samples, paste0('beta_pre[', 1:data$N_species, ']'), ylab = 'beta_pre',
                                          xticklabs = rep(NA, data$N_species),
                                          ignore_sign = T)
util$plot_disc_pushforward_quantiles_sign(new_samples, paste0('beta_pre2[', 1:data$N_species, ']'), 
                                     ylab = 'beta_pre2', xlab = 'Species', xticklabs = data$uniq_species_ids,
                                     ignore_sign = T)

par(mfrow = c(1,2), cex.lab = 0.9, cex.axis = 0.85)
util$plot_expectand_pushforward(new_samples[['mu_beta_pre2']], 30, 'mu_beta_pre2', flim = c(-0.75, 0.75))
prior <- rnorm(1e6, 0, log(1.5) / 2.57)
lines(density(prior), col = util$c_mid_teal, lwd = 2, lty = 5)
util$plot_expectand_pushforward(new_samples[['sigma_beta_pre2']], 30, 'sigma_beta_pre2', flim = c(0, 1.5))
prior <- rnorm(1e6, 0, log(1.5) / 2.57)
lines(density(prior), col = util$c_mid_teal, lwd = 2, lty = 5)

# check all parameters

params <- c(
  "mu_alpha", 
  "sigma_alpha_stand", "alpha_stand",
  "mu_beta_gdd", "sigma_beta_gdd", "beta_gdd",
  "mu_beta_pre", "sigma_beta_pre", "beta_pre",
  "mu_beta_vpd", "sigma_beta_vpd", "beta_vpd",
  "mu_beta_vpd2", "sigma_beta_vpd2", "beta_vpd2",
  "delta_clim",
  "tau_clim", # "kappa_clim_free", "kappa_clim",
  'log_kappa_clim', 'kappa_clim',
  "f_tilde_ind",
  "mu_log_rho", "sigma_log_rho", "rho_merged",
  "mu_log_gamma", "sigma_log_gamma", "log_gamma_merged", 'gamma_merged',
  "mu_log_tau_conc", "sigma_log_tau_conc", "tau_conc",
  "mu_phi", "sigma_phi", "logit_phi_sck", "phi_sck", "beta_phi_vpd", "beta_phi_pre",
  "mu_omega_conc", "sigma_omega_species", "sigma_omega_stand", "alpha_omega_species", "alpha_omega_stand",
  "mu_omega_shutdown", "sigma_omega_shutdown", "logit_omega_shutdown", "omega_shutdown",
  "thetas_idio", "tau_idio",
  # 'thetas_baseline', 'omega_thetas', 'thetas_idio',
  "sigma"
)

new_samples <- fit$draws(params)
names <- dimnames(new_samples)$variable
new_samples <- lapply(1:dim(new_samples)[3],
                      function(k){t(matrix(new_samples[1:dim(new_samples)[1],1:dim(new_samples)[2],k],
                                           nrow = dim(new_samples)[1], ncol = dim(new_samples)[2]))})
names(new_samples) <- names

for(p in params){
  message(paste0('\n\n---------\nParameter(s):',p))
  rest <- util$filter_expectands(new_samples, p, 
                                 check_arrays = TRUE)
  util$check_all_expectand_diagnostics(rest, min_ess_hat_per_chain=50)
}
