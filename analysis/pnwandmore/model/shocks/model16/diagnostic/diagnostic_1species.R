rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(cmdstanr)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)


folder <- '/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model/'

data <- readRDS(file.path(folder,'data_26mar2026_TSME_standclimate_18962024.rds'))
fit <- readRDS(file.path(folder, 'model16/fit_HSGP_26mar2026_TSME_NCPphi_NCPdeltaomega_newpriors.rds'))
gc()


# Load samples
samples <- fit$draws(c('alpha', 'beta_gdd', 'beta_vpd', 'beta_pre',
                       
                       'rho_sp', 'gamma_sp',
                       'rho_sh', 'gamma_sh',
                       'rho_ind', 'gamma_ind',
                       
                       'tau_sck',
                       'phi_sck0', 'tau_phi_sck', 'alpha_phi_sck',
                       'omega_conc_sck0', 'tau_omega_conc_sck', 'alpha_omega_conc_sck',
                       'mu_logdelta_omega_nonconc_sck', 'tau_logdelta_omega_nonconc_sck', 'logdelta_omega_nonconc_sck',
                       
                       'sigma'
))
names <- dimnames(samples)$variable
samples <- lapply(1:dim(samples)[3], 
                  function(k){t(matrix(samples[1:dim(samples)[1],1:dim(samples)[2],k], 
                                       nrow = dim(samples)[1], ncol = dim(samples)[2]))})
names(samples) <- names

# base_samples <- util$filter_expectands(samples, c("mu_logdelta_omega_nonconc_sck"), check_arrays = TRUE)
util$check_all_expectand_diagnostics(samples)

# Intercept and slopes
par(mfrow = c(2,2))
util$plot_expectand_pushforward(samples[['alpha']], 50, flim= c(-2,2), 'alpha')
prior <- rnorm(1e6, 0, log(5)/2.32)
lines(density(prior), col = 'darkblue')
util$plot_expectand_pushforward(samples[['beta_gdd']], 50, flim= c(-0.5,0.5), 'beta_gdd')
prior <- rnorm(1e6, 0, log(1.8) / 2.57)
lines(density(prior), col = 'darkblue')
util$plot_expectand_pushforward(samples[['beta_vpd']], 50, flim= c(-0.5,0.5), 'beta_vpd')
prior <- rnorm(1e6, 0, log(1.8) / 2.57)
lines(density(prior), col = 'darkblue')
util$plot_expectand_pushforward(samples[['beta_pre']], 50, flim= c(-2,2), 'beta_pre')
prior <- rnorm(1e6, 0, log(1.8) / 2.57)
lines(density(prior), col = 'darkblue')

# Gaussian processes
par(mfrow = c(3,2))
util$plot_expectand_pushforward(samples[['rho_sp']], 50, flim= c(0,100), 'rho_sp')
prior <- rlnorm(1e6, 3.55, 0.24)
lines(density(prior), col = 'darkblue')
util$plot_expectand_pushforward(samples[['gamma_sp']], 50, flim= c(0,10), 'gamma_sp')
prior <- rnorm(1e6, 0, log(10) / 2.57)
lines(density(prior), col = 'darkblue')
util$plot_expectand_pushforward(samples[['rho_sh']], 50, flim= c(1,10), 'rho_sh')
prior <- rlnorm(1e6, 0.55, 0.24)
lines(density(prior), col = 'darkblue')
util$plot_expectand_pushforward(samples[['gamma_sh']], 50, flim= c(0,10), 'gamma_sh')
prior <- rnorm(1e6, 0, log(3) / 2.57)
lines(density(prior), col = 'darkblue')
util$plot_expectand_pushforward(samples[['rho_ind']], 50, flim= c(1,5), 'rho_ind')
prior <- rlnorm(1e6, 0.80, 0.40)
lines(density(prior), col = 'darkblue')
util$plot_expectand_pushforward(samples[['gamma_ind']], 50, flim= c(0,10), 'gamma_ind')
prior <- rnorm(1e6, 0, log(3) / 2.57)
lines(density(prior), col = 'darkblue')

util$plot_pairs_by_chain(samples[['gamma_sh']], 'gamma_sh',
                         samples[['rho_sh']], 'rho_sh')

# Shocks
par(mfrow = c(4,2))
util$plot_expectand_pushforward(samples[['tau_sck']], 50, flim= c(0,6), 'tau_sck')
prior <- rnorm(1e6, 0, log(20) / 2.57)
lines(density(prior), col = 'darkblue')
plot.new()
util$plot_expectand_pushforward(samples[['phi_sck0']], 50, flim= c(0,1), 'phi_sck0')
prior <- rbeta(1e6, 2, 20)
lines(density(prior), col = 'darkblue')
util$plot_expectand_pushforward(samples[['tau_phi_sck']], 50, flim= c(0,5), 'tau_phi_sck')
prior <- rnorm(1e6, 0, 1)
lines(density(prior), col = 'darkblue')
util$plot_expectand_pushforward(samples[['omega_conc_sck0']], 50, flim= c(0,1), 'omega_conc_sck0')
prior <- rbeta(1e6, 30, 20)
lines(density(prior), col = 'darkblue')
util$plot_expectand_pushforward(samples[['tau_omega_conc_sck']], 50, flim= c(0,5), 'tau_omega_conc_sck')
prior <- rnorm(1e6, 0, 1)
lines(density(prior), col = 'darkblue')
util$plot_expectand_pushforward(samples[['mu_logdelta_omega_nonconc_sck']], 50, flim= c(-1,10), 'mu_logdelta_omega_nonconc_sck')
prior <- rnorm(1e6, log(5), log(3)/2.57)
lines(density(prior), col = 'darkblue')
util$plot_expectand_pushforward(samples[['tau_logdelta_omega_nonconc_sck']], 50, flim= c(0,2), 'tau_logdelta_omega_nonconc_sck')
prior <- rnorm(1e6, 0, 1/2.57)
lines(density(prior), col = 'darkblue')


shocks <- fit$draws(c('phi_sck', 'omega_conc_sck', 'omega_nonconc_sck'))
names <- dimnames(shocks)$variable
shocks <- lapply(1:dim(shocks)[3], 
                  function(k){t(matrix(shocks[1:dim(shocks)[1],1:dim(shocks)[2],k], 
                                       nrow = dim(shocks)[1], ncol = dim(shocks)[2]))})
names(shocks) <- names
par(mfrow = c(3,1), mar = c(4,4,1,1))
util$plot_disc_pushforward_quantiles(shocks, paste0('phi_sck[', 1:data$N_stands, ']'), display_ylim = c(0,0.5),
                                     ylab = 'p(stand concordance)')
util$plot_disc_pushforward_quantiles(shocks, paste0('omega_conc_sck[', 1:data$N_stands, ']'), display_ylim = c(0,1),
                                     ylab = 'p(tree shock | stand concordance)')
util$plot_disc_pushforward_quantiles(shocks, paste0('omega_nonconc_sck[', 1:data$N_stands, ']'), display_ylim = c(0,0.1),
                                     ylab = 'p(tree shock | no stand concordance)')




# Generate predictions!
mod_gq <- cmdstan_model(file.path(wd, 'model/stan/model16/hsgp', 'model16_HSGP_indGP_only1species_NCPphi_NCPdeltaomega_withGQ.stan'))
data_gq <- data
data_gq$grainsize <-  ceiling(data_gq$N_stands/4)
data_gq$uniq_tree_ids <- NULL
data_gq$uniq_species_ids <- NULL
data_gq$uniq_stand_ids <- NULL
data_gq$uniq_stand_lat <- NULL
data_gq$uniq_stand_lon <- NULL
fit_gq <- mod_gq$generate_quantities(fit, data = data_gq, seed = 5838293, parallel_chains = 4)
gc()
gq_samples <- fit_gq$draws(c('log_rw_pred', 'f', 'f_ind', 'delta_conc_sck', 'delta_nonconc_sck'))
gc()
names <- dimnames(gq_samples)$variable
gq_samples <- lapply(1:dim(gq_samples)[3], 
                  function(k){t(matrix(gq_samples[1:dim(gq_samples)[1],1:dim(gq_samples)[2],k], 
                                       nrow = dim(gq_samples)[1], ncol = dim(gq_samples)[2]))})
names(gq_samples) <- names

stand_samples <- fit$draws(c('f_sh'))
gc()
names <- dimnames(stand_samples)$variable
stand_samples <- lapply(1:dim(stand_samples)[3], 
                        function(k){t(matrix(stand_samples[1:dim(stand_samples)[1],1:dim(stand_samples)[2],k], 
                                             nrow = dim(stand_samples)[1], ncol = dim(stand_samples)[2]))})
names(stand_samples) <- names

pdf(file.path(wd, 'model/shocks/model16/figures', 'TSME_newpriors.pdf'), height = 10, width = 19.5)
par(mfrow = c(8,6), cex.lab = 0.8, cex.axis = 0.8, mar = c(1.5,4,0.5,1), mgp = c(1.5, 0.2, 0),tck = -0.03)
for(t in 1:data$N_trees){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- paste0("log_rw_pred[",idxs,"]")
  
  util$plot_conn_pushforward_quantiles(gq_samples, names, data$years[idxs],
                                       xlab="Year", ylab="Log ring width (per mm)", 
                                       display_ylim=c(-8, 2), display_xlim = range(data$all_years))
  points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=0.5, col="white")
  points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=0.25, col="black")
  text(x = .5, y = -6.9, labels = data$uniq_tree_ids[t], cex = 1, adj = 0)
  text(x = .5, y = -7.7, labels = paste0('Species:', data$uniq_species_ids[data$species_idxs[t]]), cex = 1, adj = 0)
  
  names <- paste0("delta_conc_sck[",idxs,"]")
  util$plot_conn_pushforward_quantiles(gq_samples, names, data$years[idxs],
                                       xlab="Year", ylab="Concordant\nshock amplitude", 
                                       display_ylim=c(-5, 1), display_xlim = range(data$all_years))
  
  names <- paste0("delta_nonconc_sck[",idxs,"]")
  util$plot_conn_pushforward_quantiles(gq_samples, names, data$years[idxs],
                                       xlab="Year", ylab="Non-concordant\nshock amplitude", 
                                       display_ylim=c(-5, 1), display_xlim = range(data$all_years))
  
  names <- paste0("f[",idxs,"]")
  util$plot_conn_pushforward_quantiles(gq_samples, names, data$years[idxs],
                                       xlab="Year", ylab="Long-term GP", 
                                       display_ylim=c(-5, 2), display_xlim = range(data$all_years))
  
  names <- paste0("f_ind[",idxs,"]")
  util$plot_conn_pushforward_quantiles(gq_samples, names, data$years[idxs],
                                       xlab="Year", ylab="Short-term GP", 
                                       display_ylim=c(-5, 2), display_xlim = range(data$all_years))
  
  s <- unique(data$stand_idxs[t])
  names <- paste0("f_sh[",s,',',1:data$N_all_years,"]")
  util$plot_conn_pushforward_quantiles(stand_samples, names, data$all_years,
                                       xlab="Year", ylab="Stand GP", 
                                       display_ylim=c(-8, 2), display_xlim = range(data$all_years))
  
  
}
dev.off()



# Diagnostics the very few divergences
samples <- fit$draws(c('mu_phi_sck', 'tau_phi_sck', 'alpha_phi_sck_tilde',
                       'mu_omega_conc_sck', 'tau_omega_conc_sck', 'alpha_omega_conc_sck',
                       'mu_logdelta_omega_nonconc_sck', 'tau_logdelta_omega_nonconc_sck', 'logdelta_omega_nonconc_sck_tilde'
))
names <- dimnames(samples)$variable
samples <- lapply(1:dim(samples)[3], 
                  function(k){t(matrix(samples[1:dim(samples)[1],1:dim(samples)[2],k], 
                                       nrow = dim(samples)[1], ncol = dim(samples)[2]))})
names(samples) <- names
util$check_all_expectand_diagnostics(samples)

diag_array <- fit$sampler_diagnostics()
div_matrix <- diag_array[, , "divergent__"]
d <- dim(div_matrix)
diagnostics <- c()
diagnostics[['divergent__']] <- t(matrix(unclass(div_matrix), nrow = d[1], ncol = d[2]))
sum(diagnostics[['divergent__']] )

util$plot_div_pairs('mu_phi_sck', 'tau_phi_sck', samples,
                    diagnostics, transforms=list('tau_phi_sck' = 1))
util$plot_div_pairs('mu_omega_conc_sck', 'tau_omega_conc_sck', samples,
                    diagnostics, transforms=list('tau_omega_conc_sck' = 1))
util$plot_div_pairs('mu_logdelta_omega_nonconc_sck', 'tau_logdelta_omega_nonconc_sck', samples,
                    diagnostics, transforms=list('tau_logdelta_omega_nonconc_sck' = 1))

util$plot_div_pairs(paste0('alpha_phi_sck_tilde[',1:data$N_stands, ']'), 'tau_phi_sck', samples,
                    diagnostics, transforms=list('tau_phi_sck' = 1))
util$plot_div_pairs(paste0('alpha_omega_conc_sck[',1:data$N_stands, ']'), 'tau_omega_conc_sck', samples,
                    diagnostics, transforms=list('tau_omega_conc_sck' = 1))
util$plot_div_pairs(paste0('logdelta_omega_nonconc_sck_tilde[',1:data$N_stands, ']'), 'tau_logdelta_omega_nonconc_sck', samples,
                    diagnostics, transforms=list('tau_logdelta_omega_nonconc_sck' = 1))
