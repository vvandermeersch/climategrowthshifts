rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(cmdstanr)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)


folder <- '/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model'

# data <- readRDS(file.path(folder,'data_19apr2026_PIPOexcept5_19202024.rds'))
# fit <- readRDS(file.path(folder, 'model17/fit_model17_HSGP_PSME_waidor_pi0_8cores.rds'))
data <- readRDS(file.path(folder,'data_22apr2026_PSME_except2_19202024.rds'))
datasets <- readRDS(file.path(folder,'datasets_22apr2026_PSME_except2_19202024.rds'))
# fit <- readRDS(file.path(folder, 'model17/fit_model17_HSGP_PSME_waidoraz_except2_pi0.rds'))
# gc()
# fit$diagnostic_summary()

# Load samples
# samples <- fit$draws(c('alpha', 'beta_gdd', 'beta_vpd', 'beta_pre',
# 
#                        'rho_sp', 'gamma_sp',
#                        'rho_sh', 'gamma_sh',
#                        'rho_ind', 'gamma_ind',
# 
#                        'tau_sck',
#                        'phi_sck0', 'tau_phi_sck', 'alpha_phi_sck',
#                        'omega_conc_sck0', 'tau_omega_conc_sck', 'alpha_omega_conc_sck',
# 
#                        # 'pi_idsc_sck',
# 
#                        'sigma'
# ))
# names <- dimnames(samples)$variable
# samples <- lapply(1:dim(samples)[3],
#                   function(k){t(matrix(samples[1:dim(samples)[1],1:dim(samples)[2],k],
#                                        nrow = dim(samples)[1], ncol = dim(samples)[2]))})
# names(samples) <- names

samples <- readRDS(file.path(folder, 'model17/basesamples_model17_HSGP_PSME_except2_pi0.rds'))
# samples <- readRDS(file.path(folder, 'model17/basesamples_model17_HSGP_19apr2026_PIPOexcept5_pi0_preXgdd.rds'))
util$check_all_expectand_diagnostics(samples)

# Intercept and slopes
par(mfrow = c(2,2))
util$plot_expectand_pushforward(samples[['alpha']], 50, flim= c(-2,2), 'alpha')
prior <- rnorm(1e6, 0, log(5)/2.32)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(samples[['beta_gdd']], 50, flim= c(-0.35,0.35), 'beta_gdd')
prior <- rnorm(1e6, 0, log(1.8) / 2.57)eue 
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(samples[['beta_vpd']], 50, flim= c(-0.35,0.35), 'beta_vpd')
prior <- rnorm(1e6, 0, log(1.8) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(samples[['beta_pre']], 50, flim= c(-0.35,0.35), 'beta_pre')
prior <- rnorm(1e6, 0, log(1.8) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(samples[['beta_pre_vpd']], 50, flim= c(-0.1,0.1), 'beta_pre_vpd')
prior <- rnorm(1e6, 0, log(1.8) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)

# Gaussian processes
par(mfrow = c(3,2))
util$plot_expectand_pushforward(samples[['rho_sp']], 50, flim= c(1,150), 'rho_sp')
prior <- rlnorm(1e6, 3.7, 0.35)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(samples[['gamma_sp']], 50, flim= c(0,2), 'gamma_sp')
prior <- rnorm(1e6, 0, log(10) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(samples[['rho_sh']], 50, flim= c(1,10), 'rho_sh')
prior <- rlnorm(1e6, 0.4, 0.3)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(samples[['gamma_sh']], 50, flim= c(0,1), 'gamma_sh')
prior <- rnorm(1e6, 0, log(3) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(samples[['rho_ind']], 50, flim= c(1,10), 'rho_ind')
prior <- rlnorm(1e6, 1.4, 0.35)
lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(samples[['gamma_ind']], 50, flim= c(0,1), 'gamma_ind')
prior <- rnorm(1e6, 0, log(3) / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)

util$plot_pairs_by_chain(samples[['gamma_sh']], 'gamma_sh',
                         samples[['tau_sck']], 'tau_sck')
util$plot_pairs_by_chain(samples[['gamma_ind']], 'gamma_ind',
                         samples[['sigma']], 'sigma')
util$plot_pairs_by_chain(samples[['rho_sp']], 'rho_sp',
                         samples[['rho_ind']], 'rho_ind')
abline(a = 0,b =1)


# Shocks
par(mfrow = c(4,2))
util$plot_expectand_pushforward(samples[['tau_sck']], 50, flim= c(0,5), 'tau_sck')
prior <- rnorm(1e6, 0, 5 / 2.57)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
plot.new()
util$plot_expectand_pushforward(samples[['phi_sck0']], 50, flim= c(0,1), 'phi_sck0')
prior <- rbeta(1e6, 2, 20)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(samples[['tau_phi_sck']], 50, flim= c(0,3), 'tau_phi_sck')
prior <- rnorm(1e6, 0, 1)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(samples[['omega_conc_sck0']], 50, flim= c(0,1), 'omega_conc_sck0')
prior <- rbeta(1e6, 12.28, 34)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
prior <- rbeta(1e6, 3, 17)
lines(density(prior), col = util$c_mid_teal, lwd = 2, lty = 1)
util$plot_expectand_pushforward(samples[['tau_omega_conc_sck']], 50, flim= c(0,3), 'tau_omega_conc_sck')
prior <- rnorm(1e6, 0, 1)
lines(density(prior), col = util$c_light_teal, lwd = 1.5, lty = 2)
util$plot_expectand_pushforward(samples[['beta_conc_pre']], 30, flim= c(-0.5,0.5), 'beta_conc_pre')
util$plot_expectand_pushforward(samples[['beta_conc_vpd']], 100, flim= c(-0.5,0.5), 'beta_conc_vpd')


states <- datasets[match(data$uniq_stand_ids, datasets$grouped_stand), 'state']
# truenames <- list(az = 'Arizona', or = 'Oregon', id = 'Idaho', wa = 'Wash.')
truenames <- list(co = 'Colorado', wy = 'Wyoming', mt = 'Montana')

# shocks <- fit$draws(c('phi_sck', 'omega_conc_sck'))
# names <- dimnames(shocks)$variable
# shocks <- lapply(1:dim(shocks)[3], 
#                   function(k){t(matrix(shocks[1:dim(shocks)[1],1:dim(shocks)[2],k], 
#                                        nrow = dim(shocks)[1], ncol = dim(shocks)[2]))})
# names(shocks) <- names
shocks <- readRDS(file.path(folder, 'model17/shocksamples_model17_HSGP_PSME_cowymt_pi0.rds'))
par(mfrow = c(2,1), mar = c(4,4,1,1), cex.axis = 0.8, cex.lab = 1)
util$plot_disc_pushforward_quantiles(shocks, paste0('phi_sck[', 1:data$N_stands, ']'), display_ylim = c(0,1),
                                     ylab = 'p(stand conc.)')
for(st in unique(states)){
  xmin <- min(which(states == st))-0.5
  xmax <- max(which(states == st))+0.5
  abline(v = c(xmax), lty = 2, col = 'grey40')
  text(x = (xmax+xmin)/2, y = 0.8, labels = truenames[st], cex = 0.8, adj = 0.5, col = "#474747")
}
util$plot_disc_pushforward_quantiles(shocks, paste0('omega_conc_sck[', 1:data$N_stands, ']'), display_ylim = c(0,1),
                                     ylab = 'p(tree sck | stand conc.)', xlab = 'Sites')
for(st in unique(states)){
  xmin <- min(which(states == st))-0.5
  xmax <- max(which(states == st))+0.5
  abline(v = c(xmax), lty = 2, col = 'grey40')
  text(x = (xmax+xmin)/2, y = 0.8, labels = truenames[st], cex = 0.8, adj = 0.5, col = "#474747")
}




# Generate predictions!
# options(cmdstanr_output_dir = file.path(folder, 'model17', 'tmp'))
# mod_gq <- cmdstan_model(file.path(wd, 'model/stan/model17/hsgp/archive', 'model17_HSGP_indGP_only1species_newpriors_pi0_withGQ.stan'))
# data_gq <- data
# data_gq$grainsize <-  ceiling(data_gq$N_stands/4)
# data_gq$uniq_tree_ids <- NULL
# data_gq$uniq_species_ids <- NULL
# data_gq$uniq_stand_ids <- NULL
# data_gq$uniq_stand_lat <- NULL
# data_gq$uniq_stand_lon <- NULL
# fit_gq <- mod_gq$generate_quantities(fit, data = data_gq, seed = 5838293, parallel_chains = 2)
# fit_gq$save_object(file = file.path(folder, 'fit_gq_model17/model17_HSGP_PIPOsub_pi0.rds'))
# gc()
# gq_samples <- fit_gq$draws(c('num_conc_sck_stdlvl', 'avg_exp_m1_delta_sck_stdlvl', 
#                              'avg_exp_m1_clim_stdlvl'))
# gc()
# names <- dimnames(gq_samples)$variable
# gq_samples <- lapply(1:dim(gq_samples)[3], 
#                   function(k){t(matrix(gq_samples[1:dim(gq_samples)[1],1:dim(gq_samples)[2],k], 
#                                        nrow = dim(gq_samples)[1], ncol = dim(gq_samples)[2]))})
# names(gq_samples) <- names
# gq_samples <- readRDS(file.path(folder, 'model17/standpred_samples_model17_HSGP_PSME_colorado_pi0.rds'))
# gc()

# stand_samples <- fit$draws(c('f_sh'))
# gc()
# names <- dimnames(stand_samples)$variable
# stand_samples <- lapply(1:dim(stand_samples)[3], 
#                         function(k){t(matrix(stand_samples[1:dim(stand_samples)[1],1:dim(stand_samples)[2],k], 
#                                              nrow = dim(stand_samples)[1], ncol = dim(stand_samples)[2]))})
# names(stand_samples) <- names


tree_samples <- readRDS(file.path(folder, 'model17/treesamples_model17_HSGP_19apr2026_PIPOexcept5_pi0_preXgdd.rds'))
gc()

pdf(file.path(wd, 'model/shocks/model17/figures', 'treesamples_PIPOexcept5_pi0_preXgdd_pi0.pdf'), height = 10, width = 16)
par(mfrow = c(8,4), cex.lab = 0.8, cex.axis = 0.8, mar = c(1.5,4,0.5,1), mgp = c(1.5, 0.2, 0),tck = -0.03)
for(t in 1:data$N_trees){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- paste0("log_rw_pred[",idxs,"]")
  
  temp <- lapply(tree_samples[names], function(m) m[1, , drop = FALSE])
  
  util$plot_conn_pushforward_quantiles(tree_samples, names, data$years[idxs],
                                       xlab="Year", ylab="Log ring width (per mm)", 
                                       display_ylim=c(-11, 2), display_xlim = range(data$all_years))
  points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=0.5, col="white")
  points(data$years[idxs], log(data$rw_obs[idxs]), pch=16, cex=0.25, col="black")
  text(x = 1920, y = -9, labels = data$uniq_tree_ids[t], cex = 0.8, adj = 0, col = "#474747")
  # text(x = 1920, y = -7.7, labels = paste0('Species:', data$uniq_species_ids[data$species_idxs[t]]), cex = 1, adj = 0)
  
  names <- paste0("delta_sck[",idxs,"]")
  util$plot_conn_pushforward_quantiles(tree_samples, names, data$years[idxs],
                                       xlab="Year", ylab="Shock amplitude", 
                                       display_ylim=c(-7, 1), display_xlim = range(data$all_years))
  
  # s <- unique(data$stand_idxs[t])
  # names <- paste0("stand_clim[",s,',',1:data$N_all_years,"]")
  # util$plot_conn_pushforward_quantiles(tree_samples, names, data$all_years,
  #                                      xlab="Year", ylab="Climate effect (stand)", 
  #                                      display_ylim=c(-1, 1), display_xlim = range(data$all_years))
  # 
  names <- paste0("f[",idxs,"]")
  util$plot_conn_pushforward_quantiles(tree_samples, names, data$years[idxs],
                                       xlab="Year", ylab="Long-term GP", 
                                       display_ylim=c(-3, 3), display_xlim = range(data$all_years))
  
  names <- paste0("f_ind[",idxs,"]")
  util$plot_conn_pushforward_quantiles(tree_samples, names, data$years[idxs],
                                       xlab="Year", ylab="Short-term GP", 
                                       display_ylim=c(-3, 3), display_xlim = range(data$all_years))
  
  
  # names <- paste0("f_sh[",s,',',1:data$N_all_years,"]")
  # util$plot_conn_pushforward_quantiles(stand_samples, names, data$all_years,
  #                                      xlab="Year", ylab="Stand GP", 
  #                                      display_ylim=c(-3, 3), display_xlim = range(data$all_years))
  
  
}
dev.off()

stand_samples <- readRDS(file.path(folder, 'model17/standsamples_model17_HSGP_PSME_waidor_pi0_16cores.rds'))
gc()

pdf(file.path(wd, 'model/shocks/model17/figures', 'standshocks_PSME_waidor_pi0.pdf'), height = 13, width = 10)
par(mfrow = c(10,2), mar = c(2.5,2.8,1,1), mgp = c(1.5,0.5,0),
    cex.lab = 0.9, cex.axis = 0.75, cex.main = 0.9)
util$plot_disc_pushforward_quantiles(shocks, paste0('phi_sck[', 1:data$N_stands, ']'), display_ylim = c(0,1),
                                     ylab = 'p(stand conc.)')
util$plot_disc_pushforward_quantiles(shocks, paste0('omega_conc_sck[', 1:data$N_stands, ']'), display_ylim = c(0,1),
                                     ylab = 'p(tree sck | stand conc.)')
# util$plot_expectand_pushforward(samples[['pi_idsc_sck']], 30, display_name = 'p(idiosync. shock)')
for(s in 1:data$N_stands){
  idxs <- seq(1 + data$N_all_years * (s-1), data$N_all_years * s, 1)
  idxs_obs <- idxs[seq(data$stand_start_years_idxs[s], data$N_stand_years[s], 1)]
  for(i in idxs){
    newname <- paste0('num_conc_sck_stdlvl_onlyobs[',i,']')
    name <- paste0('num_conc_sck_stdlvl[',i,']')
    
    if(i %in% idxs_obs){
      stand_samples[[newname]] <- stand_samples[[name]]
    }else{
      stand_samples[[newname]] <- stand_samples[[name]]+1e6 # this way, do not appear on the plot!
    }
    
  }
  names <- paste0('num_conc_sck_stdlvl_onlyobs[',idxs,']')
  util$plot_disc_pushforward_quantiles(stand_samples, names, main = 'Concordant shock (% of trees)', 
                                       ylab = paste0('Stand ', s), display_ylim = c(0, 1))
  # names <- paste0('num_idsc_sck_stdlvl[',idxs,']')
  # util$plot_disc_pushforward_quantiles(gq_samples, names, main = 'Idiosyncratic shock (% of trees)', 
  #                                      ylab = '', display_ylim = c(0, 0.5))
  # names <- paste0('num_dbl_sck_stdlvl[',idxs,']')
  # util$plot_disc_pushforward_quantiles(gq_samples, names, main = 'Both shocks (% of trees)', 
  #                                      ylab = '', display_ylim = c(0, 0.5))
}
dev.off()





# par(mfrow = c(1,1), cex.main = 0.9, mar = c(2.5,2.8,1,1))
# plot(x = gq_samples[[paste0('num_conc_sck_stdlvl[',1,']')]], 
#      y = gq_samples[[paste0('num_idsc_sck_stdlvl[',1,']')]], bty = 'n',
#      pch = 20, col = paste0(util$c_mid, 10),
#      xlab = 'Concordant shock (% of trees)',
#      ylab = 'Idiosyncratic shock (% of trees)',
#      main = 'Stand 1, year 2')


# Diagnostics the very few divergences
# samples <- fit$draws(c('mu_phi_sck', 'tau_phi_sck', 'alpha_phi_sck_tilde',
#                        'mu_omega_conc_sck', 'tau_omega_conc_sck', 'alpha_omega_conc_sck'
# ))
# names <- dimnames(samples)$variable
# samples <- lapply(1:dim(samples)[3], 
#                   function(k){t(matrix(samples[1:dim(samples)[1],1:dim(samples)[2],k], 
#                                        nrow = dim(samples)[1], ncol = dim(samples)[2]))})
# names(samples) <- names
# util$check_all_expectand_diagnostics(samples)
# 
# diag_array <- fit$sampler_diagnostics()
# div_matrix <- diag_array[, , "divergent__"]
# d <- dim(div_matrix)
# diagnostics <- c()
# diagnostics[['divergent__']] <- t(matrix(unclass(div_matrix), nrow = d[1], ncol = d[2]))
# sum(diagnostics[['divergent__']] )
# 
# util$plot_div_pairs('mu_phi_sck', 'tau_phi_sck', samples,
#                     diagnostics, transforms=list('tau_phi_sck' = 1))
# util$plot_div_pairs('mu_omega_conc_sck', 'tau_omega_conc_sck', samples,
#                     diagnostics, transforms=list('tau_omega_conc_sck' = 1))
# 
# util$plot_div_pairs(paste0('alpha_phi_sck[',1:data$N_stands, ']'), 'tau_phi_sck', samples,
#                     diagnostics, transforms=list('tau_phi_sck' = 1))
# util$plot_div_pairs(paste0('alpha_omega_conc_sck[',1:data$N_stands, ']'), 'tau_omega_conc_sck', samples,
#                     diagnostics, transforms=list('tau_omega_conc_sck' = 1))

