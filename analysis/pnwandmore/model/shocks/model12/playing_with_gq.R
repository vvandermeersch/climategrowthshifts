
library(cmdstanr)


# Parameters samples
param_samples <- fit_upd$draws(variables = c("mu_alpha","mu_gdd", "mu_vpd", "mu_pre", 
                                             "tau_alpha","tau_gdd", "tau_vpd", "tau_pre",
                                             "alpha","beta_gdd", "beta_vpd", "beta_pre",
                                             
                                             "mu_tau_sck", "tau_tau_sck", "tau_sck",
                                             "phi_sck", "omega_conc_sck",
                                             
                                             "mu_rho", "tau_rho", "rho_sp",
                                             "mu_gamma", "tau_gamma", "gamma_sp",
                                             
                                             "rho_sh",
                                             "mu_kappa", "tau_kappa", "kappa_sh",
                                             
                                             "phi_sck0", "tau_phi_sck", "alpha_phi_sck",
                                             "omega_conc_sck0", "mu_omega_conc_sck", "tau_omega_conc_sck", "alpha_omega_conc_sck",
                                             "mu_logdelta_omega_nonconc_sck", "tau_logdelta_omega_nonconc_sck", "logdelta_omega_nonconc_sck", "omega_nonconc_sck",
                                             
                                             "sigma"))
samples <- lapply(1:dim(param_samples)[3], function(k){t(matrix(param_samples[1:dim(param_samples)[1],1:dim(param_samples)[2],k], 
                                                                nrow = dim(param_samples)[1], ncol = dim(param_samples)[2]))})
names(samples) <- dimnames(param_samples)$variable
util$check_all_expectand_diagnostics(samples)

fsh_samples <- fit_upd$draws(variables = c("f_sh"))
names <- dimnames(fsh_samples)$variable
fsh_samples <- lapply(1:dim(fsh_samples)[3], function(k){t(matrix(fsh_samples[1:dim(fsh_samples)[1],1:dim(fsh_samples)[2],k], 
                                                                nrow = dim(fsh_samples)[1], ncol = dim(fsh_samples)[2]))})
names(fsh_samples) <- names

# GQ samples
mod_gq <- cmdstan_model(file.path(wd, 'model/stan', 'model13_updatedpriors_wGQ.stan'))
data_gq <- data
data_gq$grainsize <-  ceiling(data_gq$N_stands/8)
data_gq$uniq_tree_ids <- NULL
data_gq$uniq_species_ids <- NULL
data_gq$uniq_stand_ids <- NULL
data_gq$N_clades <- 1
data_gq$uniq_stand_lat <- NULL
data_gq$uniq_stand_lon <- NULL
fit_gq <- mod_gq$generate_quantities(fit_upd, data = data_gq, seed = 5838293, parallel_chains = 4)

gq_samples <- fit_gq$draws()
names <- dimnames(gq_samples)$variable
gq_samples <- lapply(1:dim(gq_samples)[3], function(k){t(matrix(gq_samples[1:dim(gq_samples)[1],1:dim(gq_samples)[2],k], 
                                                                nrow = dim(gq_samples)[1], ncol = dim(gq_samples)[2]))})
names(gq_samples) <- names

subset <- which((1-phis_q50)*omegas_nonconc_q50 > 0.1)
phis_q50[22]
omegas_conc_q50[22]
omegas_nonconc_q50[22]

phis_q50[67]
omegas_conc_q50[67]
omegas_nonconc_q50[67]

phis_q50[73]
omegas_conc_q50[73]
omegas_nonconc_q50[73]



subset <- which((phis_q50)*omegas_conc_q50 > 0.1)

tree_idxs <- which(data$stand_species_idxs == 57)
par(mfrow = c(1,3), cex.lab = 1.2)
for(t in tree_idxs){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  stand <- data$stand_idxs[t]
  names <- paste0("log_rw_pred[",idxs,"]")
  
  util$plot_conn_pushforward_quantiles(gq_samples, names, data$years[idxs]-1979,
                                       xlab="Year", ylab="Log ring width (per mm)", 
                                       display_ylim=c(-8, 2), display_xlim = c(1980, 2024)-1979)
  points(data$years[idxs]-1979, log(data$rw_obs[idxs]), pch=16, cex=1, col="white")
  points(data$years[idxs]-1979, log(data$rw_obs[idxs]), pch=16, cex=0.5, col="black")
  text(x = .5, y = -6.9, labels = data$uniq_tree_ids[t], cex = 1, adj = 0)
  text(x = .5, y = -7.7, labels = paste0('Species:', data$uniq_species_ids[data$species_idxs[t]]), cex = 1, adj = 0)
  
  start <- 1+(stand-1)*data$N_all_years
  end <- stand*data$N_all_years
  lines(1:data$N_all_years, data$pre_obs[start:end]/10, col = 'darkblue')
  lines(1:data$N_all_years, data$vpd_obs[start:end]/20, col = 'darkred')
  
  names <- paste0("sck_state[",idxs,"]")
  util$plot_disc_pushforward_quantiles(gq_samples, names, data$years[idxs],
                                       xlab="Year", ylab="Shock state", 
                                       display_ylim=c(-0.05, 1.05))
  
  names <- paste0("delta_sck[",idxs,"]")
  util$plot_conn_pushforward_quantiles(gq_samples, names, data$years[idxs]-1979,
                                       xlab="Year", ylab="Shock amplitude", 
                                       display_ylim=c(-2, 1), display_xlim = c(1980, 2024)-1979)
}


names <- paste0("f_sh[",30,',',1:data$N_all_years,"]")
util$plot_conn_pushforward_quantiles(fsh_samples, names, 1:data$N_all_years,
                                     xlab="Year", ylab="Shock amplitude", 
                                     display_ylim=c(-6, 6))


names <- paste0("beta_gdd[",1:data$N_species,"]")
util$plot_disc_pushforward_quantiles(samples, names, 1:data$N_species,
                                     xlab="Year", ylab="Shock state", 
                                     display_ylim=c(-0.1, 0.1))
