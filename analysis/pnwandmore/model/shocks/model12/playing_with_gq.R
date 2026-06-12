
library(cmdstanr)

fit_upd <- readRDS("/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model/model13/fit_29jan2025_gymnosperms_standclimate_19502024_6species_model13_updatedpriors.rds")
fit_upd$time()
fit_upd$diagnostic_summary()

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

which(omegas_conc_q50 > 0.8 & omegas_nonconc_q50 < 0.05)


subset <- which((phis_q50)*omegas_conc_q50 > 0.1)

tree_idxs <- which(data$stand_species_idxs == 29)
par(mfrow = c(3,3), cex.lab = 1.2)
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
  
  # start <- 1+(stand-1)*data$N_all_years
  # end <- stand*data$N_all_years
  # lines(1:data$N_all_years, data$pre_obs[start:end]/10, col = 'darkblue')
  # lines(1:data$N_all_years, data$vpd_obs[start:end]/20, col = 'darkred')
  
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






par(mfrow = c(1,2), cex.lab = 0.85)
plot(1, type="n", main='',
     xlim= c(0,1), xlab= bquote(omega[shock]^{concordant} ~ "(" * "individuals" * ")"),
     ylim= c(0,15), ylab="",  yaxt="n",
     bty = "n")
title(ylab="Estimated Bin\nProbabilities / Bin Width",
      mgp=c(1, 1, 0))
h <- hist(rbeta(1e6, 230,14), plot = FALSE, breaks = seq(0,1,0.01))
h$density = h$counts/sum(h$counts)*100
# plot(h, freq=FALSE, add = T, col = util$c_light_teal, border = 'white')
for(stsp in 1:data$N_stand_species){
  util$plot_expectand_pushforward(samples[[paste0('omega_conc_sck[',stsp,']')]], 50,
                                  flim = c(0,1), ylim = c(0,15),
                                  display_name = bquote(omega[shock]^{concordant}),
                                  col = '#8f2727', add = TRUE)
}
plot(1, type="n", main='',
     xlim= c(0,1), xlab= bquote(omega[shock]^{concordant} ~ "(" * "population" * ")"),
     ylim= c(0,15), ylab="",  yaxt="n",
     bty = "n")
util$plot_expectand_pushforward(samples[[paste0('omega_conc_sck0')]], 50,
                                flim = c(0,1), ylim = c(0,15),
                                display_name = bquote(omega[shock]^{concordant}),
                                col = '#278f27', add = T)

par(mfrow = c(1,2), cex.lab = 0.85)
plot(1, type="n", main='',
     xlim= c(0,1), xlab= bquote(omega[shock]^{concordant} ~ "(" * "population location" * ")"),
     ylim= c(0,20), ylab="",  yaxt="n",
     bty = "n", xaxs = "i")
title(ylab="Estimated Bin\nProbabilities / Bin Width",
      mgp=c(1, 1, 0))
h <- hist(rbeta(1e6, 230,14), plot = FALSE, breaks = seq(0,1,0.008))
h$density = h$counts/sum(h$counts)*100
plot(h, freq=FALSE, add = T, col = util$c_light_teal, border = 'white')
util$plot_expectand_pushforward(samples[[paste0('omega_conc_sck0')]], 50,
                                flim = c(0,1), ylim = c(0,20),
                                display_name = bquote(omega[shock]^{concordant}),
                                add = T)

plot(1, type="n", main='',
     xlim= c(0,1.5), xlab= bquote(tau[omega[shock]^{concordant}] ~ "(" * "population scale, logit scale" * ")"),
     ylim= c(0,10), ylab="",  yaxt="n",
     bty = "n", xaxs = "i")
title(ylab="Estimated Bin\nProbabilities / Bin Width",
      mgp=c(1, 1, 0))
h <- hist(rnorm(1e6, 0,0.3/2.57), plot = FALSE, breaks = seq(-1.5,1.5,0.008*1.5))
h$density = h$counts/sum(h$counts)*100
plot(h, freq=FALSE, add = T, col = util$c_light_teal, border = 'white')
util$plot_expectand_pushforward(samples[[paste0('tau_omega_conc_sck')]], 50,
                                flim = c(-1.5,1.5), ylim = c(0,10),
                                display_name = bquote(omega[shock]^{concordant}),
                                add = T)



par(mfrow =c(1,2))
plot(1, type="n", main='',
     xlim= c(0,6), xlab= bquote(tau[shock] ~ "(" * "population location" * ")"),
     ylim= c(0,6), ylab="",  yaxt="n",
     bty = "n", xaxs = "i")
title(ylab="Estimated Bin\nProbabilities / Bin Width",
      mgp=c(1, 1, 0))
h <- hist(rnorm(1e6, 0, log(20)/2.57), plot = FALSE, breaks = seq(-6,6,0.1))
h$density = h$counts/sum(h$counts)*100
plot(h, freq=FALSE, add = T, col = util$c_light_teal, border = 'white')
util$plot_expectand_pushforward(samples[[paste0('mu_tau_sck[1]')]], 121,
                                flim = c(-6,6), ylim = c(0,6),
                                display_name = bquote(omega[shock]^{concordant}),
                                add = T)

plot(1, type="n", main='',
     xlim= c(0,2), xlab= bquote(tau[shock] ~ "(" * "population scale" * ")"),
     ylim= c(0,6), ylab="",  yaxt="n",
     bty = "n", xaxs = "i")
title(ylab="Estimated Bin\nProbabilities / Bin Width",
      mgp=c(1, 1, 0))
h <- hist(rnorm(1e6, 0, log(2)/2.57), plot = FALSE, breaks = seq(-2,2,0.03))
h$density = h$counts/sum(h$counts)*100
plot(h, freq=FALSE, add = T, col = util$c_light_teal, border = 'white')
util$plot_expectand_pushforward(samples[[paste0('tau_tau_sck[1]')]], 134,
                                flim = c(-2,2), ylim = c(0,6),
                                display_name = bquote(omega[shock]^{concordant}),
                                add = T)





par(mfrow = c(1,1))
util$plot_expectand_pushforward(samples[[paste0('phi_sck[29]')]], 134,
                                flim = c(0,1), ylim = c(0,20),
                                display_name = bquote(phi[shock]))

util$check_all_expectand_diagnostics(util$filter_expectands(samples, 'phi_sck[29]'))
