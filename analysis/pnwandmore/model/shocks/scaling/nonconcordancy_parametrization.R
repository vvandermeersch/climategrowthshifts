rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)

# # The centered version is much faster, and no issues arise --- note that I still need to center phi!
# fit_uncent <- readRDS("/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model/model13/fit_29jan2025_gymnosperms_standclimate_19801990_5species_model13.rds")
# fit_uncent$time()
# fit_uncent$diagnostic_summary()
# 
# fit_centdelta <- readRDS("/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model/model13/fit_29jan2025_gymnosperms_standclimate_19801990_5species_model13_centdelta.rds")
# fit_centdelta$time()
# fit_centdelta$diagnostic_summary()
# 
# fit_centall <- readRDS("/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model/model13/fit_29jan2025_gymnosperms_standclimate_19801990_5species_model13_centall.rds")
# fit_centall$time()
# fit_centall$diagnostic_summary()

data <- readRDS("/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model/data_30jan2025_gymnosperms_standclimate_19502024_6species.rds")


fit_upd <- readRDS("/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model/model13/fit_29jan2025_gymnosperms_standclimate_19502024_6species_model13_updatedpriors.rds")
fit_upd$time()
fit_upd$diagnostic_summary()


# Let's check parameters now!
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

# Transform as a Rstan object (to be able to use Mike's functions)
samples <- lapply(1:dim(param_samples)[3], function(k){t(matrix(param_samples[1:dim(param_samples)[1],1:dim(param_samples)[2],k], 
                                                                nrow = dim(param_samples)[1], ncol = dim(param_samples)[2]))})
names(samples) <- dimnames(param_samples)$variable
util$check_all_expectand_diagnostics(samples)


util$plot_pairs_by_chain(samples[['mu_tau_sck[1]']], 'mu_tau_sck[1]',
                         log(samples[['tau_tau_sck[1]']]), 'log(tau_tau_sck[1])')


# Stand-specific shock probability
par(mfrow = c(1,2), cex.lab = 0.85)

plot(1, type="n", main='',
     xlim= c(0,0.7), xlab= bquote(phi[shock] ~ "(" * "individuals" * ")"),
     ylim= c(0,20), ylab="",  yaxt="n",
     bty = "n")
title(ylab="Estimated Bin\nProbabilities / Bin Width",
      mgp=c(1, 1, 0))
h <- hist(rbeta(1e6, 2,20), plot = FALSE, breaks = seq(0,1,0.01))
h$density = h$counts/sum(h$counts)*100
# plot(h, freq=FALSE, add = T, col = util$c_light_teal, border = 'white')
for(s in 1:data$N_stands){
  util$plot_expectand_pushforward(samples[[paste0('phi_sck[',s,']')]], 50,
                                  flim = c(0,0.7), ylim = c(0,20),
                                  display_name = bquote(phi[shock]),
                                  col = '#8f2727', add = TRUE)
}
plot(1, type="n", main='',
     xlim= c(0,0.7), xlab= bquote(phi[shock] ~ "(" * "population" * ")"),
     ylim= c(0,20), ylab="",  yaxt="n",
     bty = "n")
util$plot_expectand_pushforward(samples[[paste0('phi_sck0')]], 50,
                                flim = c(0,1), ylim = c(0,20),
                                display_name = bquote(phi[shock]),
                                col = '#278f27', add = T)


# Species-stand-specific shock probability given a stand in shock
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


util$plot_expectand_pushforward(samples[[paste0('mu_omega_conc_sck')]], 50,
                                flim = c(-5,5), ylim = c(0,25),
                                display_name = bquote(omega[shock]^{concordant}),
                                col = '#8f2727')
util$plot_expectand_pushforward(samples[[paste0('tau_omega_conc_sck')]], 50,
                                flim = c(0,5), ylim = c(0,25),
                                display_name = bquote(omega[shock]^{concordant}),
                                col = '#8f2727')


mu <- util$ensemble_mcmc_quantile_est(samples[[paste0('mu_omega_conc_sck')]], c(0.5))
tau <- util$ensemble_mcmc_quantile_est(samples[[paste0('tau_omega_conc_sck')]], c(0.5))
hist(rnorm(1e6, mu, tau))
hist(boot::inv.logit(rnorm(1e6, mu, tau)))
hist(boot::inv.logit(rnorm(1e6, mu, 0.25))) # tau should be much smaller than 1.7... more something like 0.25 or 0.3, my prior is not great here.
quantile(boot::inv.logit(rnorm(1e6, mu, 0.3)), c(0.025,0.975)) # everything between 0.8 and 0.92 ?


# Species-stand-specific delta
par(mfrow = c(1,2), cex.lab = 0.85)
plot(1, type="n", main='',
     xlim= c(0,1), xlab= bquote(omega[shock]^{non-concordant} ~ "(" * "individuals" * ")"),
     ylim= c(0,25), ylab="",  yaxt="n",
     bty = "n")
title(ylab="Estimated Bin\nProbabilities / Bin Width",
      mgp=c(1, 1, 0))
for(stsp in 1:data$N_stand_species){
  util$plot_expectand_pushforward(samples[[paste0('omega_nonconc_sck[',stsp,']')]], 50,
                                  flim = c(0,1), ylim = c(0,25),
                                  display_name = bquote(omega[shock]^{non-concordant}),
                                  col = '#8f2727', add = T)
}
plot(1, type="n", main='',
     xlim= c(0,1), xlab= bquote(omega[shock]^{non-concordant} ~ "(" * "population" * ")"),
     ylim= c(0,25), ylab="",  yaxt="n",
     bty = "n")
util$plot_expectand_pushforward(boot::inv.logit(samples[['mu_omega_conc_sck']] - exp(samples[['mu_logdelta_omega_nonconc_sck']])), 50,
                                flim = c(0,1), ylim = c(0,25),
                                display_name = bquote(omega[shock]^{concordant}),
                                col = '#278f27', add = T)


util$plot_expectand_pushforward(samples[['mu_logdelta_omega_nonconc_sck']], 50,
                                flim = c(0,3),
                                display_name = bquote(omega[shock]^{concordant}),
                                col = '#27278f')
prior <- rnorm(1e6, log(8), log(2)/2.57)
lines(density(prior),
      lwd = 1.5, lty = 3, col = 'grey70')

util$plot_expectand_pushforward(samples[['tau_logdelta_omega_nonconc_sck']], 50,
                                flim = c(0,1),
                                display_name = bquote(omega[shock]^{concordant}),
                                col = '#27278f')
prior <- rnorm(1e6, 0, 0.3/2.57)
lines(density(prior),
      lwd = 1.5, lty = 3, col = 'grey70')

plot(1, type="n", main='',
     xlim= c(0,5), xlab= '',
     ylim= c(0,6), ylab="",  yaxt="n",
     bty = "n")
for(stsp in 1:data$N_stand_species){
  util$plot_expectand_pushforward(samples[[paste0('logdelta_omega_nonconc_sck[',stsp,']')]], 50,
                                  flim = c(-5,5), ylim = c(0,25),
                                  display_name = bquote(omega[shock]^{non-concordant}),
                                  col = '#8f2727', add = T)
}

which(omegas_conc_q50 < 0.2)
(1-phis_q50[36])*omegas_nonconc_q50[36]
(phis_q50[36])*omegas_conc_q50[36]
util$plot_expectand_pushforward(samples[[paste0('logdelta_omega_nonconc_sck[',stsp,']')]], 50,
                                flim = c(-5,5), ylim = c(0,25),
                                display_name = bquote(omega[shock]^{non-concordant}),
                                col = '#8f2727')

subset <- which((1-phis_q50)*omegas_nonconc_q50 > 0.1)
phis_q50[22]
omegas_conc_q50[22]
omegas_nonconc_q50[22]

phis_q50[67]
omegas_conc_q50[67]
omegas_nonconc_q50[67]


par(mfrow = c(1,3), cex.lab = 0.85)
plot(1, type="n", main='',
     xlim= c(0,1), xlab= bquote(phi[shock] ~ "(" * "individuals" * ")"),
     ylim= c(0,30), ylab="",  yaxt="n",
     bty = "n")
title(ylab="Estimated Bin\nProbabilities / Bin Width",
      mgp=c(1, 1, 0))
for(stsp in subset){
  util$plot_expectand_pushforward(samples[[paste0('phi_sck[',stsp,']')]], 50,
                                  flim = c(0,1), ylim = c(0,30),
                                  display_name = bquote(omega[shock]^{non-concordant}),
                                  col = '#8f2727', add = T)
}
plot(1, type="n", main='',
     xlim= c(0,1), xlab= bquote(omega[shock]^{concordant} ~ "(" * "individuals" * ")"),
     ylim= c(0,10), ylab="",  yaxt="n",
     bty = "n")
title(ylab="Estimated Bin\nProbabilities / Bin Width",
      mgp=c(1, 1, 0))
for(stsp in subset){
  util$plot_expectand_pushforward(samples[[paste0('omega_conc_sck[',stsp,']')]], 50,
                                  flim = c(0,1), ylim = c(0,10),
                                  display_name = bquote(omega[shock]^{concordant}),
                                  col = '#8f2727', add = T)
}
plot(1, type="n", main='',
     xlim= c(0,1), xlab= bquote(omega[shock]^{non-concordant} ~ "(" * "individuals" * ")"),
     ylim= c(0,25), ylab="",  yaxt="n",
     bty = "n")
title(ylab="Estimated Bin\nProbabilities / Bin Width",
      mgp=c(1, 1, 0))
for(stsp in subset){
  util$plot_expectand_pushforward(samples[[paste0('omega_nonconc_sck[',stsp,']')]], 50,
                                  flim = c(0,1), ylim = c(0,25),
                                  display_name = bquote(omega[shock]^{non-concordant}),
                                  col = '#8f2727', add = T)
}






phis_q5 <- c()
phis_q50 <- c()
phis_q95 <- c()
omegas_conc_q5 <- c()
omegas_conc_q50 <- c()
omegas_conc_q95 <- c()
omegas_nonconc_q5 <- c()
omegas_nonconc_q50 <- c()
omegas_nonconc_q95 <- c()
for(stsp in 1:data$N_stand_species){
  s <- unique(data$stand_idxs[which(data$stand_species_idxs == stsp)])
  phis_q50 <- c(phis_q50, util$ensemble_mcmc_quantile_est(samples[[paste0('phi_sck[',s,']')]], c(0.5)))
  omegas_conc_q50 <- c(omegas_conc_q50, util$ensemble_mcmc_quantile_est(samples[[paste0('omega_conc_sck[',stsp,']')]], c(0.5)))
  omegas_nonconc_q50 <- c(omegas_nonconc_q50, util$ensemble_mcmc_quantile_est(samples[[paste0('omega_nonconc_sck[',stsp,']')]], c(0.5)))
  
  phis_q5 <- c(phis_q5, util$ensemble_mcmc_quantile_est(samples[[paste0('phi_sck[',s,']')]], c(0.05)))
  omegas_conc_q5 <- c(omegas_conc_q5, util$ensemble_mcmc_quantile_est(samples[[paste0('omega_conc_sck[',stsp,']')]], c(0.05)))
  omegas_nonconc_q5 <- c(omegas_nonconc_q5, util$ensemble_mcmc_quantile_est(samples[[paste0('omega_nonconc_sck[',stsp,']')]], c(0.05)))
  
  phis_q95 <- c(phis_q95, util$ensemble_mcmc_quantile_est(samples[[paste0('phi_sck[',s,']')]], c(0.95)))
  omegas_conc_q95 <- c(omegas_conc_q95, util$ensemble_mcmc_quantile_est(samples[[paste0('omega_conc_sck[',stsp,']')]], c(0.95)))
  omegas_nonconc_q95 <- c(omegas_nonconc_q95, util$ensemble_mcmc_quantile_est(samples[[paste0('omega_nonconc_sck[',stsp,']')]], c(0.95)))
}

par(mfrow = c(1,1))
plot(phis_q50*omegas_conc_q50, (1-phis_q50)*omegas_nonconc_q50, pch = 20)
abline(a =0, b = 1, lty = 2)

pre <- sapply(1:data$N_stand_species,function(s){
  start <- 1+(s-1)*data$N_all_years
  end <- s*data$N_all_years
  mean(data$pre_obs[start:end])
})

plot(pre, (1-phis_q50)*omegas_nonconc_q50, pch = 20)


sum(phis_q50*omegas_conc_q50 > (1-phis_q50)*omegas_nonconc_q50)
sum(phis_q50*omegas_conc_q50 < (1-phis_q50)*omegas_nonconc_q50)

for(stsp in 1:data$N_stand_species){
  segments(x0 = phis_q50[stsp]*omegas_conc_q50[stsp], 
           y0 = (1-phis_q5[stsp])*omegas_nonconc_q5[stsp], 
           y1 = (1-phis_q95[stsp])*omegas_nonconc_q95[stsp])
  segments(y0 = (1-phis_q50[stsp])*omegas_nonconc_q50[stsp], 
           x0 = (phis_q5[stsp])*omegas_conc_q5[stsp], 
           x1 = (phis_q95[stsp])*omegas_conc_q95[stsp])
}






util$plot_expectand_pushforward(samples[[paste0('mu_logdelta_omega_nonconc_sck')]], 50,
                                flim = c(0,5), ylim = c(0,25),
                                display_name = bquote(omega[shock]^{concordant}),
                                col = '#8f2727')
util$plot_expectand_pushforward(samples[[paste0('tau_logdelta_omega_nonconc_sck')]], 50,
                                flim = c(0,1), ylim = c(0,25),
                                display_name = bquote(omega[shock]^{concordant}),
                                col = '#8f2727')

mu <- util$ensemble_mcmc_quantile_est(samples[[paste0('mu_logdelta_omega_nonconc_sck')]], c(0.5))
tau <- util$ensemble_mcmc_quantile_est(samples[[paste0('tau_logdelta_omega_nonconc_sck')]], c(0.5))
hist(rnorm(1e6, mu, tau))
hist(boot::inv.logit(boot::logit(0.9) - exp(rnorm(1e6, mu, tau)))) # again, too wide right?
hist(boot::inv.logit(boot::logit(0.9) - exp(rnorm(1e6, mu, 0.1))))
hist(boot::inv.logit(boot::logit(0.9) - exp(rnorm(1e6, log(8), 0.8))), xlim = c(0,1), breaks = 50) 
round(quantile(boot::inv.logit(boot::logit(0.9) - exp(rnorm(1e6, log(8), 0.3))), c(0.025, 0.975)), 3)
hist(boot::inv.logit(boot::logit(0.9) - exp(rnorm(1e6, log(8), 0.3))), xlim = c(0,1), breaks = 50) 
