rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(terra)
library(rnaturalearth)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)


data <-readRDS(file.path(wd, 'output/model', 'data_30jan2025_gymnosperms_standclimate_19502024_6species.rds'))
fit <- readRDS(file.path(wd, 'output/model', 'fit_30jan2025_gymnosperms_standclimate_19502024_6species.rds'))

param_samples <- fit$draws(variables = c("mu_alpha","mu_gdd", "mu_vpd", "mu_pre", 
                                         "tau_alpha","tau_gdd", "tau_vpd", "tau_pre",
                                         "alpha","beta_gdd", "beta_vpd", "beta_pre",
                                         
                                         "mu_tau_sck", "tau_tau_sck", "tau_sck",
                                         "phi_sck", "omega_conc_sck",
                                         
                                         "mu_rho", "tau_rho", "rho_sp",
                                         "mu_gamma", "tau_gamma", "gamma_sp",
                                         
                                         "rho_sh",
                                         "mu_kappa", "tau_kappa", "kappa_sh",
                                         
                                         "sigma"))


# Transform as a Rstan object (to be able to use Mike's functions)
samples <- lapply(1:dim(param_samples)[3], function(k){t(matrix(param_samples[1:dim(param_samples)[1],1:dim(param_samples)[2],k], 
                                                                nrow = dim(param_samples)[1], ncol = dim(param_samples)[2]))})
names(samples) <- dimnames(param_samples)$variable
util$check_all_expectand_diagnostics(samples)


util$plot_pairs_by_chain(samples[['rho_sp[3]']], 'rho_sp[3]',
                         log(samples[['tau_rho[1]']]), 'log(tau_rho[1])')

util$plot_pairs_by_chain(samples[['rho_sp[1]']], 'rho_sp[1]',
                         log(samples[['tau_rho[1]']]), 'log(tau_rho[1])')

util$plot_pairs_by_chain(samples[['rho_sp[1]']], 'rho_sp[1]',
                         log(samples[['tau_rho[1]']]), 'log(tau_rho[1])')

util$plot_pairs_by_chain(samples[['tau_sck[3]']], 'tau_sck[3]',
                         samples[['tau_tau_sck[1]']], 'tau_tau_sck[1]')

# Species-specific shock amplitude
par(mfrow = c(1,1))
plot(1, type="n", main='',
     xlim= c(0,2), xlab= bquote(tau[shock]),
     ylim= c(0,20), ylab="",  yaxt="n",
     bty = "n")
title(ylab="Estimated Bin\nProbabilities / Bin Width",
      mgp=c(1, 1, 0))
h <- hist(abs(rnorm(1e6, 0, log(20)/2.57)), plot = FALSE, breaks = 100)
h$density = h$counts/sum(h$counts)*100
plot(h, freq=FALSE, add = T, col = 'grey90', border = 'white')
for(sp in 1:data$N_species){
  util$plot_expectand_pushforward(samples[[paste0('tau_sck[',sp,']')]], 100,
                                  flim = c(0,2), ylim = c(0,20),
                                  display_name = bquote(tau[shock]),
                                  col = '#8f2727', add = TRUE)
}

# Stand-specific shock probability
par(mfrow = c(1,2), cex.lab = 0.85)
for(i in 1:2){
  plot(1, type="n", main='',
       xlim= c(0,0.7), xlab= bquote(phi[shock]),
       ylim= c(0,20), ylab="",  yaxt="n",
       bty = "n")
  title(ylab="Estimated Bin\nProbabilities / Bin Width",
        mgp=c(1, 1, 0))
  h <- hist(rbeta(1e6, 2,20), plot = FALSE, breaks = 50)
  h$density = h$counts/sum(h$counts)*100
  plot(h, freq=FALSE, add = T, col = util$c_light_teal, border = 'white')
}
for(s in 1:data$N_stands){
  util$plot_expectand_pushforward(samples[[paste0('phi_sck[',s,']')]], 50,
                                  flim = c(0,0.7), ylim = c(0,20),
                                  display_name = bquote(phi[shock]),
                                  col = '#8f2727', add = TRUE)
}

# Species-stand-specific shock probability given a stand in shock
par(mfrow = c(1,2), cex.lab = 0.85)
for(i in 1:2){
  plot(1, type="n", main='',
       xlim= c(0.7,1), xlab= bquote(omega[shock]^{concordant}),
       ylim= c(0,25), ylab="",  yaxt="n",
       bty = "n")
  title(ylab="Estimated Bin\nProbabilities / Bin Width",
        mgp=c(1, 1, 0))
  h <- hist(rbeta(1e6, 230,14), plot = FALSE, breaks = 50)
  h$density = h$counts/sum(h$counts)*100
  plot(h, freq=FALSE, add = T, col = util$c_light_teal, border = 'white')
}
for(stsp in 1:data$N_stand_species){
  util$plot_expectand_pushforward(samples[[paste0('omega_conc_sck[',stsp,']')]], 50,
                                  flim = c(0.7,1), ylim = c(0,25),
                                  display_name = bquote(omega[shock]^{concordant}),
                                  col = '#8f2727', add = TRUE)
}

# Species-specific GDD slope
par(mfrow = c(1,1))
minf <- -1.5
maxf <- 1.5
nbins <- 200
delta <- (maxf - minf) / 100
bins <- seq(minf, maxf, delta)
plot(1, type="n", main='',
     xlim= c(-0.2,0.2), xlab= bquote(beta[GDD]),
     ylim= c(0,20), ylab="",  yaxt="n",
     bty = "n")
title(ylab="Estimated Bin\nProbabilities / Bin Width",
      mgp=c(1, 1, 0))
h <- hist(rnorm(1e6, 0, log(1.8) / 2.57), plot = FALSE, breaks = bins)
h$density = h$counts/sum(h$counts)*100
plot(h, freq=FALSE, add = T, col = 'grey90', border = 'white')
for(sp in 1:data$N_species){
  util$plot_expectand_pushforward(samples[[paste0('beta_gdd[',sp,']')]], nbins,
                                  flim = c(minf,maxf), ylim = c(0,20),
                                  display_name = bquote(beta[GDD]),
                                  col = '#8f2727', add = TRUE)
}

# Species-specific precipitation slope
par(mfrow = c(1,1))
minf <- -1.5
maxf <- 1.5
nbins <- 200
delta <- (maxf - minf) / 100
bins <- seq(minf, maxf, delta)
plot(1, type="n", main='',
     xlim= c(-0.2,0.2), xlab= bquote(beta[prec]),
     ylim= c(0,20), ylab="",  yaxt="n",
     bty = "n")
title(ylab="Estimated Bin\nProbabilities / Bin Width",
      mgp=c(1, 1, 0))
h <- hist(rnorm(1e6, 0, log(1.8) / 2.57), plot = FALSE, breaks = bins)
h$density = h$counts/sum(h$counts)*100
plot(h, freq=FALSE, add = T, col = 'grey90', border = 'white')
for(sp in 1:data$N_species){
  util$plot_expectand_pushforward(samples[[paste0('beta_pre[',sp,']')]], nbins,
                                  flim = c(minf,maxf), ylim = c(0,20),
                                  display_name = bquote(beta[prec]),
                                  col = '#8f2727', add = TRUE)
}

par(mfrow=c(1,1))
util$plot_expectand_pushforward(samples[[paste0('mu_pre[1]')]], 50,
                                flim = c(-0.1,0.1), ylim = c(0,100),
                                display_name = bquote(beta[climate]),
                                col = '#27278f')
util$plot_expectand_pushforward(samples[[paste0('mu_gdd[1]')]], 50,
                                flim = c(-0.1,0.1), ylim = c(0,100),
                                display_name = bquote(beta[prec]),
                                col = '#278f27', add = T)
util$plot_expectand_pushforward(samples[[paste0('mu_vpd[1]')]], 50,
                                flim = c(-0.1,0.1), ylim = c(0,100),
                                display_name = bquote(beta[prec]),
                                col = '#8f2727', add = T)


# Species 3 has issues
util$check_rhat(samples[['beta_gdd[3]']])
util$check_rhat(samples[['beta_vpd[3]']])
util$check_rhat(samples[['beta_pre[3]']])
util$check_rhat(samples[['alpha[3]']])
util$check_rhat(samples[['rho_sp[3]']])
util$check_rhat(samples[['kappa_sh[3]']])
util$check_rhat(samples[['tau_sck[3]']])


util$plot_pairs_by_chain(samples[['tau_sck[3]']], 'tau_sck[3]',
                         samples[['rho_sp[3]']], 'rho_sp[3]')


util$plot_pairs_by_chain(samples[['tau_sck[3]']], 'tau_sck[3]',
                         samples[['sigma']], 'sigma')

trees_sp3 <- which(data$species_idxs == 3)
stands_sp3 <- unique(data$stand_idxs[trees_sp3])
stand_species_sp3 <- unique(data$stand_species_idxs[trees_sp3])

par(mfrow = c(1,1))
plot(1, type="n", main='',
     xlim= c(0,0.7), xlab= bquote(phi[shock]),
     ylim= c(0,20), ylab="",  yaxt="n",
     bty = "n")
title(ylab="Estimated Bin\nProbabilities / Bin Width",
      mgp=c(1, 1, 0))
h <- hist(rbeta(1e6, 2,20), plot = FALSE, breaks = 50)
h$density = h$counts/sum(h$counts)*100
plot(h, freq=FALSE, add = T, col = util$c_light_teal, border = 'white')
for(s in stands_sp3){
  util$plot_expectand_pushforward(samples[[paste0('phi_sck[',s,']')]], 50,
                                  flim = c(0,0.7), ylim = c(0,20),
                                  display_name = bquote(phi[shock]),
                                  col = '#8f2727', add = TRUE)
}


plot(1, type="n", main='',
     xlim= c(0.7,1), xlab= bquote(omega[shock]^{concordant}),
     ylim= c(0,25), ylab="",  yaxt="n",
     bty = "n")
title(ylab="Estimated Bin\nProbabilities / Bin Width",
      mgp=c(1, 1, 0))
h <- hist(rbeta(1e6, 230,14), plot = FALSE, breaks = 50)
h$density = h$counts/sum(h$counts)*100
plot(h, freq=FALSE, add = T, col = util$c_light_teal, border = 'white')  
for(stsp in stand_species_sp3[12]){
  util$plot_expectand_pushforward(samples[[paste0('omega_conc_sck[',stsp,']')]], 50,
                                  flim = c(0.6,1), ylim = c(0,30),
                                  display_name = bquote(omega[shock]^{concordant}),
                                  col = '#8f2727', add = TRUE)
}

util$plot_pairs_by_chain(samples[[paste0('omega_conc_sck[',stand_species_sp3[12],']')]], 'omega_conc_sck[stand]',
                         samples[['rho_sp[3]']], 'rho_sp[3]')

trees_substand <- which(data$stand_species_idxs == stand_species_sp3[12])
stand <- unique(data$stand_idxs[trees_substand])
util$plot_expectand_pushforward(samples[[paste0('phi_sck[',stand,']')]], 50,
                                flim = c(0,1), ylim = c(0,20),
                                display_name = bquote(phi[shock]),
                                col = '#8f2727')
