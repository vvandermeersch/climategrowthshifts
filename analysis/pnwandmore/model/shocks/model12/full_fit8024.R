rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)

param_samples <- readRDS("/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model/paramsamples_29jan2025_gymnosperms_standclimate_19802924.rds")
data <- readRDS("/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/model/data_15jan2025_gymnosperms_standclimate_19802024.rds")

samples <- lapply(1:dim(param_samples)[3], function(k){t(matrix(param_samples[1:dim(param_samples)[1],1:dim(param_samples)[2],k], 
                                                                nrow = dim(param_samples)[1], ncol = dim(param_samples)[2]))})
names(samples) <- dimnames(param_samples)$variable
util$check_all_expectand_diagnostics(samples)

# Species-specific shock amplitude
par(mfrow = c(1,1))
plot(1, type="n", main='',
     xlim= c(0,8), xlab= bquote(tau[shock]),
     ylim= c(0,20), ylab="",  yaxt="n",
     bty = "n")
title(ylab="Estimated Bin\nProbabilities / Bin Width",
      mgp=c(1, 1, 0))
h <- hist(abs(rnorm(1e6, 0, log(20)/2.57)), plot = FALSE, breaks = 100)
h$density = h$counts/sum(h$counts)*100
plot(h, freq=FALSE, add = T, col = 'grey90', border = 'white')
for(sp in 1:data$N_species){
  util$plot_expectand_pushforward(samples[[paste0('tau_sck[',sp,']')]], 100,
                                  flim = c(0,8), ylim = c(0,20),
                                  display_name = bquote(tau[shock]),
                                  col = '#8f2727', add = TRUE)
}


util$plot_expectand_pushforward(samples[[paste0('phi_sck[',167,']')]], 100,
                                flim = c(0,1), ylim = c(0,15),
                                display_name = bquote(tau[shock]),
                                col = '#8f2727')
util$plot_expectand_pushforward(samples[[paste0('omega_conc_sck[',194,']')]], 100,
                                flim = c(0,1), ylim = c(0,15),
                                display_name = bquote(tau[shock]),
                                col = '#8f2727')

data$stand_idxs[which(data$stand_species_idxs == 194)]


util$plot_expectand_pushforward(samples[[paste0('rho_sp[',22,']')]], 100,
                                flim = c(0,10), ylim = c(0,5),
                                display_name = bquote(rho[species]),
                                col = '#8f2727')

util$plot_pairs_by_chain(samples[[paste0('rho_sp[',22,']')]],  bquote(rho[species]),
                         samples[[paste0('mu_rho[',1,']')]],  bquote(mu[rho]))



# Species-specific rho
par(mfrow = c(1,1))
plot(1, type="n", main='',
     xlim= c(0,30), xlab= bquote(tau[shock]),
     ylim= c(0,10), ylab="",  yaxt="n",
     bty = "n")
for(sp in 1:data$N_species){
  util$plot_expectand_pushforward(samples[[paste0('rho_sp[',sp,']')]], 100,
                                  flim = c(0,30), ylim = c(0,10),
                                  display_name = bquote(tau[shock]),
                                  col = '#8f2727', add = TRUE)
}


par(mfrow = c(1,1))
plot(1, type="n", main='',
     xlim= c(-0.2,0.2), xlab= bquote(tau[shock]),
     ylim= c(0,200), ylab="",  yaxt="n",
     bty = "n")
for(sp in 1:data$N_species){
  util$plot_expectand_pushforward(samples[[paste0('beta_pre[',sp,']')]], 100,
                                  flim = c(-0.2,0.2), ylim = c(0,200),
                                  display_name = bquote(tau[shock]),
                                  col = '#8f2727', add = TRUE)
}

par(mfrow = c(1,1))
plot(1, type="n", main='',
     xlim= c(-0.2,0.2), xlab= bquote(tau[shock]),
     ylim= c(0,200), ylab="",  yaxt="n",
     bty = "n")
for(sp in 1:data$N_species){
  util$plot_expectand_pushforward(samples[[paste0('beta_gdd[',sp,']')]], 100,
                                  flim = c(-0.2,0.2), ylim = c(0,200),
                                  display_name = bquote(tau[shock]),
                                  col = '#8f2727', add = TRUE)
}


par(mfrow = c(1,1))
plot(1, type="n", main='',
     xlim= c(-0.2,0.2), xlab= bquote(tau[shock]),
     ylim= c(0,200), ylab="",  yaxt="n",
     bty = "n")
for(sp in 1:data$N_species){
  util$plot_expectand_pushforward(samples[[paste0('beta_vpd[',sp,']')]], 100,
                                  flim = c(-0.2,0.2), ylim = c(0,200),
                                  display_name = bquote(tau[shock]),
                                  col = '#8f2727', add = TRUE)
}


# Stand-specific shock probability
stand <- which(sapply(1:data$N_stands, function(s) util$ensemble_mcmc_quantile_est(samples[[paste0('phi_sck[',s,']')]], c(0.5))) > 0.6)
substand <- unique(data$stand_species_idxs[which(data$stand_idxs == stand)])

par(mfrow = c(1,1), cex.lab = 0.85)
plot(1, type="n", main='',
       xlim= c(0,1), xlab= bquote(phi[shock]),
       ylim= c(0,20), ylab="",  yaxt="n",
       bty = "n")
title(ylab="Estimated Bin\nProbabilities / Bin Width",
        mgp=c(1, 1, 0))
h <- hist(rbeta(1e6, 2,20), plot = FALSE, breaks = 50)
h$density = h$counts/sum(h$counts)*100
for(s in 1:data$N_stands){
  util$plot_expectand_pushforward(samples[[paste0('phi_sck[',s,']')]], 50,
                                  flim = c(0,1), ylim = c(0,20),
                                  display_name = bquote(phi[shock]),
                                  col = ifelse(s == stand, 'darkorange', '#8f2727'), add = TRUE)
}
plot(h, freq=FALSE, add = T, col = util$c_light_teal, border = 'white')


par(mfrow = c(1,1), cex.lab = 0.85)
plot(1, type="n", main='',
       xlim= c(0,1), xlab= bquote(omega[shock]^{concordant}),
       ylim= c(0,40), ylab="",  yaxt="n",
       bty = "n")
title(ylab="Estimated Bin\nProbabilities / Bin Width",
        mgp=c(1, 1, 0))
h <- hist(rbeta(1e6, 230,14), plot = FALSE, breaks = 30)
h$density = h$counts/sum(h$counts)*100
for(stsp in 1:data$N_stand_species){
  util$plot_expectand_pushforward(samples[[paste0('omega_conc_sck[',stsp,']')]], 50,
                                  flim = c(0,1), ylim = c(0,40),
                                  display_name = bquote(omega[shock]^{concordant}),
                                  col = ifelse(stsp == substand, 'darkorange', '#8f2727'), add = TRUE)
}
plot(h, freq=FALSE, add = T, col = util$c_light_teal, border = 'white')


util$plot_expectand_pushforward(samples[[paste0('phi_sck[',stand,']')]], 50,
                                flim = c(0,1), ylim = c(0,20),
                                display_name = bquote(phi[shock]),
                                col = '#8f2727')

util$plot_expectand_pushforward(samples[[paste0('omega_conc_sck[',substand,']')]], 50,
                                flim = c(0,1), ylim = c(0,80),
                                display_name = bquote(omega[shock]),
                                col = '#8f2727')


phi_scks <- sapply(1:data$N_stands, function(s) util$ensemble_mcmc_quantile_est(samples[[paste0('phi_sck[',s,']')]], c(0.5)))

omega_scks <-  sapply(1:data$N_stands, function(s){
  substands <- unique(data$stand_species_idxs[which(data$stand_idxs == s)])
  mean(sapply(substands, function(sb) util$ensemble_mcmc_quantile_est(samples[[paste0('omega_conc_sck[',sb,']')]], c(0.5))))
})

plot(phi_scks ~ omega_scks, xlab = bquote(average~omega[shock]), ylab = bquote(phi[shock]),
     col = sapply(1:data$N_stands, function(s) ifelse(s == 33, 'darkorange', '#8f2727')),
     pch = 20)
