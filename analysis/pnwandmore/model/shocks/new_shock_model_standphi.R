
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"

library(rstan)
options(mc.cores = parallel::detectCores())

data <- readRDS(file = file.path(wd, 'output', 'model', 'data_15nov2025_azca_pipo2022.rds'))

# First, we start with a full non-centered parametrisation #

data$w_sck <- rep(0, data$N)
fit_nc <- stan(file=file.path(wd, 'model', 'stan', 'model8_with3predictors_nopooling_standphi.stan'),
            data=data, seed=5838293,
            chains = 4,
            warmup=1000, iter=2024, refresh=10)
# 2976.15 seconds
diagnostics <- util$extract_hmc_diagnostics(fit_nc)
util$check_all_hmc_diagnostics(diagnostics)
saveRDS(fit_nc, file = file.path(wd, 'model', 'shocks', 'output', 'pipo2022', 'fit_nc.rds'))
samples_nc <- util$extract_expectand_vals(fit_nc)

# Look at divergences?
par(mfrow=c(5, 6), cex.main = 1, mar = c(4,4,2,2))
divs <- diagnostics[['divergent__']]
C <- dim(divs)[1]
nondiv_filter <- c(sapply(1:C, function(c) divs[c,] == 0))
div_filter    <- c(sapply(1:C, function(c) divs[c,] == 1))
c_dark_trans <- c("#8F272780")
c_green_trans <- c("#00FF0080")
for(t in sample(1:data$N_trees,3)){
  tree_idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  
  for(i in tree_idxs) {
    name_x <- paste("delta_sck[", i, "]", sep='')
    name_y <- 'inner_tau_sck'
    
    plot(samples_nc[[name_x]][nondiv_filter], log(samples_nc[[name_y]][nondiv_filter]),
         col=c_dark_trans, pch=16, main=paste0("Tree ", t, ', year ', data$all_years[data$all_years_idxs[i]]),
         xlab=name_x, # xlim = c(-1,1),
         ylab=paste("log(", name_y, ")", sep=""))
    
    points(samples_nc[[name_x]][div_filter], log(samples_nc[[name_y]][div_filter]),
           col=c_green_trans, pch=16)
  }
}


# Look at divergences, only one year
par(mfrow=c(5, 6), cex.main = 1, mar = c(4.5,5,2,2), cex.lab = 1.2)
divs <- diagnostics[['divergent__']]
C <- dim(divs)[1]
nondiv_filter <- c(sapply(1:C, function(c) divs[c,] == 0))
div_filter    <- c(sapply(1:C, function(c) divs[c,] == 1))
c_dark_trans <- c("#8F272780")
c_green_trans <- c("#00FF0080")
for(t in 1:data$N_trees){
  tree_idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  tree_idxs <- tree_idxs[which(data$all_years[data$all_years_idxs[tree_idxs]] == 2002)]
  for(i in tree_idxs) {
    name_x <- paste("delta_tilde_sck[", i, "]", sep='')
    name_y <- 'inner_tau_sck'
    
    plot(samples_nc[[name_x]][nondiv_filter], log(samples_nc[[name_y]][nondiv_filter]),
         col=c_dark_trans, pch=16, main=paste0("Tree ", t, ', year ', data$all_years[data$all_years_idxs[i]]),
         xlab=name_x, # xlim = c(-1,1),
         ylab=paste("log(", name_y, ")", sep=""))
    
    points(samples_nc[[name_x]][div_filter], log(samples_nc[[name_y]][div_filter]),
           col=c_green_trans, pch=16)
  }
}


nf <- layout(
  matrix(c(0,0,1,1,0,0,
           2,3,4,5,6,7,
           8,9,10,11,12,13,
           14,15,16,17,18,19,
           20,21,22,23,24,25), 
         ncol=6, byrow=TRUE)
)

nf <- layout(
  matrix(c(0,0,1,1,0,0,
           2,3,4,5,6,7), 
         ncol=6, byrow=TRUE)
)
stand <- 1
trees_stand <- which(data$stand_idxs == stand)
years <- 1:data$N_all_years
names <- sapply(years,
                function(x) paste0('f_sh[', stand,',',x, ']'))
util$plot_conn_pushforward_quantiles(
  samples_nc, names, plot_xs = data$all_years, 
  main = paste0('Stand ', stand), 
  ylab = expression(f[short]), display_ylim = c(-4,2), 
  display_xlim = data$all_years[c(1,length(data$all_years))])
abline(v = 2002, lty = 2)
random_trees <- sample(trees_stand, 6)
for(t in random_trees){
  
  idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- sapply(idxs_t,
                  function(x) paste0('delta_sck[', x, ']'))
  util$plot_conn_pushforward_quantiles(
    samples_nc, names, plot_xs = data$years[idxs_t], 
    main = paste0('Stand ', stand, ' - ', data$uniq_species_ids[data$species_idxs[t]], ' - tree ', t), 
    ylab = "delta_tilde_sck", display_ylim = c(-4,2), 
    display_xlim = data$all_years[c(1,length(data$all_years))])
  abline(v = 1990, lty = 2)
  abline(v = 2002, lty = 2)
}  

par(mfrow=c(1, 3), cex.main = 1, mar = c(4.5,5,2,2), cex.lab = 1.2)
util$plot_expectand_pushforward(samples_nc[['omega_tree_sck']], 50,
                                flim = c(0,1), 
                                display_name=expression(omega[tree_shock]))
util$plot_expectand_pushforward(samples_nc[['phi_sck[1]']], 50,
                                flim = c(0,1), 
                                display_name=expression(phi[sck_stand1]))
util$plot_expectand_pushforward(samples_nc[['phi_sck[2]']], 50,
                                flim = c(0,1), 
                                display_name=expression(phi[sck_stand2]))




# Still non-centered, but we allow trees to have nonconcordant shocks #

data$w_sck <- rep(0, data$N)
fit_nc_nonconc <- stan(file=file.path(wd, 'model', 'stan', 'model8_with3predictors_nopooling_standphi_nonconcordantshock.stan'),
               data=data, seed=5838293,
               chains = 4,
               warmup=1000, iter=2024, refresh=10)
# 3185.37 seconds
diagnostics <- util$extract_hmc_diagnostics(fit_nc_nonconc)
util$check_all_hmc_diagnostics(diagnostics)
saveRDS(fit_nc_nonconc, file = file.path(wd, 'model', 'shocks', 'output', 'pipo2022', 'fit_nc_nonconc.rds'))
samples_nc_nonconc <- util$extract_expectand_vals(fit_nc_nonconc)


par(mfrow=c(1, 4), cex.main = 1, mar = c(4.5,5,2,2), cex.lab = 1.2)

util$plot_expectand_pushforward(samples_nc_nonconc[['phi_sck[1]']], 50,
                                flim = c(0,1), 
                                display_name=expression(phi[sck_stand1]))
util$plot_expectand_pushforward(samples_nc_nonconc[['phi_sck[2]']], 50,
                                flim = c(0,1), 
                                display_name=expression(phi[sck_stand2]))
util$plot_expectand_pushforward(samples_nc_nonconc[['omega_nonconc_sck']], 50,
                                flim = c(0,1), 
                                display_name=expression(omega[tree_nonconc_shock]))
util$plot_expectand_pushforward(samples_nc_nonconc[['omega_conc_sck']], 50,
                                flim = c(0,1), 
                                display_name=expression(omega[tree_conc_shock]))

util$plot_pairs_by_chain(samples_nc_nonconc[['phi_sck[1]']], expression(phi[sck_stand1]),
                         samples_nc_nonconc[['omega_conc_sck']], expression(omega[tree_conc_shock]))

util$plot_pairs_by_chain(samples_nc_nonconc[['phi_sck[2]']], expression(phi[sck_stand1]),
                         samples_nc_nonconc[['omega_conc_sck']], expression(omega[tree_conc_shock]))

util$plot_pairs_by_chain(samples_nc_nonconc[['omega_nonconc_sck']], expression(omega[tree_nonconc_shock]),
                         samples_nc_nonconc[['omega_conc_sck']], expression(omega[tree_conc_shock]))




# Still non-centered, but we allow trees to have nonconcordant shocks #
# Here we put more domain expertise

data$w_sck <- rep(0, data$N)
fit_nc_nonconc2 <- stan(file=file.path(wd, 'model', 'stan', 'model8_with3predictors_nopooling_standphi_nonconcordantshock.stan'),
                       data=data, seed=5838293,
                       chains = 4,
                       warmup=1000, iter=2024, refresh=10)
# 4716.6 seconds, gosh
diagnostics <- util$extract_hmc_diagnostics(fit_nc_nonconc2)
util$check_all_hmc_diagnostics(diagnostics)
saveRDS(fit_nc_nonconc2, file = file.path(wd, 'model', 'shocks', 'output', 'pipo2022', 'fit_nc_nonconc2.rds'))
samples_nc_nonconc2 <- util$extract_expectand_vals(fit_nc_nonconc2)
base_samples <- util$filter_expectands(samples_nc_nonconc2,
                                       c('alpha', 
                                         'beta_gdd', 'beta_sm', 'beta_vpd',
                                         'rho_sp', 'gamma_sp',
                                         'rho_sh', 'kappa_sh',
                                         'inner_tau_sck', 'outer_tau_sck', 'phi_sck', 
                                         'omega_conc_sck', 'omega_nonconc_sck',
                                         'sigma'),
                                       check_arrays=TRUE)
util$check_all_expectand_diagnostics(base_samples)


par(mfrow=c(1, 3), cex.main = 1, mar = c(4.5,5,2,2), cex.lab = 1.2)

util$plot_expectand_pushforward(samples_nc_nonconc2[['phi_sck[1]']], 50,
                                flim = c(0,1), 
                                display_name=expression(phi[sck_stand]))
util$plot_expectand_pushforward(samples_nc_nonconc2[['phi_sck[2]']], 50,
                                flim = c(0,1), 
                                display_name=expression(phi[sck_stand]), add = T, col=util$c_mid)
util$plot_expectand_pushforward(samples_nc_nonconc2[['omega_conc_sck']], 50,
                                flim = c(0,1), 
                                display_name=expression(omega[tree_conc_shock]))
util$plot_expectand_pushforward(samples_nc_nonconc2[['omega_nonconc_sck']], 50,
                                flim = c(0,1), 
                                display_name=expression(omega[tree_nonconc_shock]))

nf <- layout(
  matrix(c(0,0,1,1,0,0,
           2,3,4,5,6,7), 
         ncol=6, byrow=TRUE)
)
stand <- 1
trees_stand <- which(data$stand_idxs == stand)
years <- 1:data$N_all_years
names <- sapply(years,
                function(x) paste0('f_sh[', stand,',',x, ']'))
util$plot_conn_pushforward_quantiles(
  samples_nc_nonconc2, names, plot_xs = data$all_years, 
  main = paste0('Stand ', stand), 
  ylab = expression(f[short]), display_ylim = c(-4,2), 
  display_xlim = data$all_years[c(1,length(data$all_years))])
abline(v = 2002, lty = 2)
random_trees <- sample(trees_stand, 6)
for(t in random_trees){
  
  idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- sapply(idxs_t,
                  function(x) paste0('delta_sck[', x, ']'))
  util$plot_conn_pushforward_quantiles(
    samples_nc_nonconc2, names, plot_xs = data$years[idxs_t], 
    main = paste0('Stand ', stand, ' - ', data$uniq_species_ids[data$species_idxs[t]], ' - tree ', t), 
    ylab = "delta_sck", display_ylim = c(-4,2), 
    display_xlim = data$all_years[c(1,length(data$all_years))])
  abline(v = 2002, lty = 2)
}  


# Look at divergences
pdf(file.path(wd, 'model/shocks/figures', 'divergences_stanphi_nonconc_sigma.pdf'), height = 8.5, width = 11.5)
for(t in 6){
  tree_idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names_x <- paste("delta_tilde_sck[", tree_idxs, "]", sep='')
  
  util$plot_div_pairs(names_x, 'sigma', samples_nc_nonconc2, diagnostics, transforms = list('sigma' = 1, 'inner_tau_sck' = 1))
  
}
dev.off()

pdf(file.path(wd, 'model/shocks/figures', 'divergences_stanphi_nonconc_innertau.pdf'), height = 8.5, width = 11.5)
for(t in 6){
  tree_idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names_x <- paste("delta_tilde_sck[", tree_idxs, "]", sep='')
  
  util$plot_div_pairs(names_x, 'inner_tau_sck', samples_nc_nonconc2, diagnostics, transforms = list('sigma' = 1, 'inner_tau_sck' = 1))
  
}
dev.off()


util$plot_div_pairs('omega_conc_sck', 'phi_sck[2]', samples_nc_nonconc2, diagnostics)

util$plot_div_pairs(paste("delta_tilde_sck[", 238, "]", sep=''), 'sigma', samples_nc_nonconc2, diagnostics, transforms = list('sigma' = 1, 'inner_tau_sck' = 1))


names <-  paste("delta_tilde_sck[", 1:data$N, "]", sep='')
median_shocks <- sapply(names, function(x) util$ensemble_mcmc_quantile_est(samples_nc_nonconc2[[x]], 0.5))
inner_tau_median <- util$ensemble_mcmc_quantile_est(samples_nc_nonconc2[['inner_tau_sck']], 0.5)
shocks_identification <- data.frame(median_shocks = median_shocks, big = median_shocks > 5 * inner_tau_median, small = median_shocks > 3 * inner_tau_median)



# Based on previous results, we tune w_sck#
# Here we put more domain expertise

data$w_sck <- ifelse(!shocks_identification$big, 0, 1)
fit_nc_nonconc3 <- stan(file=file.path(wd, 'model', 'stan', 'model8_with3predictors_nopooling_standphi_nonconcordantshock.stan'),
                        data=data, seed=5838293,
                        chains = 4,
                        warmup=1000, iter=2024, refresh=10)
samples_nc_nonconc3 <- util$extract_expectand_vals(fit_nc_nonconc3)
saveRDS(fit_nc_nonconc3, file = file.path(wd, 'model', 'shocks', 'output', 'pipo2022', 'fit_nc_nonconc3.rds'))
diagnostics <- util$extract_hmc_diagnostics(fit_nc_nonconc3)
util$plot_div_pairs('sigma', 'inner_tau_sck', samples_nc_nonconc3, diagnostics)
util$plot_pairs_by_chain(samples_nc_nonconc3[['sigma']]^2, 'sigma', samples_nc_nonconc3[['inner_tau_sck']]^2, 'inner_tau_sck')

util$plot_div_pairs('sigma', 'inner_tau_sck', samples_nc_nonconc3, diagnostics)

pdf(file.path(wd, 'model/shocks/figures', 'divergences_stanphi_nonconc_innertau_3.pdf'), height = 8.5, width = 11.5)
for(t in 1:6){
  tree_idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names_x <- paste("delta_tilde_sck[", tree_idxs, "]", sep='')
  
  util$plot_div_pairs(names_x, 'inner_tau_sck', samples_nc_nonconc2, diagnostics, transforms = list('sigma' = 1, 'inner_tau_sck' = 1))
  
}
dev.off()
