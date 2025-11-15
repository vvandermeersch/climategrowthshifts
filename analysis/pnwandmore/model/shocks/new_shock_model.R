
library(rstan)
options(mc.cores = 10)


data$N_stand_trees <- array(data$N_stand_trees)
data$N_stand_years <- array(data$N_stand_years)
data$stand_start_years_idxs <- array(data$stand_start_years_idxs)
data$stand_trees_start_idxs <- array(data$stand_trees_start_idxs)
data$stand_trees_end_idxs <- array(data$stand_trees_end_idxs)

# data$w_sck <-  lapply(1:data$N_stands, function(s){
#   if(s %in% 1:data$N_stands){ifelse(data$all_years %in%  c(2002), 1, 0.2)}
# })

data$w_sck <- ifelse(data$all_years_idxs == which(data$all_years == 2002), 1, 0.2)
pot_issues_idxs <- c(238, 927, 975, 979, 1241, 1277, 1301, 1331, 1365, 1367, 1601)
data$w_sck[pot_issues_idxs] <- 1

fit <- stan(file=file.path(wd, 'model/stan', 'model8_with3predictors_nopooling.stan'),
            data=data, seed=5838293, 
            chains = 4,
            warmup=1000, iter=2024, refresh=10)

samples <- util$extract_expectand_vals(fit)
diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

base_samples <- util$filter_expectands(samples,
                                       c('alpha', 
                                         'beta_gdd', 'beta_sm', 'beta_vpd',
                                         'rho_sp', 'gamma_sp',
                                         # 'f_sh',
                                         'rho_sh', 'kappa_sh',
                                         # 'delta_sck', 
                                         'inner_tau_sck', 'outer_tau_sck', 'phi_sck', 'omega_tree_sck',
                                         'sigma'),
                                       check_arrays=TRUE)
util$check_all_expectand_diagnostics(base_samples)

base_samples <- util$filter_expectands(samples,
                                       c('f_sh'),
                                       check_arrays=TRUE)
util$check_all_expectand_diagnostics(base_samples)


par(mfrow = c(2,1), mar = c(4,4,2,2), cex.main = 1)
util$plot_expectand_pushforward(samples[['phi_sck']], 50,
                                flim = c(0,1), main = 'Stand-level shock',
                                display_name=expression(phi[shock]))
xs <- seq(0, 1, 0.01)
ys <- dbeta(xs, 10,500)
lines(xs, ys, lwd=1.5, col=util$c_light_teal)
util$plot_expectand_pushforward(samples[['omega_tree_sck']], 50,
                                flim = c(0,1), main = 'Tree-level shock',
                                display_name=expression(omega[tree-shock]))
xs <- seq(0, 1, 0.01)
ys <- dbeta(xs, 5,5)
lines(xs, ys, lwd=1.5, col=util$c_light_teal)

par(mfrow = c(1,1), mar = c(4,4,2,2), cex.main = 1)
util$plot_expectand_pushforward(samples[['rho_sh']], 50,
                                flim = c(0,10), main = 'Tree-level shock',
                                display_name=expression(rho[short]))

util$plot_pairs_by_chain(samples[['sigma']], expression(sigma),
                         samples[['inner_tau_sck']], expression(tau[inner]))

util$plot_pairs_by_chain(samples[['outer_tau_sck']], expression(tau[outer]),
                         samples[['inner_tau_sck']], expression(tau[inner]))

util$plot_pairs_by_chain(samples[['delta_tilde_sck[979]']], 'delta_tilde_sck[979]',
                         samples[['inner_tau_sck']], expression(tau[inner]))


pdf(file.path(wd, 'model/shocks/figures', 'divergences_outer.pdf'), height = 8.5, width = 11.5)
par(mfrow=c(5, 6), cex.main = 1, mar = c(4,4,2,2))
divs <- diagnostics[['divergent__']]
C <- dim(divs)[1]
nondiv_filter <- c(sapply(1:C, function(c) divs[c,] == 0))
div_filter    <- c(sapply(1:C, function(c) divs[c,] == 1))
c_dark_trans <- c("#8F272780")
c_green_trans <- c("#00FF0080")
for(t in 1:data$N_trees){
  tree_idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  
  for(i in tree_idxs) {
    name_x <- paste("delta_sck[", i, "]", sep='')
    name_y <- 'inner_tau_sck'
    
    plot(samples[[name_x]][nondiv_filter], log(samples[[name_y]][nondiv_filter]),
         col=c_dark_trans, pch=16, main=paste0("Tree ", t, ', year ', data$all_years[data$all_years_idxs[i]]),
         xlab=name_x, # xlim = c(-1,1),
         ylab=paste("log(", name_y, ")", sep=""))
    
    points(samples[[name_x]][div_filter], log(samples[[name_y]][div_filter]),
           col=c_green_trans, pch=16)
  }
  
}
dev.off()

name_x <- paste("delta_tilde_sck[", 300, "]", sep='')
util$plot_pairs_by_chain(samples[[name_x]], name_x, 
                         log(samples[[name_y]]), name_y)



util$plot_div_pairs('beta_sm[1]', 'delta_sck[200]', samples, diagnostics)
util$plot_div_pairs('beta_sm[1]', 'beta_sm[2]', samples, diagnostics)
util$plot_div_pairs('alpha', 'beta_gdd[2]', samples, diagnostics)
util$plot_div_pairs('alpha', 'beta_sm[1]', samples, diagnostics)
util$plot_div_pairs('alpha', 'beta_sm[2]', samples, diagnostics)

util$plot_div_pairs('omega_tree_sck', paste("delta_tilde_sck[", 100, "]", sep=''), samples, diagnostics)


par(mfrow = c(3,3), mar = c(4,5,2,2), cex.main = 1)
random_trees <- sample(1:data$N_trees, 3)
for(t in random_trees){
  idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- sapply(idxs_t,
                  function(x) paste0('delta_sck[', x, ']'))
  util$plot_conn_pushforward_quantiles(
    samples, names, plot_xs = data$years[idxs_t], 
    main = paste0('Stand ', 1, ' - ', data$uniq_species_ids[data$species_idxs[t]], ' - tree ', t), 
    ylab = expression(delta[shocks]), display_ylim = c(-12,2), 
    display_xlim = data$all_years[c(1,length(data$all_years))])
}  
for(t in random_trees){
  idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- sapply(idxs_t,
                  function(x) paste0('log_rw_pred[', x, ']'))
  util$plot_conn_pushforward_quantiles(
    samples, names, plot_xs = data$years[idxs_t], 
    main = '', 
    ylab = 'Log ring width (mm)', display_ylim = c(-12,2),
    display_xlim = data$all_years[c(1,length(data$all_years))])
  points(data$years[idxs_t], data$log_rw_obs[idxs_t], pch=16, cex=1.0, col="white")
  points(data$years[idxs_t], data$log_rw_obs[idxs_t], pch=16, cex=0.8, col="black")
}  
for(t in random_trees){
  idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- sapply(idxs_t,
                  function(x) paste0('log_rw_pred[', x, ']'))
  util$plot_conn_pushforward_quantiles(
    samples, names, plot_xs = data$years[idxs_t], 
    baseline_values = data$log_rw_obs[idxs_t],
    main = '', 
    ylab = 'Residuals', display_ylim = c(-4,4), residual = TRUE,
    display_xlim = data$all_years[c(1,length(data$all_years))])
}  

par(mfrow = c(1,1), mar = c(4,5,2,2), cex.main = 1)
s <- 1
names <- sapply(1:30,
                function(sp) paste0('f_sh[', s, ',', sp, ']'))
util$plot_realizations(samples, names, plot_xs = data$all_years[1:30], N_plots = 50,
                       main = '', ylab = 'f shortGP',
                       display_xlim = data$all_years[c(1,30)], display_ylim = c(-5,2))


par(mfrow = c(2,3), mar = c(4,5.5,2,1), cex.lab = 1.2, cex.main = 1)
util$plot_disc_pushforward_quantiles(samples, paste0('rho_sp[',1:data$N_species,']'),
                                     xticklabs = data$uniq_species_ids, ylab = expression(rho[species]))
util$plot_disc_pushforward_quantiles(samples, paste0('gamma_sp[',1:data$N_species,']'),
                                     xticklabs = data$uniq_species_ids, ylab = expression(gamma[species]))
util$plot_disc_pushforward_quantiles(samples, paste0('gamma_sck[',1:data$N_species,']'),
                                     xticklabs = data$uniq_species_ids, ylab = expression(phi[shock]))
util$plot_disc_pushforward_quantiles(samples, paste0('beta_gdd[',1:data$N_species,']'),
                                     xticklabs = data$uniq_species_ids, ylab = expression(beta[GDD]))
util$plot_disc_pushforward_quantiles(samples, paste0('beta_vpd[',1:data$N_species,']'),
                                     xticklabs = data$uniq_species_ids, ylab = expression(beta[VPD]))
util$plot_disc_pushforward_quantiles(samples, paste0('beta_sm[',1:data$N_species,']'),
                                     xticklabs = data$uniq_species_ids, ylab = expression(beta[SM]))

