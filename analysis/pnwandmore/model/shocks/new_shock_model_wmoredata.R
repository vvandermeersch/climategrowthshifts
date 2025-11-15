rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('functions/custom_functions.R', local=util)
setwd(wd)

data <- readRDS(file.path(wd, 'output/model', 'data_08nov2025_azca_2species.rds'))
# samples <- readRDS(file.path(wd, "model/shocks/output", "samples_08nov2025_azca_2species.rds"))
fit <- readRDS(file.path(wd, "model/shocks/output", "fit_08nov2025_azca_2species.rds"))

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

pdf(file.path(wd, 'model/shocks/figures', 'divergences_azca_2species.pdf'), height = 8.5, width = 11.5)
par(mfrow=c(5, 6), cex.main = 1, mar = c(4,4,2,2))
divs <- diagnostics[['divergent__']]
C <- dim(divs)[1]
nondiv_filter <- c(sapply(1:C, function(c) divs[c,] == 0))
div_filter    <- c(sapply(1:C, function(c) divs[c,] == 1))
c_dark_trans <- c("#8F272780")
c_green_trans <- c("#00FF0080")
for(t in sample(1:data$N_trees, 90)){
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


util$plot_pairs_by_chain(samples[['sigma']], expression(sigma),
                         log(samples[['inner_tau_sck']]), expression(tau[inner]))


