rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('functions/custom_functions.R', local=util)
setwd(wd)

data <- readRDS(file.path(wd, 'output/model/climatena', 'data_15oct2025_long_az618.rds'))
data$w_sck <-  lapply(1, function(x)  ifelse(data$all_years == 2002, 1, 0.2))
fit <- readRDS(file.path(wd, "model/shocks/output", "fit_15oct2025_long_az618_shocks_v1.rds"))

# Computational diagnostics
diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit)
base_samples <- util$filter_expectands(samples,
                                       c('alpha', 
                                         'beta_gdd', 'beta_winterprec', 
                                         'rho_sp', 'gamma_sp',
                                         'f_tilde_sh',
                                         'rho_sh', 'gamma_sh',
                                         'delta_tilde_sck', 'inner_tau_sck', 'outer_tau_sck', 'gamma_sck',
                                         'sigma'),
                                       check_arrays=TRUE)
util$check_all_expectand_diagnostics(base_samples)

# 
for(y in c(1,10, 101, 102)){
  util$plot_div_pairs(paste0('delta_tilde_sck[1,',y,']'), 'inner_tau_sck', 
                      samples, diagnostics, transforms = list('inner_tau_sck' = 1))
}


