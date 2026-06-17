




fit <- readRDS(file.path(wd, 'output/model/model21', 'fit_model21_HGSP_PIPO_30stands_19502022_nomu_sametau.rds'))

model_gq <- cmdstan_model('~/projects/climategrowthshifts/analysis/pnwandmore/model/stan/model21/hgsp/model21_HSGP_1species_nomu_sametau_withGQ_onlyforfind.stan')
data <- readRDS(file.path(folder,'data_18may2026_PIPO_30stands_19502022.rds'))
data$uniq_tree_ids <- NULL
data$uniq_species_ids <- NULL
data$uniq_stand_ids <- NULL
data$N_clades <- 1
data$uniq_stand_lat <- NULL
data$uniq_stand_lon <- NULL

data$grainsize <-  1
data_gq <- data
fit_gq <- model_gq$generate_quantities(fit, data = data_gq, seed = 5838293, parallel_chains = 4)
gc()



samples <- fit$draws(c('alpha',
                       'beta_gdd', 'beta_vpd', 'beta_pre',
                       
                       'delta_clim', 
                       'tau_clim',
                       
                       'rho_sp', 'gamma_sp',
                       'rho_ind', 'gamma_ind',
                       
                       'rho_sh', 'gamma_sh', 
                       'f_tilde_sh','f_sh',
                       
                       'tau_conc',
                       'phi_sck',
                       'omega_conc_sck',
                       'omega_shutdown',
                       
                       'tau_idio',
                       'thetas_idio',
                       
                       'sigma'
))
names <- dimnames(samples)$variable
samples <- lapply(1:dim(samples)[3],
                  function(k){t(matrix(samples[1:dim(samples)[1],1:dim(samples)[2],k],
                                       nrow = dim(samples)[1], ncol = dim(samples)[2]))})
names(samples) <- names
util$check_all_expectand_diagnostics(samples)


samples <- fit_gq$draws(c('f_ind'))
names <- dimnames(samples)$variable
samples <- lapply(1:dim(samples)[3],
                  function(k){t(matrix(samples[1:dim(samples)[1],1:dim(samples)[2],k],
                                       nrow = dim(samples)[1], ncol = dim(samples)[2]))})
names(samples) <- names
util$check_all_expectand_diagnostics(samples)


