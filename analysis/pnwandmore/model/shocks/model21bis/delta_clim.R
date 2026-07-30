
folder <- '/home/victor/projects/climategrowthshifts/analysis/pnwandmore/figures/model26july2026'
pdf(file.path(folder, 'delta_clim.pdf'), width = 10, height = 20)
par(mfrow = c(10,2), mar = c(4,4,1,1))
for(s in 1:data$N_stands){
  names <- paste0('delta_clim[',s,',',1:data$N_all_years,']')
  util$plot_conn_pushforward_quantiles(base_samples, names, data$all_years,
                                       display_ylim = c(-1,1))
}
dev.off()

phi_sck <- sapply(1:data$N_stands, function(s) util$ensemble_mcmc_quantile_est(base_samples[[paste0('phi_sck[',s,']')]], c(0.5)))
which(phi_sck > 0.95)


data$species_idxs[which(data$stand_idxs == 176)]
data$uniq_species_ids[unique(data$species_idxs[which(data$stand_idxs == 176)])]
data$stand_species_idxs[which(data$stand_idxs == 176)]

util$plot_pairs_by_chain(base_samples[['phi_sck[176]']], 'phi_sck[176]',
                         base_samples[['omega_conc_sck[193]']], 'omega_conc_sck[193]')
util$plot_pairs_by_chain(base_samples[['phi_sck[176]']], 'phi_sck[176]',
                         base_samples[['omega_conc_sck[197]']], 'omega_conc_sck[197]')

data$species_idxs[which(data$stand_idxs == 205)]
data$uniq_species_ids[unique(data$species_idxs[which(data$stand_idxs == 205)])]
data$stand_species_idxs[which(data$stand_idxs == 205)]

util$plot_pairs_by_chain(base_samples[['phi_sck[205]']], 'phi_sck[205]',
                         base_samples[['omega_conc_sck[228]']], 'omega_conc_sck[228]')

data$species_idxs[which(data$stand_idxs == 326)]
data$uniq_species_ids[unique(data$species_idxs[which(data$stand_idxs == 326)])]
data$stand_species_idxs[which(data$stand_idxs == 326)]

util$plot_pairs_by_chain(base_samples[['phi_sck[326]']], 'phi_sck[326]',
                         base_samples[['omega_conc_sck[361]']], 'omega_conc_sck[361]')
