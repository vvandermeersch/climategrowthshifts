rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnw"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)

# Model data
data <- readRDS(file.path(wd, 'output/model', 'data_28june2025.rds'))

# Posterior quantification - diagnostics
fit <- readRDS(file.path(wd, 'output/model', 'fit_28june2025_nopooling.rds')) # run on Margot
diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit)
base_samples <- util$filter_expectands(samples,
                                       c('alpha', 
                                         'beta_gdd', 'beta_sm', 'beta_vpd',
                                         'rho_sp', 'gamma_sp',
                                         'f_tilde_sh',
                                         'rho_sh', 'kappa_sh',
                                         'sigma'),
                                       check_arrays=TRUE)
util$check_all_expectand_diagnostics(base_samples)

base_samples <- util$filter_expectands(samples,
                                       c('log_rw_pred'),
                                       check_arrays=TRUE)
util$check_all_expectand_diagnostics(base_samples)


# Investigate species 20 and 21
for(s in 1:data$N_species){
  util$plot_pairs_by_chain(samples[[paste0('gamma_sp[',s,']')]], paste0('gamma_sp[',s,']'), samples[[paste0('rho_sp[',s,']')]], paste0('rho_sp[',s,']'))
}

for(s in 1:data$N_species){
  util$plot_pairs_by_chain(samples[[paste0('beta_gdd[',s,']')]], paste0('beta_gdd'), samples[[paste0('beta_vpd[',s,']')]], paste0('beta_vpd[',s,']'))
}

par(mfrow=c(1, 3))
util$plot_expectand_pushforward(samples[['beta_gdd']], 200,
                                flim = c(-0.15,0.15),
                                display_name="beta_gdd")
# xs <- seq(-1, 1, 0.001)
# ys <- dnorm(xs, 0, log(1.8) / 2.57)
# lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['beta_sm']], 200,
                                flim = c(-0.15,0.15),
                                display_name="beta_sm")
# xs <- seq(-1, 1, 0.01)
# ys <- dnorm(xs, 0, log(1.8) / 2.57)
# lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['beta_vpd']], 200,
                                flim = c(-0.15,0.15),
                                display_name="beta_vpd")
# xs <- seq(-1, 1, 0.01)
# ys <- dnorm(xs, 0, log(1.8) / 2.57)
# lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['rho_sh']], 20,
                                # flim = c(1,2),
                                display_name="beta_vpd")


par(mfrow=c(1, 1), mar = c(4,4,1,1))
names <- sapply(1:data$N_species,
                function(sp) paste0('rho_sp[', sp, ']'))
util$plot_disc_pushforward_quantiles(samples, names,
                                     xlab="Species",
                                     ylab="Rho")

par(mfrow=c(1, 1), mar = c(4,4,1,1))
names <- sapply(1:data$N_species,
                function(sp) paste0('gamma_sp[', sp, ']'))
util$plot_disc_pushforward_quantiles(samples, names,
                                     xlab="Species",
                                     ylab="Gamma")

par(mfrow=c(1, 1), mar = c(4,4,1,1))
names <- sapply(1:data$N_species,
                function(sp) paste0('kappa_sh[', sp, ']'))
util$plot_disc_pushforward_quantiles(samples, names,
                                     xlab="Species",
                                     ylab="Kappa")


par(mfrow=c(3, 1), mar = c(4,4,1,1))
names <- sapply(1:data$N_species,
                function(sp) paste0('beta_gdd[', sp, ']'))
util$plot_disc_pushforward_quantiles(samples, names,
                                     xlab="Species",
                                     ylab="beta_gdd")
names <- sapply(1:data$N_species,
                function(sp) paste0('beta_sm[', sp, ']'))
util$plot_disc_pushforward_quantiles(samples, names,
                                     xlab="Species",
                                     ylab="beta_sm")
names <- sapply(1:data$N_species,
                function(sp) paste0('beta_vpd[', sp, ']'))
util$plot_disc_pushforward_quantiles(samples, names,
                                     xlab="Species",
                                     ylab="beta_vpd")




par(mfrow=c(3, 2))
util$plot_expectand_pushforward(samples[[paste0('rho_sp[', 8, ']')]], 25,
                                display_name="rho sp20")
util$plot_expectand_pushforward(samples[[paste0('gamma_sp[', 8, ']')]], 25,
                                display_name="gamma sp20")
util$plot_expectand_pushforward(samples[[paste0('rho_sp[', 20, ']')]], 25,
                                display_name="rho sp20")
util$plot_expectand_pushforward(samples[[paste0('gamma_sp[', 20, ']')]], 25,
                                display_name="gamma sp20")
util$plot_expectand_pushforward(samples[[paste0('rho_sp[', 21, ']')]], 25,
                                display_name="rho sp21")
util$plot_expectand_pushforward(samples[[paste0('gamma_sp[', 21, ']')]], 25,
                                display_name="gamma sp21")


par(mfrow=c(3, 2), mar = c(4,4,1,1))
trees_sp20 <- which(data$species_idxs == 20)
for (t in  sample(trees_sp20, 18)) {
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  rw_names <- sapply(idxs,
                     function(n) paste0('log_rw_pred[', n, ']'))
  util$plot_conn_pushforward_quantiles(samples, rw_names, data$years[idxs],
                                       xlab="Year", ylab="Log Ring Width Per mm", 
                                       display_ylim=c(-2, 10), display_xlim = c(1980, 2024))
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=1.0, col="white")
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=0.8, col="black")
  abline(v=2012, lwd=1, lty=2, col="#DDDDDD")
  abline(v=2016, lwd=1, lty=2, col="#DDDDDD")
}

  
  
par(mfrow=c(3, 2), mar = c(4,4,1,1))
trees_sp20 <- which(data$species_idxs == 20)
for (t in  sample(trees_sp20, 18)) {
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  rw_names <- sapply(idxs,
                     function(n) paste0('log_rw_pred[', n, ']'))
  util$plot_realizations(samples, rw_names, data$years[idxs],N_plots=50,
                         xlab="Year", ylab="Log Ring Width Per mm", 
                         display_ylim=c(-2, 10), display_xlim = c(1980, 2024))
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=1.0, col="white")
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=0.8, col="black")
  abline(v=2012, lwd=1, lty=2, col="#DDDDDD")
  abline(v=2016, lwd=1, lty=2, col="#DDDDDD")
}

  
trees_sp20 <- which(data$species_idxs == 20)
for (t in  sample(trees_sp20, 18)) {
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  rw_names <- sapply(idxs,
                     function(n) paste0('log_rw_pred[', n, ']'))
  par(mfrow=c(2, 2), mar = c(4,4,1,1))
  for(c in 1:4){
    ss <- lapply(samples, function(s) array(s[c,], dim = c(1,1024)))
    util$plot_conn_pushforward_quantiles(ss, rw_names, data$years[idxs],
                                         xlab="Year", ylab="Log Ring Width Per mm",
                                         main = paste0('Treeid ', t),
                                         display_ylim=c(-4, 8), display_xlim = c(1980, 2024))
    points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=1.0, col="white")
    points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=0.8, col="black")
  }
  
  #abline(v=2012, lwd=1, lty=2, col="#DDDDDD")
  #abline(v=2016, lwd=1, lty=2, col="#DDDDDD")
}



t <- 3102 # 3102
idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
rw_names <- sapply(idxs,
                   function(n) paste0('log_rw_pred[', n, ']'))
par(mfrow=c(2, 2), mar = c(4,4,1,1))
gstand <- uniq_stand_ids[data$stand_idxs[t]]
stand <- unique(raw_data[raw_data$tree_id == uniq_tree_ids[t], 'stand'])
for(c in 1:4){
  ss <- lapply(samples, function(s) array(s[c,], dim = c(1,1024)))
  util$plot_conn_pushforward_quantiles(ss, rw_names, data$years[idxs],
                                       xlab="Year", ylab="Log Ring Width Per mm",
                                       main = paste0('Treeid ', t, ', stand ', 
                                                     gstand, ' (',stand,')'),
                                       display_ylim=c(-2, 2), display_xlim = c(1980, 2002))
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=1.0, col="white")
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=0.8, col="black")
  text(paste0('median(rho)=', round(median(ss[["rho_sp[20]"]]),2)), x = 1980, y = 1.85, adj = 0)
  text(paste0('median(gamma)=', round(median(ss[["gamma_sp[20]"]]),2)), x = 1980, y = 1.7, adj = 0)
  abline(v = 1994, lty = 2, col = 'grey')
}
uniq_stand_ids[data$stand_idxs[t]]
uniq_tree_ids[t]
group_dataset[group_dataset$group_dataset, 'dataset']


for(t in which(data$species_idxs == 20 & data$stand_idxs == 138)){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  rw_names <- sapply(idxs,
                     function(n) paste0('log_rw_pred[', n, ']'))
  par(mfrow=c(2, 2), mar = c(4,4,1,1))
  gstand <- uniq_stand_ids[data$stand_idxs[t]]
  stand <- unique(raw_data[raw_data$tree_id == uniq_tree_ids[t], 'stand'])
  for(c in 1:4){
    ss <- lapply(samples, function(s) array(s[c,], dim = c(1,1024)))
    util$plot_conn_pushforward_quantiles(ss, rw_names, data$years[idxs],
                                         xlab="Year", ylab="Log Ring Width Per mm",
                                         main = paste0('Treeid ', t, ', stand ', 
                                                       gstand, ' (',stand,')'),
                                         display_ylim=c(-2, 2), display_xlim = c(1980, 2002))
    points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=1.0, col="white")
    points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=0.8, col="black")
    text(paste0('median(rho)=', round(median(ss[["rho_sp[20]"]]),2)), x = 1980, y = 1.85, adj = 0)
    text(paste0('species ', uniq_species_ids[data$species_idxs[t]]), x = 2002, y = 1.85, adj = 1)
    text(paste0('median(gamma)=', round(median(ss[["gamma_sp[20]"]]),2)), x = 1980, y = 1.7, adj = 0)
    abline(v = 1994, lty = 2, col = 'grey')
  }
}

for(t in which(data$species_idxs == 8 & data$stand_idxs == 138)){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  rw_names <- sapply(idxs,
                     function(n) paste0('log_rw_pred[', n, ']'))
  par(mfrow=c(2, 2), mar = c(4,4,1,1))
  gstand <- uniq_stand_ids[data$stand_idxs[t]]
  stand <- unique(raw_data[raw_data$tree_id == uniq_tree_ids[t], 'stand'])
  for(c in 1:4){
    ss <- lapply(samples, function(s) array(s[c,], dim = c(1,1024)))
    util$plot_conn_pushforward_quantiles(ss, rw_names, data$years[idxs],
                                         xlab="Year", ylab="Log Ring Width Per mm",
                                         main = paste0('Treeid ', t, ', stand ', 
                                                       gstand, ' (',stand,')'),
                                         display_ylim=c(-2, 2), display_xlim = c(1980, 2002))
    points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=1.0, col="white")
    points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=0.8, col="black")
    text(paste0('median(rho)=', round(median(ss[["rho_sp[8]"]]),2)), x = 1980, y = 1.85, adj = 0)
    text(paste0('species ', uniq_species_ids[data$species_idxs[t]]), x = 2002, y = 1.85, adj = 1)
    text(paste0('median(gamma)=', round(median(ss[["gamma_sp[8]"]]),2)), x = 1980, y = 1.7, adj = 0)
    abline(v = 1994, lty = 2, col = 'grey')
  }
}

for(t in which(data$species_idxs == 22 & data$stand_idxs == 138)){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  rw_names <- sapply(idxs,
                     function(n) paste0('log_rw_pred[', n, ']'))
  par(mfrow=c(2, 2), mar = c(4,4,1,1))
  gstand <- uniq_stand_ids[data$stand_idxs[t]]
  stand <- unique(raw_data[raw_data$tree_id == uniq_tree_ids[t], 'stand'])
  for(c in 1:4){
    ss <- lapply(samples, function(s) array(s[c,], dim = c(1,1024)))
    util$plot_conn_pushforward_quantiles(ss, rw_names, data$years[idxs],
                                         xlab="Year", ylab="Log Ring Width Per mm",
                                         main = paste0('Tree ', t, ', stand ', 
                                                       gstand, ' (',stand,')'),
                                         display_ylim=c(-2, 2), display_xlim = c(1980, 2002))
    points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=1.0, col="white")
    points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=0.8, col="black")
    text(paste0('median(rho)=', round(median(ss[["rho_sp[22]"]]),2)), x = 1980, y = 1.85, adj = 0)
    text(paste0('species ', uniq_species_ids[data$species_idxs[t]]), x = 2002, y = 1.85, adj = 1)
    text(paste0('median(gamma)=', round(median(ss[["gamma_sp[22]"]]),2)), x = 1980, y = 1.7, adj = 0)
    abline(v = 1994, lty = 2, col = 'grey')
  }
}

stand <- 138
par(mfrow=c(2, 2), mar = c(4,4,1,1))
for(c in 1:4){
  ss <- lapply(samples, function(s) array(s[c,], dim = c(1,1024)))
  names <- sapply(1:length(data$all_years),
                  function(n) paste0('f_sh[', stand,',', n, ']'))
  util$plot_realizations(ss, names, data$all_years,N_plots=50,
                         xlab="Year",
                         display_ylim=c(-3, 3), display_xlim = c(1980, 2002))
  text(paste0('median(rho_sh)=', round(median(ss[["rho_sh"]]),2)), x = 1980, y = 2.85, adj = 0)
  # text(paste0('median(gamma_sh)=', round(median(ss[["gamma_sh"]]),2)), x = 1980, y = 2.7, adj = 0)
  abline(v = 1994, lty = 2, col = 'grey')
}


# 
unique(data$species_idxs[data$stand_idxs == 132])

t <- 2990 #2975
idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
rw_names <- sapply(idxs,
                   function(n) paste0('log_rw_pred[', n, ']'))
par(mfrow=c(2, 2), mar = c(4,4,1,1))
gstand <- uniq_stand_ids[data$stand_idxs[t]]
stand <- unique(raw_data[raw_data$tree_id == uniq_tree_ids[t], 'stand'])
for(c in 1:4){
  ss <- lapply(samples, function(s) array(s[c,], dim = c(1,1024)))
  util$plot_conn_pushforward_quantiles(ss, rw_names, data$years[idxs],
                                       xlab="Year", ylab="Log Ring Width Per mm",
                                       main = paste0('Treeid ', t, ', stand ', 
                                                     gstand, ' (',stand,')'),
                                       display_ylim=c(-2, 2), display_xlim = c(1980, 2015))
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=1.0, col="white")
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=0.8, col="black")
  text(paste0('median(rho)=', round(median(ss[["rho_sp[20]"]]),2)), x = 1980, y = 1.85, adj = 0)
  text(paste0('species ', uniq_species_ids[data$species_idxs[t]]), x = 2015, y = 1.85, adj = 1)
  text(paste0('median(gamma)=', round(median(ss[["gamma_sp[20]"]]),2)), x = 1980, y = 1.50, adj = 0)
  abline(v = 1994, lty = 2, col = 'grey')
}


stand <- 132
par(mfrow=c(2, 2), mar = c(4,4,1,1))
for(c in 1:4){
  ss <- lapply(samples, function(s) array(s[c,], dim = c(1,1024)))
  names <- sapply(1:length(data$all_years),
                  function(n) paste0('f_sh[', stand,',', n, ']'))
  util$plot_realizations(ss, names, data$all_years,N_plots=50,
                         xlab="Year",
                         display_ylim=c(-3, 3), display_xlim = c(1980, 2015))
  text(paste0('median(rho_sh)=', round(median(ss[["rho_sh"]]),2)), x = 1980, y = 2.85, adj = 0)
  # text(paste0('median(gamma_sh)=', round(median(ss[["gamma_sh"]]),2)), x = 1980, y = 2.7, adj = 0)
  abline(v = 1994, lty = 2, col = 'grey')
}






# 
unique(data$species_idxs[data$stand_idxs == 85])

t <- 2118 #2975
idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
rw_names <- sapply(idxs,
                   function(n) paste0('log_rw_pred[', n, ']'))
par(mfrow=c(2, 2), mar = c(4,4,1,1))
gstand <- uniq_stand_ids[data$stand_idxs[t]]
stand <- unique(raw_data[raw_data$tree_id == uniq_tree_ids[t], 'stand'])
for(c in 1:4){
  ss <- lapply(samples, function(s) array(s[c,], dim = c(1,1024)))
  util$plot_conn_pushforward_quantiles(ss, rw_names, data$years[idxs],
                                       xlab="Year", ylab="Log Ring Width Per mm",
                                       main = paste0('Treeid ', t, ', stand ', 
                                                     gstand, ' (',stand,')'),
                                       display_ylim=c(4, 7), display_xlim = c(1980, 1999))
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=1.0, col="white")
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=0.8, col="black")
  text(paste0('median(rho)=', round(median(ss[["rho_sp[20]"]]),2)), x = 1980, y = 6.85, adj = 0)
  text(paste0('species ', uniq_species_ids[data$species_idxs[t]]), x = 1999, y = 6.85, adj = 1)
  text(paste0('median(gamma)=', round(median(ss[["gamma_sp[20]"]]),2)), x = 1980, y = 6.60, adj = 0)
  #abline(v = 1994, lty = 2, col = 'grey')
}


stand <- 132
par(mfrow=c(2, 2), mar = c(4,4,1,1))
for(c in 1:4){
  ss <- lapply(samples, function(s) array(s[c,], dim = c(1,1024)))
  names <- sapply(1:length(data$all_years),
                  function(n) paste0('f_sh[', stand,',', n, ']'))
  util$plot_realizations(ss, names, data$all_years,N_plots=50,
                         xlab="Year",
                         display_ylim=c(-3, 3), display_xlim = c(1980, 2015))
  text(paste0('median(rho_sh)=', round(median(ss[["rho_sh"]]),2)), x = 1980, y = 2.85, adj = 0)
  # text(paste0('median(gamma_sh)=', round(median(ss[["gamma_sh"]]),2)), x = 1980, y = 2.7, adj = 0)
  abline(v = 1994, lty = 2, col = 'grey')
}



medrhos <- 
  data.frame(
    rho = sapply(1:data$N_species, function(s) median(samples[[paste0('rho_sp[',s,']')]])),
    shortname = uniq_species_ids
  )
medrhos <- merge(medrhos, sppfull[,c('shortname', 'phylo.name')])
medrhos <- medrhos[match(phy.plants.here[["tip.label"]], medrhos$phylo.name),] 

phy.plants.here[["tip.label"]] <- paste0(phy.plants.here[["tip.label"]], ' (', round(medrhos$rho, 1), ')')

rbPal <- colorRampPalette(c('#8B174DFF','#8BA749FF'))
colors <- rbPal(14)[as.numeric(cut(medrhos$rho,breaks = 14))]

par(mfrow=c(1, 1), mar = c(2,2,1,1))

plot.phylo(phy.plants.here,
           cex = 1, 
           tip.color = colors)
legend("topright",title="median(rho)",legend=seq(4,32, 4),col=rbPal(8),pch=20, cex = 1.2)


medgammas <- 
  data.frame(
    gamma = sapply(1:data$N_species, function(s) median(samples[[paste0('gamma_sp[',s,']')]])),
    shortname = uniq_species_ids
  )
medgammas <- merge(medgammas, sppfull[,c('shortname', 'phylo.name')])
medgammas <- medgammas[match(phy.plants.here[["tip.label"]], medgammas$phylo.name),] 

phy.plants.here[["tip.label"]] <- paste0(phy.plants.here[["tip.label"]], ' (', round(medgammas$gamma, 1), ')')

rbPal <- colorRampPalette(c('#8B174DFF','#8BA749FF'))
colors <- rbPal(14)[as.numeric(cut(medgammas$gamma,breaks = 14))]

par(mfrow=c(1, 1), mar = c(2,2,0,-2))

plot.phylo(phy.plants.here,
           cex = 1, 
           tip.color = colors)
legend("topright",title="median(gamma)",legend=seq(0,1, 0.2),col=rbPal(6),pch=20, cex = 1)


medkappas <- 
  data.frame(
    kappa = sapply(1:data$N_species, function(s) median(samples[[paste0('kappa_sh[',s,']')]])),
    shortname = uniq_species_ids
  )
medkappas <- merge(medkappas, sppfull[,c('shortname', 'phylo.name')])
medkappas <- medkappas[match(phy.plants.here[["tip.label"]], medkappas$phylo.name),] 

phy.plants.here[["tip.label"]] <- paste0(phy.plants.here[["tip.label"]], ' (', round(medkappas$kappa, 1), ')')

rbPal <- colorRampPalette(c('#8B174DFF','#8BA749FF'))
colors <- rbPal(14)[as.numeric(cut(medkappas$kappa,breaks = 14))]

par(mfrow=c(1, 1), mar = c(2,2,0,0))

plot.phylo(phy.plants.here,
           cex = 1, 
           tip.color = colors)


ggplot(phy.plants.here) + 
  geom_tree() +
  geom_tiplab() +
  xlim(c(0,200))


tree <- groupOTU(tree,branches)
ggtree(tree) + geom_tiplab(aes(color=group))+
  scale_color_manual(values=c(A = "#E7B800",B= "#FC4E07", C = "darkgreen"))

test <- medgammas[, c('phylo.name', 'gamma')]
names(test)[1] <- 'id'

p2 <- facet_plot(p, panel="dot", data=medgammas, geom=geom_point, aes(x=gamma), color='red3')

d1 <- data.frame(id=phy.plants.here$tip.label, val=rnorm(30, sd=3))
facet_plot(p, panel="dot", data=d1, geom=geom_point, aes(x=val), color='red3')



tree <- rtree(30)

# Make the original plot
p <- ggtree(phy.plants.here) + geom_tiplab() + xlim(c(0,200))

# generate some random values for each tip label in the data
d1 <- medrhos[, c('phylo.name', 'rho')]

# Make a second plot with the original, naming the new plot "dot", 
# using the data you just created, with a point geom.
p2 <- facet_plot(p, panel="dot", data=d1, geom=geom_point, aes(x=rho), color='red3')





pred_names <- sapply(1:data$N, function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conditional_median_quantiles(samples, pred_names, data$gdd_obs,
                                       0, 32, 1, data$log_rw_obs, 
                                       residual=FALSE,
                                       xlab="GDD (x100degC)")


