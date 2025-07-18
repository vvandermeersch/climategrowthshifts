rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)

# Model data
data <- readRDS(file.path(wd, 'output/model', 'data_03july2025.rds'))

# Posterior quantification
fit <- readRDS(file.path(wd, 'output/model', 'fit_03july2025_partialpooling_2clades_centered.rds')) # run on Margot
diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit)
base_samples <- util$filter_expectands(samples,
                                       c('alpha', 
                                         'mu_gdd', 'mu_sm', 'mu_vpd',
                                         'tau_gdd', 'tau_sm', 'tau_vpd',
                                         'beta_gdd', 'beta_sm', 'beta_vpd',
                                         'rho', 'gamma',
                                         'f_tilde_sh',
                                         'rho_sh', 'kappa_sh',
                                         'sigma'),
                                       check_arrays=TRUE)
util$check_all_expectand_diagnostics(base_samples)

# Diagnostic investigation

## species slopes
s <- 10
util$plot_pairs_by_chain(samples[[paste0('beta_gdd[',s,']')]], paste0('beta_gdd[',s,']'), samples[[paste0('beta_sm[',s,']')]], paste0('beta_sm[',s,']'))
util$plot_pairs_by_chain(samples[[paste0('beta_gdd[',s,']')]], paste0('beta_gdd[',s,']'), samples[[paste0('beta_vpd[',s,']')]], paste0('beta_vpd[',s,']'))
util$plot_pairs_by_chain(samples[[paste0('beta_vpd[',s,']')]], paste0('beta_vpd[',s,']'), samples[[paste0('beta_sm[',s,']')]], paste0('beta_sm[',s,']'))


## sigma
par(mfrow=c(1, 1), mar = c(4,4,1,1))
util$plot_expectand_pushforward(samples[[paste0('sigma')]], 25, display_name=expression(sigma))
par(mfrow=c(2, 2), mar = c(4,4,1,1))
for(c in 1:4){
  ss <- lapply(samples, function(s) array(s[c,], dim = c(1,1024)))
  util$plot_expectand_pushforward(ss[[paste0('sigma')]], 25, display_name=expression(sigma))
}
traceplot(fit, pars = c("sigma"))


## geoemtry
speciescount <- data.frame(table(data$species_idxs), clade = data$clade_idxs)
plot(Freq ~ Var1, data = speciescount, col = ifelse(speciescount$clade ==1, '#27278f', '#278f27'), type = 'b', pch = '1', xlab = 'species')
points(Freq ~ Var1, data = speciescount, col = ifelse(speciescount$clade ==1, '#27278f', '#278f27'), pch = 19, cex = 1.2)

for(c in 1:2){
  util$plot_pairs_by_chain(samples[[paste0('mu_gdd[',c,']')]], paste0('mu_gdd[',c,']'), log(samples[[paste0('tau_gdd[',c,']')]]), paste0('log(tau_gdd[',c,'])'))
  util$plot_pairs_by_chain(samples[[paste0('mu_sm[',c,']')]], paste0('mu_sm[',c,']'), log(samples[[paste0('tau_sm[',c,']')]]), paste0('log(tau_sm[',c,'])'))
  util$plot_pairs_by_chain(samples[[paste0('mu_vpd[',c,']')]], paste0('mu_vpd[',c,']'), log(samples[[paste0('tau_vpd[',c,']')]]), paste0('log(tau_vpd[',c,'])'))
}

par(mfrow=c(3, 2), mar = c(4,4,1,1))
for(c in 1:data$N_clades){
  colc <- ifelse(c ==1, '#27278f', '#278f27')
  name <- paste0("mu_gdd[", c, "]")
  namex <- ifelse(c == 1, "mu_gdd - gymno", "mu_gdd - angio")
  namey <-ifelse(c == 1, "log(tau_gdd) - gymno", "log(tau_gdd) - angio")
  plot(samples[[name]], log(samples[[paste0('tau_gdd[',c,']')]]),
       col=scales::alpha(colc, 0.2), pch=16, cex=0.8,
       xlab=namex, ylab=namey)
}
for(c in 1:data$N_clades){
  colc <- ifelse(c ==1, '#27278f', '#278f27')
  name <- paste0("mu_sm[", c, "]")
  namex <- ifelse(c == 1, "mu_sm - gymno", "mu_sm - angio")
  namey <-ifelse(c == 1, "log(tau_sm) - gymno", "log(tau_sm) - angio")
  plot(samples[[name]], log(samples[[paste0('tau_sm[',c,']')]]),
       col=scales::alpha(colc, 0.2), pch=16, cex=0.8,
       xlab=namex, ylab=namey)
}
for(c in 1:data$N_clades){
  colc <- ifelse(c ==1, '#27278f', '#278f27')
  name <- paste0("mu_vpd[", c, "]")
  namex <- ifelse(c == 1, "mu_vpd - gymno", "mu_vpd - angio")
  namey <-ifelse(c == 1, "log(tau_vpd) - gymno", "log(tau_vpd) - angio")
  plot(samples[[name]], log(samples[[paste0('tau_vpd[',c,']')]]),
       col=scales::alpha(colc, 0.2), pch=16, cex=0.8,
       xlab=namex, ylab=namey)
}

par(mfrow=c(4, 11), mar = c(4,4,1,1))
for(s in 1:data$N_species){
  c <- data$clade_idxs[s]
  colc <- ifelse(c ==1, '#27278f', '#278f27')
  name <- paste0("beta_gdd[", s, "]")
  namey <-ifelse(c == 1, "log(tau_gdd) - gymno", "log(tau_gdd) - angio")
  plot(samples[[name]], log(samples[[paste0('tau_gdd[',c,']')]]),
       col=scales::alpha(colc, 0.2), pch=16, cex=0.8,
       xlab=name, ylab=namey)
  if(speciescount[s,'Freq'] < 50){
    text(x = min(samples[[name]])+0.01, y =  max(log(samples[[paste0('tau_gdd[',c,']')]]))-0.1, 
         label = "x", col = 'darkred', cex = 1.5)
  }
}
for(s in 1:data$N_species){
  c <- data$clade_idxs[s]
  colc <- ifelse(c ==1, '#27278f', '#278f27')
  name <- paste0("beta_sm[", s, "]")
  namey <-ifelse(c == 1, "log(tau_sm) - gymno", "log(tau_sm) - angio")
  plot(samples[[name]], log(samples[[paste0('tau_sm[',c,']')]]),
       col=scales::alpha(colc, 0.2), pch=16, cex=0.8,
       xlab=name, ylab=namey)
  if(speciescount[s,'Freq'] < 50){
    text(x = min(samples[[name]])+0.01, y =  max(log(samples[[paste0('tau_sm[',c,']')]]))-0.1, 
         label = "x", col = 'darkred', cex = 1.5)
  }
}
for(s in 1:data$N_species){
  c <- data$clade_idxs[s]
  colc <- ifelse(c ==1, '#27278f', '#278f27')
  name <- paste0("beta_vpd[", s, "]")
  namey <-ifelse(c == 1, "log(tau_vpd) - gymno", "log(tau_vpd) - angio")
  plot(samples[[name]], log(samples[[paste0('tau_vpd[',c,']')]]),
       col=scales::alpha(colc, 0.2), pch=16, cex=0.8,
       xlab=name, ylab=namey)
  if(speciescount[s,'Freq'] < 50){
    text(x = min(samples[[name]])+0.01, y =  max(log(samples[[paste0('tau_vpd[',c,']')]]))-0.1, 
         label = "x", col = 'darkred', cex = 1.5)
  }
}


# Posterior inference
par(mfrow=c(1, 3), mar = c(4,4,1,1))
util$plot_expectand_pushforward(samples[['mu_gdd[1]']], 50,
                                flim = c(-0.15,0.2),
                                display_name=expression(beta[GDD]),
                                col = '#27278f')
util$plot_expectand_pushforward(samples[['mu_gdd[2]']], 50,
                                flim = c(-0.15,0.2),
                                col = '#278f27',
                                add = TRUE)
text(x = 0.115, y = 50, label = "Gymno.", col = '#27278f', cex = 1.2)
text(x = -0.08, y = 10, label = "Angio.", col = '#278f27', cex = 1.2)

util$plot_expectand_pushforward(samples[['mu_sm[1]']], 50,
                                flim = c(-0.15,0.2),
                                display_name=expression(beta[SM]),
                                col = '#27278f')
util$plot_expectand_pushforward(samples[['mu_sm[2]']], 50,
                                flim = c(-0.15,0.2),
                                col = '#278f27',
                                add = TRUE)

util$plot_expectand_pushforward(samples[['mu_vpd[1]']], 50,
                                flim = c(-0.15,0.2),
                                display_name=expression(beta[VPD]),
                                col = '#27278f')
util$plot_expectand_pushforward(samples[['mu_vpd[2]']], 50,
                                flim =c(-0.15,0.2),
                                col = '#278f27',
                                add = TRUE)

par(mfrow=c(1, 3), mar = c(4,4,1,1))
util$plot_expectand_pushforward(samples[['tau_gdd[1]']], 50,
                                flim = c(0,0.20),
                                display_name=expression(tau[GDD]),
                                col = '#27278f')
util$plot_expectand_pushforward(samples[['tau_gdd[2]']], 50,
                                flim =  c(0,0.20),
                                col = '#278f27',
                                add = TRUE)
text(x = 0.05, y = 50, label = "Gymno.", col = '#27278f', cex = 1.2)
text(x = 0.1, y = 10, label = "Angio.", col = '#278f27', cex = 1.2)

util$plot_expectand_pushforward(samples[['tau_sm[1]']], 50,
                                flim = c(0,0.20),
                                display_name=expression(tau[SM]),
                                col = '#27278f')
util$plot_expectand_pushforward(samples[['tau_sm[2]']], 50,
                                flim =  c(0,0.20),
                                col = '#278f27',
                                add = TRUE)

util$plot_expectand_pushforward(samples[['tau_vpd[1]']], 50,
                                flim =  c(0,0.20),
                                display_name=expression(tau[VPD]),
                                col = '#27278f')
util$plot_expectand_pushforward(samples[['tau_vpd[2]']], 50,
                                flim =  c(0,0.20),
                                col = '#278f27',
                                add = TRUE)



par(mfrow=c(1, 1), mar = c(4,4,1,1))
names <- sapply(1:data$N_species,
                function(sp) paste0('rho_sp[', sp, ']'))
util$plot_disc_pushforward_quantiles(samples, names,
                                     xlab="Species",
                                     ylab="Rho")


par(mfrow=c(1, 1), mar = c(4,4,1,1))
names <- sapply(1:data$N_species,
                function(sp) paste0('beta_vpd[', sp, ']'))
util$plot_disc_pushforward_quantiles(samples, names,
                                     xlab="Species",
                                     ylab="beta_vpd")

par(mfrow=c(1, 1), mar = c(4,4,1,1))
names <- sapply(1:data$N_species,
                function(sp) paste0('beta_gdd[', sp, ']'))
util$plot_disc_pushforward_quantiles(samples, names,
                                     xlab="Species",
                                     ylab="beta_gdd")



