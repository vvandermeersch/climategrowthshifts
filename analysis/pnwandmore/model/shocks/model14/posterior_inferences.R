rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(terra)
library(rnaturalearth)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)


data <-readRDS(file.path(wd, 'output/model', 'data_30jan2025_gymnosperms_standclimate_19502024_7species.rds'))
fit_conconly <- readRDS(file.path(wd, 'output/model/model14', 'fit_16feb2026_gymnosperms_standclimate_19502024_7species_model14.rds'))

# Parameter samples
base_samples <- fit_conconly$draws(
  variables =  c("mu_alpha","mu_gdd", "mu_vpd", "mu_pre", 
                 "tau_alpha","tau_gdd", "tau_vpd", "tau_pre",
                 "alpha","beta_gdd", "beta_vpd", "beta_pre",
                 
                 "mu_rho", "tau_rho", "rho_sp",
                 "mu_gamma", "tau_gamma", "gamma_sp",
                 
                 "rho_sh",
                 "mu_kappa", "tau_kappa", "kappa_sh",
                                              
                 "mu_tau_sck", "tau_tau_sck", "tau_sck",
                                              
                 "phi_sck0", "tau_phi_sck", "alpha_phi_sck",
                 "omega_conc_sck0", "tau_omega_conc_sck", "alpha_omega_conc_sck",  
                 
                 "sigma"))

names <- dimnames(base_samples)$variable
base_samples <- lapply(1:dim(base_samples)[3], function(k){t(matrix(base_samples[1:dim(base_samples)[1],1:dim(base_samples)[2],k], 
                                                                    nrow = dim(base_samples)[1], ncol = dim(base_samples)[2]))})
names(base_samples) <- names
util$check_all_expectand_diagnostics(base_samples)


base_samples2 <- fit_conconly$draws(
  variables =  c("mu_pre_phi","tau_pre_phi", "beta_pre_phi",
                 "mu_vpd_phi","tau_vpd_phi", "beta_vpd_phi"))

names <- dimnames(base_samples2)$variable
base_samples2 <- lapply(1:dim(base_samples2)[3], function(k){t(matrix(base_samples2[1:dim(base_samples2)[1],1:dim(base_samples2)[2],k], 
                                                                    nrow = dim(base_samples2)[1], ncol = dim(base_samples2)[2]))})
names(base_samples2) <- names
util$check_all_expectand_diagnostics(base_samples2)

param_samples <- fit_conconly$draws(
  variables =  c("phi_sck_baseline"))
names <- dimnames(param_samples)$variable
param_samples <- lapply(1:dim(param_samples)[3], function(k){t(matrix(param_samples[1:dim(param_samples)[1],1:dim(param_samples)[2],k], 
                                                                    nrow = dim(param_samples)[1], ncol = dim(param_samples)[2]))})
names(param_samples) <- names


util$plot_pairs_by_chain(base_samples2[["mu_pre_phi"]], "mu_pre_phi",
                         log(base_samples2[["tau_pre_phi"]]), "log(tau_pre_phi)")
util$plot_pairs_by_chain(base_samples2[["mu_vpd_phi"]], "mu_vpd_phi",
                         log(base_samples2[["tau_vpd_phi"]]), "log(tau_vpd_phi)")

util$plot_pairs_by_chain(base_samples[["mu_gdd[1]"]], "mu_gdd",
                         log(base_samples[["tau_gdd[1]"]]), "log(tau_gdd)")
util$plot_pairs_by_chain(base_samples[["mu_pre[1]"]], "mu_pre",
                         log(base_samples[["tau_pre[1]"]]), "log(tau_pre)")

util$plot_pairs_by_chain(base_samples[["mu_tau_sck[1]"]], "mu_tau_sck",
                         log(base_samples[["tau_tau_sck[1]"]]), "log(tau_tau_sck)")

util$plot_expectand_pushforward(base_samples2[["beta_pre_phi[2]"]], 20)

par(mfrow = c(1,1), mar = c(1.5,5,1,1), cex.main = 1, cex.lab = 1.2)
names <- paste0("phi_sck_baseline[",1:data$N_stands,"]")
util$plot_disc_pushforward_quantiles(param_samples, names, 
                                     xlab="Year", ylab=bquote(phi[shock ~ (stand)]^concordant),
                                     display_ylim=c(-0.05, 1.05), main = 'Only concordant shocks')

names <- paste0("beta_vpd_phi[",1:data$N_stands,"]")
util$plot_disc_pushforward_quantiles(base_samples2, names, 
                                     xlab="Year", ylab=bquote(phi[shock ~ (stand)]^concordant),
                                     display_ylim=c(-0.05, 0.1), main = 'Only concordant shocks')

names <- paste0("beta_pre[",1:data$N_species,"]")
util$plot_disc_pushforward_quantiles(base_samples, names, 
                                     xlab="Year", ylab=bquote(phi[shock ~ (stand)]^concordant),
                                     display_ylim=c(-0.05, 0.1), main = 'Only concordant shocks')

names <- paste0("tau_sck[",1:data$N_species,"]")
util$plot_disc_pushforward_quantiles(base_samples, names, 
                                     xlab="Year", ylab=bquote(phi[shock ~ (stand)]^concordant),
                                     display_ylim=c(0, 3), main = 'Only concordant shocks')

s <- 3
pre0 <- 5
pre <- seq(1,30, 1)
for(p in pre){
  newname <- paste0("phi_pre_at_",p,"")
  base_samples2[[newname]] <- boot::inv.logit(boot::logit(param_samples[[paste0("phi_sck_baseline[",s,"]")]]) +
                                                base_samples2[[paste0("beta_pre_phi[",s,"]")]]*(p-pre0))
}
par(mfrow = c(1,1), mar = c(4,4,1,1), cex.main = 1, cex.lab = 1.2)
util$plot_conn_pushforward_quantiles(base_samples2, paste0("phi_pre_at_",pre,""), pre)


s <- 3
vpd0 <- 23
vpd <- seq(1,40, 1)
for(v in vpd){
  newname <- paste0("phi_vpd_at_",v,"")
  base_samples2[[newname]] <- boot::inv.logit(boot::logit(param_samples[[paste0("phi_sck_baseline[",s,"]")]]) +
                                                base_samples2[[paste0("beta_pre_phi[",s,"]")]]*(v-vpd0))
}
par(mfrow = c(1,1), mar = c(4,4,1,1), cex.main = 1, cex.lab = 1.2)
util$plot_conn_pushforward_quantiles(base_samples2, paste0("phi_vpd_at_",vpd,""), vpd)

