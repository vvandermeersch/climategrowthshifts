rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(terra)
library(rnaturalearth)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)

data <- readRDS(file.path(wd, 'output/model', 'data_30jan2025_gymnosperms_standclimate_19502024_7species.rds'))
subdraws <- readRDS(file.path(wd, 'output/model/model14', 'subdraws_16feb2026_gymnosperms_standclimate_19502024_7species_model14.rds'))

names <- dimnames(subdraws)$variable
subdraws <- lapply(1:dim(subdraws)[3], function(k){t(matrix(subdraws[1:dim(subdraws)[1],1:dim(subdraws)[2],k], 
                                                            nrow = dim(subdraws)[1], ncol = dim(subdraws)[2]))})
names(subdraws) <- names

util$plot_expectand_pushforward(subdraws)

par(mar= c(4.5,4.5,1,1))
names <- paste0("phi_sck_baseline[",1:83,"]")
util$plot_disc_pushforward_quantiles(subdraws, names, 
                                     xlab="Year", ylab=bquote(phi[shock ~ (stand)]^concordant ~ 'at VPD = 15'),
                                     display_ylim=c(-0.05, 1.05))

names <- paste0("omega_conc_sck[",1:83,"]")
util$plot_disc_pushforward_quantiles(subdraws, names, 
                                     xlab="Year", ylab = bquote(omega[shock ~ (tree)]^concordant),
                                     display_ylim=c(-0.05, 1.05))


names <- paste0("omega_nonconc_sck[",1:83,"]")
util$plot_disc_pushforward_quantiles(subdraws, names, 
                                     xlab="Year", ylab = bquote(omega[shock ~ (tree)]^{non-concordant}),
                                     display_ylim=c(-0.05, 1.05))


s <- which(data$N_stand_years == 71)[7]

s <- 66
idxs <- seq(1+(s-1)*data$N_all_years, s*data$N_all_years, 1)
data$vpd_obs[idxs]
names <- paste0("phi_sck[", s, ",", 1:71,"]")
for(i in 1:71){
  namehere <- names[i]
  vpdhere <- data$vpd_obs[idxs[i]]
  
  subdraws[[namehere]] <- boot::inv.logit(
    boot::logit(subdraws[[paste0("phi_sck_baseline[",s,"]")]]) +
      subdraws[["beta_vpd_phi"]] * (vpdhere-15)
  )
}
order <- order(data$vpd_obs[idxs])
names <- paste0("phi_sck[", s, ",", order,"]")
util$plot_conn_pushforward_quantiles(subdraws, names, data$vpd_obs[idxs][order],
                                     xlab = 'VPD', ylab = bquote(phi[shock ~ (stand ~ '35')]^concordant))



