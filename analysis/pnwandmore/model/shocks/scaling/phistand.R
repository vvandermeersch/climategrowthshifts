rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(terra)
library(rnaturalearth)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)


data_19001940 <-readRDS(file.path(wd, 'output/model', 'data_15jan2025_PIPO_standclimate_19001940.rds'))
fit_19001940 <- readRDS(file.path(wd, 'output/model/scaling', 'fit_15jan2025_PIPO_standclimate_19001940_customgrainsize.rds'))
samples_19001940 <- fit_19001940$draws(variables = c("beta_gdd", "beta_vpd", "beta_pre", "phi_sck", "rho_sp", "rho_sh"))

data_19401980 <- readRDS(file.path(wd, 'output/model', 'data_15jan2025_PIPO_standclimate_19401980.rds'))
fit_19401980 <- readRDS(file.path(wd, 'output/model/scaling', 'fit_15jan2025_PIPO_standclimate_19401980_customgrainsize.rds'))
samples_19401980 <- fit_19401980$draws(variables = c("beta_gdd", "beta_vpd", "beta_pre", "phi_sck", "rho_sp", "rho_sh"))

data_19802024 <- readRDS(file.path(wd, 'output/model', 'data_15jan2025_PIPO_standclimate_19802024.rds'))
fit_19802024 <- readRDS(file.path(wd, 'output/model/scaling', 'fit_15jan2025_PIPO_standclimate_19802024_customgrainsize.rds'))
samples_19802024 <- fit_19802024$draws(variables = c("beta_gdd", "beta_vpd", "beta_pre", "phi_sck", "rho_sp", "rho_sh"))

stands_19001940 <- data.frame()
for(s in 1:data_19001940$N_stands){
  trees <- data_19001940$uniq_tree_ids[data_19001940$stand_idxs == s]
  original_stand_name <- unique(stringr::str_split_i(trees, '_', 1))
  if(length(original_stand_name) > 1){
    # temporary
    original_stand_name <- original_stand_name[1]
  }
  stands_19001940 <- rbind(
    stands_19001940,
    data.frame(origname = original_stand_name, id_19001940 = s)
  )
}

stands_19401980 <- data.frame()
for(s in 1:data_19401980$N_stands){
  trees <- data_19401980$uniq_tree_ids[data_19401980$stand_idxs == s]
  original_stand_name <- unique(stringr::str_split_i(trees, '_', 1))
  if(length(original_stand_name) > 1){
    # temporary
    original_stand_name <- original_stand_name[1]
  }
  stands_19401980 <- rbind(
    stands_19401980,
    data.frame(origname = original_stand_name, id_19401980 = s)
  )
}

stands_19802024 <- data.frame()
for(s in 1:data_19802024$N_stands){
  trees <- data_19802024$uniq_tree_ids[data_19802024$stand_idxs == s]
  original_stand_name <- unique(stringr::str_split_i(trees, '_', 1))
  if(length(original_stand_name) > 1){
    # temporary
    original_stand_name <- original_stand_name[1]
  }
  stands_19802024 <- rbind(
    stands_19802024,
    data.frame(origname = original_stand_name, id_19802024 = s)
  )
}

stands_corres <- merge(merge(stands_19001940, stands_19401980), stands_19802024)

par(mfrow= c(6,2), mar = c(4,1,1,1), cex.lab = 1.4, cex.axis = 0.9)
init <- 1
for(i in sample(1:nrow(stands_corres), 12)){
  util$plot_expectand_pushforward(t(matrix(samples_19001940[1:1000,1:4,paste0('phi_sck[',stands_corres[i, 'id_19001940'],']')], nrow = 1000, ncol = 4)), 50,
                                  flim = c(0,0.5), ylim = c(0,20),
                                  display_name= bquote(phi[shock] ~ (.(i))),
                                  col = '#27278f')
  util$plot_expectand_pushforward(t(matrix(samples_19401980[1:1000,1:4,paste0('phi_sck[',stands_corres[i, 'id_19401980'],']')], nrow = 1000, ncol = 4)), 50,
                                  flim = c(0,0.5), ylim = c(0,20),
                                  display_name=expression(phi[shock]),
                                  col = '#278f27', add = TRUE)
  util$plot_expectand_pushforward(t(matrix(samples_19802024[1:1000,1:4,paste0('phi_sck[',stands_corres[i, 'id_19802024'],']')], nrow = 1000, ncol = 4)), 50,
                                  flim = c(0,0.5), ylim = c(0,20),
                                  display_name=expression(phi[shock]),
                                  col = '#8f2727', add = TRUE)
  
  if(init == 1){
    legend(x = 0.35, y = 15,
           lwd = 2, 
           col = c( '#27278f', '#278f27', '#8f2727'),
           legend = c('1900-1940', '1940-1980', '1980-2024'),
           bty = 'n')
    init <- 0
  }
}


for(i in 1:nrow(stands_corres)){
  median00 <- util$ensemble_mcmc_quantile_est(t(matrix(samples_19001940[1:1000,1:4,paste0('phi_sck[',stands_corres[i, 'id_19001940'],']')], nrow = 1000, ncol = 4)), c(0.5))
  median40 <- util$ensemble_mcmc_quantile_est(t(matrix(samples_19401980[1:1000,1:4,paste0('phi_sck[',stands_corres[i, 'id_19401980'],']')], nrow = 1000, ncol = 4)), c(0.5))
  median80 <- util$ensemble_mcmc_quantile_est((matrix(samples_19802024[1:1000,1:4,paste0('phi_sck[',stands_corres[i, 'id_19802024'],']')], nrow = 1000, ncol = 4)), c(0.5))
  
  
  if(median00 > 0.4 | median40 > 0.4 | median80> 0.4){
    util$plot_expectand_pushforward(t(matrix(samples_19001940[1:1000,1:4,paste0('phi_sck[',stands_corres[i, 'id_19001940'],']')], nrow = 1000, ncol = 4)), 50,
                                    flim = c(0,1), ylim = c(0,20),
                                    display_name= bquote(phi[shock] ~ (.(i))),
                                    col = '#27278f')
    util$plot_expectand_pushforward(t(matrix(samples_19401980[1:1000,1:4,paste0('phi_sck[',stands_corres[i, 'id_19401980'],']')], nrow = 1000, ncol = 4)), 50,
                                    flim = c(0,1), ylim = c(0,20),
                                    display_name=expression(phi[shock]),
                                    col = '#278f27', add = TRUE)
    util$plot_expectand_pushforward(t(matrix(samples_19802024[1:1000,1:4,paste0('phi_sck[',stands_corres[i, 'id_19802024'],']')], nrow = 1000, ncol = 4)), 50,
                                    flim = c(0,1), ylim = c(0,20),
                                    display_name=expression(phi[shock]),
                                    col = '#8f2727', add = TRUE)
  }
}


