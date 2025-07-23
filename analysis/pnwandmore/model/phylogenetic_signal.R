rm(list = ls())
library(ggplot2)
library(phytools)
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)

data <- readRDS(file.path(wd, 'output/model', 'data_11july2025.rds'))
slopes_partial <- readRDS(file.path(wd, 'output/model/shrinkage', 'slopes_11july2025_partialpooling.rds'))
slopes_nopooling <- readRDS(file.path(wd, 'output/model/shrinkage', 'slopes_11july2025_nopooling.rds'))

phylotree <- readRDS(file.path(wd, 'output/model', 'phylotree_11july2025.rds'))

phylotree$tip.label
data$uniq_species_ids[match(phylotree$tip.label, data$uniq_species_ids)]
species_order <- match(phylotree$tip.label, data$uniq_species_ids)

slopes_df <- 
  data.frame(
    species_num = 1:data$N_species,
    species_code = data$uniq_species_ids,
    clade = ifelse(data$clade_idxs == 1, 'gymnosperms', 'angiosperms'),
    q50_beta_gdd = sapply(1:data$N_species, function(s) util$ensemble_mcmc_quantile_est(slopes_nopooling[[paste0('beta_gdd','[',s,']')]], c(0.5))),
    q50_beta_sm = sapply(1:data$N_species, function(s) util$ensemble_mcmc_quantile_est(slopes_nopooling[[paste0('beta_sm','[',s,']')]], c(0.5))),
    q50_beta_vpd = sapply(1:data$N_species, function(s) util$ensemble_mcmc_quantile_est(slopes_nopooling[[paste0('beta_vpd','[',s,']')]], c(0.5)))
    )

slopes_df <- slopes_df[species_order,]

slopes <- slopes_df$q50_beta_gdd
names(slopes) <- slopes_df$species_code
phylosig(phylotree, slopes, method = 'lambda', test= TRUE)

slopes <- slopes_df$q50_beta_sm
names(slopes) <- slopes_df$species_code
phylosig(phylotree, slopes, method = 'lambda', test= TRUE)

slopes <- slopes_df$q50_beta_vpd
names(slopes) <- slopes_df$species_code
phylosig(phylotree, slopes, method = 'lambda', test= TRUE)