rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(ggplot2)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)

# Model data
data <- readRDS(file.path(wd, 'output/model', 'data_03july2025.rds'))

# Posterior quantification - diagnostics
fit <- readRDS(file.path(wd, 'output/model', 'fit_03july2025_partialpooling_2clades_centered.rds')) 

samples <- util$extract_expectand_vals(fit)

stands <- 1:data$N_stands
fshs_df <- data.frame()

for(stand in stands){
  
  names <- sapply(1:length(data$all_years),
                  function(n) paste0('f_sh[', stand,',', n, ']'))
  names <- util$check_expectand_names(names, samples)
  
  # Extract function values
  
  fs <- t(sapply(names, function(name)
    c(samples[[name]], recursive=TRUE)))
  N <- dim(fs)[1]
  I <- dim(fs)[2]
  J <- 20
  plot_idx <- sapply(1:J, function(j) (I %/% J) * (j - 1) + 1)
  
  for (j in 1:J) {
    f <- fs[, plot_idx[j]]
    fshs_df <- rbind(fshs_df, data.frame(stand = data$uniq_stand_ids[stand], real = j, year = data$all_years, fsh = f))
  }
}

# Load datasets
datasets <- readRDS(file.path(wd, 'input', 'itrdb', 'datasets_summary_all.rds'))
datasets <- datasets[datasets$dataset != 'wa149',] # temporary 
datasets <- datasets[datasets$dataset != 'ca673',] # temporary
datasets <- datasets[datasets$dataset != 'az563',] # temporary
datasets <- datasets[datasets$dataset != 'co689',] # temporary 
datasets <- datasets[datasets$dataset != 'co690',] # temporary
datasets <- datasets[datasets$dataset != 'co691',] # temporary
datasets <- datasets[datasets$north_lat >=25.065 & datasets$north_lat <=52.925 & datasets$east_lon <= -89.025 &  datasets$east_lon >= -124.925,]
datasets <- datasets[datasets$last_year >= 1999,]


