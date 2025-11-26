rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('functions/custom_functions.R', local=util)
setwd(wd)

data <- readRDS(file.path(wd, 'output/model', 'data_07oct2025.rds'))
datasets <- readRDS(file.path(wd, 'output/model', 'datasets_07oct2025.rds'))
samples <- readRDS(file.path(wd, "output/model/samples_07oct2025_partialpooling_2clades.rds"))

species <- 'PIPO'
stands_filtered <- unique(data$stand_idxs[data$species_idxs == which(data$uniq_species_ids == species)])
stands_name <- data$uniq_stand_ids[stands_filtered]
stands_coords <- unique(datasets[datasets$grouped_stand %in% stands_name, c('grouped_stand', 'north_lat', 'east_lon')])
stands_coords$north_lat <- round(stands_coords$north_lat,2)
stands_coords$east_lon <- round(stands_coords$east_lon,2)
stands_coords <- unique(stands_coords)
stands_filtered <- stands_filtered[stands_coords[stands_coords$grouped_stand == data$uniq_stand_ids[stands_filtered],]$north_lat < 38 & 
                                     stands_coords[stands_coords$grouped_stand == data$uniq_stand_ids[stands_filtered],]$east_lon > -115]
subsets <- datasets[datasets$grouped_stand %in% c(data$uniq_stand_ids[stands_filtered]), 'dataset']
saveRDS(subsets, file.path(wd, 'output', 'subsets_datasets_Ben.rds'))

par(mfrow = c(1,1), mar = c(3,3,1,2))

plot(1, type="n", xlim = c(1980, 2020), ylim = c(-7,7), main = '', xlab = '', y ='')
for(stand in stands_filtered){
  
  trees <- which(data$stand_idxs == stand)
  observed_years <- unlist(lapply(trees, function(t){
    return(data$years[data$tree_start_idxs[t]:data$tree_end_idxs[t]])
  }))
  observed_years <- unique(observed_years)[order(unique(observed_years))]
  
  names <- sapply(which(data$all_years %in% observed_years),
                  function(n) paste0('f_tilde_sh[', stand,',', n, ']'))
  names <- util$check_expectand_names(names, samples)
  
  # Extract function values
  fs <- t(sapply(names, function(name)
    c(samples[[name]], recursive=TRUE)))
  N <- dim(fs)[1]
  I <- dim(fs)[2]
  J <- 1
  plot_idx <- sapply(1:J, function(j) (I %/% J) * (j - 1) + 1)
  
  
  plot_xs <- observed_years
  for (j in 1:J) {
    f <- fs[, plot_idx[j]]
    lines(plot_xs, f, lwd=0.5, col = "darkblue")
  }
}
abline(v=2002, lwd=1, lty=2, col="grey50")

plot(westam, xlim = c(-130,-100), ylim = c(30, 50))
points(vect(stands_coords[stands_coords$north_lat < 38 & stands_coords$east_lon > -115, ], 
            geom = c('east_lon', 'north_lat')), col = 'red')

par(mfrow = c(2,1), mar = c(2,5,1,2), cex.lab = 1.5)
plot(1, type="n", xlim = c(1980, 2020), ylim = c(10,30), main = '', xlab = '', ylab ='SM (%)')
for(stand in stands_filtered){
  
  trees <- which(data$stand_idxs == stand)
  observed_years <- unlist(lapply(trees, function(t){
    return(data$years[data$tree_start_idxs[t]:data$tree_end_idxs[t]])
  }))
  observed_years <- unique(observed_years)[order(unique(observed_years))]
  
  observed_sm <- unlist(lapply(trees, function(t){
    return(data$sm_obs[data$tree_start_idxs[t]:data$tree_end_idxs[t]])
  }))
  observed_sm <- observed_sm[as.character(observed_years)]
  
  
  plot_xs <- observed_years
  lines(plot_xs, observed_sm, lwd=0.5, col = "darkblue")
  
}
abline(v=2002, lwd=1, lty=2, col="grey50")

plot(1, type="n", xlim = c(1980, 2020), ylim = c(5,23), main = '', xlab = '', ylab ='VPD (hPa)')
for(stand in stands_filtered){
  
  trees <- which(data$stand_idxs == stand)
  observed_years <- unlist(lapply(trees, function(t){
    return(data$years[data$tree_start_idxs[t]:data$tree_end_idxs[t]])
  }))
  observed_years <- unique(observed_years)[order(unique(observed_years))]
  
  observed_vpd <- unlist(lapply(trees, function(t){
    return(data$vpd_obs[data$tree_start_idxs[t]:data$tree_end_idxs[t]])
  }))
  observed_vpd <- observed_vpd[as.character(observed_years)]
  
  
  plot_xs <- observed_years
  lines(plot_xs, observed_vpd, lwd=0.5, col = "darkblue")
  
}
abline(v=2002, lwd=1, lty=2, col="grey50")




ggplot(data = climpredictors_15oct2025_subset) +
  geom_line(aes(x = year, y = pre_amjja, group = dataset))

ggplot(data = climpredictors_15oct2025_subset) +
  geom_line(aes(x = year, y = nbd_nopre_amjja, group = dataset))

ggplot(data = climpredictors_15oct2025_subset) +
  geom_line(aes(x = year, y = vpd_mjja, group = dataset))
