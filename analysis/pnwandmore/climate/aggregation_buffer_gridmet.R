rm(list = ls());gc()
library(terra)
library(ggplot2)

wd <- "/home/victor/projects/climategrowthshifts/analysis/pnwandmore"
climdir <- "data/climate/gridmet"

datasets <- readRDS(file.path(wd, 'output/model', 'datasets_10oct2025.rds'))
years <- 1980:2020

d <- "S173"
dcoords <- vect(unique(datasets[datasets$grouped_stand %in% d, c('grouped_stand', 'north_lat','east_lon')]), geom = c('east_lon', 'north_lat'))
crs(dcoords) <- "EPSG:4326"
buffer5km <- buffer(dcoords, width = 5000)
buffer10km <- buffer(dcoords, width = 10000)
buffer20km <- buffer(dcoords, width = 20000)

pre_df <- data.frame()
for(y in years){
  cat(paste0(y, '\n'))
  
  pre <- c(rast(file.path(wd, climdir, paste0('pr_', y-1, '.nc'))), rast(file.path(wd, climdir, paste0('pr_', y, '.nc'))))
  time(pre) <- seq.Date(from = as.Date(paste0(y-1, '-01-01')), to = as.Date(paste0(y, '-12-31')), by = 1)
  
  s <- as.Date(paste0(y-1, '-11-01'))
  e <- as.Date(paste0(y, '-04-30'))
  
  onecell <- sum(unlist(extractRange(pre, dcoords, first = which(time(pre) == s), last = which(time(pre) == e), cells=FALSE, xy=FALSE, ID=FALSE)))
  
  cells5km <- sum(unlist(extractRange(pre, buffer5km, first = which(time(pre) == s), last = which(time(pre) == e), cells=FALSE, xy=FALSE, ID=FALSE, geom_fun = mean)))
  cells10km <- sum(unlist(extractRange(pre, buffer10km, first = which(time(pre) == s), last = which(time(pre) == e), cells=FALSE, xy=FALSE, ID=FALSE, geom_fun = mean)))
  cells20km <- sum(unlist(extractRange(pre, buffer20km, first = which(time(pre) == s), last = which(time(pre) == e), cells=FALSE, xy=FALSE, ID=FALSE, geom_fun = mean)))
  
  
  pre_df <- rbind(
    pre_df,
    data.frame(year = y, pre_onecell = onecell, pre_5km = cells5km, pre_10km = cells10km, pre_20km = cells20km)
  )
}

ggplot(data = pre_df) +
  geom_line(aes(x = year, y = pre_onecell)) +
  geom_line(aes(x = year, y = pre_5km), col = "#A3485A", size = 0.5) +
  geom_line(aes(x = year, y = pre_10km), col = "#842A3B", size = 0.85) + 
  geom_line(aes(x = year, y = pre_20km), col = "#662222", size = 1.2) +
  theme_classic()

clim_wldas <- readRDS("/home/victor/projects/climategrowthshifts/analysis/pnwandmore/output/climate/climpredictors_10oct2025.rds")


ggplot() +
  geom_line(data = pre_df, aes(x = year, y = pre_onecell)) +
  # geom_line(data = pre_df, aes(x = year, y = pre_5km), col = "#A3485A", size = 0.5) +
  geom_line(data = clim_wldas[clim_wldas$dataset %in% 'az565',], aes(x = year, y = pre_ndjfma), color = 'darkblue') +
  
  theme_classic()


pdsi <- rast(file.path(wd, climdir, paste0('pdsi.nc')), subds="pdsi")
time(pdsi) <- as.Date( depth(pdsi), origin = as.Date('1900-01-01'))
pdsi_df <- data.frame()
for(y in years){
  cat(paste0(y, '\n'))
  
  s <- as.Date(paste0(y, '-05-01'))
  e <- as.Date(paste0(y, '-07-31'))
  min(which(time(pdsi) %in% seq.Date(s,e,by=1)))
  
  
  onecell <- mean(unlist(extractRange(pdsi, dcoords, first =  min(which(time(pdsi) %in% seq.Date(s,e,by=1))), 
                                     last = max(which(time(pdsi) %in% seq.Date(s,e,by=1))), cells=FALSE, xy=FALSE, ID=FALSE)))
  
  pdsi_df <- rbind(
    pdsi_df,
    data.frame(year = y, pdsi_onecell = onecell)
  )
}

ggplot() +
  geom_line(data = pdsi_df, aes(x = year, y = pdsi_onecell)) +
  # geom_line(data = pre_df, aes(x = year, y = pre_5km), col = "#A3485A", size = 0.5) +
  
  theme_classic()

