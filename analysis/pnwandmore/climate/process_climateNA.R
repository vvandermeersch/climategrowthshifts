
wd <- "/home/victor/projects/climategrowthshifts/analysis/pnwandmore"
library(ClimateNAr) # this package is weird
library(data.table) # need this for ClimateNAr

source(file.path(wd, 'climate', 'climateNAr_modified.R'))

sites <- readRDS(file.path(wd, 'input', 'itrdb', 'datasets_summary_all.rds'))
sites$id1 <- 1:nrow(sites)
sites$id2 <- 1:nrow(sites)
sites <- sites[,c('id1', 'id2', 'dataset', 'north_lat', 'east_lon', 'altitude')]
names(sites) <- c('id1', 'id2', 'dataset', 'lat', 'lon', 'elev')
sites[is.na(sites$elev), 'elev'] <- 500 # need to change this!

years <- 1901:2024

ClimateNAr(
  inputFile = sites[,c('id1', 'id2', 'lat', 'lon', 'elev')],
  periodList = paste0("Year_",years,".ann"),
  varList = c("DD5", "DD1040", "NFFD", "FFP", paste0('PPT', sprintf("%02d", 1:12))),
  outDir = file.path(wd, 'data', 'climate', 'climatena/')
)

clim <- data.frame()
for(y in years){
  clim_y <- read.csv(file.path(wd, 'data', 'climate', 'climatena', paste0('clim_Year_',y,'.csv')))
  sites_y <- sites[, c('id1', 'id2', 'dataset')]
  sites_y$year <- y
  clim_y <- merge(sites_y, clim_y)
  clim <- rbind(clim, clim_y)
}
saveRDS(clim, file.path(wd, 'output', 'climate', 'climatena', 'climpredictors_11august2025.rds'))
