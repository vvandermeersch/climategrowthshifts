rm(list = ls());gc()
library(terra)
library(ggplot2)

wd <- "/home/victor/projects/climategrowthshifts/analysis/pnwandmore"
climdir <- "output/climate/rasters"

datasets <- readRDS(file.path(wd, 'output/model', 'datasets_10oct2025.rds'))
years <- 1980:2023

d <- "S210"
dcoords <- vect(unique(datasets[datasets$grouped_stand %in% d, c('grouped_stand', 'north_lat','east_lon')]), geom = c('east_lon', 'north_lat'))
crs(dcoords) <- "EPSG:4326"
buffer5km <- buffer(dcoords, width = 5000)
buffer10km <- buffer(dcoords, width = 10000)
buffer20km <- buffer(dcoords, width = 20000)

plot(sm)
points(dcoords, col = 'white', cex = 2)
points(dcoords, col = 'black', cex = 1)

sm_df <- data.frame()
for(y in years){
  
  sm <- rast(file.path(wd, climdir, paste0('sm2m_mjj_', y, '.nc')))
  
  onecell <- extract(sm, dcoords, cells=FALSE, xy=TRUE, ID=TRUE)
  cells5km <- extract(sm, buffer5km, cells=FALSE, xy=TRUE, ID=TRUE, fun = mean)
  cells10km <- extract(sm, buffer10km, cells=FALSE, xy=TRUE, ID=TRUE, fun = mean)
  cells20km <- extract(sm, buffer20km, cells=FALSE, xy=TRUE, ID=TRUE, fun = mean)
  
  
  sm_df <- rbind(
    sm_df,
    data.frame(year = y, sm_onecell = onecell$Band1, sm_5km = cells5km$Band1, sm_10km = cells10km$Band1, sm_20km = cells20km$Band1)
  )
}

ggplot(data = sm_df) +
  geom_line(aes(x = year, y = sm_onecell)) +
  geom_line(aes(x = year, y = sm_5km), col = "#A3485A", size = 0.5) +
  geom_line(aes(x = year, y = sm_10km), col = "#842A3B", size = 0.85) + 
  geom_line(aes(x = year, y = sm_20km), col = "#662222", size = 1.2) +
  theme_classic()

