rm(list = ls());gc()

library(terra)

wd <- "/home/vvandermeersch/projects/climategrowthshifts/analysis/pnw"
climdir <- "/home/vvandermeersch/data/climate"
outputdir <- 'output/climate/rasters'

for(year in 2000){
  cat(paste0(year, '\n'))
  dates <- format(seq(as.Date(paste0(year, "-01-01")), as.Date(paste0(year, "-12-31")),1), "%Y%m%d")
  
  files <- file.path(climdir,"wldas", paste0("WLDAS_NOAHMP001_DA1_",dates,".D10.nc.SUB.nc4")) 
  dat <- rast(files)
  
  # GDD over AMJJAS
  var <- "Tair_f_tavg"
  tmean <- dat[var]
  tlower <- 5+273.15
  tupper <- 35+273.15
  s <- as.Date(paste0(year, '-04-01'))
  e <- as.Date(paste0(year, '-09-30'))
  tmean_amjjas <- tmean[[time(tmean) %in% seq(s, e, 1)]]
  cat('  - applying lower threshold')
  tmean_amjjas <- ifel(tmean_amjjas < tlower, tlower, tmean_amjjas) 
  cat('  - applying upper threshold')
  tmean_amjjas <- ifel(tmean_amjjas > tupper, tupper, tmean_amjjas)
  cat('  - cumulative sum')
  gdd <- cumsum(tmean_amjjas-tlower)
  gdd <- gdd[[nlyr(gdd)]]
  writeRaster(gdd, file.path(wd, outputdir, paste0('gdd_amjjas_', year, '.nc')), overwrite = TRUE)
  
  # Average soil moisture over MJJ
  s <- as.Date(paste0(year, '-05-01'))
  e <- as.Date(paste0(year, '-07-31'))
  var <- "SoilMoi00_10cm_tavg"
  sm10 <- dat[var]
  sm10 <- mean(sm10[[time(sm10) %in% seq(s, e, 1)]])
  var <- "SoilMoi10_40cm_tavg"
  sm40 <- dat[var]
  sm40 <- mean(sm40[[time(sm40) %in% seq(s, e, 1)]])
  var <- "SoilMoi40_100cm_tavg"
  sm100 <- dat[var]
  sm100 <- mean(sm100[[time(sm100) %in% seq(s, e, 1)]])
  var <- "SoilMoi100_200cm_tavg"
  sm200 <- dat[var]
  sm200 <- mean(sm200[[time(sm200) %in% seq(s, e, 1)]])
  sm0_200 <- terra::weighted.mean(c(sm10, sm40, sm100, sm200), c(10,30,60,100))
  writeRaster(sm0_200, file.path(wd, outputdir, paste0('sm2m_mjj_', year, '.nc')), overwrite = TRUE)
  
}



