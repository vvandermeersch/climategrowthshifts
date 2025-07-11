rm(list = ls())

library(terra)

wd <- "/home/vvandermeersch/projects/climategrowthshifts/analysis/pnw"
climdir <- "/home/vvandermeersch/data/climate"

plots <- readRDS(file = file.path(wd, 'data', 'itrdb', paste0('itrdb_info.rds')))
plots <- vect(plots, geom=c("east_lon", "north_lat"))

# plotsID <- data.frame(ID = 1:16, plotname = plots$Plot, alt = plots$Elevation)

years <- seq(1980,1981,1)
vars <- list("tair" = "Tair_f_tavg", 
             "soilmoist_010" = "SoilMoi00_10cm_tavg",
             "soilmoist_1040" = "SoilMoi10_40cm_tavg",
             "soilmoist_40100" = "SoilMoi40_100cm_tavg",
             "soilmoist_100200" = "SoilMoi100_200cm_tavg")
datdf <- data.frame()

for(year in years){
  
  cat(paste0(year,"\n"))
  # files <- list.files(path = , pattern = as.character(year), full.names = TRUE)
  
  dates <- format(seq(as.Date(paste0(year, "-01-01")), as.Date(paste0(year, "-12-31")),1), "%Y%m%d")
  files <- file.path(climdir,"wldas", paste0("WLDAS_NOAHMP001_DA1_",dates,".D10.nc.SUB.nc4"))
  
  
  # if(length(files) < 365){stop("Problem with no. of files")}
  if(any(!file.exists(files))){stop("Problem with no. of files")}
  dat <- rast(files)
  
  daty <- data.frame()
  
  for(v in vars){
    newv <- names(vars)[which(vars==v)]
    cat(paste0("   ", newv,"\n"))
    
    datv <- extract(dat[v], plots)
    datv <- as.data.frame(tidyr::pivot_longer(datv, cols = -ID, names_to = "var", values_to = "value"))
    
    if(newv == "tair"){
      datv$value <- datv$value  -273.15 # K to C
    }
    
    datv$var <- newv
    daty <- rbind(daty, datv)
  }
  
  daty$date <- rep(seq(as.Date(paste0(year, "-01-01")), as.Date(paste0(year, "-12-31")),1), length(unique(daty$ID))*length(vars))
  daty$doy <- lubridate::yday(daty$date)
  datdf <- rbind(datdf, daty)
  
}

saveRDS(datdf, file = file.path(wd, 'output', 'climate', 'climvar_pnw.rds'))





