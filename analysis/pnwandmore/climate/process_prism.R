rm(list = ls());gc()
library(terra)

wd <- "/home/vvandermeersch/projects/climategrowthshifts/analysis/pnw"

plots <- readRDS(file = file.path(wd,'data/itrdb', 'datasets_summary_all.rds'))
plots <- plots[!(plots$state %in% c("can", "mexi")), ]
plots <- vect(plots, geom=c("east_lon", "north_lat"))

plotsID <- data.frame(ID = 1:dim(plots)[1], dataset = plots$dataset)

vars <- c('tmean', 'ppt', 'vpdmax')
years <- 1895:2024
months <- c('01', '02', '03', '04', '05', '06',
            '07', '08', '09', '10', '11', '12')

resalt <- '30s'
country <- 'us'

dat_allvars <- data.frame()

for(var in vars){
  cat(paste0('Processing ', var, '...\n'))
  
  climdir <- file.path('/home/vvandermeersch/data/climate/prism/800m/us', var)
  
  dat_var <- data.frame()

  for(year in years){
    
    dat <- rast(lapply(months, function(m){
      filename <- paste0(paste('prism', var, country, resalt, paste0(year, m), sep = '_'), '.tif')
      rast(file.path(climdir, filename))
    }))
    
    e <- terra::extract(dat[[1]], plots, search_radius = 2000, xy = TRUE)
    plots <- vect(e, geom = c('x','y'))
    
    dat <- extract(dat, plots)
    dat <- as.data.frame(tidyr::pivot_longer(dat, cols = -ID, names_to = "var", values_to = "value"))
    dat$year <- stringr::str_match(dat$var, "(\\d{4})(\\d{2})$")[, 2]
    dat$month <- stringr::str_match(dat$var, "(\\d{4})(\\d{2})$")[, 3]
    dat$var <- var
    
    dat_var <- rbind(dat_var, dat)

  }
  
  dat_allvars <- rbind(dat_allvars, dat_var)
  
}

dat_allvars <- merge(dat_allvars, plotsID)
saveRDS(dat_allvars, file = file.path(wd, 'output', 'climate', 'prism', 'climvar_20dec2025.rds'))





