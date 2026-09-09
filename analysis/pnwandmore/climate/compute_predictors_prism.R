rm(list = ls());gc()
wd <- "/home/vvandermeersch/projects/climategrowthshifts/analysis/pnw"

datdf <- readRDS(file = file.path(wd, 'output', 'climate', 'prism', 'climvar_20dec2025.rds'))
years <- 1896:2024

plots <- readRDS(file = file.path(wd,'data/itrdb', 'datasets_summary_all.rds'))
plots <- plots[!(plots$state %in% c("can", "mexi")), ]
plots <- vect(plots, geom=c("east_lon", "north_lat"))

plotsID <- data.frame(ID = 1:dim(plots)[1], dataset = plots$dataset)

# GDD full year
cat('GDD full year\n')
tlower <- 5
tupper <- 35
gdd_all <- data.frame()
for(year in years){
  
  months <- c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12')
  ndays <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  
  df <- datdf[datdf$var == "tmean" & datdf$year %in% year & datdf$month %in% months,]
  df$gdd <- NA
  
  for(p in unique(df$ID)){
    
    
    
    tair <- df[df$ID == p, "value"]
    tair <- ifelse(tair < tlower, tlower, ifelse(tair > tupper, tupper, tair))-tlower # apply lower and upper bounds
    
    gdd_aux <- sum(tair*ndays)
    gdd <- gdd_aux*0.982+136.05
    
    gdd_all <- rbind(gdd_all, data.frame(ID = p, year = year, gdd_all = gdd))
  }
}

# Precipitation in NDJFMA
cat('Precipitation NDJFMA\n')
ppt_ndjfma <- data.frame()
for(year in years){
  
  df <- rbind(
    datdf[datdf$var == "ppt" & datdf$year %in% c(year-1) & datdf$month %in% c('11', '12'),],
    datdf[datdf$var == "ppt" & datdf$year %in% c(year) & datdf$month %in% c('01', '02', '03', '04'),])

  for(p in unique(df$ID)){

    pre <- df[df$ID == p, "value"]
    pre <- sum(pre)

    ppt_ndjfma <- rbind(ppt_ndjfma, data.frame(ID = p, year = year, ppt_ndjfma = pre))
  }
}

# Average max. VPD in MJJA
cat('VPD MJJA\n')
vpd_mjja <- data.frame()
for(year in years){
  
  df <- datdf[datdf$var == "vpdmax" & datdf$year %in% c(year) & datdf$month %in% c('05', '06', '07', '08'),]
  
  for(p in unique(df$ID)){
    
    vpd <- df[df$ID == p, "value"]
    vpd <- mean(vpd)
    
    vpd_mjja <- rbind(vpd_mjja, data.frame(ID = p, year = year, vpd_mjja = vpd))
  }
}

climate_predictors <- merge(gdd_all, ppt_ndjfma)
climate_predictors <- merge(climate_predictors, vpd_mjja)
climate_predictors <- merge(climate_predictors, plotsID)
saveRDS(climate_predictors, file = file.path(wd, 'output', 'climate', 'prism', 'climpredictors_20dec2025.rds'))