rm(list = ls());gc()

datdf <- readRDS(file = file.path(wd, 'output', 'climate', 'prism', 'climvar_20dec2025.rds'))

# GDD in AMJJAS
cat('GDD AMJJAS\n')
tlower <- 5
tupper <- 35
gdd_amjjas <- data.frame()
for(year in unique(datdf$year){
  
  df <- datdf[datdf$var == "tmean" & lubridate::year(datdf$date) %in% year,]
  df$gdd <- NA
  
  for(p in unique(df$ID)){
    
    months <- c('04', '05', '06', '07', '08', '09')
    ndays <- c(30, 31, 30, 31, 31, 30)
    
    tair <- df[df$ID == p & df$month %in% months, "value"]
    tair <- ifelse(tair < tlower, tlower, ifelse(tair > tupper, tupper, tair))-tlower # apply lower and upper bounds
    
    gdd_aux <- sum(tair*ndays)
    gdd <- gdd_aux*0.982+136.05
    
    gdd_amjjas <- rbind(gdd_amjjas, data.frame(ID = p, year = year, gdd_amjjas = gdd))
  }
}