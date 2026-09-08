rm(list = ls());gc()
library(terra)
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
fd <- 'data/climate/nex-dcp30-cmip6'

datasets <- readRDS(file.path(wd, 'output/model', 'datasets_10july2026_24species_365stands_19502022.rds'))
datasets <- unique(datasets[,c("grouped_stand", "east_lon", "north_lat")])
plots <- vect(datasets, geom=c("east_lon", "north_lat"))
plotsID <- data.frame(ID = 1:dim(plots)[1], grouped_stand = datasets$grouped_stand)

gcms <- c('IPSL-CM6A-LR')
variants <- c('r1i1p1f1')

tmean_df <- data.frame()
gdd_df <- data.frame()
pr_df <- data.frame()
vpd_df <- data.frame()

for(m in gcms){
  cat(paste0(m, '\n'))
  
  fdm <- file.path(wd, fd, m)
  v <- variants[which(gcms == m)]
  
  # historical
  s <- 'historical'
  years <- 1951:2014
  hist_fdm <- file.path(fdm, s, v)
  
  for(y in years){
    cat(paste0('   -', y, '\n'))
    
    #----
    # tmean, all month
    var <- 'tasmin'
    file <- paste0(paste(var, 'mon', m, s, v, 'gr', y, sep = '_'), '.nc')
    rmin <- rast(file.path(hist_fdm, var, file))
    
    var <- 'tasmax'
    file <- paste0(paste(var, 'mon', m, s, v, 'gr', y, sep = '_'), '.nc')
    rmax <- rast(file.path(hist_fdm, var, file))
    
    datmin <- terra::extract(rmin, plots)
    datmax <- terra::extract(rmax, plots)
    
    datmean <- data.frame(ID = datmin$ID, (datmin[,-1]+datmax[,-1])/2)
    names(datmean)[-1] <- paste0("tasmean_", 1:12)
    
    datdf <- data.frame(gcm = m, scenario = s, year = y,
                      ID = datmean[,1], tmean_all = rowMeans(datmean[,-1]))
    tmean_df <- rbind(tmean_df, datdf)
    #----
    
    
    #----
    # GDD, all month
    ndays <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    tlower <- 5
    tupper <- 35
    
    gdd <- apply(datmean[, -1], 1, function(tasmean) {
      taux <- pmax(pmin(tasmean, tupper), tlower) - tlower # apply lower and upper bounds
      gdd_aux <- sum(taux * ndays)
      gdd <- gdd_aux*0.982+136.05
      return(gdd)
    })
    
    datdf <- data.frame(gcm = m, scenario = s, year = y,
                        ID = datmean[,1], gdd_all = gdd)
    gdd_df <- rbind(gdd_df, datdf)
    #----
    
    
    #----
    # pr in NDJFMA
    var <- 'pr'
  
    file <- paste0(paste(var, 'mon', m, s, v, 'gr', y-1, sep = '_'), '.nc')
    rp <- rast(file.path(hist_fdm, var, file)) # previous yeqr
    
    file <- paste0(paste(var, 'mon', m, s, v, 'gr', y, sep = '_'), '.nc')
    r <- rast(file.path(hist_fdm, var, file))
    
    cr <- c(subset(rp, 11:12), subset(r, 1:4))
    datpre <- terra::extract(cr, plots)
    datdf <- data.frame(gcm = m, scenario = s, year = y,
                      ID = datpre[,1], pre_ndjfma = rowSums(datpre[,-1]))
    pr_df <- rbind(pr_df, datdf)
    #----
    
    #----
    # average max. VPD in MJJA
    # I do: saturation vapor pressure at the monthly max. temperature minus actual vapor pressure
    # this is close to what PRISM does for long-term data
    var <- 'vpr'
    file <- paste0(paste(var, 'mon', m, s, v, 'gr', y, sep = '_'), '.nc')
    r <- rast(file.path(hist_fdm, var, file))
    
    datvpr <- terra::extract(r, plots)
    
    es_tmax <- 6.11 * 10^((7.5 * datmax[, -1]) / (237.3 + datmax[, -1])) # apply Magnus-Tetens to monthly tmax
    vpd_max <- es_tmax - datvpr[, -1]
    vpd_max <- vpd_max[,5:8] # mjja
    
    datdf <- data.frame(gcm = m, scenario = s, year = y,
                        ID = datvpr[,1], vpd_mjja = rowMeans(vpd_max))
    vpd_df <- rbind(vpd_df, datdf)
    #----
    
  }
}

nex_df <- merge(merge(merge(merge(tmean_df, gdd_df), pr_df), vpd_df), plotsID)
saveRDS(nex_df, file.path(wd, 'output', 'climate', 'nex', 'climprojections_08sept2026.rds'))

# par(mfrow = c(1,1), mar = c(4,4,1,1))
# plot(x = NULL, y = NULL, xlim = c(1951,2014), ylim = c(0,52))
# for(i in unique(vpd_df$ID)){
#   lines(x = vpd_df[vpd_df$ID == i, 'year'], y = vpd_df[vpd_df$ID == i, 'vpd_mjja'])
# }
# 


