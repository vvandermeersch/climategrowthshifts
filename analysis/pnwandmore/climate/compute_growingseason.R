
library(ggplot2)
library(rnaturalearth)
library(terra)
library(tidyterra)
library(lubridate)

wd <- "/home/victor/projects/climategrowthshifts/analysis/pnwandmore"
clim <- readRDS(file.path(wd, 'output/climate', 'climvar_04july2025.rds'))
datasets <- readRDS(file.path(wd, 'output/model', 'datasets_11july2025.rds'))
westam <- vect(ne_states(country = c("Canada", "Mexico", "United States of America"), returnclass = "sf"))

ggplot() +
  geom_spatvector(data = westam, color = 'white', fill = 'grey90') +
  geom_point(data = datasets[datasets$altitude >= 2500,], 
             aes(x = east_lon, y = north_lat), size = 3, color = "white") +
  geom_point(data = datasets[datasets$altitude >= 2500,], 
             aes(x = east_lon, y = north_lat), size = 1.5, color = "black") +
  coord_sf(xlim = c(-124.925, -105.025), ylim = c(25.065,52.925), expand = TRUE) +
  theme_void()

highelev_ids <- na.omit(datasets[datasets$altitude >= 2500, 'dataset'])

init <- "1980-01-01"
final <- "2023-12-31"

clim <- clim[clim$dataset %in% highelev_ids,]

# Find all the potential SOS and EOS
seos_df <- data.frame()
bg <- Sys.time()
for(s in highelev_ids){
  
  clim_subset <- clim[clim$dataset %in% s,]
  clim_subset <- clim_subset[order(clim_subset$date),]
  
  n <- 5
  # SOS
  # - thermal definition: 5 last days above 5degC
  cx <- c(0,cumsum(clim_subset[clim_subset$var %in% 'tair', 'value'] >= 5))
  thermal <- c(rep(0, n-1), (cx[(n+1):length(cx)] - cx[1:(length(cx) - n)]))
  
  # - moisture definition: 5 last days above 15%
  cx <- c(0,cumsum(clim_subset[clim_subset$var %in% 'soilmoist_1040', 'value'] >= 0.15))
  moisture <- c(rep(0, n-1), (cx[(n+1):length(cx)] - cx[1:(length(cx) - n)]))
  
  start = (thermal == 5)  & (moisture == 5)
  
  # EOS
  # - thermal critical threshold: 5 last days below 5degC
  cx <- c(0,cumsum(clim_subset[clim_subset$var %in% 'tair', 'value'] < 5))
  thermal <- c(rep(0, n-1), (cx[(n+1):length(cx)] - cx[1:(length(cx) - n)]))
  
  # - moisture critical threshold: 5 last days below 15%
  cx <- c(0,cumsum(clim_subset[clim_subset$var %in% 'soilmoist_1040', 'value'] < 0.15))
  moisture <- c(rep(0, n-1), (cx[(n+1):length(cx)] - cx[1:(length(cx) - n)]))
  
  end = thermal == 5 | moisture == 5
  
  seos_df <- rbind(seos_df, data.frame(s = s, date = as.character(seq.Date(as.Date(init), as.Date(final), "day")), start = start, end = end))
}
ed <- Sys.time()

# Find the growing intervals
grint_df <- data.frame()
for(s in highelev_ids){
  print(s)
  growing_intervals <- list()
  seos <- seos_df[seos_df$s == s,]
  
  nint <- 0
  for(d in na.omit(seos[which(seos$start), 'date'])){
    d <- as.Date(d)
    if(nint > 0){
      if(d %in% seq.Date(as.Date(growth_int[1]), as.Date(growth_int[2]), 'day')){
        next
      }else{
        end_id <- min(which(as.Date(seos[which(seos$end), 'date']) > d))
        end <- as.character(seos[which(seos$end), 'date'][end_id])
        end <- ifelse(is.na(end), as.character(final), end)
        growth_int <- c(as.character(d), as.character(end))
        growing_intervals <- append(growing_intervals, list(growth_int) )
        nint <- nint + 1
      }
    }else{
      end_id <- min(which(as.Date(seos[which(seos$end), 'date']) > d))
      end <- as.character(seos[which(seos$end), 'date'][end_id])
      end <- ifelse(is.na(end), as.character(final), end)
      growth_int <- c(as.character(d), as.character(end))
      growing_intervals <- append(growing_intervals, list(growth_int) )
      nint <- nint + 1
    }
  }
  
  if(length(growing_intervals) == 0){next}
  
  grint_df <- rbind(
    grint_df,
    data.frame(
      s = s,
      start = unlist(lapply(growing_intervals, min)),
      end = unlist(lapply(growing_intervals, max)),
      int = 1:nint
    )
  )
  
}


saveRDS(grint_df, file.path(wd, 'output/climate', 'growingseason_highelevsites.rds'))

