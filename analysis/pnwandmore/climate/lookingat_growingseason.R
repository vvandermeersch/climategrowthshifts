
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
final <- "2020-12-31"

# For one stand
clim_subset <- clim[clim$dataset %in% sample(highelev_ids, 1),]
clim_subset$year <- lubridate::year(clim_subset$date)
clim_subset <- clim_subset[clim_subset$date %in% seq.Date(as.Date(init), as.Date(final), "day"),]
ggplot() +
  geom_line(data = clim_subset[clim_subset$var == 'tair',], aes(x = doy, y = value, group = year), 
            color = "#F5CBCB", alpha = 0.1) +
  geom_line(data = clim_subset[clim_subset$var == 'soilmoist_1040',], aes(x = doy, y = value*100, group = year), 
            color = '#748DAE', alpha = 0.1) +
  stat_summary(data = clim_subset[clim_subset$var == 'tair',], aes(x = doy, y = value),
               fun=mean, colour="white", geom="line", linewidth = 2) +
  stat_summary(data = clim_subset[clim_subset$var == 'tair',], aes(x = doy, y = value),
               fun=mean, colour="#F5CBCB", geom="line",  linewidth = 1) +
  geom_hline(aes(yintercept = 5), colour="#F5CBCB", linetype = 'dashed', linewidth = 0.5) + 
  stat_summary(data = clim_subset[clim_subset$var == 'soilmoist_1040',], aes(x = doy, y = value*100),
               fun=mean, colour="white", geom="line", linewidth = 2) +
  stat_summary(data = clim_subset[clim_subset$var == 'soilmoist_1040',], aes(x = doy, y = value*100),
               fun=mean, colour="#748DAE", geom="line",  linewidth = 1) +
  geom_hline(aes(yintercept = 15), colour="#748DAE", linetype = 'dashed', linewidth = 0.5) + 
  theme_classic() +
  coord_cartesian(ylim = c(-15,40), expand = FALSE) +
  labs(y = 'degC / %', x = 'DOY', title = 'California (ca667)') +
  theme(plot.title = element_text(size = 10, hjust = 0.5))
      
            
# find the possible SOS and EOS
ncriticaldays <- 5

clim_subset$start <- NA
clim_subset$end <- NA
for(d in as.character(seq.Date(as.Date(init), as.Date(final), "day"))){
  d <- as.Date(d)
  last5days <- seq.Date(d %m-% days(ncriticaldays-1), d, "day")
  thermal <- sum(clim_subset[clim_subset$date %in% last5days & clim_subset$var %in% 'tair', 'value'] >= 5) == ncriticaldays
  moisture <- sum(clim_subset[clim_subset$date %in% last5days & clim_subset$var %in% 'soilmoist_1040', 'value'] >= 0.15) == ncriticaldays
  clim_subset[clim_subset$date %in% d, "start"] <- thermal & moisture
  
  thermal <- sum(clim_subset[clim_subset$date %in% last5days & clim_subset$var %in% 'tair', 'value'] < 5) == ncriticaldays
  moisture <- sum(clim_subset[clim_subset$date %in% last5days & clim_subset$var %in% 'soilmoist_1040', 'value'] < 0.15) == ncriticaldays
  clim_subset[clim_subset$date %in% d, "end"] <- thermal | moisture
  
}

# find the growing intervals
clim_subset <- clim_subset[order(clim_subset$date),]
growing_intervals <- list()
nint <- 0
for(d in na.omit(as.character(clim_subset[clim_subset$var == 'tair' & clim_subset$start, 'date']))){
  d <- as.Date(d)
  if(nint > 0){
    if(d %in% seq.Date(as.Date(growth_int[1]), as.Date(growth_int[2]), 'day')){
      next
    }else{
      end_id <- min(which(as.Date(clim_subset[clim_subset$var == 'tair' & clim_subset$end, 'date']) > d))
      end <- as.character(clim_subset[clim_subset$var == 'tair' & clim_subset$end, 'date'][end_id])
      end <- ifelse(is.na(end), as.character(final), end)
      growth_int <- c(as.character(d), as.character(end))
      growing_intervals <- append(growing_intervals, list(growth_int) )
      nint <- nint + 1
    }
  }else{
    end_id <- min(which(as.Date(clim_subset[clim_subset$var == 'tair' & clim_subset$end, 'date']) > d))
    end <- as.character(clim_subset[clim_subset$var == 'tair' & clim_subset$end, 'date'][end_id])
    end <- ifelse(is.na(end), as.character(final), end)
    growth_int <- c(as.character(d), as.character(end))
    growing_intervals <- append(growing_intervals, list(growth_int) )
    nint <- nint + 1
  }
}
growing_intervals_plot <- data.frame(
  start = unlist(lapply(growing_intervals, min)),
  end = unlist(lapply(growing_intervals, max)),
  int = 1:nint
)
growing_intervals_plot$startdoy <- lubridate::yday(growing_intervals_plot$start)
growing_intervals_plot$enddoy <- lubridate::yday(growing_intervals_plot$end)
growing_intervals_plot$year <- lubridate::year(growing_intervals_plot$start)

ggplot() +
  geom_line(data = clim_subset[clim_subset$var == 'tair',], aes(x = doy, y = value, group = year), 
            color = "#F5CBCB", alpha = 0.1) +
  geom_line(data = clim_subset[clim_subset$var == 'soilmoist_1040',], aes(x = doy, y = value*100, group = year), 
            color = '#748DAE', alpha = 0.1) +
  stat_summary(data = clim_subset[clim_subset$var == 'tair',], aes(x = doy, y = value),
               fun=mean, colour="white", geom="line", linewidth = 2) +
  stat_summary(data = clim_subset[clim_subset$var == 'tair',], aes(x = doy, y = value),
               fun=mean, colour="#F5CBCB", geom="line",  linewidth = 1) +
  geom_hline(aes(yintercept = 5), colour="#F5CBCB", linetype = 'dashed', linewidth = 0.5) + 
  stat_summary(data = clim_subset[clim_subset$var == 'soilmoist_1040',], aes(x = doy, y = value*100),
               fun=mean, colour="white", geom="line", linewidth = 2) +
  stat_summary(data = clim_subset[clim_subset$var == 'soilmoist_1040',], aes(x = doy, y = value*100),
               fun=mean, colour="#748DAE", geom="line",  linewidth = 1) +
  geom_hline(aes(yintercept = 15), colour="#748DAE", linetype = 'dashed', linewidth = 0.5) + 
  geom_segment(data = growing_intervals_plot, linewidth = 3, color = '#9ab784', alpha = 0.1,
               aes(x= startdoy, xend = enddoy, y = -15, yend = -15, group = paste0(int, year))) +
  theme_classic() +
  coord_cartesian(ylim = c(-20,40), expand = FALSE) +
  labs(y = 'degC / %', x = 'DOY') +
  theme(plot.title = element_text(size = 10, hjust = 0.5))
