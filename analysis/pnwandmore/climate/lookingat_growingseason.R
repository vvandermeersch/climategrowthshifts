
library(ggplot2)
library(rnaturalearth)
library(terra)
library(tidyterra)
library(lubridate)

wd <- "/home/victor/projects/climategrowthshifts/analysis/pnwandmore"
clim <- readRDS(file.path(wd, 'output/climate', 'climvar_04july2025.rds'))
datasets <- readRDS(file.path(wd, 'output/model', 'datasets_11july2025.rds'))
westam <- vect(ne_states(country = c("Canada", "Mexico", "United States of America"), returnclass = "sf"))

dataset_ids <- unique(clim$dataset) # c('mexi076', 'or099', 'ut544', 'can719', 'ca729')

ggplot() +
  geom_spatvector(data = westam, color = 'white', fill = 'grey90') +
  geom_point(data = datasets[datasets$dataset %in% c('can719', 'ca729'),], 
             aes(x = east_lon, y = north_lat), size = 3, color = "white") +
  geom_point(data = datasets[datasets$dataset %in% c('can719', 'ca729'),], 
             aes(x = east_lon, y = north_lat), size = 1.5, color = "black") +
  coord_sf(xlim = c(-124.925, -105.025), ylim = c(25.065,52.925), expand = TRUE) +
  theme_void()


clim_subset <- clim[clim$dataset %in% c('ca656'),]
clim_subset$year <- lubridate::year(clim_subset$date)
clim_subset <- clim_subset[clim_subset$date %in% seq.Date(as.Date("2016-09-01"), as.Date("2020-02-28"), "day"),]
ggplot() +
  geom_line(data = clim_subset[clim_subset$var == 'tair',], aes(x = date, y = value, group = year), 
            color = "#F5CBCB", alpha = 0.1) +
  geom_line(data = clim_subset[clim_subset$var == 'soilmoist_1040',], aes(x = date, y = value*100, group = year), 
            color = '#748DAE', alpha = 0.1) +
  stat_summary(data = clim_subset[clim_subset$var == 'tair',], aes(x = date, y = value),
               fun=mean, colour="white", geom="line", linewidth = 2) +
  stat_summary(data = clim_subset[clim_subset$var == 'tair',], aes(x = date, y = value),
               fun=mean, colour="#F5CBCB", geom="line",  linewidth = 1) +
  geom_hline(aes(yintercept = 5), colour="#F5CBCB", linetype = 'dashed', linewidth = 0.5) + 
  stat_summary(data = clim_subset[clim_subset$var == 'soilmoist_1040',], aes(x = date, y = value*100),
               fun=mean, colour="white", geom="line", linewidth = 2) +
  stat_summary(data = clim_subset[clim_subset$var == 'soilmoist_1040',], aes(x = date, y = value*100),
               fun=mean, colour="#748DAE", geom="line",  linewidth = 1) +
  geom_hline(aes(yintercept = 15), colour="#748DAE", linetype = 'dashed', linewidth = 0.5) + 
  theme_classic() +
  coord_cartesian(ylim = c(-10,40), expand = FALSE) +
  labs(y = 'degC / %', x = 'DOY', title = 'California (ca729)') +
  theme(plot.title = element_text(size = 10, hjust = 0.5))
      
            
# find the possible SOS and EOS
clim_subset$start <- NA
clim_subset$end <- NA
for(d in as.character(seq.Date(as.Date("2016-09-01"), as.Date("2019-12-31"), "day"))){
  d <- as.Date(d)
  last5days <- seq.Date(d %m-% days(4), d, "day")
  thermal <- sum(clim_subset[clim_subset$date %in% last5days & clim_subset$var %in% 'tair', 'value'] >= 5) == 5
  moisture <- sum(clim_subset[clim_subset$date %in% last5days & clim_subset$var %in% 'soilmoist_1040', 'value'] >= 0.15) == 5
  clim_subset[clim_subset$date %in% d, "start"] <- thermal & moisture
  
  thermal <- sum(clim_subset[clim_subset$date %in% last5days & clim_subset$var %in% 'tair', 'value'] < 5) == 5
  moisture <- sum(clim_subset[clim_subset$date %in% last5days & clim_subset$var %in% 'soilmoist_1040', 'value'] < 0.15) == 5
  clim_subset[clim_subset$date %in% d, "end"] <- thermal | moisture
  
}

# find the growing intervals
clim_subset <- clim_subset[order(clim_subset$date),]
growing_intervals <- list()
nint <- 0
for(d in na.omit(as.character(clim_subset[clim_subset$var == 'tair' & clim_subset$start, 'date']))){
  print(d)
  d <- as.Date(d)
  if(nint > 0){
    if(d %in% seq.Date(as.Date(growth_int[1]), as.Date(growth_int[2]), 'day')){
      next
    }else{
      end_id <- min(which(as.Date(clim_subset[clim_subset$var == 'tair' & clim_subset$end, 'date']) > d))
      end <- as.character(clim_subset[clim_subset$var == 'tair' & clim_subset$end, 'date'][end_id])
      end <- ifelse(is.na(end), as.character("2019-12-31"), end)
      growth_int <- c(as.character(d), as.character(end))
      growing_intervals <- append(growing_intervals, list(growth_int) )
      nint <- nint + 1
    }
  }else{
    end_id <- min(which(as.Date(clim_subset[clim_subset$var == 'tair' & clim_subset$end, 'date']) > d))
    end <- as.character(clim_subset[clim_subset$var == 'tair' & clim_subset$end, 'date'][end_id])
    end <- ifelse(is.na(end), "2019-12-31", end)
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





ggplot() +
  geom_line(data = clim_subset[clim_subset$var == 'tair',], aes(x = date, y = value, group = year), 
            color = "#F5CBCB", alpha = 0.1) +
  geom_line(data = clim_subset[clim_subset$var == 'soilmoist_1040',], aes(x = date, y = value*100, group = year), 
            color = '#748DAE', alpha = 0.1) +
  stat_summary(data = clim_subset[clim_subset$var == 'tair',], aes(x = date, y = value),
               fun=mean, colour="white", geom="line", linewidth = 2) +
  stat_summary(data = clim_subset[clim_subset$var == 'tair',], aes(x = date, y = value),
               fun=mean, colour="#F5CBCB", geom="line",  linewidth = 1) +
  geom_hline(aes(yintercept = 5), colour="#F5CBCB", linetype = 'dashed', linewidth = 0.5) + 
  stat_summary(data = clim_subset[clim_subset$var == 'soilmoist_1040',], aes(x = date, y = value*100),
               fun=mean, colour="white", geom="line", linewidth = 2) +
  stat_summary(data = clim_subset[clim_subset$var == 'soilmoist_1040',], aes(x = date, y = value*100),
               fun=mean, colour="#748DAE", geom="line",  linewidth = 1) +
  geom_hline(aes(yintercept = 15), colour="#748DAE", linetype = 'dashed', linewidth = 0.5) + 
  geom_segment(data = growing_intervals_plot, linewidth = 2, color = '#9ab784',
               aes(x= as.Date(start), xend = as.Date(end), y = -5, yend = -5, group = int)) +
  theme_classic() +
  coord_cartesian(ylim = c(-10,40), expand = FALSE) +
  labs(y = 'degC / %', x = 'DOY') +
  theme(plot.title = element_text(size = 10, hjust = 0.5))
