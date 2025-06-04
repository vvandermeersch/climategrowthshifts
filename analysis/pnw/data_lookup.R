rm(list = ls())
dir <- '/home/victor/projects/climategrowthshifts/analysis/pnw/data/itrdb/oregon/total_ring_width/data'

datasets <- list.files(dir, pattern = '.txt', recursive = TRUE, full.names = TRUE)
datasets_short <- list.files(dir, pattern = '.txt', recursive = TRUE, full.names = FALSE)

data_summary <- data.frame()
for(d in 1:length(datasets)){
  
  text <- read.delim(datasets[d], nrows = 120)
  subtext <- text[grep(text[,1], pattern = 'Species_Name'),]
  species_name <- unlist(strsplit(subtext, ': '))[2]
  subtext <- text[grep(text[,1], pattern = 'Tree_Species_Code'),]
  species_code <- unlist(strsplit(subtext, ': '))[2]
  
  subtext <- text[grep(text[,1], pattern = 'Northernmost_Latitude'),]
  north_lat <- as.numeric(unlist(strsplit(subtext, ': '))[2])
  subtext <- text[grep(text[,1], pattern = 'Southernmost_Latitude'),]
  south_lat <- as.numeric(unlist(strsplit(subtext, ': '))[2])
  subtext <- text[grep(text[,1], pattern = 'Easternmost_Longitude'),]
  east_lon <- as.numeric(unlist(strsplit(subtext, ': '))[2])
  subtext <- text[grep(text[,1], pattern = 'Westernmost_Longitude'),]
  west_lon <- as.numeric(unlist(strsplit(subtext, ': '))[2])
  subtext <- text[grep(text[,1], pattern = 'Elevation'),]
  altitude <- as.numeric(unlist(strsplit(subtext, ': '))[2])
  
  dataset_code <- unlist(strsplit(datasets_short[d], '/'))[1]
  
  subtext <- text[grep(text[,1], pattern = 'First_Year'),]
  first_year <- as.numeric(unlist(strsplit(subtext, ': '))[2])
  subtext <- text[grep(text[,1], pattern = 'Last_Year'),]
  last_year <- as.numeric(unlist(strsplit(subtext, ': '))[2])
  
  d_summ<- data.frame(dataset = dataset_code, species_name, species_code, first_year, last_year, north_lat, south_lat, east_lon, west_lon, altitude)
  saveRDS(d_summ, file = file.path(dir, dataset_code, paste0(dataset_code, '_info.rds')))
  
  data_summary <- rbind(data_summary,
                        d_summ)
}

data_summary$uniquesite <- ifelse(data_summary$north_lat == data_summary$south_lat & data_summary$east_lon == data_summary$west_lon, TRUE, FALSE)
saveRDS(data_summary, file = file.path('/home/victor/projects/climategrowthshifts/analysis/pnw/data/itrdb', paste0('itrdb_info.rds')))

# Quick map
library(rnaturalearth)
library(terra)
library(tidyterra)
library(ggplot2)
world <- vect(ne_countries(scale = "medium", returnclass = "sf"))
ggplot() +
  geom_spatvector(data = world, color = 'white', fill = 'grey90') +
  geom_point(data = data_summary, aes(x = east_lon, y = north_lat, color = species_code)) +
  theme_minimal() +
  theme(axis.title = element_blank()) +
  coord_sf(xlim = c(-127, -110), ylim = c(40,52))


unique(data_summary[,c('species_name', 'species_code')])
