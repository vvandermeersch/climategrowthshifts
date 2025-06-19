if (requireNamespace("rstudioapi", quietly = TRUE)) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

rm(list = ls())
library(rnaturalearth)
library(terra)
library(tidyterra)
library(ggplot2)

get_val <- function(text,p,numeric=FALSE) {
  rows <- grep(text[,1],pattern=p)
  if (length(rows)==0) {
    return(NA) #no such pattern found
  }
  subtext <- text[grep(text[,1],pattern=p),]
  split_line <- unlist(strsplit(subtext, ': '))
  if (length(split_line)<2) {
    return(NA) #format error: nothing after ": "
  }
  value <- split_line[2]
  if (numeric) {
    value <- as.numeric(value)
  }
  return(value)
}

dir <- "data/itrdb"
datasets <- list.files(dir,recursive = TRUE,pattern='.txt',full.names = TRUE)
datasets_short <- list.files(dir,pattern=".txt",recursive=TRUE,full.names=FALSE)
data_summary <- data.frame()

for(d in 1:length(datasets)){
  text <- read.delim(datasets[d], nrows = 350)
  
  species_name <- get_val(text,'Species_Name')
  species_code <- get_val(text,'Tree_Species_Code')
  
  north_lat <- get_val(text,'Northernmost_Latitude',numeric=TRUE)
  south_lat <- get_val(text,'Southernmost_Latitude',numeric=TRUE)
  east_lon <- get_val(text,'Easternmost_Longitude',numeric=TRUE)
  west_lon <- get_val(text,'Westernmost_Longitude',numeric=TRUE)
  altitude <- get_val(text,'Elevation',numeric=TRUE)
  
  state <- unlist(strsplit(datasets_short[d], '/'))[1]
  dataset_code <- unlist(strsplit(datasets_short[d], '/'))[4]
  
  first_year <- get_val(text,'First_Year',numeric=TRUE)
  last_year <- get_val(text,'Last_Year',numeric=TRUE)
  d_summ<- data.frame(state = state, dataset = dataset_code, species_name, species_code, first_year, last_year, north_lat, south_lat, east_lon, west_lon, altitude)
  saveRDS(d_summ, file = file.path(str_remove(datasets[d],basename(datasets[d])), paste0(dataset_code, '_info.rds')))
  data_summary <- rbind(data_summary, d_summ)
}

data_summary$uniquesite <- ifelse(data_summary$north_lat == data_summary$south_lat & data_summary$east_lon == data_summary$west_lon, TRUE, FALSE)
saveRDS(data_summary, file = file.path('data/itrdb', paste0('itrdb_info.rds')))


# Quick map
world <- vect(ne_countries(scale = "medium", returnclass = "sf"))
quick_map <- ggplot() +
  geom_spatvector(data = world, color = 'white', fill = 'grey90') +
  geom_point(data = data_summary, aes(x = east_lon, y = north_lat, color = species_code)) +
  theme_minimal() +
  theme(axis.title = element_blank()) +
  coord_sf(xlim = c(-127, -110), ylim = c(40,52))


unique(data_summary[,c('species_name', 'species_code')])
# save the plot in the folder
ggsave(plot=quick_map,filename="quick_map.pdf",path="data/itrdb",width=6,height=4,dpi=300)
