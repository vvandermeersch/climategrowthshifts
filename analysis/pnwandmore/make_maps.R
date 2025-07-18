rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"

library(rnaturalearth)
library(ggplot2)
library(terra)
library(tidyterra)
library(patchwork)
world <- vect(ne_countries(scale = "medium", returnclass = "sf"))
westam <- vect(ne_states(country = c("Canada", "Mexico", "United States of America"), returnclass = "sf"))

# Load treering data
datasets <- readRDS(file.path(wd, 'input', 'itrdb', 'datasets_summary_all.rds'))

# Temporary - missing stands due to land mask
datasets <- datasets[datasets$dataset != 'wa149',] # temporary 
datasets <- datasets[datasets$dataset != 'ca673',] # temporary
datasets <- datasets[datasets$dataset != 'az563',] # temporary

# Temporary - missing species 
datasets <- datasets[datasets$dataset != 'co689',] # temporary 
datasets <- datasets[datasets$dataset != 'co690',] # temporary
datasets <- datasets[datasets$dataset != 'co691',] # temporary

# datasets <- datasets[datasets$dataset != 'can673',] # temporary

# toremove <- c("az589","ca681","ca682","ca683","ca684","ca685","ca736","mexi060","mexi063","mt130","ut576")
# datasets <- datasets[!(datasets$dataset %in% toremove),]

# Keep only datasets within WLDAS extent
datasets <- datasets[datasets$north_lat >=25.065 & datasets$north_lat <=52.925 & datasets$east_lon <= -89.025 &  datasets$east_lon >= -124.925,]

datasets <- datasets[datasets$last_year >= 1999,] # at least 20 years of observations

# Same species
datasets[datasets$species_code == 'ABBI', c('species_name', 'species_code')] <- 
  unique(datasets[datasets$species_code == 'ABLA', c('species_name', 'species_code')])

# Get shapefile name
temp <- str_split_fixed(datasets$species_name, " ", 3)
datasets$species_name_simplified <- paste(temp[,1], temp[,2], sep=" ")
species_treerings <- data.frame(species = unique(datasets$species_name_simplified))

species_little <- read.csv(file.path(wd, 'data', 'USTreeAtlas', 'Little_datatable.csv'))
species_little <- species_little[,c('Latin.Name', 'SHP..')]
names(species_little) <- c('species', 'shp')

species_treerings <- merge(species_treerings, species_little,  all.x = TRUE)

# make species maps
for(i in 1:nrow(species_treerings)){
  sp <- species_treerings[i, 'species']
  shpname <- species_treerings[i, 'shp']
  rangemap <- aggregate(vect(file.path(wd, 'data', 'USTreeAtlas', 'shp', shpname, paste0(shpname, '.shp'))))
  crs(rangemap) <- 'EPSG:4267'
  plotmap_ext <- ggplot() +
    geom_spatvector(data = westam, color = 'white', fill = 'grey90') +
    geom_spatvector(data = rangemap, color = 'white', fill = '#51ba69') +
    annotate('rect', xmin = -124.925, xmax = -89.025, ymin = 25.065, ymax = 52.925,
             fill = NA, color = 'grey30', linewidth = 0.5) +
    ggplot2::theme_minimal() +
    theme(axis.title = ggplot2::element_blank()) +
    ggplot2::coord_sf(xlim = c(-160, -87), ylim = c(20,70), expand = FALSE)
  ggsave(plotmap_ext, file = file.path(wd, 'figures', 'rangemaps', paste0( tolower(gsub(" ", "_", sp)), '_extended.pdf')),
         height = 8, width = 8)
  
  totsurf <- expanse(rangemap, unit = 'km')
  wldassurf <- expanse(crop(rangemap, ext(-124.925, -89.025, 25.065,52.925)), unit = 'km')
  
  records <- vect(unique(datasets[datasets$species_name_simplified == sp, c('north_lat', 'east_lon')]), 
                  geom = c('east_lon', 'north_lat'))
  crs(records) <- 'EPSG:4267'
  plotmap <- ggplot() +
    geom_spatvector(data = westam, color = 'white', fill = 'grey90') +
    geom_spatvector(data = rangemap, color = 'white', fill = '#51ba69') +
    geom_spatvector(data = records, color = 'white', fill = 'grey30', shape = 21) +
    ggplot2::theme_minimal() +
    theme(axis.title = ggplot2::element_blank()) +
    ggplot2::coord_sf(xlim = c(-124.925, -89.025), ylim = c(25.065,52.925), expand = FALSE) +
    annotate('rect', xmin = -100, xmax = -90, ymin = 50, ymax = 52,
             fill = 'white', color = 'grey30', linewidth = 0.5, alpha = 0.7) +
    annotate('text', x = -99.5, y = 51.5, label = paste0('No. of records: ', dim(records)[1]),
             color = 'grey30', hjust = 0) +
    annotate('text', x = -99.5, y = 50.7, label = paste0('Range in WLDAS: ',  round(wldassurf/totsurf*100,1),'%'),
             color = 'grey30', hjust = 0) 
  ggsave(plotmap, file = file.path(wd, 'figures', 'rangemaps', paste0( tolower(gsub(" ", "_", sp)), '.pdf')),
         height = 8, width = 8)
  
}

# make Angio/Gymno maps
c('PLRA', 'QUDG', 'QULO', 'QUGA', 'PPFR', 'PPTR', 'PPDE')
angiosperms <- c('Platanus racemosa', 'Populus deltoides', 'Populus fremontii', 'Populus tremuloides',
                 'Quercus douglasii', 'Quercus gambelii', 'Quercus lobata')


angiospermsrange <- lapply(which(species_treerings$species %in% angiosperms),
       function(i) vect(file.path(wd, 'data', 'USTreeAtlas', 'shp', 
                                  species_treerings[i, 'shp'], paste0(species_treerings[i, 'shp'], '.shp'))))
angiospermsrange <- aggregate(vect(angiospermsrange))
crs(angiospermsrange) <- 'EPSG:4267'
records <- vect(unique(datasets[datasets$species_name_simplified %in% angiosperms, 
                                c('north_lat', 'east_lon')]), geom = c('east_lon', 'north_lat'))
crs(records) <- 'EPSG:4267'
angiomap <- ggplot() +
  geom_spatvector(data = westam, color = 'white', fill = 'grey90') +
  geom_spatvector(data = angiospermsrange, color = 'white', fill = '#73c787') +
  geom_spatvector(data = westam, color = 'white', fill = NA) +
  geom_spatvector(data = records, color = 'white', fill = 'grey30', shape = 21) +
  ggplot2::theme_minimal() +
  theme(axis.title = ggplot2::element_blank()) +
  ggplot2::coord_sf(xlim = c(-124.925, -89.025), ylim = c(25.065,52.925), expand = FALSE) +
  annotate('rect', xmin = -97, xmax = -90, ymin = 50.2, ymax = 52,
           fill = 'white', color = 'grey30', linewidth = 0.5, alpha = 0.7) +
  annotate('text', x = -96.5, y = 51.5, label = 'Angiosperms',
           color = 'grey30', hjust = 0) +
  annotate('text', x = -96.5, y = 50.7, label = paste0('No. of sites: ', dim(records)[1]),
           color = 'grey30', hjust = 0)

gymnospermssrange <- lapply(which(!(species_treerings$species %in% angiosperms)),
                           function(i) vect(file.path(wd, 'data', 'USTreeAtlas', 'shp', 
                                                      species_treerings[i, 'shp'], paste0(species_treerings[i, 'shp'], '.shp'))))
gymnospermssrange <- aggregate(vect(gymnospermssrange))
crs(gymnospermssrange) <- 'EPSG:4267'
records <- vect(unique(datasets[!(datasets$species_name_simplified %in% angiosperms), 
                                c('north_lat', 'east_lon')]), geom = c('east_lon', 'north_lat'))
crs(records) <- 'EPSG:4267'
gymnomap <- ggplot() +
  geom_spatvector(data = westam, color = 'white', fill = 'grey90') +
  geom_spatvector(data = gymnospermssrange, color = 'white', fill = '#73c787') +
  geom_spatvector(data = westam, color = 'white', fill = NA) +
  geom_spatvector(data = records, color = 'white', fill = 'grey30', shape = 21) +
  ggplot2::theme_minimal() +
  theme(axis.title = ggplot2::element_blank()) +
  ggplot2::coord_sf(xlim = c(-124.925, -89.025), ylim = c(25.065,52.925), expand = FALSE) +
  annotate('rect', xmin = -97, xmax = -90, ymin = 50.2, ymax = 52,
           fill = 'white', color = 'grey30', linewidth = 0.5, alpha = 0.7) +
  annotate('text', x = -96.5, y = 51.5, label = 'Gymnosperms',
           color = 'grey30', hjust = 0) +
  annotate('text', x = -96.5, y = 50.7, label = paste0('No. of sites: ', dim(records)[1]),
           color = 'grey30', hjust = 0)

clademaps <- gymnomap + angiomap
ggsave(clademaps, file = file.path(wd, 'figures', 'rangemaps', 'gymno_angio.pdf'),
       height = 8, width = 16)
