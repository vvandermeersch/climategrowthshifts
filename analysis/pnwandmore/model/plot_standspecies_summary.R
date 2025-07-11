
library(rnaturalearth)
library(ggplot2)
library(terra)
library(tidyterra)
world <- vect(ne_countries(scale = "medium", returnclass = "sf"))
westam <- vect(ne_states(country = c("Canada", "Mexico", "United States of America"), returnclass = "sf"))


# Gather stands if they have the same latitude and longitude (rounded to 1e-2 degree, ie ~1.1km)
group_keys <- interaction(
  round(datasets$north_lat,2),
  round(datasets$south_lat,2),
  round(datasets$east_lon,2),
  round(datasets$west_lon,2),
  drop = TRUE
)
datasets$grouped_stand <- paste0("S", as.integer(factor(group_keys)))
datasets$east_lon <- round(datasets$east_lon, 2)
datasets$north_lat <- round(datasets$north_lat, 2)

# datasets_woAngiosperms <- datasets[!(datasets$species_code %in% todrop),]

datasetscount <- aggregate(dataset ~ grouped_stand + east_lon + north_lat, datasets, function(x) length(unique(x)))
sum(datasetscount$dataset)
speciescount <- aggregate(species_code ~ grouped_stand + east_lon + north_lat, datasets, function(x) length(unique(x)))
names(speciescount)[which(names(speciescount)=='species_code')] <- 'count'

ggplot() +
  geom_spatvector(data = westam, color = 'white', fill = 'grey90') +
  geom_point(data = speciescount, aes(x = east_lon, y = north_lat, size = count), col = "black", fill = NA, shape = 1) +
  ggplot2::theme_minimal() +
  theme(axis.title = ggplot2::element_blank()) +
  ggplot2::coord_sf(xlim = c(-124.925, -89.025), ylim = c(25.065,52.925))

standcount <- aggregate(grouped_stand ~ species_code, datasets, function(x) length(unique(x)))
names(standcount)[which(names(standcount)=='grouped_stand')] <- 'count'
names(standcount)[which(names(standcount)=='species_code')] <- 'shortname'
standcount <- merge(standcount, sppfull[,c('shortname', 'phylo.name')])
names(standcount)[names(standcount)=='phylo.name'] <- 'label'

source(file.path(wd, 'getphylo.R'))
standcount <- standcount[match(phy.plants.here$tip.label, standcount$label),]
phy.plants.here$tip.label <- paste0(standcount$label, ' (n=',standcount$count,')')
plotTree(phy.plants.here)
## reset margins to default
par(mar=c(5.1,4.1,4.1,2.1))