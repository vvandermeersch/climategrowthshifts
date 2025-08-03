wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(rnaturalearth)
library(ggplot2)
library(terra)
library(tidyterra)
library(geodata)
world <- vect(ne_countries(scale = "medium", returnclass = "sf"))
westam <- vect(ne_states(country = c("Canada", "Mexico", "United States of America"), returnclass = "sf"))

datasets <- readRDS(file.path(wd, 'output/model', 'datasets_11july2025.rds'))

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
  geom_point(data = speciescount[speciescount$grouped_stand == 'S469',], aes(x = east_lon, y = north_lat), col = "red", fill = 'red', shape = 16,size = 2) +
  ggplot2::theme_minimal() +
  theme(axis.title = ggplot2::element_blank()) +
  ggplot2::coord_sf(xlim = c(-124.925, -89.025), ylim = c(25.065,52.925)) +
  labs(size = '# species')

standcount <- aggregate(grouped_stand ~ species_code, datasets, function(x) length(unique(x)))
names(standcount)[which(names(standcount)=='grouped_stand')] <- 'count'
names(standcount)[which(names(standcount)=='species_code')] <- 'shortname'
standcount <- merge(standcount, sppfull[,c('shortname', 'phylo.name')])
names(standcount)[names(standcount)=='phylo.name'] <- 'label'

# source(file.path(wd, 'getphylo.R'))
# standcount <- standcount[match(phy.plants.here$tip.label, standcount$label),]
# phy.plants.here$tip.label <- paste0(standcount$label, ' (n=',standcount$count,')')
# plotTree(phy.plants.here)
# ## reset margins to default
# par(mar=c(5.1,4.1,4.1,2.1))


# alt_N00W120 <- rast(file.path(wd, 'data/srtm', 'N00W120', 'cut_n00w120.tif'))
# alt_N30W120 <- rast(file.path(wd, 'data/srtm', 'N30W120', 'cut_n30w120.tif'))
# alt_N30W150 <- rast(file.path(wd, 'data/srtm', 'N30W150', 'cut_n30w150.tif'))
# alt <- terra::sprc(alt_N00W120, alt_N30W120, alt_N30W150)
# alt <- mosaic(alt)
# alt_agg <- aggregate(alt, 50, mean)
# slope <- terrain(alt_agg, "slope", unit = "radians")
# aspect <- terrain(alt_agg, "aspect", unit = "radians")
# hill <- shade(slope, aspect, 40, 270)
# writeRaster(alt_agg, file.path(wd, 'data/srtm', 'alt_agg.tif'))
# writeRaster(hill, file.path(wd, 'data/srtm', 'hill.tif'))


alt_agg <- rast(file.path(wd, 'data/srtm', 'alt_agg.tif'))
hill <- rast(file.path(wd, 'data/srtm', 'hill.tif'))
# hill_mask <- hill
# hill_mask <- ifel(alt_agg< 1000, NA, hill_mask)
# hill <- mask(hill, hill_mask)

grad_hypso <- hypso.colors(10, "wiki-2.0_hypso")
pal_greys <- hcl.colors(10, "Grays")
ggplot() +
  geom_spatraster(data = hill) +
  scale_fill_gradientn(colors = pal_greys, na.value = NA) +
  ggnewscale::new_scale_fill() +
  geom_spatraster(data = alt_agg,  alpha = 0.6) +
  scale_fill_gradientn(colours = grad_hypso, na.value = NA) +
  geom_spatvector(data = westam, color = 'white', fill = NA) +
  geom_point(data = speciescount, aes(x = east_lon, y = north_lat), col = "black", fill = 'white', shape = 21, size = 2, stroke = 1.5) +
  # geom_point(data = speciescount[speciescount$grouped_stand == 'S469',], aes(x = east_lon, y = north_lat), col = "red", fill = 'red', shape = 16,size = 2) +
  ggplot2::theme_void() +
  theme(axis.title = ggplot2::element_blank()) +
  ggplot2::coord_sf(xlim = c(-129.925, -99.025), ylim = c(25.065,52.925), expand = TRUE) +
  labs(size = '# species') +
  theme(legend.position = 'none')

geom_spatraster(data = alt_N00W120,  alpha = 0.4) +grad_hypso <- hypso.colors2(10, "dem_poster")
autoplot(alt_N00W120) +
  scale_fill_gradientn(colours = grad_hypso, na.value = NA) +
  ggplot2::coord_sf(xlim = c(-124.925, -89.025), ylim = c(25.065,52.925)) 


index <- hill %>%
  mutate(index_col = rescale(shades, to = c(1, length(pal_greys)))) %>%
  mutate(index_col = round(index_col)) %>%
  pull(index_col)


# Get cols
vector_cols <- pal_greys[index]

# Need to avoid resampling
# and dont use aes

hill_plot <- ggplot() +
  geom_spatraster(
    data = hill, fill = vector_cols, maxcell = Inf,
    alpha = 1
  )

hill_plot

