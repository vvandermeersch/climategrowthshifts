library(terra)
library(leaflet)

wd <- "/home/victor/projects/climategrowthshifts/analysis/mountrainier/data"

temp <- rast(file.path(wd, "climate/wldas", "WLDAS_NOAHMP001_DA1_19790303.D10.nc.SUB.nc4"))
plots <- read.csv(file.path(wd, "treerings_ailene", "tree_plot_climate_temp.csv"))
plots <- vect(unique(plots[,c("Longitude", "Latitude")]), geom=c("Longitude", "Latitude"))

plet(temp$SWdown_f_tavg) |> 
  points(plots, col = "white", cex = 2) |> 
  points(plots, col = "black", cex = 0.7)

-122.3,46.6,-121.7,47.2


leaflet() %>% addTiles() %>%
  addRectangles(
    lng1=-122, lat1=46.65,
    lng2=-121.45, lat2=47.05,
    fillColor = "transparent"
  )
