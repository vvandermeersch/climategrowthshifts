rm(list=ls())
wd <- "~/projects/climategrowthshifts/analysis/mountrainier"

ailene_gdd <- read.csv(file.path(wd, 'data/treerings_ailene', 'tree_plot_climate_grdd5.csv'))
ailene_swe <- read.csv(file.path(wd, 'data/treerings_ailene', 'tree_plot_climate_swe.csv'))
ailene_snowdur <- read.csv(file.path(wd, 'data/treerings_ailene', 'tree_plot_climate_snowdur.csv'))

ailene_gdd <- ailene_gdd[c('year', 'AV06')]
names(ailene_gdd) <- c('year', 'gddb5 (degC)')
ailene_swe <- ailene_swe[c('Year', 'AV06')]
names(ailene_swe) <- c('year', 'swe (mm)')
# ailene_snowdur <- ailene_snowdur[c('Year', 'AV06')]
# names(ailene_snowdur) <- c('year', 'snowdur (days)')

av06_climate <- merge(ailene_gdd, ailene_swe)

write.csv(av06_climate, file = file.path(wd, 'input/climate/av06_climate.csv'))
