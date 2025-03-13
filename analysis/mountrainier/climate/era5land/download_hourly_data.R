library(reticulate)

wd <- "~/projects/climategrowthshifts/analysis/mountrainier"

datafolder <- file.path(wd, 'data', 'climate', 'era5land')
setwd(datafolder) # folder where files will be downloaded

# Python config
use_virtualenv(file.path(wd, 'climate', 'era5land', "bcookenv"), required=TRUE)
cdsapi <- import("cdsapi")
source_python(file.path(wd, 'climate', 'era5land', "download_hourly_era5land.py")) # custom script

#Run script
for(year in as.character(1950:1980)){
  download_era5land(year = year)
}

