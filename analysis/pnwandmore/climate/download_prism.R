

server <- 'https://data.prism.oregonstate.edu/time_series'
country <- 'us'
network <- 'an'
res <- '800m'
resalt <- '30s'
temp <- 'monthly'

vars <- c('tmean', 'ppt', 'vpdmax')
years <- 1895:2024
months <- c('01', '02', '03', '04', '05', '06',
            '07', '08', '09', '10', '11', '12')


for(var in vars){
  destdir <- file.path('/home/vvandermeersch/data/climate/prism/800m/us', var)
  dir.create(destdir)
  
  for(year in years){
    
    for(month in months){
      
      filename <- paste0(paste('prism', var, country, resalt, paste0(year, month), sep = '_'), '.zip')
      download.file(file.path(server, country, network, res, var, temp, year, filename), destfile = file.path(destdir, filename))
      
    }
  }
}

for(var in vars){
  destdir <- file.path('/home/vvandermeersch/data/climate/prism/800m/us', var)
  for(year in years){
    for(month in months){
      
      zipname <- paste0(paste('prism', var, country, resalt, paste0(year, month), sep = '_'), '.zip')
      filename <- paste0(paste('prism', var, country, resalt, paste0(year, month), sep = '_'), '.tif')
      utils::unzip(file.path(destdir, zipname), exdir = destdir, files = filename)
      
    }
  }
}

prism_tdmean_us_30s_1895.tif

