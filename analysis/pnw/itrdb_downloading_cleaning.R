
# Download ITRDB data
download_until_success <- function(url, destfile, ..., maxcount = 5) {
  count <- 0
  repeat{
    Sys.sleep(0.5)
    try(download.file(url, destfile, ...))
    count <- count + 1
    if (file.exists(destfile) || count >= maxcount){break}
  }
}
raw_dir <- "/home/victor/projects/climategrowthshifts/analysis/pnw/data/itrdb_victor/usa"
server <- "https://www.ncei.noaa.gov/pub/data/paleo/treering/measurements/northamerica/usa"
filenames <- unlist(xml2::as_list(httr::content(httr::GET(server))))
filenames <- filenames[grepl(filenames, pattern = '*.rwl$') & !grepl(filenames, pattern = '*-noaa.rwl$')]
for (file in filenames) {
  download_until_success(file.path(server, file), file.path(raw_dir, file))
}

# Pre-process ITRDB data
ringwidth_series <- data.frame()
years <- c(1980:2025)
states <- c('az', 'ca', 'id', 'nv', 'or', 'ut', 'wa')
minn_years <- 20 # minimum timeseries length

raw_dir <- "/home/victor/projects/climategrowthshifts/analysis/pnw/data/itrdb_victor/usa"
rwl_files <- list.files(raw_dir, pattern = '.rwl') 

for(f in rwl_files){
  state <- sub("(\\D*).*", "\\1", f)
  if(!(state %in% states)){next}
  
  dataset <- sub("\\.[^.]*$", "", f)
  # there is something wrong with this dataset, and it doesnot span our period of interest
  if(dataset %in% c('or042x')){next}
  
  rwdat <- dplR::read.rwl(file.path(raw_dir, f), format = 'tucson', verbose = FALSE)
  # tryCatch({m   <-  dplR::read.rwl(file.path(raw_dir, f), format = 'tucson', verbose = FALSE)}, warning=function(w) print(f))
  idrows <- which(rownames(rwdat) %in% as.character(years))
  if(length(idrows) >= minn_years){
    #print(f)
    rwdat <- rwdat[idrows,]
    rwdat$year <- rownames(rwdat)
    rwdat_long <- reshape(rwdat, varying = names(rwdat)[names(rwdat) != "year"],
                          v.names = "rw_mm", timevar = "core_id",
                          times = names(rwdat)[names(rwdat) != "year"],
                          direction = "long")
    rwdat_long$tree_id <- sub("^(.*\\d)[^\\d]*$", "\\1", rwdat_long$core_id)
    rwdat_long$dataset <- dataset
    rwdat_long <- rwdat_long[, c('dataset', "year", 'tree_id', "core_id", "rw_mm")]
    
    # datasets with errors in the end-of-marker (should be -9999)
    if(dataset %in% c('ca684', 'ca681', 'ca682', 'ca683', 'ca685', 'ca686',
                      'nv524', 'nv525', 'nv526', 'nv523', 'nv527')){
      rwdat_long$rw_mm <- rwdat_long$rw_mm/10 
    }
    
    ringwidth_series <- rbind(ringwidth_series, rwdat_long)
  }else{next}
}





