rm(list = ls())
wd <- '/home/victor/projects/climategrowthshifts/analysis/pnw'

# ------------------- #
# Download ITRDB data #
# ------------------- #
if(FALSE){
  download_until_success <- function(url, destfile, ..., maxcount = 5) {
    count <- 0
    repeat{
      Sys.sleep(0.5)
      try(download.file(url, destfile, ...))
      count <- count + 1
      if (file.exists(destfile) || count >= maxcount){break}
    }
  }
  states <- c('az', 'ca', 'id', 'nv', 'or', 'ut', 'wa',
              'can', 'cana')
  raw_dir <- file.path(wd, "data", "itrdb", "usa")
  server <- "https://www.ncei.noaa.gov/pub/data/paleo/treering/measurements/northamerica/usa"
  filenames <- unlist(xml2::as_list(httr::content(httr::GET(server))))
  filenames <- filenames[grepl(filenames, pattern = '*.rwl$') & !grepl(filenames, pattern = '*-noaa.rwl$')]
  for (file in filenames) {
    id <- gsub("\\.rwl", "", file)
    if(id == "or015b"){next}
    state <- sub("(\\D*).*", "\\1", id)
    if(!(state %in% states)){
      next
    }
    type <- regmatches(id, regexec("[0-9]+(.*)", id))[[1]][2]
    if(type == ''){
      download_until_success(file.path(server, file), file.path(raw_dir, file))
    }else if(type %in% c('l', 'e', 'i', 'n', 't', 'x', 'd', 'p')){
      # we don't want late/early wood width or density
      next
    }else{
      stop(paste0(file))
    }
  }
  raw_dir <- file.path(wd, "data", "itrdb", "canada")
  server <- "https://www.ncei.noaa.gov/pub/data/paleo/treering/measurements/northamerica/canada"
  filenames <- unlist(xml2::as_list(httr::content(httr::GET(server))))
  filenames <- filenames[grepl(filenames, pattern = '*.rwl$') & !grepl(filenames, pattern = '*-noaa.rwl$')]
  for (file in filenames) {
    id <- gsub("\\.rwl", "", file)
    if(id == "cana" | sub("(.*\\d)[^\\d]*$", "\\1", id) == 'cana168'){next}
    state <- sub("(\\D*).*", "\\1", id)
    if(!(state %in% states)){
      next
    }
    type <- regmatches(id, regexec("[0-9]+(.*)", id))[[1]][2]
    if(type == '' | type == 'rw'){
      download_until_success(file.path(server, file), file.path(raw_dir, file))
    }else if(type %in% c('l', 'e', 'i', 'n', 't', 'x', 'd', 'p',
                         'rlw', 'rew', 'rd', 'ld', 'mxd', 'mnd', 'lw',
                         'ew', 'ed', 'bm', 'ba')){
      # we don't want late/early wood width or density
      next
    }else{
      stop(paste0(file))
    }
  }
}

# ---------------------- #
# Pre-process ITRDB data #
# ---------------------- #
ringwidth_series <- data.frame()
datasets_summary <- data.frame()
years <- c(1980:2025)
min_nyears <- 20 # minimum timeseries length

# US datasets
raw_dir <- file.path(wd, "data", "itrdb", "usa")
rwl_files <- list.files(raw_dir, pattern = '.rwl') 
for(f in rwl_files){
  state <- sub("(\\D*).*", "\\1", f)
  if(!(state %in% states)){next}
  
  dataset <- sub("\\.[^.]*$", "", f)
  # empty datasets, or other issues
  # not critical because they don't span our perido of interest (I think?)
  if(dataset %in% c('or042x', 'ut528', 'ca726', 'az580', 'ca581')){next}
  
  rwdat <- dplR::read.rwl(file.path(raw_dir, f), format = 'tucson', verbose = FALSE)
  # tryCatch({m   <-  dplR::read.rwl(file.path(raw_dir, f), format = 'tucson', verbose = FALSE)}, warning=function(w) print(f))
  idrows <- which(rownames(rwdat) %in% as.character(years))
  if(length(idrows) >= min_nyears){
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
    
    # datasets with obvious errors in the end-of-series marker (should be -9999)
    if(dataset %in% c('ca684', 'ca681', 'ca682', 'ca683', 'ca685', 'ca686',
                      'nv524', 'nv525', 'nv526', 'nv523', 'nv527', 'or090')){
      rwdat_long$rw_mm <- rwdat_long$rw_mm/10 # all the cores have wrong end-of-series markers
    }
    # if(dataset %in% c('or090')){
    #   coresid <- c('LCF16A', 'LCF30B', 'LCF37B') # three cores obviously have wrong end-of-series markers
    #   rwdat_long[rwdat_long$core_id %in% coresid, 'rw_mm'] <- rwdat_long[rwdat_long$core_id %in% coresid, 'rw_mm']/10
    # }
    
    ringwidth_series <- rbind(ringwidth_series, rwdat_long)
  }else{next}
}
saveRDS(ringwidth_series, file.path(wd, 'input', 'itrdb', 'ringwidth_series_usonly.rds'))




# Canada datasets
raw_dir <- "/home/victor/projects/climategrowthshifts/analysis/pnw/data/itrdb_victor/canada"
rwl_files <- list.files(raw_dir, pattern = '.rwl') 
for(f in rwl_files){
  state <- sub("(\\D*).*", "\\1", f)
  if(!(state %in% states)){next}
  
  dataset <- sub("\\.[^.]*$", "", f)
  # empty datasets, or other issues
  # not critical because they don't span our perido of interest
  if(dataset %in% c('cana097p', 'cana362')){next}
  
  rwdat <- dplR::read.rwl(file.path(raw_dir, f), format = 'tucson', verbose = FALSE)
  # tryCatch({m   <-  dplR::read.rwl(file.path(raw_dir, f), format = 'tucson', verbose = FALSE)}, warning=function(w) print(f))
  idrows <- which(rownames(rwdat) %in% as.character(years))
  if(length(idrows) >= min_nyears){
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
    
    # datasets with obvious errors in the end-of-series marker (should be -9999)
    if(dataset %in% c('ca684', 'ca681', 'ca682', 'ca683', 'ca685', 'ca686',
                      'nv524', 'nv525', 'nv526', 'nv523', 'nv527', 'or090')){
      rwdat_long$rw_mm <- rwdat_long$rw_mm/10 # all the cores have wrong end-of-series markers
    }
    # if(dataset %in% c('or090')){
    #   coresid <- c('LCF16A', 'LCF30B', 'LCF37B') # three cores obviously have wrong end-of-series markers
    #   rwdat_long[rwdat_long$core_id %in% coresid, 'rw_mm'] <- rwdat_long[rwdat_long$core_id %in% coresid, 'rw_mm']/10
    # }
    
    ringwidth_series <- rbind(ringwidth_series, rwdat_long)
  }else{next}
}
