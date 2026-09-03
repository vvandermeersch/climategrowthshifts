rm(list = ls())
wd <- '/home/victor/projects/climategrowthshifts/analysis/pnwandmore'

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
              'wy', 'co', 'mt', 'nm',
              'can', 'cana', 'mexi', 'mex') 
  raw_dir <- file.path(wd, "data", "itrdb", "usa")
  server <- "https://www.ncei.noaa.gov/pub/data/paleo/treering/measurements/northamerica/usa"
  filenames <- unlist(xml2::as_list(httr::content(httr::GET(server))))
  filenames <- filenames[grepl(filenames, pattern = '*.rwl$') & !grepl(filenames, pattern = '*-noaa.rwl$')]
  for (file in filenames) {
    id <- gsub("\\.rwl", "", file)
    if(id == "or015b" | file == 'co.rwl'){next}
    state <- sub("(\\D*).*", "\\1", id)
    if(!(state %in% states)){
      next
    }
    type <- regmatches(id, regexec("[0-9]+(.*)", id))[[1]][2]
    if(type == '' | type == 'w'){
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
  # raw_dir <- file.path(wd, "data", "itrdb", "mexico")
  # server <- "https://www.ncei.noaa.gov/pub/data/paleo/treering/measurements/northamerica/mexico"
  # filenames <- unlist(xml2::as_list(httr::content(httr::GET(server))))
  # filenames <- filenames[grepl(filenames, pattern = '*.rwl$') & !grepl(filenames, pattern = '*-noaa.rwl$')]
  # for (file in filenames) {
  #   id <- gsub("\\.rwl", "", file)
  #   if(id == "cana" | sub("(.*\\d)[^\\d]*$", "\\1", id) == 'cana168'){next}
  #   state <- sub("(\\D*).*", "\\1", id)
  #   if(!(state %in% states)){
  #     next
  #   }
  #   type <- regmatches(id, regexec("[0-9]+(.*)", id))[[1]][2]
  #   if(type == '' | type == 'rw'){
  #     download_until_success(file.path(server, file), file.path(raw_dir, file))
  #   }else if(type %in% c('l', 'e', 'i', 'n', 't', 'x', 'd', 'p',
  #                        'rlw', 'rew', 'rd', 'ld', 'mxd', 'mnd', 'lw',
  #                        'ew', 'ed', 'bm', 'ba')){
  #     # we don't want late/early wood width or density
  #     next
  #   }else{
  #     stop(paste0(file))
  #   }
  # }
  # For Mexico, we directly used the API, to easily avoid tropical tree species
  # see this url: https://www.ncei.noaa.gov/access/paleo-search/study/search.json?dataPublisher=NOAA&minLat=25&locations=Continent%3ENorth%20America%3EMexico&cvWhats=physical%20property%3Ewidth%3Etotal%20ring%20width&headersOnly=true
}

# ---------------------- #
# Pre-process ITRDB data #
# ---------------------- #
years <- c(1896:2025)
min_nyears <- 50 # minimum timeseries length
# US datasets
ringwidth_series <- data.frame()
manualchecks <- c()
match_idxs_all <- c()
raw_dir <- file.path(wd, "data", "itrdb", "usa")
rwl_files <- list.files(raw_dir, pattern = '.rwl') 
for(f in rwl_files){
  print(f)
  state <- sub("(\\D*).*", "\\1", f)
  if(!(state %in% states)){next}
  
  dataset <- sub("\\.[^.]*$", "", f)
  # empty datasets, or other issues
  # not critical because they don't span our perido of interest (I think?)
  if(dataset %in% c('or042x', 'ut528', 'ca726', 'az580', 'ca581')){next}
  
  if(dataset %in% c('mt139', 'co660', 'nm549', 'nm555', 'nm560', 'nm587', 
                    'co083', 'co566', 'co661', 'co673', 'az512', 'ca721',
                    'co570', 'co572', 'co586', 'co657', 'co664', 'co676', 
                    'co677', 'co692', 'wy038', 'mt115', 'mt116',
                    'ca717', 'co653', 'co654', 'mt149', 'mt157', 'nm606', 
                    'wy028', 'ca603', 'ca605', 'ca606', 'ca643', 'ca670', 
                    'co679', 'co693', 'id020', 'mt117', 'mt121', 'mt124', 
                    'nm583', 'nv523', 'nv524', 'nv525', 'nv526', 'nv527', 
                    'or047', 'or096', 'or097', 'ut525', 'ut533', 'ut553',
                    'wa135', 'wy064', 'nm583')){ 
    # some lower/upper case thing... but still look like different trees
    # or other different reasons
    manualchecks <- c(manualchecks, f)
    next
  }
  
  rwdat <- dplR::read.rwl(file.path(raw_dir, f), format = 'tucson', verbose = FALSE)
  
  if(dataset == 'wy028'){
    colnames(rwdat) <- gsub(" ", "", colnames(rwdat))
  }
  
  idrows <- which(rownames(rwdat) %in% as.character(years))
  
  if(length(idrows) < min_nyears){
    next
  }
  
  if(any(grepl("-", colnames(rwdat)))){ # these datasets are not handled by autoread.ids
    manualchecks <- c(manualchecks, f)
    next
  }else{
    
    if(dataset %in% c('ca534', 'ca535', 'co569', 'ut524', 'az036', 'co578w',
                      'az573', 'co658')){
      idxs <- dplR::read.ids(rwdat, fix.typos = FALSE, stc = c(3, 2, 1))
    }else{
      idxs <- dplR::autoread.ids(rwdat, fix.typos = FALSE)
    }
    idxs$original <- colnames(rwdat)[1:nrow(idxs)] 
  }
  
  if(max(idxs$core) > 3){
    message(' - Dataset: ', dataset)
    message('Avg. #cores per tree: ', mean(table(idxs$tree)))
    cat('Trees IDs:', sample(colnames(rwdat)[1:nrow(idxs)],5), '\n')
    manualchecks <- c(manualchecks, f)
    next
  }
  match_idxs <-
    data.frame(
      dataset = dataset,
      original_id = colnames(rwdat),
      tree = idxs$tree, 
      core = idxs$core
    )
  match_idxs_all <- rbind(match_idxs_all, match_idxs)
  
  colnames(rwdat)[1:nrow(idxs)] <- paste0('tree', idxs$tree, '_core', idxs$core)
  
  if(mean(table(idxs$tree)) >= 3){
    message(' - Dataset: ', dataset)
    message('#unique trees: ', length(unique(idxs$tree)))
    message('Avg. #cores per tree: ', mean(table(idxs$tree)))
  }else if(mean(table(idxs$tree)) > 2){
    # cat('#unique trees:', length(unique(idxs$tree)), '\n')
    # cat('Avg. #cores per tree:', mean(table(idxs$tree)), '\n')
  }
  
  keep <- colSums(!is.na(rwdat[idrows, ])) >= 50
  rwdat <- rwdat[, keep]
  idrows <- which(rownames(rwdat) %in% as.character(years))
  
  # tryCatch({m   <-  dplR::read.rwl(file.path(raw_dir, f), format = 'tucson', verbose = FALSE)}, warning=function(w) print(f))
  idrows <- which(rownames(rwdat) %in% as.character(years))
  if(length(idrows) >= min_nyears){
    rwdat <- rwdat[idrows,]
    rwdat$year <- rownames(rwdat)
    rwdat_long <- reshape(rwdat, varying = names(rwdat)[names(rwdat) != "year"],
                          v.names = "rw_mm", timevar = "full_id",
                          times = names(rwdat)[names(rwdat) != "year"],
                          direction = "long")

    rwdat_long$tree_id <- sub(".*tree(\\d+)_core.*", "\\1", rwdat_long$full_id)
    rwdat_long$core_id <- sub(".*core(\\d+).*", "\\1", rwdat_long$full_id)

    rwdat_long$dataset <- dataset
    rwdat_long <- rwdat_long[, c('dataset', "year", 'full_id', 'tree_id', "core_id", "rw_mm")]

    # datasets with obvious errors in the end-of-series marker (should be -9999)
    if(dataset %in% c('ca684', 'ca681', 'ca682', 'ca683', 'ca685', 'ca686', 'ca598',
                      'nv524', 'nv525', 'nv526', 'nv523', 'nv527', 'ca600', 'ca601')){
      rwdat_long$rw_mm <- rwdat_long$rw_mm/10 # all the cores have wrong end-of-series markers
    }
    if(dataset %in% c('or090')){
      # coresid <- c('LCF16A', 'LCF30B', 'LCF37B') # three cores obviously have wrong end-of-series markers?
      coresid <- c('tree16_core1') # at least one
      rwdat_long[rwdat_long$full_id %in% coresid, 'rw_mm'] <- rwdat_long[rwdat_long$full_id %in% coresid, 'rw_mm']/10
    }
    if(dataset %in% c('mt165')){
      stop()
      # from D. Schook: you should start this tree core series in 1939 (12/08/25 email)
      coresid <- c('LY38A','LY38b' )
      rwdat_long[rwdat_long$core_id %in% coresid & rwdat_long$year < 1939, 'rw_mm'] <- NA
    }

    ringwidth_series <- rbind(ringwidth_series, rwdat_long)
  }else{next}
}
# write.csv(match_idxs_all, '~/match_idxs_all.csv')

# Deal with manual checks -- NOT PUSHED, need cleaning
# source('manual_checks_final.R')
# write.csv(match_idxs_all, '~/match_idxs_all.csv')
# saveRDS(ringwidth_series, file.path(wd, 'input', 'itrdb', 'ringwidth_series_usonly_from1896_update19Apr2026.rds'))




# Canada datasets
# ringwidth_series <- data.frame()
# raw_dir <- file.path(wd, "data", "itrdb", "canada")
# rwl_files <- list.files(raw_dir, pattern = '.rwl') 
# for(f in rwl_files){
#   state <- sub("(\\D*).*", "\\1", f)
#   if(!(state %in% states)){next}
#   
#   dataset <- sub("\\.[^.]*$", "", f)
#   # empty datasets, or other issues
#   # not critical because they don't span our period of interest
#   if(dataset %in% c('cana097p', 'cana362')){next}
#   
#   rwdat <- dplR::read.rwl(file.path(raw_dir, f), format = 'tucson', verbose = FALSE)
#   # tryCatch({m   <-  dplR::read.rwl(file.path(raw_dir, f), format = 'tucson', verbose = FALSE)}, warning=function(w) print(f))
#   idrows <- which(rownames(rwdat) %in% as.character(years))
#   if(length(idrows) >= min_nyears){
#     #print(f)
#     rwdat <- rwdat[idrows,]
#     rwdat$year <- rownames(rwdat)
#     rwdat_long <- reshape(rwdat, varying = names(rwdat)[names(rwdat) != "year"],
#                           v.names = "rw_mm", timevar = "core_id",
#                           times = names(rwdat)[names(rwdat) != "year"],
#                           direction = "long")
#     rwdat_long$tree_id <- sub("^(.*\\d)[^\\d]*$", "\\1", rwdat_long$core_id)
#     rwdat_long$dataset <- dataset
#     rwdat_long <- rwdat_long[, c('dataset', "year", 'tree_id', "core_id", "rw_mm")]
#     
#     ringwidth_series <- rbind(ringwidth_series, rwdat_long)
#   }else{next}
# }
# saveRDS(ringwidth_series, file.path(wd, 'input', 'itrdb', 'ringwidth_series_canadaonly_from1902.rds'))
# ringwidth_series_all <- rbind(ringwidth_series_all, ringwidth_series)
# 
# # Mexico datasets
# ringwidth_series <- data.frame()
# raw_dir <- file.path(wd, "data", "itrdb", "mexico")
# rwl_files <- list.files(raw_dir, pattern = '*\\.rwl$') 
# for(f in rwl_files){
#   state <- sub("(\\D*).*", "\\1", f)
#   if(!(state %in% states)){next}
#   
#   dataset <- sub("\\.[^.]*$", "", f)
#   
#   rwdat <- dplR::read.rwl(file.path(raw_dir, f), format = 'tucson', verbose = FALSE)
#   # tryCatch({m   <-  dplR::read.rwl(file.path(raw_dir, f), format = 'tucson', verbose = FALSE)}, warning=function(w) print(f))
#   idrows <- which(rownames(rwdat) %in% as.character(years))
#   if(length(idrows) >= min_nyears){
#     #print(f)
#     rwdat <- rwdat[idrows,]
#     rwdat$year <- rownames(rwdat)
#     rwdat_long <- reshape(rwdat, varying = names(rwdat)[names(rwdat) != "year"],
#                           v.names = "rw_mm", timevar = "core_id",
#                           times = names(rwdat)[names(rwdat) != "year"],
#                           direction = "long")
#     rwdat_long$tree_id <- sub("^(.*\\d)[^\\d]*$", "\\1", rwdat_long$core_id)
#     rwdat_long$dataset <- dataset
#     rwdat_long <- rwdat_long[, c('dataset', "year", 'tree_id', "core_id", "rw_mm")]
#     
#     # obivous problem in this core
#     if(dataset %in% c('mexi108')){
#       rwdat_long[rwdat_long$core_id == 'SPS06A' & rwdat_long$year == 1981, 'rw_mm'] <- NA
#     }
#     
#     # datasets with obvious errors in the end-of-series marker (should be -9999)
#     # if(dataset %in% c('mex134', 'mex118')){
#     #   rwdat_long$rw_mm <- rwdat_long$rw_mm/10 # all the cores have wrong end-of-series markers
#     # }
#     
#     ringwidth_series <- rbind(ringwidth_series, rwdat_long)
#   }else{next}
# }
# saveRDS(ringwidth_series, file.path(wd, 'input', 'itrdb', 'ringwidth_series_mexicoonly_from1902.rds'))
# 
# ringwidth_series_all <- rbind(readRDS(file = file.path(wd, 'input', 'itrdb', 'ringwidth_series_usonly_from1902.rds')),
#                           readRDS(file = file.path(wd, 'input', 'itrdb', 'ringwidth_series_canadaonly_from1902.rds')),
#                           readRDS(file = file.path(wd, 'input', 'itrdb', 'ringwidth_series_mexicoonly_from1902.rds')))
# saveRDS(ringwidth_series_all, file = file.path(wd, 'input', 'itrdb', 'ringwidth_series_all_from1902.rds'))



