rm(list = ls())
wd <- '/home/victor/projects/climategrowthshifts/analysis/pnw'
options(timeout=120)

# ----------------------------------- #
# Retrieve metadata from ITRDB stands #
# ----------------------------------- #
download_until_success <- function(url, destfile, ..., maxcount = 5) {
  count <- 0
  repeat{
    Sys.sleep(0.5)
    try(download.file(url, destfile, ...))
    count <- count + 1
    if (file.exists(destfile) || count >= maxcount){break}
  }
}
get_val <- function(text,p,numeric=FALSE) {
  rows <- grep(text[,1],pattern=p)
  if (length(rows)==0) {
    return(NA) #no such pattern found
  }
  subtext <- text[grep(text[,1],pattern=p),]
  split_line <- unlist(strsplit(subtext, ': '))
  if (length(split_line)<2) {
    return(NA) #format error: nothing after ": "
  }
  value <- split_line[2]
  if (numeric) {
    value <- as.numeric(value)
  }
  return(value)
}

server <- "https://www.ncei.noaa.gov/pub/data/paleo/treering/measurements/northamerica/usa"
filenames <- unlist(xml2::as_list(httr::content(httr::GET(server))))
filenames <- filenames[grepl(filenames, pattern = '*.txt$')]
datasets <- sub("-.*", "", filenames)
ringwidth_series <- readRDS(file.path(wd, 'input', 'itrdb', 'ringwidth_series_usonly.rds'))
datasets_filtered <- unique(ringwidth_series$dataset)
data_summary <- data.frame()
for(f in filenames[which(datasets %in% datasets_filtered)]){
  
  download_until_success(file.path(server, f), file.path(tempdir(), f))
  text <- read.delim(file.path(tempdir(), f), nrows = 350)
  file.remove(file.path(tempdir(), f))
  
  species_name <- get_val(text,'Species_Name')
  species_code <- get_val(text,'Tree_Species_Code')
  
  north_lat <- get_val(text,'Northernmost_Latitude',numeric=TRUE)
  south_lat <- get_val(text,'Southernmost_Latitude',numeric=TRUE)
  east_lon <- get_val(text,'Easternmost_Longitude',numeric=TRUE)
  west_lon <- get_val(text,'Westernmost_Longitude',numeric=TRUE)
  altitude <- get_val(text,'Elevation',numeric=TRUE)
  
  state <- sub("(\\D*).*", "\\1", f)
  dataset <- sub("\\.[^.]*$", "", f)
  
  first_year <- get_val(text,'First_Year',numeric=TRUE)
  last_year <- get_val(text,'Last_Year',numeric=TRUE)
  
  if(is.na(first_year)){
    first_year <- get_val(text,'Earliest_Year',numeric=TRUE)
  }
  if(is.na(last_year)){
    last_year <- get_val(text,'Most_Recent_Year',numeric=TRUE)
  }
  
  time_unit <- get_val(text,'Time_Unit',numeric=FALSE)
  
  subtext <- text[grep(text[,1],pattern='millimeter'),][1]
  split_line <- unlist(strsplit(subtext, ','))
  ringwidth_unit <- split_line[grep(split_line,pattern='millimeter')]
  ringwidth_unit_notes <- sub(".*are ", "", subtext)
  if(length(ringwidth_unit)==0){ringwidth_unit <- 'unknown'}
  if(length(ringwidth_unit_notes)==0){ringwidth_unit_notes <- 'none'}
  
  d_summ<- data.frame(state = state, dataset = dataset, species_name, species_code, first_year, last_year, ringwidth_unit, time_unit, 
                      north_lat, south_lat, east_lon, west_lon, altitude)
  data_summary <- rbind(data_summary, d_summ)
}
saveRDS(data_summary, file = file.path('input', 'itrdb', 'datasets_summary_usonly.rds'))