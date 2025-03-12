rm(list=ls())
wd <- "~/projects/climategrowthshifts/analysis/mountrainier"


# LONGMIRE #
#----------#

# 1909 - 1978
USC00456894 <- read.csv(file.path(wd, 'data', 'climate', 'USC00456894.csv'), header = FALSE,
                        col.names = c('id', 'date', 'var', 'value', 'x1', 'x2', 'x3', 'x4')) 
# 1978 - today
USC00454764 <- read.csv(file.path(wd, 'data', 'climate', 'USC00454764.csv'), header = FALSE,
                        col.names = c('id', 'date', 'var', 'value', 'x1', 'x2', 'x3', 'x4'))

climate_data <- rbind(USC00456894, USC00454764)
climate_data <- climate_data[c('id', 'date', 'var', 'value')]
climate_data$date <- as.Date(as.character(climate_data$date),format="%Y%m%d")
climate_data$year <- as.numeric(format(climate_data$date, "%Y"))

snowdepth <- climate_data[climate_data$var == 'SNWD',]
snowdepth <- snowdepth[order(snowdepth$date), ]
snowdepth$nobs_year <- ave(seq_len(nrow(snowdepth)), snowdepth$year, FUN = length) 
check <- unique(snowdepth[c('year', 'nobs_year')])
plot(check$nobs_year ~ check$year)

# Missing values
snowdepth <- snowdepth[snowdepth$year >= 1930,] # before 1931, data are not great
missingdates <- seq(as.Date('1930-01-01'), as.Date('2024-12-31'), 1) # missing days
missingdates <- missingdates[!(missingdates %in% snowdepth$date) & format(missingdates, "%Y") %in% snowdepth$year]
hist(as.numeric(format(missingdates, "%m")))
missingsnowdepth <- data.frame(
  id = 'missingdata', date = missingdates,
  var = 'SNOW', value = NA, year = format(missingdates, "%Y")
)
snowdepth <- rbind(snowdepth[c('id', 'date', 'var', 'value', 'year')], missingsnowdepth)
snowdepth <- snowdepth[order(snowdepth$date), ]
sum(is.na(snowdepth[snowdepth$year > 1930, 'value']))
names(snowdepth) <- c('id', 'date', 'var', 'snowdepth', 'year')

# Load snowfall
snowfall <- climate_data[climate_data$var == 'SNOW',]
snowfall <- snowfall[order(snowfall$date), ]
snowfall$nobs_year <- ave(seq_len(nrow(snowfall)), snowfall$year, FUN = length) 
check <- unique(snowfall[c('year', 'nobs_year')])
plot(check$nobs_year ~ check$year)
# Fill missing values (=0?)
snowfall <- snowfall[snowfall$year >= 1930,]
missingdates <- seq(as.Date('1930-01-01'), as.Date('2024-12-31'), 1) # missing days
missingdates <- missingdates[!(missingdates %in% snowfall$date) & format(missingdates, "%Y") %in% snowfall$year]
missingsnowfall <- data.frame(
  id = 'missingdata', date = missingdates,
  var = 'SNOW', value = 0, year = format(missingdates, "%Y")
)
snowfall <- rbind(snowfall[c('id', 'date', 'var', 'value', 'year')], missingsnowfall)
snowfall <- snowfall[order(snowfall$date), ]
names(snowfall) <- c('id', 'date', 'var', 'snowfall', 'year')

# Complete missing snowdepth data
snow <- merge(snowdepth[c('date', 'year', 'snowdepth')], snowfall[c('date', 'year', 'snowfall')])
snow$fill <- FALSE
ids_missing <- which(is.na(snow[, 'snowdepth']) & snow[, 'year'] > 1930)
length(ids_missing)
for(i in ids_missing){
  print(i)
  
  k <- 1
  while((sum(!is.na(snow[c((i-k):(i-1)),'snowdepth'])) < 1 | sum(!is.na(snow[c((i+1):(i+k)),'snowdepth'])) < 1) & k < 6){
    k <- k + 1
  }
  
  if(sum(!is.na(snow[c((i-k):(i-1)),'snowdepth'])) > 0 & sum(!is.na(snow[c((i+1):(i+k)),'snowdepth'])) > 0){
    
    prev <- snow[c((i-k):(i-1)),]
    prev <- prev[!is.na(prev$snowdepth),]
    prevclosest <- prev[prev$date == max(prev$date),]
    
    after <- snow[c((i+1):(i+k)),]
    after <- after[!is.na(after$snowdepth),]
    afterclosest <- after[after$date == min(after$date),]
    
    if(prevclosest$snowdepth == afterclosest$snowdepth){
      snow[i, 'snowdepth'] <- prevclosest$snowdepth
      snow[i, 'fill'] <- TRUE
    }
    
  }
  
}
length(which(is.na(snow[, 'snowdepth']) & snow[, 'year'] > 1930))
snow[366,]
