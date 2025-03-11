
wd <- "~/projects/climategrowthshifts/analysis/mountrainier"

# 1909 - 1978
USC00456894 <- read.csv(file.path(wd, 'data', 'climate', 'USC00456894.csv'), header = FALSE,
                        col.names = c('id', 'date', 'var', 'value', 'x1', 'x2', 'x3', 'x4')) 
# 1978 - today
USC00454764 <- read.csv(file.path(wd, 'data', 'climate', 'USC00454764.csv'), header = FALSE,
                        col.names = c('id', 'date', 'var', 'value', 'x1', 'x2', 'x3', 'x4'))

climate_data <- rbind(USC00456894, USC00454764)
climate_data <- climate_data[c('id', 'date', 'var', 'value')]
climate_data$date <- as.Date(as.character(climate_data$date),format="%Y%m%d")
climate_data$year <- format(climate_data$date, "%Y")

snowfall <- climate_data[climate_data$var == 'SNOW',]
snowfall <- snowfall[order(snowfall$date), ]
snowfall$nobs_year <- ave(seq_len(nrow(snowfall)), snowfall$year, FUN = length) 
check <- unique(snowfall[c('year', 'nobs_year')])
plot(check$nobs_year ~ check$year)

# Fill missing values (=0?)
snowfall <- snowfall[snowfall$nobs_year >= 290 | snowfall$year == 2025,]
missingdates <- seq(as.Date('1914-01-01'), as.Date('2024-12-31'), 1) # missing days
missingdates <- missingdates[!(missingdates %in% snowfall$date) & format(missingdates, "%Y") %in% snowfall$year]
missingsnowfall <- data.frame(
  id = 'missingdata', date = missingdates,
  var = 'SNOW', value = 0, year = format(missingdates, "%Y")
)
snowfall <- rbind(snowfall[c('id', 'date', 'var', 'value', 'year')], missingsnowfall)
snowfall <- snowfall[order(snowfall$date), ]

# Sum of snowfall in the growing season (sept - sept)
snowfall$growyear_id <- cumsum(format(snowfall$date, "%d-%m") == "01-09")
snowfall <- snowfall[snowfall$growyear_id > 0,]

yearly_snowfall <- aggregate(snowfall$value, by=list(growyear = snowfall$growyear_id), FUN=sum)
yearly_snowfall$year <- unique(snowfall$year)[-1]

pdf(file.path(wd, 'figures', 'yearlysnowfall_longmire.pdf'), width = 8, height = 6.5)
par(mfrow=c(1, 1) , mar = c(5,5,1,2))
plot(yearly_snowfall$x/1000 ~ yearly_snowfall$year,
     xlab = 'Year', ylab = 'Total snowfall (m)')
dev.off()

