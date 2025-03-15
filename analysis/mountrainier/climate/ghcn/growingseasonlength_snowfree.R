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
  
  k <- 1
  while((sum(!is.na(snow[c((i-k):(i-1)),'snowdepth'])) < 1 | sum(!is.na(snow[c((i+1):(i+k)),'snowdepth'])) < 1) & k < 10){
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
    }else if(prevclosest$snowdepth > 0 & afterclosest$snowdepth == 0){
      prev <- snow[c((i-2):(i-1)),]
      prevclose <- prev[prev$date %in% c(max(prev$date)-1, max(prev$date)),]
      slope <- (prevclose[2,'snowdepth']-prevclose[1,'snowdepth'])/as.numeric(prevclose[2,'date']-prevclose[1,'date'])
      snow[i, 'snowdepth'] <- max(prevclosest$snowdepth+slope*as.numeric(snow[i, 'date']-prevclosest$date), 0)
      snow[i, 'fill'] <- TRUE
    }else if(prevclosest$snowdepth == 0 & prevclosest$snowdepth > 0){
      
      afterclose <- after[after$date %in% c(min(after$date), min(after$date)+1),]
      slope <- (afterclose[2,'snowdepth']-afterclose[1,'snowdepth'])/as.numeric(afterclose[2,'date']-afterclose[1,'date'])
      snow[i, 'snowdepth'] <- max(afterclosest$snowdepth+slope*as.numeric(snow[i, 'date']-afterclosest$date), 0)
      snow[i, 'fill'] <- TRUE
      
    }else{
      
      slope <- (afterclosest$snowdepth-prevclosest$snowdepth)/as.numeric((afterclosest$date-prevclosest$date))
      snow[i, 'snowdepth'] <- max(prevclosest$snowdepth+slope*as.numeric(snow[i, 'date']-prevclosest$date), 0)
      snow[i, 'fill'] <- TRUE
      
      
    }
    
  }
  
}

pdf(file.path(wd, 'figures', 'gslength_snowfree_longmire.pdf'), width = 8, height = 20)
par(mfrow = c(5,2))
snow_gsl <- data.frame()
gsl_snowfree <- data.frame()
for(y in as.character(1931:2010)){
  
  df <- snow[snow$date >= as.Date(paste0(as.numeric(y)-1, '-10-01')) & snow$date < as.Date(paste0(y, '-10-01')),]
 
  df$snowfreecount <- NA
  snowfreecount <- zoo::rollapply(df$snowdepth == 0, 14, sum, na.rm = TRUE, align = 'left')
  df[1:length(snowfreecount), 'snowfreecount'] <- snowfreecount
  snowfreedays <- df[df$date >=  as.Date(paste0(y, '-03-01')) & df$snowfreecount == 14,]
  sos <-  min(snowfreedays$date, na.rm = TRUE)
  eos <- as.Date(paste0(y, '-09-30'))
  
  df$cumsnowfall <- NA
  wdw <- 3
  cumsnowfall <- zoo::rollapply(df$snowfall, wdw, sum, na.rm = TRUE, align = 'right')
  df[wdw:(length(cumsnowfall)+wdw-1), 'cumsnowfall'] <- cumsnowfall
  
  plot.new()
  plot.window(ylim =  c(0, max(df$snowdepth, na.rm = TRUE)), xlim =  range(df$date))
  axis.Date(1, format="%d-%m", cex.axis = 0.7)
  axis(2, las = 2, cex.axis = 0.7, tck=-0.02, labels=TRUE)
  
  abline(v=as.Date(paste0(y, '-03-01')), col="grey80", lty = 'dashed')
  rect(xleft = sos, xright = eos, ybottom = 0, ytop = max(df$snowdepth, na.rm = TRUE)*1.01, 
       col = '#f0f5df', border = NA)
  
  newsos <- NA
  df$newsnowfreecount <- NA
  if(length(df[is.na(df$snowdepth) & df$date >=  as.Date(paste0(y, '-03-01')) & df$snowfall > 0, 'snowfall']) > 0){
    
    newsnowfreecount <- zoo::rollapply(df[is.na(df$snowdepth) & df$date >=  as.Date(paste0(y, '-03-15')), 'cumsnowfall'] == 0, 
                                       7, sum, na.rm = TRUE, align = 'left', fill = NA)
    df[is.na(df$snowdepth) & df$date >=  as.Date(paste0(y, '-03-15')), 'newsnowfreecount'] <- newsnowfreecount
    newsnowfreedays <- df[df$date >=  as.Date(paste0(y, '-03-15')) & df$newsnowfreecount %in% 7,]
    newsos <- min(newsnowfreedays$date, na.rm = TRUE)
    
    rect(xleft = newsos, xright = sos-1, ybottom = 0, ytop = max(df$snowdepth, na.rm = TRUE)*1.01, 
         col = '#f2f9fb', border = NA)
    
  }
  
  title(ylab = "Snowdepth", cex.lab = 0.8)
  title(xlab = "Month", cex.lab = 0.8, line = 2.2)
  title(main = y, cex.main = 0.8)
  lines(df$snowdepth ~ df$date, col = 'white', lwd = 3)
  lines(df$snowdepth ~ df$date, col = 'grey40', lwd = 1)
  points(df[df$fill, 'snowdepth'] ~ df[df$fill, 'date'],
         col = 'darkred', cex = 0.6, pch = 20)
  
  lines(df[is.na(df$snowdepth) & df$date >=  as.Date(paste0(y, '-03-15')), 'cumsnowfall'] ~
          df[is.na(df$snowdepth) & df$date >=  as.Date(paste0(y, '-03-15')), 'date'],
        col = 'white', lwd = 3)
  lines(df[is.na(df$snowdepth) & df$date >=  as.Date(paste0(y, '-03-15')), 'cumsnowfall'] ~
          df[is.na(df$snowdepth) & df$date >=  as.Date(paste0(y, '-03-15')), 'date'],
        col = '#6bbad1', lwd = 1)
  
  points(df[is.na(df$snowdepth) & df$date >=  as.Date(paste0(y, '-03-15')) & df$snowfall > 0, 'snowfall'] ~
           df[is.na(df$snowdepth) & df$date >=  as.Date(paste0(y, '-03-15')) & df$snowfall > 0, 'date'],
         col = 'darkblue', cex = 0.6, pch = 8)
  
  text(x = (eos-sos)/2+sos , y = max(df$snowdepth, na.rm = TRUE)/2, labels = paste0(eos-sos, ' days'), col = '#9ab739')
  
  if(length(df[is.na(df$snowdepth) & df$date >=  as.Date(paste0(y, '-03-15')) & df$snowfall > 0, 'snowfall']) > 0){
    
    text(x = (eos-newsos)/2+newsos-10 , y = max(df$snowdepth, na.rm = TRUE)/1.5, labels = paste0(eos-newsos, ' days'), col = '#6bbad1')
    
  }
  
  snow_gsl <- rbind(snow_gsl, df)
  
  gsl_year <- data.frame(year = y, sos, eos, newsos, gsl = as.numeric(eos-sos), newgsl = as.numeric(eos-newsos))
  gsl_snowfree <- rbind(gsl_snowfree, gsl_year)
  
}
dev.off()

gsl_snowfree$gsl_merge <- apply(gsl_snowfree[c('gsl', 'newgsl')],1,max,na.rm=TRUE)
snowfreegsl_longmire <- gsl_snowfree[c('year', 'gsl_merge')]
names(snowfreegsl_longmire) <- c('year', 'gsl (days)')

write.csv(snowfreegsl_longmire, file = file.path(wd, 'output/climate/snowfreegsl_longmire.csv'))

# Load temperature
source(file.path(wd, 'climate', 'ghcn', 'gdd_insnowfreegsl.R'))

snow_gsl_temp <- merge(snow_gsl, temp[c('year', 'date', 'tmean')])


pdf(file.path(wd, 'figures', 'gslength_snowfree_temp5_longmire.pdf'), width = 8, height = 20)
par(mfrow = c(5,2))
gsl_snowfree <- data.frame()
for(y in as.character(1931:2010)){
  
  df <- snow_gsl_temp[snow_gsl_temp$date >= as.Date(paste0(as.numeric(y)-1, '-10-01')) & snow_gsl_temp$date < as.Date(paste0(y, '-10-01')),]
  
  df$snowfreecount <- NA
  snowfreecount <- zoo::rollapply(df$snowdepth == 0, 14, sum, na.rm = TRUE, align = 'left')
  df[1:length(snowfreecount), 'snowfreecount'] <- snowfreecount
  snowfreedays <- df[df$date >=  as.Date(paste0(y, '-03-01')) & df$snowfreecount == 14,]
  
  df$tempcount <- zoo::rollapply(df$tmean >= 5, 7, sum, na.rm = TRUE, align = 'right', fill = NA)
  temp5days <- df[df$date >=  as.Date(paste0(y, '-03-01')) & df$tempcount >= 3,]
  
  sosdays <- df[df$date >=  as.Date(paste0(y, '-03-01')) & df$snowfreecount == 14 & df$tempcount >= 3,]
  
  sos <-  min(sosdays$date, na.rm = TRUE)
  eos <- as.Date(paste0(y, '-09-30'))
  
  df$cumsnowfall <- NA
  wdw <- 3
  cumsnowfall <- zoo::rollapply(df$snowfall, wdw, sum, na.rm = TRUE, align = 'right')
  df[wdw:(length(cumsnowfall)+wdw-1), 'cumsnowfall'] <- cumsnowfall
  
  plot.new()
  plot.window(ylim =  c(0, max(df$snowdepth, na.rm = TRUE)), xlim =  range(df$date))
  axis.Date(1, format="%d-%m", cex.axis = 0.7)
  axis(2, las = 2, cex.axis = 0.7, tck=-0.02, labels=TRUE)
  
  abline(v=as.Date(paste0(y, '-03-01')), col="grey80", lty = 'dashed')
  rect(xleft = sos, xright = eos, ybottom = 0, ytop = max(df$snowdepth, na.rm = TRUE)*1.01, 
       col = '#f0f5df', border = NA)
  
  newsos <- NA
  df$newsnowfreecount <- NA
  if(length(df[is.na(df$snowdepth) & df$date >=  as.Date(paste0(y, '-03-01')) & df$snowfall > 0, 'snowfall']) > 0){
    
    newsnowfreecount <- zoo::rollapply(df[is.na(df$snowdepth) & df$date >=  as.Date(paste0(y, '-03-01')), 'cumsnowfall'] == 0, 
                                       14, sum, na.rm = TRUE, align = 'left', fill = NA)
    df[is.na(df$snowdepth) & df$date >=  as.Date(paste0(y, '-03-01')), 'newsnowfreecount'] <- newsnowfreecount
    newsosdays <- df[df$date >=  as.Date(paste0(y, '-03-01')) & df$newsnowfreecount %in% 14 & df$tempcount >= 3,]
    
  }
  
  lines(df[is.na(df$snowdepth) & df$date >=  as.Date(paste0(y, '-03-01')), 'cumsnowfall'] ~
          df[is.na(df$snowdepth) & df$date >=  as.Date(paste0(y, '-03-01')), 'date'],
        col = 'white', lwd = 3)
  lines(df[is.na(df$snowdepth) & df$date >=  as.Date(paste0(y, '-03-01')), 'cumsnowfall'] ~
          df[is.na(df$snowdepth) & df$date >=  as.Date(paste0(y, '-03-01')), 'date'],
        col = '#6bbad1', lwd = 1)
  
  title(ylab = "Snowdepth", cex.lab = 0.8)
  title(xlab = "Month", cex.lab = 0.8, line = 2.2)
  title(main = y, cex.main = 0.8)
  lines(df$snowdepth ~ df$date, col = 'white', lwd = 3)
  lines(df$snowdepth ~ df$date, col = 'grey40', lwd = 1)
  
  abline(v= min(temp5days$date), col="orange", lty = 'dashed')
  
 
  
  text(x = (eos-sos)/2+sos , y = max(df$snowdepth, na.rm = TRUE)/2, labels = paste0(eos-sos, ' days'), col = '#9ab739')
  
  gsl_year <- data.frame(year = y, sos, eos, newsos, gsl = as.numeric(eos-sos), newgsl = as.numeric(eos-newsos))
  gsl_snowfree <- rbind(gsl_snowfree, gsl_year)
  
}
dev.off()

snowfreegsl_longmire <- gsl_snowfree[c('year', 'gsl')]
names(snowfreegsl_longmire) <- c('year', 'gsl (days)')

write.csv(snowfreegsl_longmire, file = file.path(wd, 'output/climate/snowfreegsl_addtempcriteria_longmire.csv'))




# PARADISE #
#----------#

# 1916 - today
USC00456898 <- read.csv(file.path(wd, 'data', 'climate', 'USC00456898.csv'), header = FALSE,
                        col.names = c('id', 'date', 'var', 'value', 'x1', 'x2', 'x3', 'x4'))

climate_data <- USC00456898
climate_data <- climate_data[c('id', 'date', 'var', 'value')]
climate_data$date <- as.Date(as.character(climate_data$date),format="%Y%m%d")
climate_data$year <- as.numeric(format(climate_data$date, "%Y"))

snowdepth <- climate_data[climate_data$var == 'SNWD',]
snowdepth <- snowdepth[order(snowdepth$date), ]
snowdepth$nobs_year <- ave(seq_len(nrow(snowdepth)), snowdepth$year, FUN = length) 
check <- unique(snowdepth[c('year', 'nobs_year')])
plot(check$nobs_year ~ check$year)

# Missing values
snowdepth <- snowdepth[snowdepth$year >= 1960,] # before 1960, data are not great
missingdates <- seq(as.Date('1960-01-01'), as.Date('2024-12-31'), 1) # missing days
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
snowfall <- snowfall[snowfall$year >= 1960,]
missingdates <- seq(as.Date('1960-01-01'), as.Date('2024-12-31'), 1) # missing days
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
  
  k <- 1
  while((sum(!is.na(snow[c((i-k):(i-1)),'snowdepth'])) < 1 | sum(!is.na(snow[c((i+1):(i+k)),'snowdepth'])) < 1) & k < 30){
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
    }else if(prevclosest$snowdepth > 0 & afterclosest$snowdepth == 0){
      prev <- snow[c((i-2):(i-1)),]
      prevclose <- prev[prev$date %in% c(max(prev$date)-1, max(prev$date)),]
      slope <- (prevclose[2,'snowdepth']-prevclose[1,'snowdepth'])/as.numeric(prevclose[2,'date']-prevclose[1,'date'])
      snow[i, 'snowdepth'] <- max(prevclosest$snowdepth+slope*as.numeric(snow[i, 'date']-prevclosest$date), 0)
      snow[i, 'fill'] <- TRUE
    }else if(prevclosest$snowdepth == 0 & prevclosest$snowdepth > 0){
      
      afterclose <- after[after$date %in% c(min(after$date), min(after$date)+1),]
      slope <- (afterclose[2,'snowdepth']-afterclose[1,'snowdepth'])/as.numeric(afterclose[2,'date']-afterclose[1,'date'])
      snow[i, 'snowdepth'] <- max(afterclosest$snowdepth+slope*as.numeric(snow[i, 'date']-afterclosest$date), 0)
      snow[i, 'fill'] <- TRUE
      
    }else{
      
      slope <- (afterclosest$snowdepth-prevclosest$snowdepth)/as.numeric((afterclosest$date-prevclosest$date))
      snow[i, 'snowdepth'] <- max(prevclosest$snowdepth+slope*as.numeric(snow[i, 'date']-prevclosest$date), 0)
      snow[i, 'fill'] <- TRUE
      
      
    }
    
  }
  
}

pdf(file.path(wd, 'figures', 'gslength_snowfree_paradise.pdf'), width = 8, height = 20)
par(mfrow = c(5,2))
snow_gsl <- data.frame()
gsl_snowfree <- data.frame()
for(y in as.character(1961:2010)){
  
  df <- snow[snow$date >= as.Date(paste0(as.numeric(y)-1, '-10-01')) & snow$date < as.Date(paste0(y, '-10-01')),]
  
  df$snowfreecount <- NA
  snowfreecount <- zoo::rollapply(df$snowdepth == 0, 14, sum, na.rm = TRUE, align = 'left')
  df[1:length(snowfreecount), 'snowfreecount'] <- snowfreecount
  snowfreedays <- df[df$date >=  as.Date(paste0(y, '-03-01')) & df$snowfreecount == 14,]
  sos <-  min(snowfreedays$date, na.rm = TRUE)
  eos <- as.Date(paste0(y, '-09-30'))
  
  df$cumsnowfall <- NA
  wdw <- 3
  cumsnowfall <- zoo::rollapply(df$snowfall, wdw, sum, na.rm = TRUE, align = 'right')
  df[wdw:(length(cumsnowfall)+wdw-1), 'cumsnowfall'] <- cumsnowfall
  
  plot.new()
  plot.window(ylim =  c(0, max(df$snowdepth, na.rm = TRUE)), xlim =  range(df$date))
  axis.Date(1, format="%d-%m", cex.axis = 0.7)
  axis(2, las = 2, cex.axis = 0.7, tck=-0.02, labels=TRUE)
  
  abline(v=as.Date(paste0(y, '-03-01')), col="grey80", lty = 'dashed')
  rect(xleft = sos, xright = eos, ybottom = 0, ytop = max(df$snowdepth, na.rm = TRUE)*1.01, 
       col = '#f0f5df', border = NA)
  
  newsos <- NA
  df$newsnowfreecount <- NA
  if(length(df[is.na(df$snowdepth) & df$date >=  as.Date(paste0(y, '-03-01')) & df$snowfall > 0, 'snowfall']) > 0){
    
    newsnowfreecount <- zoo::rollapply(df[is.na(df$snowdepth) & df$date >=  as.Date(paste0(y, '-03-15')), 'cumsnowfall'] == 0, 
                                       7, sum, na.rm = TRUE, align = 'left', fill = NA)
    df[is.na(df$snowdepth) & df$date >=  as.Date(paste0(y, '-03-15')), 'newsnowfreecount'] <- newsnowfreecount
    newsnowfreedays <- df[df$date >=  as.Date(paste0(y, '-03-15')) & df$newsnowfreecount %in% 7,]
    newsos <- min(newsnowfreedays$date, na.rm = TRUE)
    
    rect(xleft = newsos, xright = sos-1, ybottom = 0, ytop = max(df$snowdepth, na.rm = TRUE)*1.01, 
         col = '#f2f9fb', border = NA)
    
  }
  
  title(ylab = "Snowdepth", cex.lab = 0.8)
  title(xlab = "Month", cex.lab = 0.8, line = 2.2)
  title(main = y, cex.main = 0.8)
  lines(df$snowdepth ~ df$date, col = 'white', lwd = 3)
  lines(df$snowdepth ~ df$date, col = 'grey40', lwd = 1)
  points(df[df$fill, 'snowdepth'] ~ df[df$fill, 'date'],
         col = 'darkred', cex = 0.6, pch = 20)
  
  lines(df[is.na(df$snowdepth) & df$date >=  as.Date(paste0(y, '-03-15')), 'cumsnowfall'] ~
          df[is.na(df$snowdepth) & df$date >=  as.Date(paste0(y, '-03-15')), 'date'],
        col = 'white', lwd = 3)
  lines(df[is.na(df$snowdepth) & df$date >=  as.Date(paste0(y, '-03-15')), 'cumsnowfall'] ~
          df[is.na(df$snowdepth) & df$date >=  as.Date(paste0(y, '-03-15')), 'date'],
        col = '#6bbad1', lwd = 1)
  
  points(df[is.na(df$snowdepth) & df$date >=  as.Date(paste0(y, '-03-15')) & df$snowfall > 0, 'snowfall'] ~
           df[is.na(df$snowdepth) & df$date >=  as.Date(paste0(y, '-03-15')) & df$snowfall > 0, 'date'],
         col = 'darkblue', cex = 0.6, pch = 8)
  
  text(x = (eos-sos)/2+sos , y = max(df$snowdepth, na.rm = TRUE)/2, labels = paste0(eos-sos, ' days'), col = '#9ab739')
  
  if(length(df[is.na(df$snowdepth) & df$date >=  as.Date(paste0(y, '-03-15')) & df$snowfall > 0, 'snowfall']) > 0){
    
    text(x = (eos-newsos)/2+newsos-10 , y = max(df$snowdepth, na.rm = TRUE)/1.5, labels = paste0(eos-newsos, ' days'), col = '#6bbad1')
    
  }
  
  snow_gsl <- rbind(snow_gsl, df)
  
  gsl_year <- data.frame(year = y, sos, eos, newsos, gsl = as.numeric(eos-sos), newgsl = as.numeric(eos-newsos))
  gsl_snowfree <- rbind(gsl_snowfree, gsl_year)
  
}
dev.off()

plot.new()
plot.window(ylim =  c(0, 150), xlim =  range(as.numeric(gsl_snowfree$year)))
axis(1, cex.axis = 0.7)
axis(2, las = 2, cex.axis = 0.7, tck=-0.02, labels=TRUE)
lines(gsl_snowfree$gsl ~ as.numeric(gsl_snowfree$year), lwd = 0.5)
points(gsl_snowfree$gsl ~ as.numeric(gsl_snowfree$year), cex = 0.5, pch = 20)
abline(v=1974, col="grey80", lty = 'dashed')

gsl_snowfree$gsl_merge <- apply(gsl_snowfree[c('gsl', 'newgsl')],1,max,na.rm=TRUE)
snowfreegsl_paradise <- gsl_snowfree[c('year', 'gsl_merge')]
names(snowfreegsl_paradise) <- c('year', 'gsl (days)')

write.csv(snowfreegsl_paradise, file = file.path(wd, 'output/climate/snowfreegsl_paradise.csv'))
