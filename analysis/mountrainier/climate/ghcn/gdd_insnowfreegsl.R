


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

tmax <- climate_data[climate_data$var == 'TMAX', c('year', 'date', 'value')]
names(tmax) <- c('year','date', 'tmax')
tmin <- climate_data[climate_data$var == 'TMIN', c('year', 'date', 'value')]
names(tmin) <- c('year','date', 'tmin')
temp <- merge(tmax, tmin)
temp$tmean <- (temp$tmin+temp$tmax)/2
temp$tmean <- temp$tmean/10 # tenths of degrees C

temp <- temp[order(temp$date), ]
temp$nobs_year <- ave(seq_len(nrow(temp)), temp$year, FUN = length) 
check <- unique(temp[c('year', 'nobs_year')])
plot(check$nobs_year ~ check$year)

temp <- temp[temp$year >= 1930,]
missingdates <- seq(as.Date('1930-01-01'), as.Date('2024-12-31'), 1) # missing days
missingdates <- missingdates[!(missingdates %in% temp$date) & format(missingdates, "%Y") %in% temp$year]
missingtemp <- data.frame(
  year = format(missingdates, "%Y"), date = missingdates,
  tmax = NA, tmin = NA, tmean = NA
)
temp <- rbind(temp[c('year', 'date', 'tmax', 'tmin', 'tmean')], missingtemp)
temp <- temp[order(temp$date), ]


pdf(file.path(wd, 'figures', 'tmean.pdf'), width = 8, height = 20)
par(mfrow = c(5,2))
for(y in as.character(1931:2010)){
  
  df <- temp[temp$date >= as.Date(paste0(as.numeric(y)-1, '-10-01')) & temp$date < as.Date(paste0(y, '-10-01')),]
  
  plot.new()
  plot.window(ylim =  c(min(df$tmean, na.rm = TRUE), max(df$tmean, na.rm = TRUE)), xlim =  range(df$date))
  axis.Date(1, format="%d-%m", cex.axis = 0.7)
  axis(2, las = 2, cex.axis = 0.7, tck=-0.02, labels=TRUE)
  
  sos <- min(gsl_snowfree[gsl_snowfree$year == y,c('sos')], gsl_snowfree[gsl_snowfree$year == y,c('newsos')], na.rm = T)
  eos <- gsl_snowfree[gsl_snowfree$year == y,c('eos')]
  rect(xleft = sos, xright = eos, ybottom = min(df$tmean, na.rm = TRUE), ytop = max(df$tmean, na.rm = TRUE), 
       col = '#f0f5df', border = NA)
  
  abline(v=as.Date(paste0(y, '-03-15')), col="grey80", lty = 'dashed')
  abline(h = 0, col="grey80", lty = 'dashed')
  
  title(ylab = "Average temp.", cex.lab = 0.8)
  title(xlab = "Month", cex.lab = 0.8, line = 2.2)
  title(main = y, cex.main = 0.8)
  lines(df$tmean ~ df$date, col = 'white', lwd = 2.7)
  lines(df$tmean ~ df$date, col = '#00468b', lwd = 0.7)
  
  countmissing <- nrow(df[is.na(df$tmean) & df$date < as.Date(paste0(y, '-10-01')) & df$date >=  sos,])
  text(x = (eos-sos)/2+sos , y = -1, labels = paste0(countmissing, ' missing days'), col = 'darkred')
  
}
dev.off()
