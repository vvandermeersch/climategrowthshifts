

# Snowfree growing season
snowfree_threshold <- 7
growingseason <- data.frame()
for(year in unique(lubridate::year(datdf$date))){
  df <- datdf[datdf$var == "swe" & lubridate::year(datdf$date) %in% year,]
  df$snowfreecount <- NA
  for(p in unique(df$ID)){
    snowfreecount <- zoo::rollapply(df[df$ID == p, "value"]==0, snowfree_threshold, sum, na.rm = TRUE, align = 'left', fill = NA)
    df[df$ID == p, 'snowfreecount'] <- snowfreecount
    snowfreedays <- df[df$ID == p & df$date >=  as.Date(paste0(year, '-03-01')) & df$snowfreecount == snowfree_threshold,]
    sos <-  min(snowfreedays$date, na.rm = TRUE)
    eos <- as.Date(paste0(year, '-09-30'))
    growingseason <- rbind(growingseason, data.frame(ID = p, year = year, sos = sos, eos = eos))
  }
}

limits <- data.frame(minsos = as.Date(paste0(unique(lubridate::year(datdf$date)), '-03-01')), 
                     maxeos = as.Date(paste0(unique(lubridate::year(datdf$date)), '-09-30')))
ggplot() +
  geom_segment(data = growingseason, aes(x = sos, y = ID, xend = eos, yend = ID, color = ID), linewidth = 2) +
  geom_vline(data = limits, aes(xintercept = minsos), linewidth = 0.8, color = "white") +
  geom_vline(data = limits, aes(xintercept = minsos), linetype = "dashed", linewidth = 0.3) +
  geom_vline(data = limits, aes(xintercept = maxeos), linewidth = 0.8, color = "white") +
  geom_vline(data = limits, aes(xintercept = maxeos), linetype = "dashed", linewidth = 0.3) +
  scale_color_gradientn(colors = kippenberger) + 
  theme_bw() +
  scale_x_date(expand = c(0,0)) +
  theme(legend.position = 'none',
        panel.grid = element_blank()) +
  labs(x = "DOY", y = "")


growingseason$gslength <- growingseason$eos - growingseason$sos

growingseason <-
  growingseason %>%
  group_by(ID) %>%
  mutate(meanlength = mean(gslength)) %>%
  ungroup()

ggplot() +
  geom_line(data = growingseason, aes(x = as.Date(paste(year, 01, 01, sep = "-")), y = gslength-meanlength, color = ID, group = ID), linewidth = 0.5) +
  scale_color_gradientn(colors = kippenberger) + 
  theme_bw() +
  scale_x_date(expand = c(0,0), date_labels = "%Y", breaks =  as.Date(paste(seq(1980,2020,5), 01, 01, sep = "-"))) +
  theme(legend.position = 'none',
        panel.grid = element_blank()) +
  labs(x = "Year", y = "GSL (anomaly, days)")

# GDD in growing season
tlower <- 5
tupper <- 35
growingdegrees <- data.frame()
for(year in unique(lubridate::year(datdf$date))){
  df <- datdf[datdf$var == "tair" & lubridate::year(datdf$date) %in% year,]
  df$gdd <- NA
  for(p in unique(df$ID)){
    
    gs <- growingseason[growingseason$ID == p & growingseason$year == year,]
    
    tair <- df[df$ID == p & df$date %in% seq(gs$sos, gs$eos, 1), "value"]
    tair <- ifelse(tair< tlower, tlower, ifelse(tair > tupper, tupper, tair)) # apply lower and upper bounds
    gdd <- cumsum(tair-tlower)
    growingdegrees <- rbind(growingdegrees, data.frame(ID = p, year = year, gdd_ings = max(gdd)))
  }
}

growingdegrees <-
  growingdegrees %>%
  group_by(ID) %>%
  mutate(meangdd = mean(gdd_ings)) %>%
  ungroup()

ggplot() +
  geom_line(data = growingdegrees, aes(x = as.Date(paste(year, 01, 01, sep = "-")), y = gdd_ings-meangdd, color = ID, group = ID), linewidth = 0.5) +
  scale_color_gradientn(colors = kippenberger) + 
  theme_bw() +
  scale_x_date(expand = c(0,0), date_labels = "%Y", breaks =  as.Date(paste(seq(1980,2020,5), 01, 01, sep = "-"))) +
  theme(legend.position = 'none',
        panel.grid = element_blank()) +
  labs(x = "Year", y = "GDD in GS (anomaly, °C)")

# Soil moisture in growing season
soilmoisture <- data.frame()
for(year in unique(lubridate::year(datdf$date))){
  
  df010 <- datdf[datdf$var == "soilmoist_010" & lubridate::year(datdf$date) %in% year,]
  df1040 <- datdf[datdf$var == "soilmoist_1040" & lubridate::year(datdf$date) %in% year,]
  df40100 <- datdf[datdf$var == "soilmoist_40100" & lubridate::year(datdf$date) %in% year,]
  df100200 <- datdf[datdf$var == "soilmoist_100200" & lubridate::year(datdf$date) %in% year,]
  
  for(p in unique(df$ID)){
    
    gs <- growingseason[growingseason$ID == p & growingseason$year == year,]
    
    sm010 <- mean(df010[df010$ID == p & df010$date %in% seq(gs$sos, gs$eos, 1), "value"])
    sm1040 <- mean(df1040[df1040$ID == p & df1040$date %in% seq(gs$sos, gs$eos, 1), "value"])
    sm40100 <- mean(df40100[df40100$ID == p & df40100$date %in% seq(gs$sos, gs$eos, 1), "value"])
    sm100200 <- mean(df100200[df100200$ID == p & df100200$date %in% seq(gs$sos, gs$eos, 1), "value"])
    sm <- weighted.mean(c(sm010, sm1040, sm40100, sm100200), c(10,30,60,100))

    soilmoisture <- rbind(soilmoisture, data.frame(ID = p, year = year, soilmoist_ings = sm))
  }
}

soilmoisture <-
  soilmoisture %>%
  group_by(ID) %>%
  mutate(meansm = mean(soilmoist_ings)) %>%
  ungroup()

ggplot() +
  geom_line(data = soilmoisture, aes(x = as.Date(paste(year, 01, 01, sep = "-")), y = soilmoist_ings-meansm, color = ID, group = ID), linewidth = 0.5) +
  scale_color_gradientn(colors = kippenberger) + 
  theme_bw() +
  scale_x_date(expand = c(0,0), date_labels = "%Y", breaks =  as.Date(paste(seq(1980,2020,5), 01, 01, sep = "-"))) +
  theme(legend.position = 'none',
        panel.grid = element_blank()) +
  labs(x = "Year", y = "Soil moisture in GS (anomaly, m3.m-3)")

# Soil moisture in JJA
soilmoisturejja <- data.frame()
for(year in unique(lubridate::year(datdf$date))){
  
  df010 <- datdf[datdf$var == "soilmoist_010" & lubridate::year(datdf$date) %in% year,]
  df1040 <- datdf[datdf$var == "soilmoist_1040" & lubridate::year(datdf$date) %in% year,]
  df40100 <- datdf[datdf$var == "soilmoist_40100" & lubridate::year(datdf$date) %in% year,]
  df100200 <- datdf[datdf$var == "soilmoist_100200" & lubridate::year(datdf$date) %in% year,]
  
  for(p in unique(df$ID)){
    
    s <- as.Date(paste0(year, '-06-01'))
    e <- as.Date(paste0(year, '-08-31'))
    
    sm010 <- mean(df010[df010$ID == p & df010$date %in% seq(s, e, 1), "value"])
    sm1040 <- mean(df1040[df1040$ID == p & df1040$date %in% seq(s, e, 1), "value"])
    sm40100 <- mean(df40100[df40100$ID == p & df40100$date %in% seq(s, e, 1), "value"])
    sm100200 <- mean(df100200[df100200$ID == p & df100200$date %in% seq(s, e, 1), "value"])
    sm <- weighted.mean(c(sm010, sm1040, sm40100, sm100200), c(10,30,60,100))
    
    soilmoisturejja <- rbind(soilmoisturejja, data.frame(ID = p, year = year, soilmoist_jja = sm))
  }
}

soilmoisturejja <-
  soilmoisturejja %>%
  group_by(ID) %>%
  mutate(meansmjja = mean(soilmoist_jja)) %>%
  ungroup()

ggplot() +
  geom_line(data = soilmoisturejja, aes(x = as.Date(paste(year, 01, 01, sep = "-")), y = soilmoist_jja-meansmjja, color = ID, group = ID), linewidth = 0.5) +
  scale_color_gradientn(colors = kippenberger) + 
  theme_bw() +
  scale_x_date(expand = c(0,0), date_labels = "%Y", breaks =  as.Date(paste(seq(1980,2020,5), 01, 01, sep = "-"))) +
  theme(legend.position = 'none',
        panel.grid = element_blank()) +
  labs(x = "Year", y = "Soil moisture in GS (anomaly, m3.m-3)")

climate_predictors <- merge(merge(merge(merge(growingseason, growingdegrees), soilmoisture), soilmoisturejja), plotsID)
saveRDS(climate_predictors, file = file.path(wd, "climate", "processed_predictors", "wldas_climpredictors.rds"))
