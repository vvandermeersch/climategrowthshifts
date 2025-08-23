# Started 8 June 2025 
# by Mao
# housekeeping
rm(list=ls()) 
library(ggplot2)
library(tidyr)
library(dplyr)

setwd("C:/PhD/Project/climategrowthshifts/analysis/mountrainier/data/climate")
longmirestation <- read.csv('USC00456894.csv', header = F)
longmirestation1 <- read.csv('USC00454764.csv', header = F)
# Full longmire data
longmirestation <- rbind(longmirestation, longmirestation1)
# only need date and data
longmirestation <- longmirestation[,2:4]
colnames(longmirestation)<- c("date","type","data")
longmirestation <- longmirestation[longmirestation$type %in% c("TMAX", "TMIN"), ]
# rearrage the format
tab <- xtabs(data ~ date + type, data = longmirestation)
longmirestation1 <- as.data.frame(tab)
longmirestation1 <- reshape(longmirestation1,
                   timevar = "type",
                   idvar = "date",
                   direction = "wide")
# Clean column names
names(longmirestation1) <- sub("Freq\\.", "", names(longmirestation1))
# Adjust the real temp
longmirestation1$TMAX <- longmirestation1$TMAX/10
longmirestation1$TMIN <- longmirestation1$TMIN/10
# Adjust the date
longmirestation1$date <- as.Date(as.character(longmirestation1$date), format = "%Y%m%d")
longmirestation1$year <- as.numeric(format(longmirestation1$date, "%Y"))
longmirestation1$month <- as.numeric(format(longmirestation1$date, "%m"))
longmirestation1$doy <- as.numeric(format(longmirestation1$date, "%j"))

# There are some missing data in this dataset, and my easy solution for now is generating a full set of year and DOY, merge it with your existing data, and fill any missing days using long-term averages.
# Remember to consider the leap years!
is_leap_year <- function(y) {
  (y %% 4 == 0 & y %% 100 != 0) | (y %% 400 == 0)
}

day_counts <- table(longmirestation1$year)
years <- as.numeric(names(day_counts))
expected_days <- ifelse(is_leap_year(years), 366, 365)
missing_years <- years[day_counts < expected_days]
missing_years
# Lots of missing data actually...
years <- sort(unique(longmirestation1$year))

doy_full <- 1:366

# Create full year+doy
full_grid <- expand.grid(year = years, doy = doy_full)

# Remove 366 for non-leap years
full_grid <- full_grid[!(full_grid$doy == 366 & !is_leap_year(full_grid$year)), ]
longmirestation1_full <- merge(full_grid, longmirestation1, by = c("year", "doy"), all.x = TRUE)

# Calculate the climate average
clim_avg <- aggregate(cbind(TMAX, TMIN) ~ doy, data = longmirestation1, FUN = mean, na.rm = TRUE)

longmirestation_full <- merge(longmirestation1_full, clim_avg, by = "doy", all.x = TRUE, suffixes = c("", "_clim"))

# Fill missing data using averages, for only the missing data
longmirestation_full$TMAX[is.na(longmirestation_full$TMAX)] <- longmirestation_full$TMAX_clim[is.na(longmirestation_full$TMAX)]
longmirestation_full$TMIN[is.na(longmirestation_full$TMIN)] <- longmirestation_full$TMIN_clim[is.na(longmirestation_full$TMIN)]

# Reconstruct date
longmirestation_full$date <- as.Date(ISOdate(longmirestation_full$year, 1, 1)) + longmirestation_full$doy - 1
longmirestation_full$date <- as.Date(longmirestation_full$date, format = "%Y%m%d")
longmirestation_full$month <- as.numeric(format(longmirestation_full$date, "%m"))

# Calculate GDD, I use base temp 5 here for now
Tbase <- 5

longmirestation_full$MeanTemp <- (longmirestation_full$TMAX+longmirestation_full$TMIN) / 2
longmirestation_full$GDD_daily <- pmax(0, longmirestation_full$MeanTemp-Tbase)

longmirestation_full$year <- as.numeric(format(longmirestation_full$date, "%Y"))
longmirestation_full$month <- as.numeric(format(longmirestation_full$date, "%m"))

# Monthly GDD
GDD_monthly_st <- aggregate(GDD_daily ~ year + month, data = longmirestation_full, sum, na.rm = TRUE)

# Yearly GDD
GDD_yearly_st <- aggregate(GDD_daily ~ year, data = longmirestation_full, sum, na.rm = TRUE)
colnames(GDD_yearly_st) <- c("year","GDD_yearly_st")

###############################################################################################################ClimateNA#####################################################################################################################
runclimateNA <- FALSE
# Code to download data from climateNA
if(runclimateNA){
library("ClimateNAr")
periodList <- 1901:2024
outDir <- 'C:/PhD/Project/model/climateNA'

clm <- ClimateNAr(inputFile='C:/PhD/Project/model/rainier.csv', periodList,varList=c('Tmax01','Tmax02','Tmax03','Tmax04','Tmax05','Tmax06','Tmax07','Tmax08','Tmax09','Tmax10','Tmax11','Tmax12','Tmin01','Tmin02','Tmin03','Tmin04','Tmin05','Tmin06','Tmin07','Tmin08','Tmin09','Tmin10','Tmin11','Tmin12'),outDir)

setwd("C:/PhD/Project/model/climateNA")
files <- list.files(pattern = "\\.csv$")
all_data <- list()
for (file in files) {
  # Read the file
  df <- read.csv(file, stringsAsFactors = FALSE)
  
  # Extract year from the filename (adjust if needed)
  year <- as.numeric(gsub("\\D", "", file))  # This keeps all digits in filename (e.g., 2005)
  
  # Get Tmax and Tmin column names
  tmax_cols <- paste0("Tmax", sprintf("%02d", 1:12))
  tmin_cols <- paste0("Tmin", sprintf("%02d", 1:12))
  
  # Process each station (row)
  for (i in 1:nrow(df)) {
    station <- df[i, "id1"]
    tmax <- as.numeric(df[i, tmax_cols])
    tmin <- as.numeric(df[i, tmin_cols])
    month <- 1:12
    
    tidy_row <- data.frame(
      station = rep(station, 12),
      year = rep(year, 12),
      month = month,
      TMAX = tmax,
      TMIN = tmin,
      stringsAsFactors = FALSE
    )
    
    all_data[[length(all_data) + 1]] <- tidy_row
  }
}

final_df <- do.call(rbind, all_data)
longmire_df <- final_df[final_df$station == "longmire", ]
paradise_df <- final_df[final_df$station == "paradise", ]

write.csv(longmire_df, "longmire_climateNA.csv")
write.csv(paradise_df, "paradise_climateNA.csv")

clm <- ClimateNAr(inputFile='C:/PhD/Project/model/rainier.csv', periodList,varList= 'DD5',outDir = 'C:/PhD/Project/climategrowthshifts/analysis/mountrainier/data/climate/climateNA/single')

setwd("C:/PhD/Project/climategrowthshifts/analysis/mountrainier/data/climate/climateNA/station")
files <- list.files(pattern = "\\.csv$")

all_data <- list()

for (file in files) {
  df <- read.csv(file, stringsAsFactors = FALSE)
  
  # Extract year from filename
  year <- as.numeric(gsub("\\D", "", file))
  
  # Add year column
  df$year <- year
  
  # Reorder columns: id1, year, then the rest
  other_cols <- setdiff(names(df), c("id1", "year"))
  df <- df[, c("id1", "year", other_cols)]
  
  # Store the processed data
  all_data[[length(all_data) + 1]] <- df
}

# Combine into one data frame
combined_df <- do.call(rbind, all_data)
combined_df <- combined_df[, -c(3:6)]

longmire_DD5 <- combined_df[combined_df$id1 == "longmire", ]
paradise_DD5 <- combined_df[combined_df$id1 == "paradise", ]
write.csv(longmire_DD5,"C:/PhD/Project/climategrowthshifts/analysis/mountrainier/data/climate/longmire_DD5.csv")
write.csv(paradise_DD5,"C:/PhD/Project/climategrowthshifts/analysis/mountrainier/data/climate/paradise_DD5.csv")
}


longmire_df <- read.csv("C:/PhD/Project/climategrowthshifts/analysis/mountrainier/data/climate/longmire_climateNA.csv", header = T)
longmire_DD5 <- read.csv("C:/PhD/Project/climategrowthshifts/analysis/mountrainier/data/climate/longmire_DD5.CSV", header = T)

# Calculate GDD based on climateNA data
Tbase <- 5

longmire_df$MeanTemp <- (longmire_df$TMAX+longmire_df$TMIN) / 2
longmire_df$GDD_monthly_average <- pmax(0, longmire_df$MeanTemp-Tbase)

# This is just one number, and I have to multiply it by the total days of each month
days_in_month <- c(31, 28, 31, 30, 31, 30,
                   31, 31, 30, 31, 30, 31)
longmire_df$GDD_monthly <- NA

for (i in 1:nrow(longmire_df)) {
  yr <- longmire_df$year[i]
  mo <- longmire_df$month[i]
  
  # Adjust days for leap year if month is Feb
  if (mo == 2 && is_leap_year(yr)) {
    days <- 29
  } else {
    days <- days_in_month[mo]
  }
  
  # Calculate daily average temperature
  Tavg <- (longmire_df$TMAX[i] + longmire_df$TMIN[i]) / 2
  
  # Compute GDD for the month
  gdd <- (Tavg - Tbase) * days
  
  # Ensure GDD is not negative
  longmire_df$GDD_monthly[i] <- ifelse(gdd > 0, gdd, 0)
}

# Get yearly GDD
GDD_yearly_longmire_climateNA <- aggregate(GDD_monthly ~ year, data = longmire_df, sum, na.rm = TRUE)
colnames(GDD_yearly_longmire_climateNA) <- c("year","GDD_yearly_climateNA")

# Merge yearly GDD from weather station data and from climateNA so we can compare
longmire_GDD <- merge(GDD_yearly_st, GDD_yearly_longmire_climateNA, by = "year", all.x = TRUE)
longmire_GDD <- merge(longmire_GDD, longmire_DD5, by = "year", all.x = TRUE)
longmire_GDD <- longmire_GDD[, -c(4:5)]
#write.csv(longmire_GDD, "longmire_GDD.csv")

# Plotting

par(mfrow = c(2, 1))
longmire_long <- pivot_longer(longmire_GDD, cols = -year, names_to = "source", values_to = "GDD5")
non_missing_years <- longmire_long$year %>% unique() %>% setdiff(missing_years)
ggplot(longmire_long, aes(x = year, y = GDD5, color = source)) +
  geom_line(size = 1) +
  geom_point() +
  # Add vertical dashed lines for non-missing years
  geom_vline(data = data.frame(year = non_missing_years),
             aes(xintercept = year),
             linetype = "dashed", color = "gray50") +
  labs(x = "Year", y = "GDD5 (Longmire)") +
  theme_minimal()
# Plot GDD from weather station versus from climateNA
longmire_GDD$MY <- ifelse(longmire_GDD$year %in% missing_years, "missing", "full")

# Plot
ggplot(longmire_GDD, aes(x = GDD_yearly_st, y = DD5, color = MY)) + scale_color_manual(values = c("missing" = "red", "full" = "black")) + geom_point() + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +  
  theme_classic() + xlim(1200, 2100) + ylim(1200, 2100) +
  labs(
    x = "GDD from weather station",
    y = "GDD from climateNA",
    color = "Missing years",
    title = "Longmire"
  )
# Paradise
# housekeeping
rm(list=ls()) 

paradisestation <- read.csv('USC00456898.csv', header = F)
# only need date and data
paradisestation <- paradisestation[,2:4]
colnames(paradisestation)<- c("date","type","data")
paradisestation <- paradisestation[paradisestation$type %in% c("TMAX", "TMIN"), ]
# rearrage the format
tab <- xtabs(data ~ date + type, data = paradisestation)
paradisestation1 <- as.data.frame(tab)
paradisestation1 <- reshape(paradisestation1,
                            timevar = "type",
                            idvar = "date",
                            direction = "wide")
# Clean column names
names(paradisestation1) <- sub("Freq\\.", "", names(paradisestation1))
# Adjust the real temp
paradisestation1$TMAX <- paradisestation1$TMAX/10
paradisestation1$TMIN <- paradisestation1$TMIN/10
# Adjust the date
paradisestation1$date <- as.Date(as.character(paradisestation1$date), format = "%Y%m%d")
paradisestation1$year <- as.numeric(format(paradisestation1$date, "%Y"))
paradisestation1$month <- as.numeric(format(paradisestation1$date, "%m"))
paradisestation1$doy <- as.numeric(format(paradisestation1$date, "%j"))

#leap year
is_leap_year <- function(y) {
  (y %% 4 == 0 & y %% 100 != 0) | (y %% 400 == 0)
}

day_counts <- table(paradisestation1$year)
years <- as.numeric(names(day_counts))
expected_days <- ifelse(is_leap_year(years), 366, 365)
missing_years <- years[day_counts < expected_days]
missing_years
# Lots of missing data actually...
years <- sort(unique(paradisestation1$year))

doy_full <- 1:366

# Create full year+doy
full_grid <- expand.grid(year = years, doy = doy_full)

# Remove 366 for non-leap years
full_grid <- full_grid[!(full_grid$doy == 366 & !is_leap_year(full_grid$year)), ]
paradisestation1_full <- merge(full_grid, paradisestation1, by = c("year", "doy"), all.x = TRUE)

# Calculate the climate average
clim_avg <- aggregate(cbind(TMAX, TMIN) ~ doy, data = paradisestation1, FUN = mean, na.rm = TRUE)

paradisestation_full <- merge(paradisestation1_full, clim_avg, by = "doy", all.x = TRUE, suffixes = c("", "_clim"))

# Fill missing data using averages, for only the missing data
paradisestation_full$TMAX[is.na(paradisestation_full$TMAX)] <- paradisestation_full$TMAX_clim[is.na(paradisestation_full$TMAX)]
paradisestation_full$TMIN[is.na(paradisestation_full$TMIN)] <- paradisestation_full$TMIN_clim[is.na(paradisestation_full$TMIN)]

# Reconstruct date
paradisestation_full$date <- as.Date(ISOdate(paradisestation_full$year, 1, 1)) + paradisestation_full$doy - 1
paradisestation_full$date <- as.Date(paradisestation_full$date, format = "%Y%m%d")
paradisestation_full$month <- as.numeric(format(paradisestation_full$date, "%m"))

# Calculate GDD, I use base temp 5 here for now
Tbase <- 5

paradisestation_full$MeanTemp <- (paradisestation_full$TMAX+paradisestation_full$TMIN) / 2
paradisestation_full$GDD_daily <- pmax(0, paradisestation_full$MeanTemp-Tbase)

paradisestation_full$year <- as.numeric(format(paradisestation_full$date, "%Y"))
paradisestation_full$month <- as.numeric(format(paradisestation_full$date, "%m"))

# Monthly GDD
GDD_monthly_st <- aggregate(GDD_daily ~ year + month, data = paradisestation_full, sum, na.rm = TRUE)

# Yearly GDD
GDD_yearly_st <- aggregate(GDD_daily ~ year, data = paradisestation_full, sum, na.rm = TRUE)
colnames(GDD_yearly_st) <- c("year","GDD_yearly_st")

###############################################################################################################ClimateNA#####################################################################################################################
paradise_df <- read.csv("paradise_climateNA.csv", header = T)
paradise_DD5 <- read.csv("C:/PhD/Project/climategrowthshifts/analysis/mountrainier/data/climate/paradise_DD5.CSV", header = T)

# Calculate GDD based on climateNA data
Tbase <- 5

paradise_df$MeanTemp <- (paradise_df$TMAX+paradise_df$TMIN) / 2
paradise_df$GDD_monthly_average <- pmax(0, paradise_df$MeanTemp-Tbase)

# This is just one number, and I have to multiply it by the total days of each month
days_in_month <- c(31, 28, 31, 30, 31, 30,
                   31, 31, 30, 31, 30, 31)
paradise_df$GDD_monthly <- NA

for (i in 1:nrow(paradise_df)) {
  yr <- paradise_df$year[i]
  mo <- paradise_df$month[i]
  
  # Adjust days for leap year if month is Feb
  if (mo == 2 && is_leap_year(yr)) {
    days <- 29
  } else {
    days <- days_in_month[mo]
  }
  
  # Calculate daily average temperature
  Tavg <- (paradise_df$TMAX[i] + paradise_df$TMIN[i]) / 2
  
  # Compute GDD for the month
  gdd <- (Tavg - Tbase) * days
  
  # Ensure GDD is not negative
  paradise_df$GDD_monthly[i] <- ifelse(gdd > 0, gdd, 0)
}

# Get yearly GDD
GDD_yearly_paradise_climateNA <- aggregate(GDD_monthly ~ year, data = paradise_df, sum, na.rm = TRUE)
colnames(GDD_yearly_paradise_climateNA) <- c("year","GDD_yearly_climateNA")

# Merge yearly GDD from weather station data and from climateNA so we can compare
paradise_GDD <- merge(GDD_yearly_st, GDD_yearly_paradise_climateNA, by = "year", all.x = TRUE)
paradise_GDD <- merge(paradise_GDD, paradise_DD5, by = "year", all.x = TRUE)
paradise_GDD <- paradise_GDD[, -4]

write.csv(paradise_GDD, "paradise_GDD.csv")
#Plotting for paradise

paradise_GDD$MY <- ifelse(paradise_GDD$year %in% missing_years, "missing", "full")

# Plot
ggplot(paradise_GDD, aes(x = GDD_yearly_st, y = DD5, color = MY)) + scale_color_manual(values = c("missing" = "red", "full" = "black")) + geom_point() + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +  
  theme_classic() + xlim(500, 1200) + ylim(500, 1200) +
  labs(
    x = "GDD from weather station",
    y = "GDD from climateNA",
    color = "Missing years",
    title = "Paradise"
  )

## Precipitation
rm(list=ls()) 
setwd("C:/PhD/Project/climategrowthshifts/analysis/mountrainier/data/climate")
longmirestation <- read.csv('USC00456894.csv', header = F)
longmirestation1 <- read.csv('USC00454764.csv', header = F)
# Full longmire data
longmirestation <- rbind(longmirestation, longmirestation1)
# only need date and data
longmirestation <- longmirestation[,2:4]
colnames(longmirestation)<- c("date","type","data")
longmirestation <- longmirestation[longmirestation$type %in% "PRCP", ]
# rearrage the format
tab <- xtabs(data ~ date + type, data = longmirestation)
longmirestation1 <- as.data.frame(tab)
longmirestation1 <- reshape(longmirestation1,
                            timevar = "type",
                            idvar = "date",
                            direction = "wide")
# Clean column names
names(longmirestation1) <- sub("Freq\\.", "", names(longmirestation1))

# Adjust the date
longmirestation1$date <- as.Date(as.character(longmirestation1$date), format = "%Y%m%d")
longmirestation1$year <- as.numeric(format(longmirestation1$date, "%Y"))
longmirestation1$month <- as.numeric(format(longmirestation1$date, "%m"))
longmirestation1$doy <- as.numeric(format(longmirestation1$date, "%j"))


#leap year
is_leap_year <- function(y) {
  (y %% 4 == 0 & y %% 100 != 0) | (y %% 400 == 0)
}

day_counts <- table(longmirestation1$year)
years <- as.numeric(names(day_counts))
expected_days <- ifelse(is_leap_year(years), 366, 365)
missing_years <- years[day_counts < expected_days]
missing_years

doy_full <- 1:366

# Create full year+doy
full_grid <- expand.grid(year = years, doy = doy_full)

# Remove 366 for non-leap years
full_grid <- full_grid[!(full_grid$doy == 366 & !is_leap_year(full_grid$year)), ]
longmirestation1_full <- merge(full_grid, longmirestation1, by = c("year", "doy"), all.x = TRUE)

# Calculate the climate average
clim_avg <- aggregate(PRCP ~ doy, data = longmirestation1, FUN = mean, na.rm = TRUE)

longmirestation_full <- merge(longmirestation1_full, clim_avg, by = "doy", all.x = TRUE, suffixes = c("", "_clim"))

# Fill missing data using averages, for only the missing data
longmirestation_full$PRCP[is.na(longmirestation_full$PRCP)] <- longmirestation_full$PRCP_clim[is.na(longmirestation_full$PRCP)]

# Reconstruct date
longmirestation_full$date <- as.Date(ISOdate(longmirestation_full$year, 1, 1)) + longmirestation_full$doy - 1
longmirestation_full$date <- as.Date(longmirestation_full$date, format = "%Y%m%d")
longmirestation_full$month <- as.numeric(format(longmirestation_full$date, "%m"))
longmirestation_full$PRCP <- longmirestation_full$PRCP/10

# Get yearly precipitation from station
PRCP_yearly_longmire_station <- aggregate(PRCP ~ year, data = longmirestation_full, sum, na.rm = TRUE)
colnames(PRCP_yearly_longmire_station) <- c("year","PRCP_yearly_station")

# Get yearly precipitation from climateNA

paradise_climateNA <- read.csv("C:/PhD/Project/climategrowthshifts/analysis/mountrainier/data/climate/paradise_climateNA.csv", header = T)

PRCP_yearly_longmire_climateNA <- aggregate(PPT ~ year, data = paradise_climateNA, sum, na.rm = TRUE)

paradise_PRCP <- merge(PRCP_yearly_longmire_station, PRCP_yearly_longmire_climateNA, by = "year", all.x = TRUE)

write.csv(paradise_PRCP, "paradise_PRCP.csv")
#Plotting for paradise

paradise_PRCP$MY <- ifelse(paradise_PRCP$year %in% missing_years, "missing", "full")

# Plot
ggplot(paradise_PRCP, aes(x = PRCP_yearly_station, y = PPT, color = MY)) + scale_color_manual(values = c("missing" = "red", "full" = "black")) + geom_point() + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +  
  theme_classic() + xlim(1000, 3000) + ylim(1000, 3000) +
  labs(
    x = "PRCP from weather station",
    y = "PRCP from climateNA",
    color = "Missing years",
    title = "Paradise"
  )
setwd("C:/PhD/Project/climategrowthshifts/analysis/mountrainier/data/climate")
longmirestation <- read.csv('USC00456894.csv', header = F)
longmirestation1 <- read.csv('USC00454764.csv', header = F)
# Full longmire data
longmirestation <- rbind(longmirestation, longmirestation1)
# only need date and data
longmirestation <- longmirestation[,2:4]
colnames(longmirestation)<- c("date","type","data")
longmirestation <- longmirestation[longmirestation$type %in% "PRCP", ]
# rearrage the format
tab <- xtabs(data ~ date + type, data = longmirestation)
longmirestation1 <- as.data.frame(tab)
longmirestation1 <- reshape(longmirestation1,
                            timevar = "type",
                            idvar = "date",
                            direction = "wide")
# Clean column names
names(longmirestation1) <- sub("Freq\\.", "", names(longmirestation1))

# Adjust the date
longmirestation1$date <- as.Date(as.character(longmirestation1$date), format = "%Y%m%d")
longmirestation1$year <- as.numeric(format(longmirestation1$date, "%Y"))
longmirestation1$month <- as.numeric(format(longmirestation1$date, "%m"))
longmirestation1$doy <- as.numeric(format(longmirestation1$date, "%j"))


#leap year
is_leap_year <- function(y) {
  (y %% 4 == 0 & y %% 100 != 0) | (y %% 400 == 0)
}

day_counts <- table(longmirestation1$year)
years <- as.numeric(names(day_counts))
expected_days <- ifelse(is_leap_year(years), 366, 365)
missing_years <- years[day_counts < expected_days]
missing_years

doy_full <- 1:366

# Create full year+doy
full_grid <- expand.grid(year = years, doy = doy_full)

# Remove 366 for non-leap years
full_grid <- full_grid[!(full_grid$doy == 366 & !is_leap_year(full_grid$year)), ]
longmirestation1_full <- merge(full_grid, longmirestation1, by = c("year", "doy"), all.x = TRUE)

# Calculate the climate average
clim_avg <- aggregate(PRCP ~ doy, data = longmirestation1, FUN = mean, na.rm = TRUE)

longmirestation_full <- merge(longmirestation1_full, clim_avg, by = "doy", all.x = TRUE, suffixes = c("", "_clim"))

# Fill missing data using averages, for only the missing data
longmirestation_full$PRCP[is.na(longmirestation_full$PRCP)] <- longmirestation_full$PRCP_clim[is.na(longmirestation_full$PRCP)]

# Reconstruct date
longmirestation_full$date <- as.Date(ISOdate(longmirestation_full$year, 1, 1)) + longmirestation_full$doy - 1
longmirestation_full$date <- as.Date(longmirestation_full$date, format = "%Y%m%d")
longmirestation_full$month <- as.numeric(format(longmirestation_full$date, "%m"))
longmirestation_full$PRCP <- longmirestation_full$PRCP/10

# Get yearly precipitation from station
PRCP_yearly_longmire_station <- aggregate(PRCP ~ year, data = longmirestation_full, sum, na.rm = TRUE)
colnames(PRCP_yearly_longmire_station) <- c("year","PRCP_yearly_station")

# Get yearly precipitation from climateNA

longmire_climateNA <- read.csv("C:/PhD/Project/climategrowthshifts/analysis/mountrainier/data/climate/longmire_climateNA.csv", header = T)

PRCP_yearly_longmire_climateNA <- aggregate(PPT ~ year, data = longmire_climateNA, sum, na.rm = TRUE)

longmire_PRCP <- merge(PRCP_yearly_longmire_station, PRCP_yearly_longmire_climateNA, by = "year", all.x = TRUE)

write.csv(longmire_PRCP, "longmire_PRCP.csv")
#Plotting for paradise

longmire_PRCP$MY <- ifelse(longmire_PRCP$year %in% missing_years, "missing", "full")

# Plot
ggplot(longmire_PRCP, aes(x = PRCP_yearly_station, y = PPT, color = MY)) + scale_color_manual(values = c("missing" = "red", "full" = "black")) + geom_point() + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +  
  theme_classic() + xlim(1000, 3000) + ylim(1000, 3000) +
  labs(
    x = "PRCP from weather station",
    y = "PRCP from climateNA",
    color = "Missing years",
    title = "Longmire"
  )
# paradise
paradisestation <- read.csv('USC00456898.csv', header = F)
# only need date and data
paradisestation <- paradisestation[,2:4]
colnames(paradisestation)<- c("date","type","data")
paradisestation <- paradisestation[paradisestation$type %in% "PRCP", ]
# rearrage the format
tab <- xtabs(data ~ date + type, data = paradisestation)
paradisestation1 <- as.data.frame(tab)
paradisestation1 <- reshape(paradisestation1,
                            timevar = "type",
                            idvar = "date",
                            direction = "wide")
# Clean column names
names(paradisestation1) <- sub("Freq\\.", "", names(paradisestation1))

# Adjust the date
paradisestation1$date <- as.Date(as.character(paradisestation1$date), format = "%Y%m%d")
paradisestation1$year <- as.numeric(format(paradisestation1$date, "%Y"))
paradisestation1$month <- as.numeric(format(paradisestation1$date, "%m"))
paradisestation1$doy <- as.numeric(format(paradisestation1$date, "%j"))


#leap year
is_leap_year <- function(y) {
  (y %% 4 == 0 & y %% 100 != 0) | (y %% 400 == 0)
}

day_counts <- table(paradisestation1$year)
years <- as.numeric(names(day_counts))
expected_days <- ifelse(is_leap_year(years), 366, 365)
missing_years <- years[day_counts < expected_days]
missing_years

doy_full <- 1:366

# Create full year+doy
full_grid <- expand.grid(year = years, doy = doy_full)

# Remove 366 for non-leap years
full_grid <- full_grid[!(full_grid$doy == 366 & !is_leap_year(full_grid$year)), ]
paradisestation1_full <- merge(full_grid, paradisestation1, by = c("year", "doy"), all.x = TRUE)

# Calculate the climate average
clim_avg <- aggregate(PRCP ~ doy, data = paradisestation1, FUN = mean, na.rm = TRUE)

paradisestation1_full <- merge(paradisestation1_full, clim_avg, by = "doy", all.x = TRUE, suffixes = c("", "_clim"))

# Fill missing data using averages, for only the missing data
paradisestation1_full$PRCP[is.na(paradisestation1_full$PRCP)] <- paradisestation1_full$PRCP_clim[is.na(paradisestation1_full$PRCP)]

# Reconstruct date
paradisestation1_full$date <- as.Date(ISOdate(paradisestation1_full$year, 1, 1)) + paradisestation1_full$doy - 1
paradisestation1_full$date <- as.Date(paradisestation1_full$date, format = "%Y%m%d")
paradisestation1_full$month <- as.numeric(format(paradisestation1_full$date, "%m"))
paradisestation1_full$PRCP <- paradisestation1_full$PRCP/10

# Get yearly precipitation from station
PRCP_yearly_paradise_station <- aggregate(PRCP ~ year, data = paradisestation1_full, sum, na.rm = TRUE)
colnames(PRCP_yearly_paradise_station) <- c("year","PRCP_yearly_station")

# Get yearly precipitation from climateNA

paradise_climateNA <- read.csv("C:/PhD/Project/climategrowthshifts/analysis/mountrainier/data/climate/paradise_climateNA.csv", header = T)

PRCP_yearly_paradise_climateNA <- aggregate(PPT ~ year, data = paradise_climateNA, sum, na.rm = TRUE)

paradise_PRCP <- merge(PRCP_yearly_paradise_station, PRCP_yearly_paradise_climateNA, by = "year", all.x = TRUE)

write.csv(paradise_PRCP, "paradise_PRCP.csv")
#Plotting for paradise

paradise_PRCP$MY <- ifelse(paradise_PRCP$year %in% missing_years, "missing", "full")

# Plot
ggplot(paradise_PRCP, aes(x = PRCP_yearly_station, y = PPT, color = MY)) + scale_color_manual(values = c("missing" = "red", "full" = "black")) + geom_point() + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +  
  theme_classic() + xlim(2000, 4000) + ylim(2000, 4000) +
  labs(
    x = "PRCP from weather station",
    y = "PRCP from climateNA",
    color = "Missing years",
    title = "Paradise"
  )
