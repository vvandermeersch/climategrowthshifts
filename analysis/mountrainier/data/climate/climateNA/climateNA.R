# Started 21 June 2025 
# by Mao

library("ClimateNAr")

setwd("C:/PhD/Project/climategrowthshifts/analysis/mountrainier/data/climate/climateNA")
inputFile <- 'MORAstands.csv'
periodList <- 1901:2024
varList <- c('MAT','DD5','NFFD','FFP','PAS','Tmax01','Tmax02','Tmax03','Tmax04','Tmax05','Tmax06','Tmax07','Tmax08','Tmax09','Tmax10','Tmax11','Tmax12','Tmin01','Tmin02','Tmin03','Tmin04','Tmin05','Tmin06','Tmin07','Tmin08','Tmin09','Tmin10','Tmin11','Tmin12','Tave01','Tave02','Tave03','Tave04','Tave05','Tave06','Tave07','Tave08','Tave09','Tave10','Tave11','Tave12','PPT01','PPT02','PPT03','PPT04','PPT05','PPT06','PPT07','PPT08','PPT09','PPT10','PPT11','PPT12')
outDir <- 'C:/PhD/Project/climategrowthshifts/analysis/mountrainier/data/climate/climateNA/single/'
clm <- ClimateNAr(inputFile, periodList,varList,outDir)

setwd("C:/PhD/Project/climategrowthshifts/analysis/mountrainier/data/climate/climateNA/single")
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
names(combined_df)[names(combined_df) == "id1"] <- "stand_id"
write.csv(combined_df,"C:/PhD/Project/climategrowthshifts/analysis/mountrainier/data/climate/climateNAMORAStands.csv")
