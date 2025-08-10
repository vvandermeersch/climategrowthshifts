##load some libraries
library(data.table)
library(lubridate)
library(dplyr)
options(stringsAsFactors = FALSE) ##disable strings as factors
options(scipen=999) ##disable scientific notation

##This script calculates O-18 source water from climate data (air temperature and humidity) and measured cellulose O-18 data 

## The code requires two input files to run, "WY1992to2019_LochVale_Climate.csv" and "treeiso.csv"

## The file "WY1992to2019_LochVale_Climate.csv" is hourly climate data for the Andrews Creek Weather Station and can be downloaded from:
##  Akie, G.A., McDermott, W.R., Sexstone, G.A., and Clow D.W., 2020, Climatological data for the Loch Vale 
##  watershed in Rocky Mountain National Park, Colorado, water years 1992-2019: U.S. Geological Survey data release, 
##  https://doi.org/10.5066/P92ULNAG.

## The file "treeiso.csv" contains the measured cellulose O-18 values and is available in this data release.

## To run the code put the R script and 2 input files in same directory and set working directory to 'source file location'

## read in climate data file obtained from the data release listed above, make sure the file name matches the script
wxdata<-read.csv("WY1992to2019_LochVale_Climate.csv")
wxdata$Date<-as_datetime(wxdata$timestamp)
wxdata$month<-as.numeric(month(wxdata$Date))
wxdata$year<-as.numeric(year(wxdata$Date))
wxdata$time<-as.numeric(hour(wxdata$Date))

#filter data based on growing season (June, Jul, Aug, Sept) and daylight hours (8 am to 7 pm) for Andrews Creek weather station
wxand<-wxdata%>%filter(station_name=="Andrews Creek weather station")%>%filter(month>5&month<10)%>%filter(time>7&time<20)

#calculate average annual growing season values for Air Temperature and Relative Humidiy (RH)
andsumm<-wxand%>%group_by(station_name,year)%>%summarize(aveT_a=mean(T_air,na.rm=TRUE),countT_a=sum(!is.na(T_air))/1486*100,RH_a=mean(RH,na.rm=TRUE),countRH_a=sum(!is.na(RH))/1486*100)%>%data.frame()

##### calculate source water

#read in tree isotope values from the file treeiso.csv available in this data release, make sure the file name matches the script
isodata<-read.csv("treeiso.csv",colClasses="numeric") 
anddata<-andsumm[3:25,]

#Loop used to calculate source water 18-O for each tree following equation 2 in the manuscript

allsource<-data.frame(matrix(ncol=11,nrow=23)) #initialize empty dataframe
for (k in 2:9){
  Sourcewater<-data.frame()
for (i in 1:length(isodata$year)){
  leaftemp<-anddata$aveT_a[i]+3+273.15
  Epsilon<-(-7.685 + (6.7123*(1000/(leaftemp))) - (1.6664*(1000000/leaftemp^2)) + (0.35041*(1000000000/leaftemp^3)))
  Dewpt<-(243.04*(log(anddata$RH_a[i]*.01)+((17.625*anddata$aveT_a[i])/(243.04+anddata$aveT_a[i])))/(17.625-log(anddata$RH_a[i]*.01)+((17.625*anddata$aveT_a[i])/(243.04+anddata$aveT_a[i]))))
  RHleaf<-(100*(exp((17.625*Dewpt)/(243.04+Dewpt))/exp((17.625*(anddata$aveT_a[i]+5))/(243.04+(anddata$aveT_a[i]+5))))/100)
  result<-(isodata[i,k]-((1-.4)*(1-RHleaf)*(Epsilon + 28)) - 27)
  year<-isodata$year[i]
  AirT<-anddata$aveT_a[i]
  RH<-anddata$RH_a[i]
  Sourcewater<-rbind(data.frame(year,AirT,RH,result),Sourcewater)
}
  allsource[,1]<-Sourcewater[,1]
  allsource[,2]<-Sourcewater[,2]
  allsource[,3]<-Sourcewater[,3]
  allsource[,k+2]<-Sourcewater[,4]
}

##assign column names to the script output and write results to .csv file to directory containing the source code

names(allsource)<-c("year","AirT","RH", "AMS06B", "AMS07A", "AMC01B", "AMC08A", "AMC09A", "Spring", "Creek" , "ALL")   
write.csv(allsource,"sourcewater.csv",row.names=FALSE)
