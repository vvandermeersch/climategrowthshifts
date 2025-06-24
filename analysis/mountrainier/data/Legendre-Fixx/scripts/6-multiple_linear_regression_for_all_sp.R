##=================================================================================================================================================================================
##Organizing climate and RWI data for each species/site and testing the chosen linear model. File with coefficients, p-values, and standard error is saved for each species/site
##=================================================================================================================================================================================

#species/site are coded individually since the length of the chronologies vary
#PARA_CALNOO
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Climate WMA data")
climate<- read.csv("PARA_CALNOO_clim_1901-2015MSYT.csv", header=TRUE) #change filename for species/site
#load in climate variables from site/species file
DD5_sp <- climate$DD5_sp
DD5_sm <- climate$DD5_sm
RH_sp <- climate$RH_sp
PAS <- climate$PAS
prev_DD5_sm <- c(NA,DD5_sm)
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#remove the first value of all of the climate variables to remove the NA at the beginning of prev_GDD_sm
DD5_sp <- DD5_sp[-1]
DD5_sm <- DD5_sm[-1]
RH_sp <- RH_sp[-1]
PAS <- PAS[-1]
prev_DD5_sm <- prev_DD5_sm[-1]
#now the climate variable time series go from 1902-2015

#RWI data
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Master Chronologies 2017")
master <- read.csv(file="PARA_CALNOO_master.csv",header=TRUE)
master_trim <- master[-c(1:22),] #trim file to go from 1902-2015
rwi <- master_trim$RWI
rwi <- rwi[1:(length(rwi)-1)]#trim the value for 2015 if there is an NA (the multiple linear regression will not take an NA)
#If rwi vector is trimmed, then climate variable vectors should be trimmed as well
DD5_sp <- DD5_sp[1:(length(DD5_sp)-1)]
DD5_sm <- DD5_sm[1:(length(DD5_sm)-1)]
RH_sp <- RH_sp[1:(length(RH_sp)-1)]
PAS <- PAS[1:(length(PAS)-1)]
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#standardizing climate variables
DD5_sp<-(DD5_sp-mean(DD5_sp))/sd(DD5_sp)
DD5_sm<-(DD5_sm-mean(DD5_sm))/sd(DD5_sm)
RH_sp<-(RH_sp-mean(RH_sp))/sd(RH_sp)
PAS<-(PAS-mean(PAS))/sd(PAS)
prev_DD5_sm<-(prev_DD5_sm-mean(prev_DD5_sm))/sd(prev_DD5_sm)

test1 <- lm(rwi~DD5_sm+DD5_sp+RH_sp+PAS+prev_DD5_sm)
summary(test1)

betas <- summary(test1)$coefficients[,c(1,2,4)] #gives estimate (correlation coeff.) std. error, and p-value
betas

#saving the results as a dataframe
ParaCALNOO_growth_vars <- as.data.frame(betas)
ParaCALNOO_growth_vars <- ParaCALNOO_growth_vars[-1,]
ParaCALNOO_growth_vars[,4] <- "PARA"
ParaCALNOO_growth_vars[,5] <- "CALNOO"
#write the created dataframe to a csv file 
write.csv(ParaCALNOO_growth_vars, file ="C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/MLR_data/PARA_CALNOO_growth_all_vars_std.csv", row.names=TRUE) 

#SPRY CALNOO
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Climate WMA data")
climate<- read.csv("SPRY_CALNOO_clim_1901-2015MSYT.csv", header=TRUE) #change filename for species/site
DD5_sp <- climate$DD5_sp
DD5_sm <- climate$DD5_sm
RH_sp <- climate$RH_sp
PAS <- climate$PAS
prev_DD5_sm <- c(NA,DD5_sm)
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#remove the first value of all of the climate variables to remove the NA at the beginning of prev_GDD_sm
DD5_sp <- DD5_sp[-1]
DD5_sm <- DD5_sm[-1]
RH_sp <- RH_sp[-1]
PAS <- PAS[-1]
prev_DD5_sm <- prev_DD5_sm[-1]
#now the climate variable time series go from 1902-2015

#RWI
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Master Chronologies 2017")
master <- read.csv(file="SPRY_CALNOO_master.csv",header=TRUE)
master_trim <- master[-c(1:22),] #trim file to go from 1902-2015
rwi <- master_trim$RWI
rwi <- rwi[1:(length(rwi)-1)]#trim the value for 2015 if there is an NA (the multiple linear regression will not take an NA)
#If rwi vector is trimmed, climate variable vectors must be trimmed as well
DD5_sp <- DD5_sp[1:(length(DD5_sp)-1)]
DD5_sm <- DD5_sm[1:(length(DD5_sm)-1)]
RH_sp <- RH_sp[1:(length(RH_sp)-1)]
PAS <- PAS[1:(length(PAS)-1)]
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#standardizing climate variables 
DD5_sp<-(DD5_sp-mean(DD5_sp))/sd(DD5_sp)
DD5_sm<-(DD5_sm-mean(DD5_sm))/sd(DD5_sm)
RH_sp<-(RH_sp-mean(RH_sp))/sd(RH_sp)
PAS<-(PAS-mean(PAS))/sd(PAS)
prev_DD5_sm<-(prev_DD5_sm-mean(prev_DD5_sm))/sd(prev_DD5_sm)

test1 <- lm(rwi~DD5_sm+DD5_sp+RH_sp+PAS+prev_DD5_sm)
summary(test1)

betas <- summary(test1)$coefficients[,c(1,2,4)] #gives estimate (correlation coeff.), std. error, and p-value
betas

#saving the results as a dataframe 
SpryCALNOO_growth_vars <- as.data.frame(betas)
SpryCALNOO_growth_vars <- SpryCALNOO_growth_vars[-1,]
SpryCALNOO_growth_vars[,4] <- "SPRY"
SpryCALNOO_growth_vars[,5] <- "CALNOO"
#write the created dataframe to a csv file 
write.csv(SpryCALNOO_growth_vars, file ="C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/MLR_data/SPRY_CALNOO_growth_all_vars_std.csv", row.names=TRUE) 

#SUNR CALNOO
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Climate WMA data")
climate<- read.csv("SUNR_CALNOO_clim_1901-2015MSYT.csv", header=TRUE) #change filename for species/site
DD5_sp <- climate$DD5_sp
DD5_sm <- climate$DD5_sm
RH_sp <- climate$RH_sp
PAS <- climate$PAS
prev_DD5_sm <- c(NA,DD5_sm)
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#remove the first value of all of the climate variables to remove the NA at the beginning of prev_GDD_sm
DD5_sp <- DD5_sp[-1]
DD5_sm <- DD5_sm[-1]
RH_sp <- RH_sp[-1]
PAS <- PAS[-1]
prev_DD5_sm <- prev_DD5_sm[-1]
#now the climate variable time series go from 1902-2015

#RWI
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Master Chronologies 2017")
master <- read.csv(file="SUNR_CALNOO_master.csv",header=TRUE)
master_trim <- master[-c(1:22),] #trim file to go from 1902-2015
rwi <- master_trim$RWI
rwi <- rwi[1:(length(rwi)-1)]#trim the value for 2015 if there is an NA (the multiple linear regression will not take an NA)
#If rwi vector is trimmed, then climate variable vectors must be trimmed as well
DD5_sp <- DD5_sp[1:(length(DD5_sp)-1)]
DD5_sm <- DD5_sm[1:(length(DD5_sm)-1)]
RH_sp <- RH_sp[1:(length(RH_sp)-1)]
PAS <- PAS[1:(length(PAS)-1)]
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#standardizing climate variables
DD5_sp<-(DD5_sp-mean(DD5_sp))/sd(DD5_sp)
DD5_sm<-(DD5_sm-mean(DD5_sm))/sd(DD5_sm)
RH_sp<-(RH_sp-mean(RH_sp))/sd(RH_sp)
PAS<-(PAS-mean(PAS))/sd(PAS)
prev_DD5_sm<-(prev_DD5_sm-mean(prev_DD5_sm))/sd(prev_DD5_sm)

test1 <- lm(rwi~DD5_sm+DD5_sp+RH_sp+PAS+prev_DD5_sm) #use for all species
summary(test1)

betas <- summary(test1)$coefficients[,c(1,2,4)] #gives estimate (correlation coeff.) std. error, and p-value
betas

#saving the results as a dataframe
SunrCALNOO_growth_vars <- as.data.frame(betas)
SunrCALNOO_growth_vars <- SunrCALNOO_growth_vars[-1,]
SunrCALNOO_growth_vars[,4] <- "SUNR"
SunrCALNOO_growth_vars[,5] <- "CALNOO"
#write the created dataframe to a csv file 
write.csv(SunrCALNOO_growth_vars, file ="C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/MLR_data/SUNR_CALNOO_growth_all_vars_std.csv", row.names=TRUE) 

#PARA TSUMER
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Climate WMA data")
climate<- read.csv("PARA_TSUMER_clim_1901-2015MSYT.csv", header=TRUE) #change filename for species/site
DD5_sp <- climate$DD5_sp
DD5_sm <- climate$DD5_sm
RH_sp <- climate$RH_sp
PAS <- climate$PAS
prev_DD5_sm <- c(NA,DD5_sm)
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#remove the first value of all of the climate variables to remove the NA at the beginning of prev_GDD_sm
DD5_sp <- DD5_sp[-1]
DD5_sm <- DD5_sm[-1]
RH_sp <- RH_sp[-1]
PAS <- PAS[-1]
prev_DD5_sm <- prev_DD5_sm[-1]
#now the climate variable time series go from 1902-2015

#RWI
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Master Chronologies 2017")
master <- read.csv(file="PARA_TSUMER_master.csv",header=TRUE)
master_trim <- master[-c(1:22),] #trim file to go from 1902-2015
rwi <- master_trim$RWI
rwi <- rwi[1:(length(rwi)-1)]#trim the value for 2015 if there is an NA (the multiple linear regression will not take an NA)
#If I trim the rwi vector, then I need to trim the climate variable vectors as well
DD5_sp <- DD5_sp[1:(length(DD5_sp)-1)]
DD5_sm <- DD5_sm[1:(length(DD5_sm)-1)]
RH_sp <- RH_sp[1:(length(RH_sp)-1)]
PAS <- PAS[1:(length(PAS)-1)]
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#standardizing climate variables 
DD5_sp<-(DD5_sp-mean(DD5_sp))/sd(DD5_sp)
DD5_sm<-(DD5_sm-mean(DD5_sm))/sd(DD5_sm)
RH_sp<-(RH_sp-mean(RH_sp))/sd(RH_sp)
PAS<-(PAS-mean(PAS))/sd(PAS)
prev_DD5_sm<-(prev_DD5_sm-mean(prev_DD5_sm))/sd(prev_DD5_sm)

test1 <- lm(rwi~DD5_sm+DD5_sp+RH_sp+PAS+prev_DD5_sm) #use for all species
summary(test1)

betas <- summary(test1)$coefficients[,c(1,2,4)] #gives estimate (correlation coeff.) std. error, and p-value
betas

#saving the results as a dataframe 
ParaTSUMER_growth_vars <- as.data.frame(betas)
ParaTSUMER_growth_vars <- ParaTSUMER_growth_vars[-1,]
ParaTSUMER_growth_vars[,4] <- "PARA"
ParaTSUMER_growth_vars[,5] <- "TSUMER"
#write the created dataframe to a csv file 
write.csv(ParaTSUMER_growth_vars, file ="C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/MLR_data/PARA_TSUMER_growth_all_vars_std.csv", row.names=TRUE) 

#SPRY TSUMER
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Climate WMA data")
climate<- read.csv("SPRY_TSUMER_clim_1901-2015MSYT.csv", header=TRUE) #change filename for species/site
DD5_sp <- climate$DD5_sp
DD5_sm <- climate$DD5_sm
RH_sp <- climate$RH_sp
PAS <- climate$PAS
prev_DD5_sm <- c(NA,DD5_sm)
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#remove the first value of all of the climate variables to remove the NA at the beginning of prev_GDD_sm
DD5_sp <- DD5_sp[-1]
DD5_sm <- DD5_sm[-1]
RH_sp <- RH_sp[-1]
PAS <- PAS[-1]
prev_DD5_sm <- prev_DD5_sm[-1]
#now the climate variable time series go from 1902-2015

#RWI
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Master Chronologies 2017")
master <- read.csv(file="SPRY_TSUMER_master.csv",header=TRUE)
master_trim <- master[-c(1:22),] #trim file to go from 1902-2015
rwi <- master_trim$RWI
rwi <- rwi[1:(length(rwi)-1)]#trim the value for 2015 if there is an NA (the multiple linear regression will not take an NA)
#If I trim the rwi vector, then I need to trim the climate variable vectors as well
DD5_sp <- DD5_sp[1:(length(DD5_sp)-1)]
DD5_sm <- DD5_sm[1:(length(DD5_sm)-1)]
RH_sp <- RH_sp[1:(length(RH_sp)-1)]
PAS <- PAS[1:(length(PAS)-1)]
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#standardizing climate variables
DD5_sp<-(DD5_sp-mean(DD5_sp))/sd(DD5_sp)
DD5_sm<-(DD5_sm-mean(DD5_sm))/sd(DD5_sm)
RH_sp<-(RH_sp-mean(RH_sp))/sd(RH_sp)
PAS<-(PAS-mean(PAS))/sd(PAS)
prev_DD5_sm<-(prev_DD5_sm-mean(prev_DD5_sm))/sd(prev_DD5_sm)

test1 <- lm(rwi~DD5_sm+DD5_sp+RH_sp+PAS+prev_DD5_sm) #use for all species
summary(test1)

betas <- summary(test1)$coefficients[,c(1,2,4)] #gives estimate (correlation coeff.) std. error, and p-value
betas

#saving the results as a dataframe
SpryTSUMER_growth_vars <- as.data.frame(betas)
SpryTSUMER_growth_vars <- SpryTSUMER_growth_vars[-1,]
SpryTSUMER_growth_vars[,4] <- "SPRY"
SpryTSUMER_growth_vars[,5] <- "TSUMER"
#write the created dataframe to a csv file --at the end I will read all of these in and combine them
write.csv(SpryTSUMER_growth_vars, file ="C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/MLR_data/SPRY_TSUMER_growth_all_vars_std.csv", row.names=TRUE) 

#SUNR TSUMER
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Climate WMA data")
climate<- read.csv("SUNR_TSUMER_clim_1901-2015MSYT.csv", header=TRUE) #change filename for species/site
DD5_sp <- climate$DD5_sp
DD5_sm <- climate$DD5_sm
RH_sp <- climate$RH_sp
PAS <- climate$PAS
prev_DD5_sm <- c(NA,DD5_sm)
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#remove the first value of all of the climate variables to remove the NA at the beginning of prev_GDD_sm
DD5_sp <- DD5_sp[-1]
DD5_sm <- DD5_sm[-1]
RH_sp <- RH_sp[-1]
PAS <- PAS[-1]
prev_DD5_sm <- prev_DD5_sm[-1]
#now the climate variable time series go from 1902-2015

#RWI
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Master Chronologies 2017")
master <- read.csv(file="SUNR_TSUMER_master.csv",header=TRUE)
master_trim <- master[-c(1:22),] #trim file to go from 1902-2015
rwi <- master_trim$RWI#no NA value for the year of 2015 in this chronology

#standardizing climate variables 
DD5_sp<-(DD5_sp-mean(DD5_sp))/sd(DD5_sp)
DD5_sm<-(DD5_sm-mean(DD5_sm))/sd(DD5_sm)
RH_sp<-(RH_sp-mean(RH_sp))/sd(RH_sp)
PAS<-(PAS-mean(PAS))/sd(PAS)
prev_DD5_sm<-(prev_DD5_sm-mean(prev_DD5_sm))/sd(prev_DD5_sm)

test1 <- lm(rwi~DD5_sm+DD5_sp+RH_sp+PAS+prev_DD5_sm) #use for all species
summary(test1)

betas <- summary(test1)$coefficients[,c(1,2,4)] #gives estimate (correlation coeff.) std. error, and p-value
betas

#saving the results as a dataframe 
SunrTSUMER_growth_vars <- as.data.frame(betas)
SunrTSUMER_growth_vars <- SunrTSUMER_growth_vars[-1,]
SunrTSUMER_growth_vars[,4] <- "SUNR"
SunrTSUMER_growth_vars[,5] <- "TSUMER"
#write the created dataframe to a csv file
write.csv(SunrTSUMER_growth_vars, file ="C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/MLR_data/SUNR_TSUMER_growth_all_vars_std.csv", row.names=TRUE) 

## PARA ABIAMA
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Climate WMA data")
climate<- read.csv("PARA_ABIAMA_clim_1901-2015MSYT.csv", header=TRUE) #change filename for species/site
DD5_sp <- climate$DD5_sp
DD5_sm <- climate$DD5_sm
RH_sp <- climate$RH_sp
PAS <- climate$PAS
prev_DD5_sm <- c(NA,DD5_sm)
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#remove the first value of all of the climate variables to remove the NA at the beginning of prev_GDD_sm
DD5_sp <- DD5_sp[-1]
DD5_sm <- DD5_sm[-1]
RH_sp <- RH_sp[-1]
PAS <- PAS[-1]
prev_DD5_sm <- prev_DD5_sm[-1]
#now the climate variable time series go from 1902-2015

#RWI data
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Master Chronologies 2017")
master <- read.csv(file="PARA_ABIAMA_master.csv",header=TRUE)
master_trim <- master[-c(1:22),] #trim file to go from 1902-2015
rwi <- master_trim$RWI
rwi <- rwi[1:(length(rwi)-1)]#trim the value for 2015 if there is an NA (the multiple linear regression will not take an NA)
#If rwi vector is trimmed, then climate variable vectors should be trimmed as well
DD5_sp <- DD5_sp[1:(length(DD5_sp)-1)]
DD5_sm <- DD5_sm[1:(length(DD5_sm)-1)]
RH_sp <- RH_sp[1:(length(RH_sp)-1)]
PAS <- PAS[1:(length(PAS)-1)]
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#standardizing climate variables
DD5_sp<-(DD5_sp-mean(DD5_sp))/sd(DD5_sp)
DD5_sm<-(DD5_sm-mean(DD5_sm))/sd(DD5_sm)
RH_sp<-(RH_sp-mean(RH_sp))/sd(RH_sp)
PAS<-(PAS-mean(PAS))/sd(PAS)
prev_DD5_sm<-(prev_DD5_sm-mean(prev_DD5_sm))/sd(prev_DD5_sm)

test1 <- lm(rwi~DD5_sm+DD5_sp+RH_sp+PAS+prev_DD5_sm)
summary(test1)

betas <- summary(test1)$coefficients[,c(1,2,4)] #gives estimate (correlation coeff.) std. error, and p-value
betas

#saving the results as a dataframe
ParaABIAMA_growth_vars <- as.data.frame(betas)
ParaABIAMA_growth_vars <- ParaABIAMA_growth_vars[-1,]
ParaABIAMA_growth_vars[,4] <- "PARA"
ParaABIAMA_growth_vars[,5] <- "ABIAMA"
#write the created dataframe to a csv file 
write.csv(ParaABIAMA_growth_vars, file ="C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/MLR_data/PARA_ABIAMA_growth_all_vars_std.csv", row.names=TRUE)
#SPRY ABIAMA
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Climate WMA data")
climate<- read.csv("SPRY_ABIAMA_clim_1901-2015MSYT.csv", header=TRUE) #change filename for species/site
DD5_sp <- climate$DD5_sp
DD5_sm <- climate$DD5_sm
RH_sp <- climate$RH_sp
PAS <- climate$PAS
prev_DD5_sm <- c(NA,DD5_sm)
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#remove the first value of all of the climate variables to remove the NA at the beginning of prev_GDD_sm
DD5_sp <- DD5_sp[-1]
DD5_sm <- DD5_sm[-1]
RH_sp <- RH_sp[-1]
PAS <- PAS[-1]
prev_DD5_sm <- prev_DD5_sm[-1]
#now the climate variable time series go from 1902-2015

#RWI
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Master Chronologies 2017")
master <- read.csv(file="SPRY_ABIAMA_master.csv",header=TRUE)
master_trim <- master[-c(1:22),] #trim file to go from 1902-2015
rwi <- master_trim$RWI
rwi <- rwi[1:(length(rwi)-1)]#trim the value for 2015 if there is an NA (the multiple linear regression will not take an NA)
#If rwi vector is trimmed, climate variable vectors must be trimmed as well
DD5_sp <- DD5_sp[1:(length(DD5_sp)-1)]
DD5_sm <- DD5_sm[1:(length(DD5_sm)-1)]
RH_sp <- RH_sp[1:(length(RH_sp)-1)]
PAS <- PAS[1:(length(PAS)-1)]
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#standardizing climate variables 
DD5_sp<-(DD5_sp-mean(DD5_sp))/sd(DD5_sp)
DD5_sm<-(DD5_sm-mean(DD5_sm))/sd(DD5_sm)
RH_sp<-(RH_sp-mean(RH_sp))/sd(RH_sp)
PAS<-(PAS-mean(PAS))/sd(PAS)
prev_DD5_sm<-(prev_DD5_sm-mean(prev_DD5_sm))/sd(prev_DD5_sm)

test1 <- lm(rwi~DD5_sm+DD5_sp+RH_sp+PAS+prev_DD5_sm)
summary(test1)

betas <- summary(test1)$coefficients[,c(1,2,4)] #gives estimate (correlation coeff.), std. error, and p-value
betas

#saving the results as a dataframe 
SpryABIAMA_growth_vars <- as.data.frame(betas)
SpryABIAMA_growth_vars <- SpryABIAMA_growth_vars[-1,]
SpryABIAMA_growth_vars[,4] <- "SPRY"
SpryABIAMA_growth_vars[,5] <- "ABIAMA"
#write the created dataframe to a csv file 
write.csv(SpryABIAMA_growth_vars, file ="C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/MLR_data/SPRY_ABIAMA_growth_all_vars_std.csv", row.names=TRUE) 

#SUNR ABIAMA
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Climate WMA data")
climate<- read.csv("SUNR_ABIAMA_clim_1901-2015MSYT.csv", header=TRUE) #change filename for species/site
DD5_sp <- climate$DD5_sp
DD5_sm <- climate$DD5_sm
RH_sp <- climate$RH_sp
PAS <- climate$PAS
prev_DD5_sm <- c(NA,DD5_sm)
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#remove the first value of all of the climate variables to remove the NA at the beginning of prev_GDD_sm
DD5_sp <- DD5_sp[-1]
DD5_sm <- DD5_sm[-1]
RH_sp <- RH_sp[-1]
PAS <- PAS[-1]
prev_DD5_sm <- prev_DD5_sm[-1]
#now the climate variable time series go from 1902-2015

#RWI
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Master Chronologies 2017")
master <- read.csv(file="SUNR_ABIAMA_master.csv",header=TRUE)
master_trim <- master[-c(1:22),] #trim file to go from 1902-2015
rwi <- master_trim$RWI
rwi <- rwi[1:(length(rwi)-1)]#trim the value for 2015 if there is an NA (the multiple linear regression will not take an NA)
#If rwi vector is trimmed, then climate variable vectors must be trimmed as well
DD5_sp <- DD5_sp[1:(length(DD5_sp)-1)]
DD5_sm <- DD5_sm[1:(length(DD5_sm)-1)]
RH_sp <- RH_sp[1:(length(RH_sp)-1)]
PAS <- PAS[1:(length(PAS)-1)]
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#standardizing climate variables
DD5_sp<-(DD5_sp-mean(DD5_sp))/sd(DD5_sp)
DD5_sm<-(DD5_sm-mean(DD5_sm))/sd(DD5_sm)
RH_sp<-(RH_sp-mean(RH_sp))/sd(RH_sp)
PAS<-(PAS-mean(PAS))/sd(PAS)
prev_DD5_sm<-(prev_DD5_sm-mean(prev_DD5_sm))/sd(prev_DD5_sm)

test1 <- lm(rwi~DD5_sm+DD5_sp+RH_sp+PAS+prev_DD5_sm) #use for all species
summary(test1)

betas <- summary(test1)$coefficients[,c(1,2,4)] #gives estimate (correlation coeff.) std. error, and p-value
betas

#saving the results as a dataframe
SunrABIAMA_growth_vars <- as.data.frame(betas)
SunrABIAMA_growth_vars <- SunrABIAMA_growth_vars[-1,]
SunrABIAMA_growth_vars[,4] <- "SUNR"
SunrABIAMA_growth_vars[,5] <- "ABIAMA"
#write the created dataframe to a csv file 
write.csv(SunrABIAMA_growth_vars, file ="C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/MLR_data/SUNR_ABIAMA_growth_all_vars_std.csv", row.names=TRUE)

##PARA ABILAS
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Climate WMA data")
climate<- read.csv("PARA_ABILAS_clim_1901-2015MSYT.csv", header=TRUE) #change filename for species/site
DD5_sp <- climate$DD5_sp
DD5_sm <- climate$DD5_sm
RH_sp <- climate$RH_sp
PAS <- climate$PAS
prev_DD5_sm <- c(NA,DD5_sm)
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#remove the first value of all of the climate variables to remove the NA at the beginning of prev_GDD_sm
DD5_sp <- DD5_sp[-1]
DD5_sm <- DD5_sm[-1]
RH_sp <- RH_sp[-1]
PAS <- PAS[-1]
prev_DD5_sm <- prev_DD5_sm[-1]
#now the climate variable time series go from 1902-2015

#RWI data
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Master Chronologies 2017")
master <- read.csv(file="PARA_ABILAS_master.csv",header=TRUE)
master_trim <- master[-c(1:22),] #trim file to go from 1902-2015
rwi <- master_trim$RWI
#If rwi vector is trimmed, then climate variable vectors should be trimmed as well
DD5_sp <- DD5_sp[1:(length(DD5_sp)-4)]
DD5_sm <- DD5_sm[1:(length(DD5_sm)-4)]
RH_sp <- RH_sp[1:(length(RH_sp)-4)]
PAS <- PAS[1:(length(PAS)-4)]
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-4)]

#standardizing climate variables
DD5_sp<-(DD5_sp-mean(DD5_sp))/sd(DD5_sp)
DD5_sm<-(DD5_sm-mean(DD5_sm))/sd(DD5_sm)
RH_sp<-(RH_sp-mean(RH_sp))/sd(RH_sp)
PAS<-(PAS-mean(PAS))/sd(PAS)
prev_DD5_sm<-(prev_DD5_sm-mean(prev_DD5_sm))/sd(prev_DD5_sm)

test1 <- lm(rwi~DD5_sm+DD5_sp+RH_sp+PAS+prev_DD5_sm)
summary(test1)

betas <- summary(test1)$coefficients[,c(1,2,4)] #gives estimate (correlation coeff.) std. error, and p-value
betas

#saving the results as a dataframe
ParaABILAS_growth_vars <- as.data.frame(betas)
ParaABILAS_growth_vars <- ParaABILAS_growth_vars[-1,]
ParaABILAS_growth_vars[,4] <- "PARA"
ParaABILAS_growth_vars[,5] <- "ABILAS"
#write the created dataframe to a csv file 
write.csv(ParaABILAS_growth_vars, file ="C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/MLR_data/PARA_ABILAS_growth_all_vars_std.csv", row.names=TRUE)
#SPRY ABILAS
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Climate WMA data")
climate<- read.csv("SPRY_ABILAS_clim_1901-2015MSYT.csv", header=TRUE) #change filename for species/site
DD5_sp <- climate$DD5_sp
DD5_sm <- climate$DD5_sm
RH_sp <- climate$RH_sp
PAS <- climate$PAS
prev_DD5_sm <- c(NA,DD5_sm)
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#remove the first value of all of the climate variables to remove the NA at the beginning of prev_GDD_sm
DD5_sp <- DD5_sp[-1]
DD5_sm <- DD5_sm[-1]
RH_sp <- RH_sp[-1]
PAS <- PAS[-1]
prev_DD5_sm <- prev_DD5_sm[-1]
#now the climate variable time series go from 1902-2015

#RWI
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Master Chronologies 2017")
master <- read.csv(file="SPRY_ABILAS_master.csv",header=TRUE)
master_trim <- master[-c(1:22),] #trim file to go from 1902-2015
rwi <- master_trim$RWI

#standardizing climate variables 
DD5_sp<-(DD5_sp-mean(DD5_sp))/sd(DD5_sp)
DD5_sm<-(DD5_sm-mean(DD5_sm))/sd(DD5_sm)
RH_sp<-(RH_sp-mean(RH_sp))/sd(RH_sp)
PAS<-(PAS-mean(PAS))/sd(PAS)
prev_DD5_sm<-(prev_DD5_sm-mean(prev_DD5_sm))/sd(prev_DD5_sm)

test1 <- lm(rwi~DD5_sm+DD5_sp+RH_sp+PAS+prev_DD5_sm)
summary(test1)

betas <- summary(test1)$coefficients[,c(1,2,4)] #gives estimate (correlation coeff.), std. error, and p-value
betas

#saving the results as a dataframe 
SpryABILAS_growth_vars <- as.data.frame(betas)
SpryABILAS_growth_vars <- SpryABILAS_growth_vars[-1,]
SpryABILAS_growth_vars[,4] <- "SPRY"
SpryABILAS_growth_vars[,5] <- "ABILAS"
#write the created dataframe to a csv file 
write.csv(SpryABILAS_growth_vars, file ="C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/MLR_data/SPRY_ABILAS_growth_all_vars_std.csv", row.names=TRUE) 

#SUNR ABILAS
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Climate WMA data")
climate<- read.csv("SUNR_ABILAS_clim_1901-2015MSYT.csv", header=TRUE) #change filename for species/site
DD5_sp <- climate$DD5_sp
DD5_sm <- climate$DD5_sm
RH_sp <- climate$RH_sp
PAS <- climate$PAS
prev_DD5_sm <- c(NA,DD5_sm)
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#remove the first value of all of the climate variables to remove the NA at the beginning of prev_GDD_sm
DD5_sp <- DD5_sp[-1]
DD5_sm <- DD5_sm[-1]
RH_sp <- RH_sp[-1]
PAS <- PAS[-1]
prev_DD5_sm <- prev_DD5_sm[-1]
#now the climate variable time series go from 1902-2015

#RWI
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Master Chronologies 2017")
master <- read.csv(file="SUNR_ABILAS_master.csv",header=TRUE)
master_trim <- master[-c(1:22),] #trim file to go from 1902-2015
rwi <- master_trim$RWI

#standardizing climate variables
DD5_sp<-(DD5_sp-mean(DD5_sp))/sd(DD5_sp)
DD5_sm<-(DD5_sm-mean(DD5_sm))/sd(DD5_sm)
RH_sp<-(RH_sp-mean(RH_sp))/sd(RH_sp)
PAS<-(PAS-mean(PAS))/sd(PAS)
prev_DD5_sm<-(prev_DD5_sm-mean(prev_DD5_sm))/sd(prev_DD5_sm)

test1 <- lm(rwi~DD5_sm+DD5_sp+RH_sp+PAS+prev_DD5_sm) #use for all species
summary(test1)

betas <- summary(test1)$coefficients[,c(1,2,4)] #gives estimate (correlation coeff.) std. error, and p-value
betas

#saving the results as a dataframe
SunrABILAS_growth_vars <- as.data.frame(betas)
SunrABILAS_growth_vars <- SunrABILAS_growth_vars[-1,]
SunrABILAS_growth_vars[,4] <- "SUNR"
SunrABILAS_growth_vars[,5] <- "ABILAS"
#write the created dataframe to a csv file 
write.csv(SunrABILAS_growth_vars, file ="C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/MLR_data/SUNR_ABILAS_growth_all_vars_std.csv", row.names=TRUE)



