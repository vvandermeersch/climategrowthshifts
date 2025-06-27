###===========================================================================================================================================================================================
##Organizing climate variables and RWI data for each species at every site and then performing the dredge() function to determine the importance of the climate vars selected on tree growth
##============================================================================================================================================================================================

#PARA ABIAMA
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Climate WMA data")
climate <- read.csv("PARA_ABIAMA_clim_1901-2015MSYT.csv", header=TRUE) #change filename for site/species

#load in climate variables from species/site file
DD5_sp <- climate$DD5_sp
DD5_sm <- climate$DD5_sm
RH_sp <- climate$RH_sp
PAS <- climate$PAS
prev_DD5_sm <- c(NA,DD5_sm)
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#removes the first value of all of the climate variables to remove NA at the beginning of prev_GDD_sm
DD5_sp <- DD5_sp[-1]
DD5_sm <- DD5_sm[-1]
RH_sp <- RH_sp[-1]
PAS <- PAS[-1]
prev_DD5_sm <- prev_DD5_sm[-1]

#now the climate variable time series go from 1902-2015

#RWI
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Master Chronologies 2017")
master <- read.csv(file="PARA_ABIAMA_master.csv",header=TRUE)
master_trim <- master[-c(1:22),] #shortens chronology from extending to 1880 to 1902
rwi <- master_trim$RWI
rwi <- rwi[1:(length(rwi)-1)]#trim the value for 2015 if the RWI value is an NA
#If rwi is trimmed, climate variables must be trimmed as well #all values should be the same length
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

require(MuMIn)
mod <- lm(rwi~DD5_sm+DD5_sp+RH_sp+PAS+prev_DD5_sm) #fit full model
options(na.action = "na.fail")
moddredge <- dredge(mod, extra = "R^2") # fits all nested models, creates a table of their coefficients and AICs, and R2s
#modavgs <- model.avg(moddredge, subset= delta <=2) #takes all models with deltaAIC<2 from the best, and averages their coefficients to get the 'best model estimate'
modavgs <- model.avg(moddredge, cumsum(weight) <= .95)
modavgs$importance

#SPRY ABIAMA
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Climate WMA data")
climate <- read.csv("SPRY_ABIAMA_clim_1901-2015MSYT.csv", header=TRUE) #change filename for site/species

#load in climate variables from species/site file
DD5_sp <- climate$DD5_sp
DD5_sm <- climate$DD5_sm
RH_sp <- climate$RH_sp
PAS <- climate$PAS
prev_DD5_sm <- c(NA,DD5_sm)
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#removes the first value of all of the climate variables so that I can get rid of the NA at the beginning of prev_GDD_sm
DD5_sp <- DD5_sp[-1]
DD5_sm <- DD5_sm[-1]
RH_sp <- RH_sp[-1]
PAS <- PAS[-1]
prev_DD5_sm <- prev_DD5_sm[-1]

#now the climate variable time series go from 1902-2015

#RWI
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Master Chronologies 2017")
master <- read.csv(file="SPRY_ABIAMA_master.csv",header=TRUE)
master_trim <- master[-c(1:22),] #shortens chronology from extending to 1880 to 1902
rwi <- master_trim$RWI
rwi <- rwi[1:(length(rwi)-1)]#trim the value for 2015 if there is an NA
#If rwi is trimmed, climate variables must be trimmed as well
DD5_sp <- DD5_sp[1:(length(DD5_sp)-1)]
DD5_sm <- DD5_sm[1:(length(DD5_sm)-1)]
RH_sp <- RH_sp[1:(length(RH_sp)-1)]
PAS <- PAS[1:(length(PAS)-1)]
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#all values should be the same length

#standardizing climate variables
DD5_sp<-(DD5_sp-mean(DD5_sp))/sd(DD5_sp)
DD5_sm<-(DD5_sm-mean(DD5_sm))/sd(DD5_sm)
RH_sp<-(RH_sp-mean(RH_sp))/sd(RH_sp)
PAS<-(PAS-mean(PAS))/sd(PAS)
prev_DD5_sm<-(prev_DD5_sm-mean(prev_DD5_sm))/sd(prev_DD5_sm)

require(MuMIn)
mod <- lm(rwi~DD5_sm+DD5_sp+RH_sp+PAS+prev_DD5_sm) #fit full model
options(na.action = "na.fail")
moddredge <- dredge(mod, extra = "R^2") # fits all nested models, creates a table of their coefficients and AICs, and R2s
#modavgs <- model.avg(moddredge, subset= delta <=2) #takes all models with deltaAIC<2 from the best, and averages their coefficients to get the 'best model estimate'
modavgs <- model.avg(moddredge, cumsum(weight) <= .95)
modavgs$importance

#SUNR ABIAMA
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Climate WMA data")
climate <- read.csv("SUNR_ABIAMA_clim_1901-2015MSYT.csv", header=TRUE) #change filename for site/species

#load in climate variables from site/species file
DD5_sp <- climate$DD5_sp
DD5_sm <- climate$DD5_sm
RH_sp <- climate$RH_sp
PAS <- climate$PAS
prev_DD5_sm <- c(NA,DD5_sm)
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#removes the first value of all of the climate variables so that I can get rid of the NA at the beginning of prev_GDD_sm
DD5_sp <- DD5_sp[-1]
DD5_sm <- DD5_sm[-1]
RH_sp <- RH_sp[-1]
PAS <- PAS[-1]
prev_DD5_sm <- prev_DD5_sm[-1]

#now the climate variable time series go from 1902-2015

#RWI
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Master Chronologies 2017")
master <- read.csv(file="SUNR_ABIAMA_master.csv",header=TRUE)
master_trim <- master[-c(1:22),] #shortens chronology from extending to 1880 to 1902
rwi <- master_trim$RWI
rwi <- rwi[1:(length(rwi)-1)]#trim the value for 2015 if there is an NA
#If rwi is trimmed, climate variables must be trimmed as well
DD5_sp <- DD5_sp[1:(length(DD5_sp)-1)]
DD5_sm <- DD5_sm[1:(length(DD5_sm)-1)]
RH_sp <- RH_sp[1:(length(RH_sp)-1)]
PAS <- PAS[1:(length(PAS)-1)]
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#all values should be the same length

#standardizing climate variables
DD5_sp<-(DD5_sp-mean(DD5_sp))/sd(DD5_sp)
DD5_sm<-(DD5_sm-mean(DD5_sm))/sd(DD5_sm)
RH_sp<-(RH_sp-mean(RH_sp))/sd(RH_sp)
PAS<-(PAS-mean(PAS))/sd(PAS)
prev_DD5_sm<-(prev_DD5_sm-mean(prev_DD5_sm))/sd(prev_DD5_sm)

require(MuMIn)
mod <- lm(rwi~DD5_sm+DD5_sp+RH_sp+PAS+prev_DD5_sm) #fit full model
options(na.action = "na.fail")
moddredge <- dredge(mod, extra = "R^2") # fits all nested models, creates a table of their coefficients and AICs, and R2s
#modavgs <- model.avg(moddredge, subset= delta <=2) #takes all models with deltaAIC<2 from the best, and averages their coefficients to get the 'best model estimate'
modavgs <- model.avg(moddredge, cumsum(weight) <= .95)
modavgs$importance

#PARA ABILAS
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Climate WMA data")
climate <- read.csv("PARA_ABILAS2_clim_1901-2015MSYT.csv", header=TRUE) #change filename for site/species

#load in climate variables from site/species file
DD5_sp <- climate$DD5_sp
DD5_sm <- climate$DD5_sm
RH_sp <- climate$RH_sp
PAS <- climate$PAS
prev_DD5_sm <- c(NA,DD5_sm)
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#removes the first value of all of the climate variables so that I can get rid of the NA at the beginning of prev_GDD_sm
DD5_sp <- DD5_sp[-1]
DD5_sm <- DD5_sm[-1]
RH_sp <- RH_sp[-1]
PAS <- PAS[-1]
prev_DD5_sm <- prev_DD5_sm[-1]

#now the climate variable time series go from 1902-2015

#RWI
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Master Chronologies 2017")
master <- read.csv(file="PARA_ABILAS2_master.csv",header=TRUE)
master_trim <- master[-c(1:22),] #shortens chronology from extending to 1880 to 1902
rwi <- master_trim$RWI
DD5_sp <- DD5_sp[1:(length(DD5_sp)-4)]
DD5_sm <- DD5_sm[1:(length(DD5_sm)-4)]
RH_sp <- RH_sp[1:(length(RH_sp)-4)]
PAS <- PAS[1:(length(PAS)-4)]
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-4)]
#all values should now be the same length


#standardizing climate variables
DD5_sp<-(DD5_sp-mean(DD5_sp))/sd(DD5_sp)
DD5_sm<-(DD5_sm-mean(DD5_sm))/sd(DD5_sm)
RH_sp<-(RH_sp-mean(RH_sp))/sd(RH_sp)
PAS<-(PAS-mean(PAS))/sd(PAS)
prev_DD5_sm<-(prev_DD5_sm-mean(prev_DD5_sm))/sd(prev_DD5_sm)

require(MuMIn)
mod <- lm(rwi~DD5_sm+DD5_sp+RH_sp+PAS+prev_DD5_sm) #fit full model
options(na.action = "na.fail")
moddredge <- dredge(mod, extra = "R^2") # fits all nested models, creates a table of their coefficients and AICs, and R2s
#modavgs <- model.avg(moddredge, subset= delta <=2) #takes all models with deltaAIC<2 from the best, and averages their coefficients to get the 'best model estimate'
modavgs <- model.avg(moddredge, cumsum(weight) <= .95)
modavgs$importance

#SPRY ABILAS
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Climate WMA data")
climate <- read.csv("SPRY_ABILAS_clim_1901-2015MSYT.csv", header=TRUE) #change filename for site/species

#load in climate variables from site/species file
DD5_sp <- climate$DD5_sp
DD5_sm <- climate$DD5_sm
RH_sp <- climate$RH_sp
PAS <- climate$PAS
prev_DD5_sm <- c(NA,DD5_sm)
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#removes the first value of all of the climate variables so that I can get rid of the NA at the beginning of prev_GDD_sm
DD5_sp <- DD5_sp[-1]
DD5_sm <- DD5_sm[-1]
RH_sp <- RH_sp[-1]
PAS <- PAS[-1]
prev_DD5_sm <- prev_DD5_sm[-1]

#now the climate variable time series go from 1902-2015

#RWI
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Master Chronologies 2017")
master <- read.csv(file="SPRY_ABILAS_master.csv",header=TRUE)
master_trim <- master[-c(1:22),] #shortens chronology from extending to 1880 to 1902
rwi <- master_trim$RWI

#all values should be the same length

#standardizing climate variables
DD5_sp<-(DD5_sp-mean(DD5_sp))/sd(DD5_sp)
DD5_sm<-(DD5_sm-mean(DD5_sm))/sd(DD5_sm)
RH_sp<-(RH_sp-mean(RH_sp))/sd(RH_sp)
PAS<-(PAS-mean(PAS))/sd(PAS)
prev_DD5_sm<-(prev_DD5_sm-mean(prev_DD5_sm))/sd(prev_DD5_sm)

require(MuMIn)
mod <- lm(rwi~DD5_sm+DD5_sp+RH_sp+PAS+prev_DD5_sm) #fit full model
options(na.action = "na.fail")
moddredge <- dredge(mod, extra = "R^2") # fits all nested models, creates a table of their coefficients and AICs, and R2s
#modavgs <- model.avg(moddredge, subset= delta <=2) #takes all models with deltaAIC<2 from the best, and averages their coefficients to get the 'best model estimate'
modavgs <- model.avg(moddredge, cumsum(weight) <= .95)
modavgs$importance

#SUNR ABILAS
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Climate WMA data")
climate <- read.csv("SUNR_ABILAS_clim_1901-2015MSYT.csv", header=TRUE) #change filename for site/species

#load in climate variables from site/species file
DD5_sp <- climate$DD5_sp
DD5_sm <- climate$DD5_sm
RH_sp <- climate$RH_sp
PAS <- climate$PAS
prev_DD5_sm <- c(NA,DD5_sm)
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#removes the first value of all of the climate variables so that I can get rid of the NA at the beginning of prev_GDD_sm
DD5_sp <- DD5_sp[-1]
DD5_sm <- DD5_sm[-1]
RH_sp <- RH_sp[-1]
PAS <- PAS[-1]
prev_DD5_sm <- prev_DD5_sm[-1]

#now the climate variable time series go from 1902-2015

#RWI
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Master Chronologies 2017")
master <- read.csv(file="SUNR_ABILAS_master.csv",header=TRUE)
master_trim <- master[-c(1:22),] #shortens chronology from extending to 1880 to 1902
rwi <- master_trim$RWI

#all values should be the same length

#standardizing climate variables
DD5_sp<-(DD5_sp-mean(DD5_sp))/sd(DD5_sp)
DD5_sm<-(DD5_sm-mean(DD5_sm))/sd(DD5_sm)
RH_sp<-(RH_sp-mean(RH_sp))/sd(RH_sp)
PAS<-(PAS-mean(PAS))/sd(PAS)
prev_DD5_sm<-(prev_DD5_sm-mean(prev_DD5_sm))/sd(prev_DD5_sm)

require(MuMIn)
mod <- lm(rwi~DD5_sm+DD5_sp+RH_sp+PAS+prev_DD5_sm) #fit full model
options(na.action = "na.fail")
moddredge <- dredge(mod, extra = "R^2") # fits all nested models, creates a table of their coefficients and AICs, and R2s
#modavgs <- model.avg(moddredge, subset= delta <=2) #takes all models with deltaAIC<2 from the best, and averages their coefficients to get the 'best model estimate'
modavgs <- model.avg(moddredge, cumsum(weight) <= .95)
modavgs$importance

#PARA TSUMER
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Climate WMA data")
climate <- read.csv("PARA_TSUMER_clim_1901-2015MSYT.csv", header=TRUE) #change filename for site/species

#load in climate variables from site/species file
DD5_sp <- climate$DD5_sp
DD5_sm <- climate$DD5_sm
RH_sp <- climate$RH_sp
PAS <- climate$PAS
prev_DD5_sm <- c(NA,DD5_sm)
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#removes the first value of all of the climate variables so that I can get rid of the NA at the beginning of prev_GDD_sm
DD5_sp <- DD5_sp[-1]
DD5_sm <- DD5_sm[-1]
RH_sp <- RH_sp[-1]
PAS <- PAS[-1]
prev_DD5_sm <- prev_DD5_sm[-1]

#now the climate variable time series go from 1902-2015

#RWI
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Master Chronologies 2017")
master <- read.csv(file="PARA_TSUMER_master.csv",header=TRUE)
master_trim <- master[-c(1:22),] #shortens chronology from extending to 1880 to 1902
rwi <- master_trim$RWI
rwi <- rwi[1:(length(rwi)-1)]#trim the value for 2015 if there is an NA
#If rwi is trimmed, climate variables must be trimmed as well
DD5_sp <- DD5_sp[1:(length(DD5_sp)-1)]
DD5_sm <- DD5_sm[1:(length(DD5_sm)-1)]
RH_sp <- RH_sp[1:(length(RH_sp)-1)]
PAS <- PAS[1:(length(PAS)-1)]
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]
#all values should be the same length

#standardizing climate variables
DD5_sp<-(DD5_sp-mean(DD5_sp))/sd(DD5_sp)
DD5_sm<-(DD5_sm-mean(DD5_sm))/sd(DD5_sm)
RH_sp<-(RH_sp-mean(RH_sp))/sd(RH_sp)
PAS<-(PAS-mean(PAS))/sd(PAS)
prev_DD5_sm<-(prev_DD5_sm-mean(prev_DD5_sm))/sd(prev_DD5_sm)

require(MuMIn)
mod <- lm(rwi~DD5_sm+DD5_sp+RH_sp+PAS+prev_DD5_sm) #fit full model
options(na.action = "na.fail")
moddredge <- dredge(mod, extra = "R^2") # fits all nested models, creates a table of their coefficients and AICs, and R2s
#modavgs <- model.avg(moddredge, subset= delta <=2) #takes all models with deltaAIC<2 from the best, and averages their coefficients to get the 'best model estimate'
modavgs <- model.avg(moddredge, cumsum(weight) <= .95)
modavgs$importance

#SPRY TSUMER
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Climate WMA data")
climate <- read.csv("SPRY_TSUMER_clim_1901-2015MSYT.csv", header=TRUE) #change filename for site/species

#load in climate variables from site/species file
DD5_sp <- climate$DD5_sp
DD5_sm <- climate$DD5_sm
RH_sp <- climate$RH_sp
PAS <- climate$PAS
prev_DD5_sm <- c(NA,DD5_sm)
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#removes the first value of all of the climate variables so that I can get rid of the NA at the beginning of prev_GDD_sm
DD5_sp <- DD5_sp[-1]
DD5_sm <- DD5_sm[-1]
RH_sp <- RH_sp[-1]
PAS <- PAS[-1]
prev_DD5_sm <- prev_DD5_sm[-1]

#now the climate variable time series go from 1902-2015

#RWI
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Master Chronologies 2017")
master <- read.csv(file="SPRY_TSUMER_master.csv",header=TRUE)
master_trim <- master[-c(1:22),] #shortens chronology from extending to 1880 to 1902
rwi <- master_trim$RWI
rwi <- rwi[1:(length(rwi)-1)]#trim the value for 2015 if there is an NA
#If rwi is trimmed, climate variables must be trimmed as well
DD5_sp <- DD5_sp[1:(length(DD5_sp)-1)]
DD5_sm <- DD5_sm[1:(length(DD5_sm)-1)]
RH_sp <- RH_sp[1:(length(RH_sp)-1)]
PAS <- PAS[1:(length(PAS)-1)]
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]
#all values should be the same length

#standardizing climate variables
DD5_sp<-(DD5_sp-mean(DD5_sp))/sd(DD5_sp)
DD5_sm<-(DD5_sm-mean(DD5_sm))/sd(DD5_sm)
RH_sp<-(RH_sp-mean(RH_sp))/sd(RH_sp)
PAS<-(PAS-mean(PAS))/sd(PAS)
prev_DD5_sm<-(prev_DD5_sm-mean(prev_DD5_sm))/sd(prev_DD5_sm)

require(MuMIn)
mod <- lm(rwi~DD5_sm+DD5_sp+RH_sp+PAS+prev_DD5_sm) #fit full model
options(na.action = "na.fail")
moddredge <- dredge(mod, extra = "R^2") # fits all nested models, creates a table of their coefficients and AICs, and R2s
#modavgs <- model.avg(moddredge, subset= delta <=2) #takes all models with deltaAIC<2 from the best, and averages their coefficients to get the 'best model estimate'
modavgs <- model.avg(moddredge, cumsum(weight) <= .95)
modavgs$importance

#SUNR TSUMER
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Climate WMA data")
climate <- read.csv("SUNR_TSUMER_clim_1901-2015MSYT.csv", header=TRUE) #change filename for site/species

#load in climate variables from site/species file
DD5_sp <- climate$DD5_sp
DD5_sm <- climate$DD5_sm
RH_sp <- climate$RH_sp
PAS <- climate$PAS
prev_DD5_sm <- c(NA,DD5_sm)
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#removes the first value of all of the climate variables so that I can get rid of the NA at the beginning of prev_GDD_sm
DD5_sp <- DD5_sp[-1]
DD5_sm <- DD5_sm[-1]
RH_sp <- RH_sp[-1]
PAS <- PAS[-1]
prev_DD5_sm <- prev_DD5_sm[-1]

#now the climate variable time series go from 1902-2015

#RWI
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Master Chronologies 2017")
master <- read.csv(file="SUNR_TSUMER_master.csv",header=TRUE)
master_trim <- master[-c(1:22),] #shortens chronology from extending to 1880 to 1902
rwi <- master_trim$RWI

#standardizing climate variables
DD5_sp<-(DD5_sp-mean(DD5_sp))/sd(DD5_sp)
DD5_sm<-(DD5_sm-mean(DD5_sm))/sd(DD5_sm)
RH_sp<-(RH_sp-mean(RH_sp))/sd(RH_sp)
PAS<-(PAS-mean(PAS))/sd(PAS)
prev_DD5_sm<-(prev_DD5_sm-mean(prev_DD5_sm))/sd(prev_DD5_sm)

require(MuMIn)
mod <- lm(rwi~DD5_sm+DD5_sp+RH_sp+PAS+prev_DD5_sm) #fit full model
options(na.action = "na.fail")
moddredge <- dredge(mod, extra = "R^2") # fits all nested models, creates a table of their coefficients and AICs, and R2s
#modavgs <- model.avg(moddredge, subset= delta <=2) #takes all models with deltaAIC<2 from the best, and averages their coefficients to get the 'best model estimate'
modavgs <- model.avg(moddredge, cumsum(weight) <= .95)
modavgs$importance

#PARA XANNOO
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Climate WMA data")
climate <- read.csv("PARA_CALNOO_clim_1901-2015MSYT.csv", header=TRUE) #change filename for site/species

#load in climate variables from site/species file
DD5_sp <- climate$DD5_sp
DD5_sm <- climate$DD5_sm
RH_sp <- climate$RH_sp
PAS <- climate$PAS
prev_DD5_sm <- c(NA,DD5_sm)
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#removes the first value of all of the climate variables so that I can get rid of the NA at the beginning of prev_GDD_sm
DD5_sp <- DD5_sp[-1]
DD5_sm <- DD5_sm[-1]
RH_sp <- RH_sp[-1]
PAS <- PAS[-1]
prev_DD5_sm <- prev_DD5_sm[-1]

#now the climate variable time series go from 1902-2015

#RWI
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Master Chronologies 2017")
master <- read.csv(file="PARA_CALNOO_master.csv",header=TRUE)
master_trim <- master[-c(1:22),] #shortens chronology from extending to 1880 to 1902
rwi <- master_trim$RWI
rwi <- rwi[1:(length(rwi)-1)]#trim the value for 2015 if there is an NA
#If rwi is trimmed, climate variables must be trimmed as well
DD5_sp <- DD5_sp[1:(length(DD5_sp)-1)]
DD5_sm <- DD5_sm[1:(length(DD5_sm)-1)]
RH_sp <- RH_sp[1:(length(RH_sp)-1)]
PAS <- PAS[1:(length(PAS)-1)]
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]
#all values should be the same length

#standardizing climate variables
DD5_sp<-(DD5_sp-mean(DD5_sp))/sd(DD5_sp)
DD5_sm<-(DD5_sm-mean(DD5_sm))/sd(DD5_sm)
RH_sp<-(RH_sp-mean(RH_sp))/sd(RH_sp)
PAS<-(PAS-mean(PAS))/sd(PAS)
prev_DD5_sm<-(prev_DD5_sm-mean(prev_DD5_sm))/sd(prev_DD5_sm)

require(MuMIn)
mod <- lm(rwi~DD5_sm+DD5_sp+RH_sp+PAS+prev_DD5_sm) #fit full model
options(na.action = "na.fail")
moddredge <- dredge(mod, extra = "R^2") # fits all nested models, creates a table of their coefficients and AICs, and R2s
#modavgs <- model.avg(moddredge, subset= delta <=2) #takes all models with deltaAIC<2 from the best, and averages their coefficients to get the 'best model estimate'
modavgs <- model.avg(moddredge, cumsum(weight) <= .95)
modavgs$importance

#SPRY XANNOO
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Climate WMA data")
climate <- read.csv("SPRY_CALNOO_clim_1901-2015MSYT.csv", header=TRUE) #change filename for site/species

#load in climate variables from site/species file
DD5_sp <- climate$DD5_sp
DD5_sm <- climate$DD5_sm
RH_sp <- climate$RH_sp
PAS <- climate$PAS
prev_DD5_sm <- c(NA,DD5_sm)
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#removes the first value of all of the climate variables so that I can get rid of the NA at the beginning of prev_GDD_sm
DD5_sp <- DD5_sp[-1]
DD5_sm <- DD5_sm[-1]
RH_sp <- RH_sp[-1]
PAS <- PAS[-1]
prev_DD5_sm <- prev_DD5_sm[-1]

#now the climate variable time series go from 1902-2015

#RWI
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Master Chronologies 2017")
master <- read.csv(file="SPRY_CALNOO_master.csv",header=TRUE)
master_trim <- master[-c(1:22),] #shortens chronology from extending to 1880 to 1902
rwi <- master_trim$RWI

#standardizing climate variables
DD5_sp<-(DD5_sp-mean(DD5_sp))/sd(DD5_sp)
DD5_sm<-(DD5_sm-mean(DD5_sm))/sd(DD5_sm)
RH_sp<-(RH_sp-mean(RH_sp))/sd(RH_sp)
PAS<-(PAS-mean(PAS))/sd(PAS)
prev_DD5_sm<-(prev_DD5_sm-mean(prev_DD5_sm))/sd(prev_DD5_sm)

require(MuMIn)
mod <- lm(rwi~DD5_sm+DD5_sp+RH_sp+PAS+prev_DD5_sm) #fit full model
options(na.action = "na.fail")
moddredge <- dredge(mod, extra = "R^2") # fits all nested models, creates a table of their coefficients and AICs, and R2s
#modavgs <- model.avg(moddredge, subset= delta <=2) #takes all models with deltaAIC<2 from the best, and averages their coefficients to get the 'best model estimate'
modavgs <- model.avg(moddredge, cumsum(weight) <= .95)
modavgs$importance

#SUNR XANNOO
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Climate WMA data")
climate <- read.csv("SUNR_CALNOO_clim_1901-2015MSYT.csv", header=TRUE) #change filename for site/species

#load in climate variables from site/species file
DD5_sp <- climate$DD5_sp
DD5_sm <- climate$DD5_sm
RH_sp <- climate$RH_sp
PAS <- climate$PAS
prev_DD5_sm <- c(NA,DD5_sm)
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]

#removes the first value of all of the climate variables so that I can get rid of the NA at the beginning of prev_GDD_sm
DD5_sp <- DD5_sp[-1]
DD5_sm <- DD5_sm[-1]
RH_sp <- RH_sp[-1]
PAS <- PAS[-1]
prev_DD5_sm <- prev_DD5_sm[-1]

#now the climate variable time series go from 1902-2015

#RWI
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Master Chronologies 2017")
master <- read.csv(file="SUNR_CALNOO_master.csv",header=TRUE)
master_trim <- master[-c(1:22),] #shortens chronology from extending to 1880 to 1902
rwi <- master_trim$RWI
rwi <- rwi[1:(length(rwi)-1)]#trim the value for 2015 if there is an NA
#If rwi is trimmed, climate variables must be trimmed as well
DD5_sp <- DD5_sp[1:(length(DD5_sp)-1)]
DD5_sm <- DD5_sm[1:(length(DD5_sm)-1)]
RH_sp <- RH_sp[1:(length(RH_sp)-1)]
PAS <- PAS[1:(length(PAS)-1)]
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]
#all values should be the same length

#standardizing climate variables
DD5_sp<-(DD5_sp-mean(DD5_sp))/sd(DD5_sp)
DD5_sm<-(DD5_sm-mean(DD5_sm))/sd(DD5_sm)
RH_sp<-(RH_sp-mean(RH_sp))/sd(RH_sp)
PAS<-(PAS-mean(PAS))/sd(PAS)
prev_DD5_sm<-(prev_DD5_sm-mean(prev_DD5_sm))/sd(prev_DD5_sm)

require(MuMIn)
mod <- lm(rwi~DD5_sm+DD5_sp+RH_sp+PAS+prev_DD5_sm) #fit full model
options(na.action = "na.fail")
moddredge <- dredge(mod, extra = "R^2") # fits all nested models, creates a table of their coefficients and AICs, and R2s
#modavgs <- model.avg(moddredge, subset= delta <=2) #takes all models with deltaAIC<2 from the best, and averages their coefficients to get the 'best model estimate'
modavgs <- model.avg(moddredge, cumsum(weight) <= .95)
modavgs$importance