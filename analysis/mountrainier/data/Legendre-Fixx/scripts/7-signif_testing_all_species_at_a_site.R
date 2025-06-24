###===================================================================================
#Organizing datasets to compare for sensitivity differences for all species at a site
#=======================================================================================

##PARA
setwd("C:/Users/Myesa/Desktop/e-output_from_WNA")
climate <- read.csv("e-PARA_ABILAS_clim_1901-2015MSYT.csv", header=TRUE) 
library(graphics)
library(utils)

climate$sp <- "ABILAS"
climate$site <- "PARA"

DD5_sp <- climate$DD5_sp
DD5_sm <- climate$DD5_sm
RH_sp <- climate$RH_sp
PAS <- climate$PAS
prev_DD5_sm <- c(NA,DD5_sm)
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]
sp <- climate$sp
site <- climate$site

#remove the first value of all of the climate variables to remove the NA at the beginning of prev_GDD_sm
DD5_sp <- DD5_sp[-1]
DD5_sm <- DD5_sm[-1]
RH_sp <- RH_sp[-1]
PAS <- PAS[-1]
prev_DD5_sm <- prev_DD5_sm[-1]
sp <- sp[-1]
site <- site[-1]
#now the climate variable time series go from 1902-2015

#RWI
setwd("C:/Users/Myesa/Desktop/d-master_chronologies")
master <- read.csv(file="d-PARA_ABILAS_master.csv",header=TRUE)
master_trim <- master[-c(1:22),]
rwi <- master_trim$RWI
#If rwi vector is shorter than 2015, then climate variable vectors should be trimmed
DD5_sp <- DD5_sp[1:(length(DD5_sp)-4)]
DD5_sm <- DD5_sm[1:(length(DD5_sm)-4)]
RH_sp <- RH_sp[1:(length(RH_sp)-4)]
PAS <- PAS[1:(length(PAS)-4)]
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-4)]
sp <- sp[1:(length(sp)-4)]
site <- site[1:(length(site)-4)]

#standardizing climate variables
DD5_sp<-(DD5_sp-mean(DD5_sp))/sd(DD5_sp)
DD5_sm<-(DD5_sm-mean(DD5_sm))/sd(DD5_sm)
RH_sp<-(RH_sp-mean(RH_sp))/sd(RH_sp)
PAS<-(PAS-mean(PAS))/sd(PAS)
prev_DD5_sm<-(prev_DD5_sm-mean(prev_DD5_sm))/sd(prev_DD5_sm)

#PARA ABIAMA
setwd("C:/Users/Myesa/Desktop/e-output_from_WNA")
climate2 <- read.csv("e-PARA_ABIAMA_clim_1901-2015MSYT.csv", header=TRUE) 
climate2$sp <- "ABIAMA"
climate2$site <- "PARA"

DD5_sp2 <- climate2$DD5_sp
DD5_sm2 <- climate2$DD5_sm
RH_sp2 <- climate2$RH_sp
PAS2 <- climate2$PAS
prev_DD5_sm2 <- c(NA,DD5_sm2)
prev_DD5_sm2 <- prev_DD5_sm2[1:(length(prev_DD5_sm2)-1)]
sp2 <-climate2$sp
site2 <- climate2$site

#remove the first value of all of the climate variables to remove the NA at the beginning of prev_GDD_sm
DD5_sp2 <- DD5_sp2[-1]
DD5_sm2 <- DD5_sm2[-1]
RH_sp2 <- RH_sp2[-1]
PAS2 <- PAS2[-1]
prev_DD5_sm2 <- prev_DD5_sm2[-1]
sp2 <- sp2[-1]
site2 <- site2[-1]
#now the climate variable time series go from 1902-2015

#RWI
setwd("C:/Users/Myesa/Desktop/d-master_chronologies")
master2 <- read.csv(file="d-PARA_ABIAMA_master.csv",header=TRUE)
master_trim2 <- master2[-c(1:22),]
rwi2 <- master_trim2$RWI
rwi2 <- rwi2[1:(length(rwi2)-1)]

DD5_sp2 <- DD5_sp2[1:(length(DD5_sp2)-1)]
DD5_sm2 <- DD5_sm2[1:(length(DD5_sm2)-1)]
RH_sp2 <- RH_sp2[1:(length(RH_sp2)-1)]
PAS2 <- PAS2[1:(length(PAS2)-1)]
prev_DD5_sm2 <- prev_DD5_sm2[1:(length(prev_DD5_sm2)-1)]
sp2 <- sp2[1:(length(sp2)-1)]
site2 <- site2[1:(length(site2)-1)]

#standardizing climate variables
DD5_sp2<-(DD5_sp2-mean(DD5_sp2))/sd(DD5_sp2)
DD5_sm2<-(DD5_sm2-mean(DD5_sm2))/sd(DD5_sm2)
RH_sp2<-(RH_sp2-mean(RH_sp2))/sd(RH_sp2)
PAS2<-(PAS2-mean(PAS2))/sd(PAS2)
prev_DD5_sm2<-(prev_DD5_sm2-mean(prev_DD5_sm2))/sd(prev_DD5_sm2)

##PARA TSUMER
setwd("C:/Users/Myesa/Desktop/e-output_from_WNA")
climate3 <- read.csv("e-PARA_TSUMER_clim_1901-2015MSYT.csv", header=TRUE) 
climate3$sp <- "TSUMER"
climate3$site <- "PARA"

DD5_sp3 <- climate3$DD5_sp
DD5_sm3 <- climate3$DD5_sm
RH_sp3 <- climate3$RH_sp
PAS3 <- climate3$PAS
prev_DD5_sm3 <- c(NA,DD5_sm3)
prev_DD5_sm3 <- prev_DD5_sm3[1:(length(prev_DD5_sm3)-1)]
sp3 <-climate3$sp
site3 <-climate3$site

#remove the first value of all of the climate variables to remove the NA at the beginning of prev_GDD_sm
DD5_sp3 <- DD5_sp3[-1]
DD5_sm3 <- DD5_sm3[-1]
RH_sp3 <- RH_sp3[-1]
PAS3 <- PAS3[-1]
prev_DD5_sm3 <- prev_DD5_sm3[-1]
sp3 <- sp3[-1]
site3 <- site3[-1]
#now the climate variable time series go from 1902-2015

#RWI
setwd("C:/Users/Myesa/Desktop/d-master_chronologies")
master3 <- read.csv(file="d-PARA_TSUMER_master.csv",header=TRUE)
master_trim3 <- master3[-c(1:22),] 
rwi3 <- master_trim3$RWI
rwi3 <- rwi3[1:(length(rwi3)-1)]#trim the value for 2015 if there is an NA (the multiple linear regression will not take an NA)
DD5_sp3 <- DD5_sp3[1:(length(DD5_sp3)-1)]
DD5_sm3 <- DD5_sm3[1:(length(DD5_sm3)-1)]
RH_sp3 <- RH_sp3[1:(length(RH_sp3)-1)]
PAS3 <- PAS3[1:(length(PAS3)-1)]
prev_DD5_sm3 <- prev_DD5_sm3[1:(length(prev_DD5_sm3)-1)]
sp3 <- sp3[1:(length(sp3)-1)]
site3 <- site3[1:(length(site3)-1)]

#standardizing climate variables
DD5_sp3<-(DD5_sp3-mean(DD5_sp3))/sd(DD5_sp3)
DD5_sm3<-(DD5_sm3-mean(DD5_sm3))/sd(DD5_sm3)
RH_sp3<-(RH_sp3-mean(RH_sp3))/sd(RH_sp3)
PAS3<-(PAS3-mean(PAS3))/sd(PAS3)
prev_DD5_sm3<-(prev_DD5_sm3-mean(prev_DD5_sm3))/sd(prev_DD5_sm3)

#PARA CALNOO
setwd("C:/Users/Myesa/Desktop/e-output_from_WNA")
climate4 <- read.csv("e-PARA_CALNOO_clim_1901-2015MSYT.csv", header=TRUE)
climate4$sp <- "CALNOO"
climate4$site <- "PARA"

DD5_sp4 <- climate4$DD5_sp
DD5_sm4 <- climate4$DD5_sm
RH_sp4 <- climate4$RH_sp
PAS4 <- climate4$PAS
prev_DD5_sm4 <- c(NA,DD5_sm4)
prev_DD5_sm4 <- prev_DD5_sm4[1:(length(prev_DD5_sm4)-1)]
sp4 <-climate4$sp
site4 <-climate4$site

DD5_sp4 <- DD5_sp4[-1]
DD5_sm4 <- DD5_sm4[-1]
RH_sp4 <- RH_sp4[-1]
PAS4 <- PAS4[-1]
prev_DD5_sm4 <- prev_DD5_sm4[-1]
sp4 <- sp4[-1]
site4 <- site4[-1]
#now the climate variable time series go from 1902-2015
setwd("C:/Users/Myesa/Desktop/d-master_chronologies")
master4 <- read.csv(file="d-PARA_CALNOO_master.csv",header=TRUE)
master_trim4 <- master4[-c(1:22),] 
#should trim file to go from 1902-2015 
rwi4 <- master_trim4$RWI
rwi4 <- rwi4[1:(length(rwi4)-1)]#trim the value for 2015 if there is an NA (the multiple linear regression will not take an NA)
DD5_sp4 <- DD5_sp4[1:(length(DD5_sp4)-1)]
DD5_sm4 <- DD5_sm4[1:(length(DD5_sm4)-1)]
RH_sp4 <- RH_sp4[1:(length(RH_sp4)-1)]
PAS4 <- PAS4[1:(length(PAS4)-1)]
prev_DD5_sm4 <- prev_DD5_sm4[1:(length(prev_DD5_sm4)-1)]
sp4 <- sp4[1:(length(sp4)-1)]
site4 <- site4[1:(length(site4)-1)]

#standardizing climate variables
DD5_sp4<-(DD5_sp4-mean(DD5_sp4))/sd(DD5_sp4)
DD5_sm4<-(DD5_sm4-mean(DD5_sm4))/sd(DD5_sm4)
RH_sp4<-(RH_sp4-mean(RH_sp4))/sd(RH_sp4)
PAS4<-(PAS4-mean(PAS4))/sd(PAS4)
prev_DD5_sm4<-(prev_DD5_sm4-mean(prev_DD5_sm4))/sd(prev_DD5_sm4)


#append all together
DD5_sp<-append(DD5_sp,DD5_sp2)
DD5_sp<-append(DD5_sp,DD5_sp3)
DD5_sp<-append(DD5_sp,DD5_sp4)

DD5_sm<-append(DD5_sm,DD5_sm2)
DD5_sm<-append(DD5_sm,DD5_sm3)
DD5_sm<-append(DD5_sm,DD5_sm4)

RH_sp<-append(RH_sp,RH_sp2)
RH_sp<-append(RH_sp,RH_sp3)
RH_sp<-append(RH_sp,RH_sp4)

PAS<-append(PAS,PAS2)
PAS<-append(PAS,PAS3)
PAS<-append(PAS,PAS4)

prev_DD5_sm<-append(prev_DD5_sm,prev_DD5_sm2)
prev_DD5_sm<-append(prev_DD5_sm,prev_DD5_sm3)
prev_DD5_sm<-append(prev_DD5_sm,prev_DD5_sm4)

rwi<-append(rwi,rwi2)
rwi<-append(rwi,rwi3)
rwi<-append(rwi,rwi4)

sp<-append(sp,sp2)
sp<-append(sp,sp3)
sp<-append(sp,sp4)

site<-append(site,site2)
site<-append(site,site3)
site<-append(site,site4)

#this dataset has climate variables, rw, site and species id for all four species at PARA (S-side)
alldata_PARA <- data.frame(DD5_sp,DD5_sm,RH_sp,PAS,prev_DD5_sm,rwi,sp,site)

#Making a SPRY (NW-side dataset)
setwd("C:/Users/Myesa/Desktop/e-output_from_WNA")
climate<- read.csv("e-SPRY_ABILAS_clim_1901-2015MSYT.csv", header=TRUE)
climate$sp <- "ABILAS"
climate$site <- "SPRY"

DD5_sp <- climate$DD5_sp
DD5_sm <- climate$DD5_sm
RH_sp <- climate$RH_sp
PAS <- climate$PAS
prev_DD5_sm <- c(NA,DD5_sm)
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]
sp <-climate$sp
site <-climate$site

DD5_sp <- DD5_sp[-1]
DD5_sm <- DD5_sm[-1]
RH_sp <- RH_sp[-1]
PAS <- PAS[-1]
prev_DD5_sm <- prev_DD5_sm[-1]
sp <- sp[-1]
site <- site[-1]
#now the climate variable time series go from 1902-2015

setwd("C:/Users/Myesa/Desktop/d-master_chronologies")
master <- read.csv(file="d-SPRY_ABILAS_master.csv",header=TRUE)
master_trim <- master[-c(1:22),] 
rwi <- master_trim$RWI

#standardizing climate variables 
DD5_sp<-(DD5_sp-mean(DD5_sp))/sd(DD5_sp)
DD5_sm<-(DD5_sm-mean(DD5_sm))/sd(DD5_sm)
RH_sp<-(RH_sp-mean(RH_sp))/sd(RH_sp)
PAS<-(PAS-mean(PAS))/sd(PAS)
prev_DD5_sm<-(prev_DD5_sm-mean(prev_DD5_sm))/sd(prev_DD5_sm)

#SPRY ABIAMA
setwd("C:/Users/Myesa/Desktop/e-output_from_WNA")
climate2 <- read.csv("e-SPRY_ABIAMA_clim_1901-2015MSYT.csv", header=TRUE)

climate2$sp <- "ABIAMA"
climate2$site <- "SPRY"

DD5_sp2 <- climate2$DD5_sp
DD5_sm2 <- climate2$DD5_sm
RH_sp2 <- climate2$RH_sp
PAS2 <- climate2$PAS
prev_DD5_sm2 <- c(NA,DD5_sm2)
prev_DD5_sm2 <- prev_DD5_sm2[1:(length(prev_DD5_sm2)-1)]
sp2 <-climate2$sp
site2 <- climate2$site

DD5_sp2 <- DD5_sp2[-1]
DD5_sm2 <- DD5_sm2[-1]
RH_sp2 <- RH_sp2[-1]
PAS2 <- PAS2[-1]
prev_DD5_sm2 <- prev_DD5_sm2[-1]
sp2 <- sp2[-1]
site2 <- site2[-1]
#now the climate variable time series go from 1902-2015
setwd("C:/Users/Myesa/Desktop/d-master_chronologies")
master2 <- read.csv(file="d-SPRY_ABIAMA_master.csv",header=TRUE)
master_trim2 <- master2[-c(1:22),] 
#should trim file to go from 1902-2015 
rwi2 <- master_trim2$RWI
rwi2 <- rwi2[1:(length(rwi2)-1)]
DD5_sp2 <- DD5_sp2[1:(length(DD5_sp2)-1)]
DD5_sm2 <- DD5_sm2[1:(length(DD5_sm2)-1)]
RH_sp2 <- RH_sp2[1:(length(RH_sp2)-1)]
PAS2 <- PAS2[1:(length(PAS2)-1)]
prev_DD5_sm2 <- prev_DD5_sm2[1:(length(prev_DD5_sm2)-1)]
sp2 <- sp2[1:(length(sp2)-1)]
site2 <- site2[1:(length(site2)-1)]

#standardizing climate variables 
DD5_sp2<-(DD5_sp2-mean(DD5_sp2))/sd(DD5_sp2)
DD5_sm2<-(DD5_sm2-mean(DD5_sm2))/sd(DD5_sm2)
RH_sp2<-(RH_sp2-mean(RH_sp2))/sd(RH_sp2)
PAS2<-(PAS2-mean(PAS2))/sd(PAS2)
prev_DD5_sm2<-(prev_DD5_sm2-mean(prev_DD5_sm2))/sd(prev_DD5_sm2)

##SPRY TSUMER
setwd("C:/Users/Myesa/Desktop/e-output_from_WNA")
climate3 <- read.csv("e-SPRY_TSUMER_clim_1901-2015MSYT.csv", header=TRUE)
climate3$sp <- "TSUMER"
climate3$site <- "SPRY"

DD5_sp3 <- climate3$DD5_sp
DD5_sm3 <- climate3$DD5_sm
RH_sp3 <- climate3$RH_sp
PAS3 <- climate3$PAS
prev_DD5_sm3 <- c(NA,DD5_sm3)
prev_DD5_sm3 <- prev_DD5_sm3[1:(length(prev_DD5_sm3)-1)]
sp3 <-climate3$sp
site3 <-climate3$site

DD5_sp3 <- DD5_sp3[-1]
DD5_sm3 <- DD5_sm3[-1]
RH_sp3 <- RH_sp3[-1]
PAS3 <- PAS3[-1]
prev_DD5_sm3 <- prev_DD5_sm3[-1]
sp3 <- sp3[-1]
site3 <- site3[-1]
#now the climate variable time series go from 1902-2015

setwd("C:/Users/Myesa/Desktop/d-master_chronologies")
master3 <- read.csv(file="d-SPRY_TSUMER_master.csv",header=TRUE)
master_trim3 <- master3[-c(1:22),] 
rwi3 <- master_trim3$RWI
rwi3 <- rwi3[1:(length(rwi3)-1)]#trim the value for 2015 if there is an NA (the multiple linear regression will not take an NA)
DD5_sp3 <- DD5_sp3[1:(length(DD5_sp3)-1)]
DD5_sm3 <- DD5_sm3[1:(length(DD5_sm3)-1)]
RH_sp3 <- RH_sp3[1:(length(RH_sp3)-1)]
PAS3 <- PAS3[1:(length(PAS3)-1)]
prev_DD5_sm3 <- prev_DD5_sm3[1:(length(prev_DD5_sm3)-1)]
sp3 <- sp3[1:(length(sp3)-1)]
site3 <- site3[1:(length(site3)-1)]
#standardizing climate variables 
DD5_sp3<-(DD5_sp3-mean(DD5_sp3))/sd(DD5_sp3)
DD5_sm3<-(DD5_sm3-mean(DD5_sm3))/sd(DD5_sm3)
RH_sp3<-(RH_sp3-mean(RH_sp3))/sd(RH_sp3)
PAS3<-(PAS3-mean(PAS3))/sd(PAS3)
prev_DD5_sm3<-(prev_DD5_sm3-mean(prev_DD5_sm3))/sd(prev_DD5_sm3)

#SPRY CALNOO
setwd("C:/Users/Myesa/Desktop/e-output_from_WNA")
climate4 <- read.csv("e-SPRY_CALNOO_clim_1901-2015MSYT.csv", header=TRUE) 
climate4$sp <- "CALNOO"
climate4$site <- "SPRY"

DD5_sp4 <- climate4$DD5_sp
DD5_sm4 <- climate4$DD5_sm
RH_sp4 <- climate4$RH_sp
PAS4 <- climate4$PAS
prev_DD5_sm4 <- c(NA,DD5_sm4)
prev_DD5_sm4 <- prev_DD5_sm4[1:(length(prev_DD5_sm4)-1)]
sp4 <-climate4$sp
site4 <-climate4$site

DD5_sp4 <- DD5_sp4[-1]
DD5_sm4 <- DD5_sm4[-1]
RH_sp4 <- RH_sp4[-1]
PAS4 <- PAS4[-1]
prev_DD5_sm4 <- prev_DD5_sm4[-1]
sp4 <- sp4[-1]
site4 <- site4[-1]
#now the climate variable time series go from 1902-2015

setwd("C:/Users/Myesa/Desktop/d-master_chronologies")
master4 <- read.csv(file="d-SPRY_CALNOO_master.csv",header=TRUE)
master_trim4 <- master4[-c(1:22),] 
rwi4 <- master_trim4$RWI

#standardizing climate variables 
DD5_sp4<-(DD5_sp4-mean(DD5_sp4))/sd(DD5_sp4)
DD5_sm4<-(DD5_sm4-mean(DD5_sm4))/sd(DD5_sm4)
RH_sp4<-(RH_sp4-mean(RH_sp4))/sd(RH_sp4)
PAS4<-(PAS4-mean(PAS4))/sd(PAS4)
prev_DD5_sm4<-(prev_DD5_sm4-mean(prev_DD5_sm4))/sd(prev_DD5_sm4)

#append all together
DD5_sp<-append(DD5_sp,DD5_sp2)
DD5_sp<-append(DD5_sp,DD5_sp3)
DD5_sp<-append(DD5_sp,DD5_sp4)

DD5_sm<-append(DD5_sm,DD5_sm2)
DD5_sm<-append(DD5_sm,DD5_sm3)
DD5_sm<-append(DD5_sm,DD5_sm4)

RH_sp<-append(RH_sp,RH_sp2)
RH_sp<-append(RH_sp,RH_sp3)
RH_sp<-append(RH_sp,RH_sp4)

PAS<-append(PAS,PAS2)
PAS<-append(PAS,PAS3)
PAS<-append(PAS,PAS4)

prev_DD5_sm<-append(prev_DD5_sm,prev_DD5_sm2)
prev_DD5_sm<-append(prev_DD5_sm,prev_DD5_sm3)
prev_DD5_sm<-append(prev_DD5_sm,prev_DD5_sm4)

rwi<-append(rwi,rwi2)
rwi<-append(rwi,rwi3)
rwi<-append(rwi,rwi4)

sp<-append(sp,sp2)
sp<-append(sp,sp3)
sp<-append(sp,sp4)

site<-append(site,site2)
site<-append(site,site3)
site<-append(site,site4)

#this dataset has climate variables, rw, site and species id for all four species at SPRY (NW-side)
alldata_SPRY <- data.frame(DD5_sp,DD5_sm,RH_sp,PAS,prev_DD5_sm,rwi,sp,site)

#Making a SUNR (E-side dataset)

setwd("C:/Users/Myesa/Desktop/e-output_from_WNA")
climate <- read.csv("e-SUNR_ABILAS_clim_1901-2015MSYT.csv", header=TRUE) 
climate$sp <- "ABILAS"
climate$site <- "SUNR"

#load in climate variables from site/species file
DD5_sp <- climate$DD5_sp
DD5_sm <- climate$DD5_sm
RH_sp <- climate$RH_sp
PAS <- climate$PAS
prev_DD5_sm <- c(NA,DD5_sm)
prev_DD5_sm <- prev_DD5_sm[1:(length(prev_DD5_sm)-1)]
sp <-climate$sp
site <-climate$site

DD5_sp <- DD5_sp[-1]
DD5_sm <- DD5_sm[-1]
RH_sp <- RH_sp[-1]
PAS <- PAS[-1]
prev_DD5_sm <- prev_DD5_sm[-1]
sp <- sp[-1]
site <- site[-1]
#now the climate variable time series go from 1902-2015

setwd("C:/Users/Myesa/Desktop/d-master_chronologies")
master <- read.csv(file="d-SUNR_ABILAS_master.csv",header=TRUE)
master_trim <- master[-c(1:22),] 
rwi <- master_trim$RWI

#standardizing climate variables 
DD5_sp<-(DD5_sp-mean(DD5_sp))/sd(DD5_sp)
DD5_sm<-(DD5_sm-mean(DD5_sm))/sd(DD5_sm)
RH_sp<-(RH_sp-mean(RH_sp))/sd(RH_sp)
PAS<-(PAS-mean(PAS))/sd(PAS)
prev_DD5_sm<-(prev_DD5_sm-mean(prev_DD5_sm))/sd(prev_DD5_sm)

#SUNR ABIAMA
setwd("C:/Users/Myesa/Desktop/e-output_from_WNA")
climate2 <- read.csv("e-SUNR_ABIAMA_clim_1901-2015MSYT.csv", header=TRUE)

climate2$sp <- "ABIAMA"
climate2$site <- "SUNR"

#load in climate variables from site/species file
DD5_sp2 <- climate2$DD5_sp
DD5_sm2 <- climate2$DD5_sm
RH_sp2 <- climate2$RH_sp
PAS2 <- climate2$PAS
prev_DD5_sm2 <- c(NA,DD5_sm2)
prev_DD5_sm2 <- prev_DD5_sm2[1:(length(prev_DD5_sm2)-1)]
sp2 <-climate2$sp
site2 <-climate2$site

DD5_sp2 <- DD5_sp2[-1]
DD5_sm2 <- DD5_sm2[-1]
RH_sp2 <- RH_sp2[-1]
PAS2 <- PAS2[-1]
prev_DD5_sm2 <- prev_DD5_sm2[-1]
sp2 <- sp2[-1]
site2 <- site2[-1]
#now the climate variable time series go from 1902-2015

setwd("C:/Users/Myesa/Desktop/d-master_chronologies")
master2 <- read.csv(file="d-SUNR_ABIAMA_master.csv",header=TRUE)
master_trim2 <- master2[-c(1:22),] 
rwi2 <- master_trim2$RWI
rwi2 <- rwi2[1:(length(rwi2)-1)]

DD5_sp2 <- DD5_sp2[1:(length(DD5_sp2)-1)]
DD5_sm2 <- DD5_sm2[1:(length(DD5_sm2)-1)]
RH_sp2 <- RH_sp2[1:(length(RH_sp2)-1)]
PAS2 <- PAS2[1:(length(PAS2)-1)]
prev_DD5_sm2 <- prev_DD5_sm2[1:(length(prev_DD5_sm2)-1)]
sp2 <- sp2[1:(length(sp2)-1)]
site2 <- site2[1:(length(site2)-1)]
#standardizing climate variables 
DD5_sp2<-(DD5_sp2-mean(DD5_sp2))/sd(DD5_sp2)
DD5_sm2<-(DD5_sm2-mean(DD5_sm2))/sd(DD5_sm2)
RH_sp2<-(RH_sp2-mean(RH_sp2))/sd(RH_sp2)
PAS2<-(PAS2-mean(PAS2))/sd(PAS2)
prev_DD5_sm2<-(prev_DD5_sm2-mean(prev_DD5_sm2))/sd(prev_DD5_sm2)

##SUNR TSUMER
setwd("C:/Users/Myesa/Desktop/e-output_from_WNA")
climate3 <- read.csv("e-SUNR_TSUMER_clim_1901-2015MSYT.csv", header=TRUE) 
climate3$sp <- "TSUMER"
climate3$site <- "SUNR"

DD5_sp3 <- climate3$DD5_sp
DD5_sm3 <- climate3$DD5_sm
RH_sp3 <- climate3$RH_sp
PAS3 <- climate3$PAS
prev_DD5_sm3 <- c(NA,DD5_sm3)
prev_DD5_sm3 <- prev_DD5_sm3[1:(length(prev_DD5_sm3)-1)]
sp3 <-climate3$sp
site3 <-climate3$site

DD5_sp3 <- DD5_sp3[-1]
DD5_sm3 <- DD5_sm3[-1]
RH_sp3 <- RH_sp3[-1]
PAS3 <- PAS3[-1]
prev_DD5_sm3 <- prev_DD5_sm3[-1]
sp3 <- sp3[-1]
site3 <- site3[-1]
#now the climate variable time series go from 1902-2015

setwd("C:/Users/Myesa/Desktop/d-master_chronologies")
master3 <- read.csv(file="d-SUNR_TSUMER_master.csv",header=TRUE)
master_trim3 <- master3[-c(1:22),] 
rwi3 <- master_trim3$RWI

#standardizing climate variables 
DD5_sp3<-(DD5_sp3-mean(DD5_sp3))/sd(DD5_sp3)
DD5_sm3<-(DD5_sm3-mean(DD5_sm3))/sd(DD5_sm3)
RH_sp3<-(RH_sp3-mean(RH_sp3))/sd(RH_sp3)
PAS3<-(PAS3-mean(PAS3))/sd(PAS3)
prev_DD5_sm3<-(prev_DD5_sm3-mean(prev_DD5_sm3))/sd(prev_DD5_sm3)

#SUNR CALNOO
setwd("C:/Users/Myesa/Desktop/e-output_from_WNA")
climate4 <- read.csv("e-SUNR_CALNOO_clim_1901-2015MSYT.csv", header=TRUE) 
climate4$sp <- "CALNOO"
climate4$site <- "SUNR"

DD5_sp4 <- climate4$DD5_sp
DD5_sm4 <- climate4$DD5_sm
RH_sp4 <- climate4$RH_sp
PAS4 <- climate4$PAS
prev_DD5_sm4 <- c(NA,DD5_sm4)
prev_DD5_sm4 <- prev_DD5_sm4[1:(length(prev_DD5_sm4)-1)]
sp4 <-climate4$sp
site4 <-climate4$site

DD5_sp4 <- DD5_sp4[-1]
DD5_sm4 <- DD5_sm4[-1]
RH_sp4 <- RH_sp4[-1]
PAS4 <- PAS4[-1]
prev_DD5_sm4 <- prev_DD5_sm4[-1]
sp4 <- sp4[-1]
site4 <- site4[-1]
#now the climate variable time series go from 1902-2015

setwd("C:/Users/Myesa/Desktop/d-master_chronologies")
master4 <- read.csv(file="d-SUNR_CALNOO_master.csv",header=TRUE)
master_trim4 <- master4[-c(1:22),] 
rwi4 <- master_trim4$RWI
rwi4 <- rwi4[1:(length(rwi4)-1)]#trim the value for 2015 if there is an NA (the multiple linear regression will not take an NA)

DD5_sp4 <- DD5_sp4[1:(length(DD5_sp4)-1)]
DD5_sm4 <- DD5_sm4[1:(length(DD5_sm4)-1)]
RH_sp4 <- RH_sp4[1:(length(RH_sp4)-1)]
PAS4 <- PAS4[1:(length(PAS4)-1)]
prev_DD5_sm4 <- prev_DD5_sm4[1:(length(prev_DD5_sm4)-1)]
sp4 <- sp4[1:(length(sp4)-1)]
site4 <- site4[1:(length(site4)-1)]
#standardizing climate variables 
DD5_sp4<-(DD5_sp4-mean(DD5_sp4))/sd(DD5_sp4)
DD5_sm4<-(DD5_sm4-mean(DD5_sm4))/sd(DD5_sm4)
RH_sp4<-(RH_sp4-mean(RH_sp4))/sd(RH_sp4)
PAS4<-(PAS4-mean(PAS4))/sd(PAS4)
prev_DD5_sm4<-(prev_DD5_sm4-mean(prev_DD5_sm4))/sd(prev_DD5_sm4)


#append all together
DD5_sp<-append(DD5_sp,DD5_sp2)
DD5_sp<-append(DD5_sp,DD5_sp3)
DD5_sp<-append(DD5_sp,DD5_sp4)

DD5_sm<-append(DD5_sm,DD5_sm2)
DD5_sm<-append(DD5_sm,DD5_sm3)
DD5_sm<-append(DD5_sm,DD5_sm4)

RH_sp<-append(RH_sp,RH_sp2)
RH_sp<-append(RH_sp,RH_sp3)
RH_sp<-append(RH_sp,RH_sp4)

PAS<-append(PAS,PAS2)
PAS<-append(PAS,PAS3)
PAS<-append(PAS,PAS4)

prev_DD5_sm<-append(prev_DD5_sm,prev_DD5_sm2)
prev_DD5_sm<-append(prev_DD5_sm,prev_DD5_sm3)
prev_DD5_sm<-append(prev_DD5_sm,prev_DD5_sm4)

rwi<-append(rwi,rwi2)
rwi<-append(rwi,rwi3)
rwi<-append(rwi,rwi4)

sp<-append(sp,sp2)
sp<-append(sp,sp3)
sp<-append(sp,sp4)

site<-append(site,site2)
site<-append(site,site3)
site<-append(site,site4)

alldata_SUNR <- data.frame(DD5_sp,DD5_sm,RH_sp,PAS,prev_DD5_sm,rwi,sp,site)
#=====================================================================================================
#Testing for significance between ABIAMA growth and the rest of the species on the south side (PARA)
#=====================================================================================================
DD5_sp <- alldata_PARA$DD5_sp
DD5_sm <- alldata_PARA$DD5_sm
RH_sp <- alldata_PARA$RH_sp
PAS <- alldata_PARA$PAS
prev_DD5_sm <- alldata_PARA$prev_DD5_sm
rwi <- alldata_PARA$rwi
sp <- alldata_PARA$sp

testPARA <- lm(rwi~DD5_sm*sp+DD5_sp*sp+RH_sp*sp+PAS*sp+prev_DD5_sm*sp-sp)
summary(testPARA)

betas <- summary(testPARA)$coefficients[,c(1:4)] #gives estimate (correlation coeff.) std. error, t-value,and p-value
betas
#saving the results as a dataframe
PARA_signif_diff <- as.data.frame(betas)
write.csv(PARA_signif_diff, file ="C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/MLR_data/sigdif_comparisons/PARA_between_ABIAMA_and_all_species.csv", row.names=TRUE) 



#======================================================================================
####testing for across species for differences in sensitivity at the NW side (SPRY)
#======================================================================================
#Testing for significance between ABIAMA growth and the rest of the species on the NW side (SPRY)
DD5_sp <- alldata_SPRY$DD5_sp
DD5_sm <- alldata_SPRY$DD5_sm
RH_sp <- alldata_SPRY$RH_sp
PAS <- alldata_SPRY$PAS
prev_DD5_sm <- alldata_SPRY$prev_DD5_sm
rwi <- alldata_SPRY$rwi
sp <- alldata_SPRY$sp

testSPRY <- lm(rwi~DD5_sm*sp+DD5_sp*sp+RH_sp*sp+PAS*sp+prev_DD5_sm*sp-sp)
summary(testSPRY)

betas <- summary(testSPRY)$coefficients[,c(1:4)] #gives estimate (correlation coeff.) std. error, t-value,and p-value
betas
#saving the results as a dataframe
SPRY_signif_diff <- as.data.frame(betas)
write.csv(SPRY_signif_diff, file ="C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/MLR_data/sigdif_comparisons/SPRY_between_ABIAMA_and_all_species.csv", row.names=TRUE) 


#======================================================================================
####testing for across species for differences in sensitivity at the East side (SUNR)
#======================================================================================

#Testing for significance between ABIAMA growth and the rest of the species on the East side (SUNR)
DD5_sp <- alldata_SUNR$DD5_sp
DD5_sm <- alldata_SUNR$DD5_sm
RH_sp <- alldata_SUNR$RH_sp
PAS <- alldata_SUNR$PAS
prev_DD5_sm <- alldata_SUNR$prev_DD5_sm
rwi <- alldata_SUNR$rwi
sp <- alldata_SUNR$sp

testSUNR <- lm(rwi~DD5_sm*sp+DD5_sp*sp+RH_sp*sp+PAS*sp+prev_DD5_sm*sp-sp)
summary(testSUNR)
betas <- summary(testSUNR)$coefficients[,c(1:4)] #gives estimate (correlation coeff.) std. error, t-values, and p-value
betas
#saving the results as a dataframe
SUNR_signif_diff <- as.data.frame(betas)
write.csv(SUNR_signif_diff, file ="C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/MLR_data/sigdif_comparisons/SUNR_between_ABIAMA_and_all_species.csv", row.names=TRUE) 











