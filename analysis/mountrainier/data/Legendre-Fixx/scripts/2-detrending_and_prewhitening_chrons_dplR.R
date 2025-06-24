#=======================================================
## 			detrending 
#=======================================================
#first, set working directory
library(dplR)
# looking at the trends for each core:

for(i in 1:ncol(spstdat.rwl.df)){windows()
  detrend(rwl=spstdat.rwl.df[,c(i,i+1)], method="Spline", nyrs=67,make.plot=T, pos.slope=T)}


## Creating a detrended rw data frame
rw.det <- detrend(rwl=spstdat.rwl.df, method="Spline", nyrs=67, pos.slope=T)
write.csv(rw.det, file ="C:/Users/Western Hemlock/Desktop/WinDENDRO/Myesa/Detrended CSV files 2017/PARA_ABILAS_detrended_67yr.csv") #Change the output filename here

#=======================================================
## 			prewhittening
#=======================================================


# set wd to source location
setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO/Myesa/Detrended CSV files 2017")
library(dplR)

RWdet <- read.csv(file="SUNR_CALNOO_detrended_67yr.csv",header=TRUE) #read in detrended file

RWdet[,1] <- NULL # removing year core was collected column (which is year "Yr")
Yr <- 2015#year of last ring; for me: 2011 for PARA ABILAS, and 2015 for the rest. Double check year of last ring
newnrow=nrow(RWdet)-1
rownames(RWdet)<-seq(from=(Yr), to= (Yr-newnrow)) # this renames the dataframe with the year (must be in this format for following code to run)


######## AR1 detrending function called if prewhiten =T
# takes a single series, and turns it into resid(AR1) + mean / mean
# this shortens the time series by 1 on each side.
det.ar <- function(y, ARorder=1){
  good.y <- which(!is.na(y))
  if (any(diff(good.y) != 1)){
    stop("series has internal NAs")
  }
  if(is.array(y) == FALSE){
    y <- t(y)
  }
  y2 <- y[good.y]
  det <- rep(NA, times=length(y))
  Ar <- ar(y2, order.max = ARorder)
  ydet <- (Ar$resid + Ar$x.mean)/Ar$x.mean
  det[good.y] <- ydet
  return(det)
}

#using the apply() function: quick opperations applied to a matrix, vector, or array
#apply(variable,margin, function)
prew <- apply(RWdet,2,det.ar) ##margin=1 apply by row, margin=2 apply by column, margin=3 apply to each element in matrix

sink("prew",append=F, split=T)


Yr <- 2015#year of last ring; for me: 2011 for PARA ABILAS, and 2015 for the rest. Double check year of last ring
newnrow=nrow(prew)-1
rownames(prew)<-seq(from=(Yr), to= (Yr-newnrow)) # this adds the year column to the file

write.csv(prew, file ="C:/Users/Western Hemlock/Desktop/WinDENDRO/Myesa/Prewhitened CSV files 2017/SUNR_CALNOO_prew.csv") #Change the output filename here

####plot the lagged autocorrelation for a single core

acf(RWdet$Q43, type="correlation", na.action=na.pass)

