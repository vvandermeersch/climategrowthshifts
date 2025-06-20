##=================================================================================================================
##### Making a master chronology for each species at every site after prewhitening the tree chronologies#######
##=================================================================================================================

setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO/Myesa/Prewhitened CSV files 2017")
library(dplR)

prewh <- read.csv(file="PARA_ABIAMA_prew.csv",header=TRUE)

prewh[,1] <- NULL # removing year core was collected column(which is year "Yr")
Yr <- 2015#year of last ring; for me: 2011 for PARA ABILAS and 2015 for all others, doublecheck year of last ring 
newnrow=nrow(prewh)-1
rownames(prewh)<-seq(from=(Yr), to= (Yr-newnrow)) #this renames the column in R with the years, I found that the functions work better when in this format, because when the csv is read in, the years are incorporated as part of the dataframe 

#apply() the tbrm function across rows to get a master column 
master <-apply(prewh,1,tbrm) #margin=1 applies across rows
sink("master",append=F, split=T)

write.csv(master, file ="C:/Users/Western Hemlock/Desktop/WinDENDRO/Myesa/Master Chronologies 2017/SUNR_CALNOO_master1.csv") #Change the output filename here

setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO/Myesa/Master Chronologies 2017")
#read the csv in again to lable the RW column RWI (it defaults to "X")
to_label <- read.csv(file="SUNR_CALNOO_master1.csv",header=TRUE) #change filename to the file just written
colnames(to_label) <- c("Year", "RWI")
write.csv(to_label, file ="C:/Users/Western Hemlock/Desktop/WinDENDRO/Myesa/Master Chronologies 2017/SUNR_CALNOO_master.csv", row.names=FALSE) #Change the output filename here

#=========================================
#ploting a master chronology for a figure
#=========================================

setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO/Myesa/Master Chronologies 2017")
library(dplR)
master_plt <- read.csv(file="SUNR_TSUMER_master.csv",header=TRUE)
plot(master_plt$RWI~master_plt$Year, type="l")
