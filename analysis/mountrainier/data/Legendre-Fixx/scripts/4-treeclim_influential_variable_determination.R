##====================================================================================================
## package treeclim, using the function dcc() to idendify variables with a strong influence on growth
##====================================================================================================

setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Master Chronologies 2017")
RWI <- read.csv("SUNR_TSUMER_master.csv", header=TRUE)

setwd("C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/Climate WMA data")
climate <- read.csv("SUNR_TSUMER_clim_1901-2015MSYT.csv", header=TRUE)

#load libraries
library(treeclim)
library(graphics)
library(utils)

climate <- climate[,-c(2:6)]
climate <- climate[,c(1:169)] #matrix with only the monthly data

Tmax <- climate[,grep("^Tmax",names(climate))]#selects all columns beginning with only Tmax  
Tmin <- climate[,grep("^Tmin",names(climate))]#create a data frame for each climate variable  
Tave <- climate[,grep("^Tave",names(climate))]#each column is a monthly parameter for each variable 
PPT <- climate[,grep("^PPT",names(climate))]  
RAD <- climate[,grep("^Rad",names(climate))]  
DD_0 <- climate[,grep("^DD_0",names(climate))] 
DD5 <- climate[,grep("^DD5",names(climate))]  
DD18 <- climate[,grep("^DD18",names(climate))]
DD_18 <- climate[,grep("^DD_18",names(climate))]  
NFFD <- climate[,grep("^NFFD",names(climate))] 
PAS <- climate[,grep("^PAS",names(climate))] 
Eref <- climate[,grep("^Eref",names(climate))] 
CMD <- climate[,grep("^CMD",names(climate))] 
RH <- climate[,grep("^RH",names(climate))]  
Year <- climate$Year

Tmax <- cbind(Year,Tmax)
Tmin <- cbind(Year,Tmin)
Tave <- cbind(Year,Tave)
PPT <- cbind(Year,PPT)
RAD <- cbind(Year,RAD)
DD_0 <- cbind(Year,DD_0) 
DD5 <- cbind(Year,DD5)
DD18 <- cbind(Year,DD18) 
DD_18 <- cbind(Year,DD_18)  
NFFD <- cbind(Year,NFFD)
PAS <- cbind(Year,PAS)
Eref <- cbind(Year,Eref)
CMD <- cbind(Year,CMD) 
RH <- cbind(Year,RH)

colnames(Tmax) <- c("Year", "JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")
colnames(Tmin) <- c("Year", "JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")
colnames(Tave) <- c("Year", "JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")
colnames(PPT) <- c("Year", "JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")
colnames(RAD) <- c("Year", "JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")
colnames(DD_0) <- c("Year", "JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")
colnames(DD5) <- c("Year", "JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")
colnames(DD18) <- c("Year", "JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")
colnames(DD_18) <- c("Year", "JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")
colnames(NFFD) <- c("Year", "JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")
colnames(PAS) <- c("Year", "JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")
colnames(Eref) <- c("Year", "JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")
colnames(CMD) <- c("Year", "JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")
colnames(RH) <- c("Year", "JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")

colnames(RWI) <- c("","RWI")

RWI[,1] <- NULL # removing year column
rownames(RWI)<-seq(from=1880, to=2015) # this renames the rownames w/ the years, must use the appropriate values for the chronology

##=========Important=======
#Do not use the 5 lines of code below if you are continuing adding to the sum_results dataframe. Instead, use the 5 lines of code just below that are commented out

##climate variables from previous May to current October are tested (selection=-5:10)
##For the variables DD_0, DD18, Eref and CMD have to adjust the selection range manually;month columns with too many zeros will cause an error in treeclim
results <- dcc(RWI,DD_0, selection= -5:-6, method="correlation", dynamic="static", var_names= "DD_0")
results
plot(results)
sum_results <- summary(results)
sum_results <- data.frame(sum_results)
#=========================Important==================================
##For repeating this script for the 2nd and 3rd sites, do not create a new results file.
#Add on to the results of the same species
#Use the code below
#temp <- dcc(RWI,DD_0, selection= -5:-6, method="correlation", dynamic="static", var_names= "DD_0")
#plot(temp)
#sum_temp <- summary(temp)
#sum_temp <- data.frame(sum_temp)
#sum_results <- rbind(sum_results,sum_temp)
#======================================================================================================

temp <- dcc(RWI,DD_0, selection = -9:-10, method="correlation", dynamic="static", var_names= "DD_0") 
plot(temp)
sum_temp <- summary(temp)
sum_temp <- data.frame(sum_temp)
sum_results <- rbind(sum_results,sum_temp)

temp <- dcc(RWI,DD_0, selection = -9:6, method="correlation", dynamic="static", var_names= "DD_0") 
plot(temp)
sum_temp <- summary(temp)
sum_temp <- data.frame(sum_temp)
sum_results <- rbind(sum_results,sum_temp)

temp <- dcc(RWI,DD_0, selection = 9:10, method="correlation", dynamic="static", var_names= "DD_0") 
plot(temp)
sum_temp <- summary(temp)
sum_temp <- data.frame(sum_temp)
sum_results <- rbind(sum_results,sum_temp)
#===================End of DD_0 correlation=========================================
temp <- dcc(RWI,DD18, selection = -5:-10, method="correlation", dynamic="static", var_names= "DD18") 
plot(temp)
sum_temp <- summary(temp)
sum_temp <- data.frame(sum_temp)
sum_results <- rbind(sum_results,sum_temp)

temp <- dcc(RWI,DD18, selection = 5:10, method="correlation", dynamic="static", var_names= "DD18") 
plot(temp)
sum_temp <- summary(temp)
sum_temp <- data.frame(sum_temp)
sum_results <- rbind(sum_results,sum_temp)

#===================End of DD18 correlation=========================================
temp <- dcc(RWI,Eref, selection = -5:-11, method="correlation", dynamic="static", var_names= "Eref") 
plot(temp)
sum_temp <- summary(temp)
sum_temp <- data.frame(sum_temp)
sum_results <- rbind(sum_results,sum_temp)

temp <- dcc(RWI,Eref, selection = 3:10, method="correlation", dynamic="static", var_names= "Eref") 
plot(temp)
sum_temp <- summary(temp)
sum_temp <- data.frame(sum_temp)
sum_results <- rbind(sum_results,sum_temp)

#==================End of Eref correlation=====================================================
temp <- dcc(RWI,CMD, selection = -5:-9, method="correlation", dynamic="static", var_names= "CMD") 
plot(temp)
sum_temp <- summary(temp)
sum_temp <- data.frame(sum_temp)
sum_results <- rbind(sum_results,sum_temp)

temp <- dcc(RWI,CMD, selection = 5:9, method="correlation", dynamic="static", var_names= "CMD") 
plot(temp)
sum_temp <- summary(temp)
sum_temp <- data.frame(sum_temp)
sum_results <- rbind(sum_results,sum_temp)
#==================End of CMD correlation=====================================================
temp <- dcc(RWI,Tmax, selection = -5:10, method="correlation", dynamic="static", var_names= "Tmax")
plot(temp)
sum_temp <- summary(temp)
sum_temp <- data.frame(sum_temp)
sum_results <- rbind(sum_results,sum_temp)

temp <- dcc(RWI,Tmin, selection = -5:10, method="correlation", dynamic="static", var_names= "Tmin") 
plot(temp)
sum_temp <- summary(temp)
sum_temp <- data.frame(sum_temp)
sum_results <- rbind(sum_results,sum_temp)

temp <- dcc(RWI,Tave, selection = -5:10, method="correlation", dynamic="static", var_names= "Tave") 
plot(temp)
sum_temp <- summary(temp)
sum_temp <- data.frame(sum_temp)
sum_results <- rbind(sum_results,sum_temp)

temp <- dcc(RWI,PPT, selection = -5:10, method="correlation", dynamic="static", var_names= "PPT") 
plot(temp)
sum_temp <- summary(temp)
sum_temp <- data.frame(sum_temp)
sum_results <- rbind(sum_results,sum_temp)

temp <- dcc(RWI,RAD, selection = -5:10, method="correlation", dynamic="static", var_names= "RAD") 
plot(temp)
sum_temp <- summary(temp)
sum_temp <- data.frame(sum_temp)
sum_results <- rbind(sum_results,sum_temp)

temp <- dcc(RWI,DD5, selection = -5:10, method="correlation", dynamic="static", var_names= "DD5") 
plot(temp)
sum_temp <- summary(temp)
sum_temp <- data.frame(sum_temp)
sum_results <- rbind(sum_results,sum_temp)

temp <- dcc(RWI,DD_18, selection = -5:10, method="correlation", dynamic="static", var_names= "DD_18") 
plot(temp)
sum_temp <- summary(temp)
sum_temp <- data.frame(sum_temp)
sum_results <- rbind(sum_results,sum_temp)


temp <- dcc(RWI,NFFD, selection = -5:10, method="correlation", dynamic="static", var_names= "NFFD") 
plot(temp)
sum_temp <- summary(temp)
sum_temp <- data.frame(sum_temp)
sum_results <- rbind(sum_results,sum_temp)

temp <- dcc(RWI,PAS, selection = -5:10, method="correlation", dynamic="static", var_names= "PAS") 
plot(temp)
sum_temp <- summary(temp)
sum_temp <- data.frame(sum_temp)
sum_results <- rbind(sum_results,sum_temp)

temp <- dcc(RWI,RH, selection = -5:10, method="correlation", dynamic="static", var_names= "RH") 
plot(temp)
sum_temp <- summary(temp)
sum_temp <- data.frame(sum_temp)
sum_results <- rbind(sum_results,sum_temp)

#now, the above code should be repeated for the same species at the other two sites

#######============================================================================
##only write .csv after results from all sites of a species are put into the dataframe sum_results
write.csv(sum_results, file ="C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/MLR_data/Correl Output from treeclim/PARA_SPRY_SUNR_TSUMER_correl_all_variables.csv", row.names=TRUE) #Change the output filename here

#returns subset of the dataframe that meets the conditions
signif_results<-subset(sum_results,sum_results$significant == TRUE)
#sort the dataframe by coef value and use head() and tail() to select the variables that have the most influence
attach(signif_results)
signif_results_sorted <- signif_results[order(coef),]
first_30 <- head(signif_results_sorted, 30)
last_30 <- tail(signif_results_sorted, 30)

#Combine the first and last rows to make into a condensed df for export
most_sig_results <- rbind(first_30,last_30)
#summary(most_sig_results$varname) #this shows how many times a variable was most significant
#rownames(most_sig_results)  #This shows the names of the variables that were significant for this species
write.csv(most_sig_results, file ="C:/Users/Western Hemlock/Desktop/WinDENDRO_do_not_delete/Myesa/MLR_data/Correl Output from treeclim/PARA_SPRY_SUNR_TSUMER_correl_most_all_variables.csv", row.names=TRUE) #Change the output filename here

