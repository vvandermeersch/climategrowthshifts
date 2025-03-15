##Code Yields a table called "dat" of tree ring data averaged between two cores per tree (for trees with 2 cores collected, if only one core- it is a single core). Each row is a different tree. "dat.df" is the same data, but in dataframe form."rngdat" contains ONLY the tree ring data
##Read in data
setwd("~/Dropbox/Documents/Work/UW/Research/Mount Rainier/2010ManuscriptTreeGrowthClimate/Analyses for 2010 Manuscript/LMEANALYSES.2010")
rawdat1<-read.csv("SouthSideCoresX.csv", header=TRUE)
#load lme4 and languageR packages
#first average cores from same trees, get rid of 1st and last year; go back 1850
rawdat<-rawdat1[,c(1:3,5:6,9:169)] #only columns we need
tags<-unique(rawdat[,3])
inddat<-matrix(NA, nrow=length(tags), ncol=5)
dimnames(inddat)<-list(c(), c("stand","species","tag","DBH","lastyear"))
rngdat<-c()

#Code to get averages from 2 trees: NOTE THIS WILL ONLY WORK WITH 2007 AND 2008 CORES MIXED
for(i in 1:length(tags)){
	tmpdat<-rawdat[rawdat[,3]==tags[i],]
	lastring<-min(tmpdat[,5])-1
	if(lastring==2007){firstyr<-8}
	if(lastring==2008){firstyr<-7}
	lastyr<-c()
	for(j in 1:dim(tmpdat)[1]){
		rngs<-tmpdat[j,6:dim(tmpdat)[2]]
		tmp<-which(is.na(rngs)==TRUE)
		tmp<-tmp[tmp>3]
		if(length(tmp)==0){mintmp<-161}
		if(length(tmp)>0){mintmp<-min(tmp)}
		lastyr<-c(lastyr,mintmp+4)
	}
	lastyr<-min(lastyr)
	avgrng<-colMeans(tmpdat[,firstyr:lastyr])
	for(j in 1:4){
		inddat[i,j]<-as.character(tmpdat[1,j])}
	inddat[i,5]<-lastring
	if(lastring==2007){
		treerng<-c(NA,avgrng,rep(NA, times=dim(tmpdat)[2]-lastyr))}
	if(lastring==2008){
		treerng<-c(avgrng,rep(NA, times=dim(tmpdat)[2]-lastyr))}
	rngdat<-rbind(rngdat,treerng)
}

dimnames(rngdat)<-list(c(), seq(2008,1849))
dat<-cbind(inddat,rngdat[,1:159]) #merge data together
dat.df<-as.data.frame(cbind(inddat,rngdat[,1:159]))
rownames(rngdat)<-dat[,3]

#to detrend all rings at once using detrend() function
rwidat<-detrend(t(rngdat),method="Spline",nyrs=100, f=0.5)
dat.rwi<-cbind(inddat,t(rwidat))
dat.rwi<-dat.rwi[!dat.rwi[,3]=="F38",]#because tree F38 is doing weird things! has some really high rwis!
dat.rwi.df<-data.frame(inddat,t(rwidat))
dat.rwi.df<-dat.rwi.df[!dat.rwi.df[,3]=="F38",]#because tree F38 is doing weird things! has some really high rwis!
###could instead detrend all ring widths one tree at a time- does this make a difference?
#allrwi<-c()
#for (i in 1:dim(rngdat)[1]){
#	tree<-rngdat[i,]
#rwiseries<-detrend.series(tree,method="Spline",nyrs=100, #f=0.5,make.plot = FALSE)
#allrwi<-rbind(allrwi,rwiseries)
#}
#all.rwi<-cbind(inddat,allrwi)
#rownames(all.rwi)<-dat[,3]
#dat.rwi<-all.rwi
#dat.rwi<-dat.rwi[!dat.rwi[,3]=="F38",]#because tree F38 is doing weird things! has some really high rwis!
#dat.rwi.df<-data.frame(dat.rwi)
#check to see if results from detrend.series and detrend are correlated
#plot(allrwi[15,],rwidat[,15])#they area for just about every tree.

#summarize climate data: up to 2007 growing season
climate<-read.csv("Clim.csv", header=TRUE)
minyr<-min(climate[,1])+1; maxyr<-2007 #WE SHOULD GET MORE RECENT DATA
MAT<-c(); TP<-c(); Tsnow<-c(); GP<-c(); DP<-c(); GT<-c(); DT<-c()
for(i in minyr:maxyr){
	tmpclim<-climate[climate[,1]==i|climate[,1]==i-1,]
	climyr<-tmpclim[10:21,]
	tmpTP<-sum(climyr[,3]); TP<-c(TP,tmpTP)
	tmpMAT<-mean(climyr[,4]); MAT<-c(MAT,tmpMAT)
	tmpTsnow<-sum(climyr[,5]); Tsnow<-c(Tsnow,tmpTsnow)
	tmpGP<-sum(climyr[8:12,3]); GP<-c(GP,tmpGP)
	tmpDP<-sum(climyr[2:6,3]); DP<-c(DP,tmpDP)
	tmpGT<-mean(climyr[8:12,4]); GT<-c(GT,tmpGT)
	tmpDT<-mean(climyr[2:6,4]); DT<-c(DT,tmpDT)
}
clim<-cbind(seq(minyr,maxyr),MAT,TP,Tsnow,GP,DP,GT,DT)
dimnames(clim)<-list(c(), c("year","MAT","TPrecip","Tsnow","GPrecip","DPrecip","GTemp","DTemp"))
MAT<-clim[,2]
TP<-clim[,3]
GP<-clim[,5]
DP<-clim[,6]
GT<-clim[,7]
DT<-clim[,8]
#Next correlate with tree growth
#standardize climatic variables, if desired
MAT2<-(MAT-mean(MAT))/sd(MAT)
TP2<-(TP-mean(TP))/sd(TP)
Tsnow2<-(Tsnow-mean(Tsnow))/sd(Tsnow)
GP2<-(GP-mean(GP))/sd(GP)
DP2<-(DP-mean(DP))/sd(DP)
GT2<-(GT-mean(GT))/sd(GT)
DT2<-(DT-mean(DT))/sd(DT)
climst<-cbind(seq(minyr,maxyr),MAT2,TP2,Tsnow2,GP2,DP2,GT2,DT2)
dimnames(climst)<-list(c(), c("year","MATST","TPrecipST","TsnowST","GPrecipST","DPrecipST","GTempST","DTempST"))
#reorder clim to match tree ring order
climst<-climst[order(climst[,1], decreasing=TRUE),]
head(climst)
MATST<-climst[,2]
TPST<-climst[,3]
GPST<-climst[,5]
DPST<-climst[,6]
GTST<-climst[,7]
DTST<-climst[,8]


#Kevin's climate data estimates for each stand, which are based on PRISM estimates
standtemp<-read.csv("tree_plot_climate_temp.csv", header=TRUE)
standprecip<-read.csv("tree_plot_climate_precip.csv", header=TRUE)
unique(standprecip[,"Plot"])
dimnames(standprecip)
stands<-standtemp[,1]
month<-rep(1:12,length.out=1204)
year<-rep(1909:2009, each=12,length.out=1204)
monthyear<-cbind(year,month)
sttemp<-standtemp[,5:1208]
sttemp<-t(sttemp)
stdtemp<-cbind(monthyear,sttemp)
stdtemp<-stdtemp[61:1186,]#because snow data only goes back to 1914 and tree ring data ends at 2007 for most of the trees
colnames(stdtemp)<-c("year","month","AB08","AE10","AG05","AM16","AO03","AR07","AV02","AV06","AV14","AX15","PARA","PP17","TA01","TB13","TO04","TO11")
stprecip<-standprecip[,5:1208]
stprecip<-t(stprecip)
stdprecip<-cbind(monthyear,stprecip)
stdprecip<-stdprecip[61:1186,]
colnames(stdprecip)<-c("year","month","AB08","AE10","AG05","AM16","AO03","AR07","AV02","AV06","AV14","AX15","PARA","PP17","TA01","TB13","TO04","TO11")
row.names(stdtemp)<-NULL
row.names(stdprecip)<-NULL
sites<-c("PARA","AE10","AR07","AM16","AX15","AV06","AG05","TB13","TO04")
stdswe<-read.csv("tree_plot_climate_swe.csv", header=TRUE)
stdsnowdur<-read.csv("tree_plot_climate_snowdur.csv", header=TRUE)
stdgrdd5<-read.csv("tree_plot_climate_grdd5.csv", header=TRUE)
dim(stdsnowdur)
dim(stdswe)
dim(stdtemp)

#create a for loop to get all climate variables for each stand and examine correlations between them
for (i in 1:length(sites)){
	sitetemp<-cbind(stdtemp[,1:2],stdtemp[,dimnames(stdtemp)[[2]]==sites[i]])
	mntemp<-tapply(sitetemp[,3],sitetemp[,1],mean)
	sitetempgr<-sitetemp[sitetemp[,2]>5&sitetemp[,2]<11,]
	growtemp<-tapply(sitetempgr[,3],sitetempgr[,1],mean)
	sitetempdr<-sitetemp[sitetemp[,2]<5,]
	sitetempdr2<-sitetemp[sitetemp[,2]==12,]
	sitetempdr2[,1]<-sitetempdr2[,1]+1
	sitetempdr<-rbind(sitetempdr,sitetempdr2)
	dortemp<-tapply(sitetempdr[,3],sitetempdr[,1],mean)
	#now precip
	siteprecip<-cbind(stdprecip[,1:2],stdprecip[,dimnames(stdprecip)[[2]]==sites[i]])
	totprecip<-tapply(siteprecip[,3],siteprecip[,1],sum)
	siteprecipgr<-siteprecip[siteprecip[,2]>5&siteprecip[,2]<11,]
	growprecip<-tapply(siteprecipgr[,3],siteprecipgr[,1],sum)
	siteprecipdr<-siteprecip[siteprecip[,2]<5,]
	siteprecipdr2<-siteprecip[siteprecip[,2]==12,]
	siteprecipdr2[,1]<-siteprecipdr2[,1]+1
	siteprecipdr<-rbind(siteprecipdr,siteprecipdr2)
	dorprecip<-tapply(siteprecipdr[,3],siteprecipdr[,1],sum)
	#now snowdata&growingdegreedays(5c threshold)
	siteswe<-stdswe[,dimnames(stdswe)[[2]]==sites[i]]
	names(siteswe)<-stdswe[,1]
	swe<-siteswe[1:94]
	sitesnowdur<-stdsnowdur[,dimnames(stdsnowdur)[[2]]==sites[i]]
	names(sitesnowdur)<-stdsnowdur[,1]
	snowdur<-sitesnowdur[1:94]
	sitegrowdd<-stdgrdd5[,dimnames(stdgrdd5)[[2]]==sites[i]]
	names(sitegrowdd)<-stdgrdd5[,1]
	growdd<-sitegrowdd[1:94]
	#create matrix with all climate variables clim using cbind
	climsite<-cbind(mntemp,growtemp,dortemp,totprecip,growprecip,dorprecip,swe,snowdur,growdd)
	climsite<-climsite[order(dimnames(climsite)[[1]],decreasing=TRUE),]#to reorder climsite so that it is in the same order as the tree ring data!
	#use command pairs to look at data; cor to get correlations
	quartz(width=10, height=7)
	pairs(climsite,main=sites[i])
	print(sites[i])
	print(cor(climsite))#print and save correlations between climate variables
	#standardize site climatic variables, if desired
STMAT2<-(mntemp-mean(mntemp))/sd(mntemp)
STTP2<-(totprecip-mean(totprecip))/sd(totprecip)
STGP2<-(growprecip-mean(growprecip))/sd(growprecip)
STDP2<-(dorprecip-mean(dorprecip))/sd(dorprecip)
STGT2<-(growtemp-mean(growtemp))/sd(growtemp)
STDT2<-(dortemp-mean(dortemp))/sd(dortemp)
STSWE2<-(swe-mean(swe))/sd(swe)
STSNDUR2<-(snowdur-mean(snowdur))/sd(snowdur)
STGROWDD2<-(growdd-mean(growdd))/sd(growdd)
climsitestd<-cbind(seq(1914,2007),STMAT2,STTP2,STGP2,STDP2,STGT2,STDT2,STSWE2,STSNDUR2,STGROWDD2)
dimnames(climsitestd)<-list(c(), c("year","MATST","TPrecipST","GPrecipST","DPrecipST","GTempST","DTempST","SWEST","SnowDurST","GrowDDST"))
#reorder clim
climsitestd<-climsitestd[order(climsitestd[,1], decreasing=TRUE),]
head(climsitestd)
	#now save correlations and climate data for each site
	if(sites[i]=="PARA"){corclimPARA<-cor(climsite)}
	if(sites[i]=="AE10"){corclimAE10<-cor(climsite)}
	if(sites[i]=="AR07"){corclimAR07<-cor(climsite)}
	if(sites[i]=="AV06"){corclimAV06<-cor(climsite)}
	if(sites[i]=="AM16"){corclimAM16<-cor(climsite)}
	if(sites[i]=="AX15"){corclimAX15<-cor(climsite)}
	if(sites[i]=="AG05"){corclimAG05<-cor(climsite)}
	if(sites[i]=="TB13"){corclimTB13<-cor(climsite)}
	if(sites[i]=="TO04"){corclimTO04<-cor(climsite)}
	#create matrices for each stand's climate data(unstand & stand)
	if(sites[i]=="PARA"){climstdPARA<-climsitestd}
	if(sites[i]=="AE10"){climstdAE10<-climsitestd}
	if(sites[i]=="AR07"){climstdAR07<-climsitestd}
	if(sites[i]=="AV06"){climstdAV06<-climsitestd}
	if(sites[i]=="AM16"){climstdAM16<-climsitestd}
	if(sites[i]=="AX15"){climstdAX15<-climsitestd}
	if(sites[i]=="AG05"){climstdAG05<-climsitestd}
	if(sites[i]=="TB13"){climstdTB13<-climsitestd}
	if(sites[i]=="TO04"){climstdTO04<-climsitestd}
	
	if(sites[i]=="PARA"){climPARA<-climsite}
	if(sites[i]=="AE10"){climAE10<-climsite}
	if(sites[i]=="AR07"){climAR07<-climsite}
	if(sites[i]=="AV06"){climAV06<-climsite}
	if(sites[i]=="AM16"){climAM16<-climsite}
	if(sites[i]=="AX15"){climAX15<-climsite}
	if(sites[i]=="AG05"){climAG05<-climsite}
	if(sites[i]=="TB13"){climTB13<-climsite}
	if(sites[i]=="TO04"){climTO04<-climsite}
	}
#Next correlate with tree growth
MATST<-climsitestd[,2]
TPST<-climsitestd[,3]
GPST<-climsitestd[,4]
DPST<-climsitestd[,5]
GTST<-climsitestd[,6]
DTST<-climsitestd[,7]

SNDST<-climsitestd[,9]
GDDST<-climsitestd[,10]

##Now look at residuals, etc of rwi vs MAT, TP, and both together.
#stands<-unique(dat.rwi.df[,1])
#stands<-c#("TO04","TB13","AG05","AV06","AX15","AM16","AR07","AE10","PARA")
#elev<-c(704,851,950,1064,1091,1197,1454,1460,1603)
#r.elev<-elev[order(elev, decreasing=TRUE)]
#reorder stands by elevation
#stands<-stands[order(elev, decreasing=TRUE)]
#species<-unique(dat.rwi.df[,2])
#Look at boxcox and qq plots of individual tree-MAT and PPT models
#for (i in 1: length(species)){
#	tmp<-dat.rwi.df[dat.rwi.df[,2]==species[i],1:dim(dat.rwi.df)[2]]
#	for(j in 1:length(stands)){
# 		tmp2<-tmp[tmp[,1]==stands[j],]
# 		if(dim(tmp2)[1]==0){next}
# 		quartz(width=16,height=9)
 #		numtrees<-dim(tmp2)[[1]]
 #		par(mfcol=c(10,numtrees), omi=c(0.1,0.1,0.1,0.1), mai=c(0.1,0.1,0.1,0.1))
# 		for (k in 1:dim(tmp2)[[1]]){
# 			indrwi<-as.numeric(tmp2[k,7:104])
# 			hist(indrwi,main=paste(tmp2[k,3],tmp2[k,2]), col="gray")#to look at distribution of growth-looks pretty normal
 #			temptest1<-lm(indrwi~MATST)
 #			temptest2<-lm(indrwi~TPST)
 #			temptestint<-lm(indrwi~MATST+TPST)
#			plot(MAT,indrwi, main=paste("MAT",tmp2[k,3],tmp2[k,2]))
#			abline(temptest1)
#			plot(temptest1$fitted,temptest1$resid)
 #			qqnorm(temptest1$resid)
#			qqline(temptest1$resid)
#			boxcox(temptest1, plotit=T)
#			plot(TP,indrwi, main=paste("TP"))
#			abline(temptest2)
#			plot(temptest2$fitted,temptest2$resid)
# 			qqnorm(temptest2$resid)
#			qqline(temptest2$resid)
#			boxcox(temptest2, plotit=T)
#			boxcox(temptestint, plotit=T)
# 	}}}


####Now analyze tree growth and climate relationships, using lmer
stands<-c("TO04","TB13","AG05","AV06","AX15","AM16","AR07","AE10","PARA")
elev<-c(704,851,950,1064,1091,1197,1454,1460,1603)
r.elev<-elev[order(elev, decreasing=TRUE)]
#reorder stands by elevation
stands<-stands[order(elev, decreasing=TRUE)]
species<-unique(dat.rwi[,2])
#this code is really slow, so generally better to run one species at a time
for (i in 1: length(species)){
	spdat<-dat.rwi[dat.rwi[,2]==species[i],1:dim(dat.rwi)[2]]
	allbestmod<-c()
	allsites<-c()
	allsimplebestmod<-c()
	for(j in 1:length(stands)){
 		spstanddat<-spdat[spdat[,1]==stands[j],]
 		if(dim(spstanddat)[1]==0){next}
		standdat<-c() 
		for(k in 1:dim(spstanddat)[1]){
			inddat<-spstanddat[k,6:165]
			yrs<-seq(2008,1849)
			inds<-rep(k, times=length(yrs))
			tmp<-cbind(yrs,inds,inddat)
			tmp2<-tmp[2:95,] #start at 2007, end at 1910 (we have ppt and mat data for this
			if(stands[j]=="PARA"){climst<-climstdPARA}
			if(stands[j]=="AE10"){climst<-climstdAE10}
			if(stands[j]=="AR07"){climst<-climstdAR07}
			if(stands[j]=="AV06"){climst<-climstdAV06}
			if(stands[j]=="AM16"){climst<-climstdAM16}
			if(stands[j]=="AX15"){climst<-climstdAX15}
			if(stands[j]=="AG05"){climst<-climstdAG05}
			if(stands[j]=="TB13"){climst<-climstdTB13}
			if(stands[j]=="TO04"){climst<-climstdTO04}
			yrs<-climst[,1]
			MAT<-climst[,2]
			PPT<-climst[,3]
			GPT<-climst[,4]
			DPT<-climst[,5]
			GST<-climst[,6]
			DST<-climst[,7]
			SWE<-climst[,8]
			SNDR<-climst[,9]
			GDD<-climst[,10]
			#put data into standdat (for later analysis)
			tmp3<-cbind(yrs,MAT,GST,DST,PPT,GPT,DPT,SWE,SNDR,GDD,tmp2)
			tmp4<-tmp3[is.na(tmp3[,11])==F,]
			standdat<-rbind(standdat,tmp4)
			}
			yrs<-as.factor(standdat[,1])
			MAT<-as.numeric(standdat[,2])
			GST<-as.numeric(standdat[,3])
			DST<-as.numeric(standdat[,4])
			PPT<-as.numeric(standdat[,5])
			GPT<-as.numeric(standdat[,6])
			DPT<-as.numeric(standdat[,7])
			SWE<-as.numeric(standdat[,8])
			SNDR<-as.numeric(standdat[,9])
			GDD<-as.numeric(standdat[,10])
			ind<-as.factor(standdat[,12])
			rwi<-as.numeric(standdat[,13])
		#first plot things so you can see relationships
		#quartz()
		#xyplot(rwi ~ MAT | ind,groups = ind,pch=16,type =c("p","r"),main=paste(species[i],stands[j]),lwd=2)
		#quartz()
		#xyplot(rwi ~ PPT | ind,groups = ind,pch=16,type =c("p","r"),main=paste(species[i],stands[j]),lwd=2)
		#quartz()
		#xyplot(rwi ~ SWE | ind,groups = ind,pch=16,type =c("p","r"),main=paste(species[i],stands[j]),lwd=2)
		#quartz()
		#xyplot(rwi ~ GDD | ind,groups = ind,pch=16,type =c("p","r"),main=paste(species[i],stands[j]),lwd=2)
		#now fit detrended rwi as a function of climate
		test0<-lmer(rwi~1+(0+1|ind)+(1|yrs),control=list(maxIter=5000),REML=FALSE)
		test1<-lmer(rwi~MAT+(0+MAT|ind)+(1|yrs), control=list(maxIter=5000),REML=FALSE)# 
		test2<-lmer(rwi~PPT+(0+PPT|ind)+(1|yrs), control=list(maxIter=5000),REML=FALSE)
		test3<-lmer(rwi~MAT+PPT+(0+MAT|ind)+(0+PPT|ind)+(1|yrs), control=list(maxIter=5000),REML=FALSE)
		test4<-lmer(rwi~MAT*PPT+(0+MAT|ind)+(0+PPT|ind)+(0+MAT:PPT|ind)+(1|yrs), control=list(maxIter=5000),REML=FALSE)
		test5<-lmer(rwi~GST+(0+GST|ind)+(1|yrs),control=list(maxIter=5000),REML=FALSE)#
		test6<-lmer(rwi~DST+(0+DST|ind)+(1|yrs), control=list(maxIter=5000),REML=FALSE)#
		test7<-lmer(rwi~GPT+(0+GPT|ind)+(1|yrs), control=list(maxIter=5000),REML=FALSE)#
		test8<-lmer(rwi~DPT+(0+DPT|ind)+(1|yrs), control=list(maxIter=5000),REML=FALSE)#
		test9<-lmer(rwi~GST+DST+(0+GST|ind)+(0+DST|ind)+(1|yrs),control=list(maxIter=5000),REML=FALSE)
		test10<-lmer(rwi~GPT+DPT+(0+GPT|ind)+(0+DPT|ind)+(1|yrs),control=list(maxIter=5000),REML=FALSE)
		test11<-lmer(rwi~GST+DST+GPT+DPT+(0+GST|ind)+(0+DST|ind)+(0+GPT|ind)+(0+DPT|ind)+(1|yrs),control=list(maxIter=5000),REML=FALSE)
		test12<-lmer(rwi~GST+GPT+(0+GST|ind)+(0+GPT|ind)+(1|yrs),control=list(maxIter=5000),REML=FALSE)
		test13<-lmer(rwi~GST*GPT+(0+GST|ind)+(0+GPT|ind)+(0+GST:GPT|ind)+(1|yrs),control=list(maxIter=5000),REML=FALSE)
		test14<-lmer(rwi~DST+DPT+(0+DST|ind)+(0+DPT|ind)+(1|yrs),control=list(maxIter=5000),REML=FALSE)
		test15<-lmer(rwi~DST*DPT+(0+DST|ind)+(0+DPT|ind)+(0+DST:DPT|ind)+(1|yrs),control=list(maxIter=5000),REML=FALSE)
		test16<-lmer(rwi~SWE+(0+SWE|ind)+(1|yrs),control=list(maxIter=5000),REML=FALSE)
		test17<-lmer(rwi~SNDR+(0+SNDR|ind)+(1|yrs),control=list(maxIter=5000),REML=FALSE)
		test18<-lmer(rwi~GDD+(0+GDD|ind)+(1|yrs),control=list(maxIter=5000),REML=FALSE)
		test19<-lmer(rwi~SWE+GST+(0+SWE|ind)+(0+GST|ind)+(1|yrs),control=list(maxIter=5000),REML=FALSE)
		test20<-lmer(rwi~SNDR+GST+(0+SNDR|ind)+(0+GST|ind)+(1|yrs),control=list(maxIter=5000),REML=FALSE)
		test21<-lmer(rwi~SWE+GPT+(0+SWE|ind)+(0+GPT|ind)+(1|yrs),control=list(maxIter=5000),REML=FALSE)
		test22<-lmer(rwi~SNDR+GPT+(0+SNDR|ind)+(0+GPT|ind)+(1|yrs),control=list(maxIter=5000),REML=FALSE)
		test23<-lmer(rwi~SWE+GST*GPT+(0+SWE|ind)+(0+GST|ind)+(0+GPT|ind)+(0+GST:GPT|ind)+(1|yrs),control=list(maxIter=5000),REML=FALSE)
		test24<-lmer(rwi~SWE+GST+GPT+(0+SWE|ind)+(0+GST|ind)+(0+GPT|ind)+(1|yrs),control=list(maxIter=5000),REML=FALSE)
		test25<-lmer(rwi~SNDR+GST*GPT+(0+SNDR|ind)+(0+GST|ind)+(0+GPT|ind)+(0+GST:GPT|ind)+(1|yrs),control=list(maxIter=5000),REML=FALSE)
		test26<-lmer(rwi~SNDR+GST+GPT+(0+SNDR|ind)+(0+GST|ind)+(0+GPT|ind)+(1|yrs),control=list(maxIter=5000),REML=FALSE)
		test27<-lmer(rwi~GDD+DST+GPT+DPT+(0+GDD|ind)+(0+DST|ind)+(0+GPT|ind)+(0+DPT|ind)+(1|yrs),control=list(maxIter=5000),REML=FALSE)
		test28<-lmer(rwi~GDD+GPT+DPT+(0+GDD|ind)+(0+GPT|ind)+(0+DPT|ind)+(1|yrs),control=list(maxIter=5000),REML=FALSE)
		test29<-lmer(rwi~GDD+GPT+(0+GDD|ind)+(0+GPT|ind)+(1|yrs),control=list(maxIter=5000),REML=FALSE)
		test30<-lmer(rwi~GDD*GPT+(0+GDD|ind)+(0+GPT|ind)+(0+GDD:GPT|ind)+(1|yrs),control=list(maxIter=5000),REML=FALSE)
		test31<-lmer(rwi~GDD+SWE+(0+GDD|ind)+(0+SWE|ind)+(1|yrs),control=list(maxIter=5000),REML=FALSE)
		print(species[i]); print(stands[j])
		print("rwi tests")
		#print(summary(test0))
		#print(summary(test1))
		#print(pvals.fnc(test1))
		#print(summary(test2))
		#print(summary(test5))
		#print(pvals.fnc(test5))
		#print(summary(test6))
		#print(pvals.fnc(test6))
		#print(summary(test7))
		#print(summary(test8))
		#print(summary(test9))
		#print(pvals.fnc(test9))
		#print(summary(test14))
		#print(summary(test16))
		#print(pvals.fnc(test16))
		#print(summary(test17))
		#print(pvals.fnc(test17))
		#print(summary(test18))
		#print(pvals.fnc(test18))
		#print(summary(test19))
		#print(pvals.fnc(test19))
		#print(summary(test20))
		#print(pvals.fnc(test20))
		#print(summary(test24))
		#print(pvals.fnc(test24))
		#print(summary(test27))
		#print(pvals.fnc(test27))
		#print(summary(test28))
		#print(pvals.fnc(test28))
		#print(summary(test29))
		#print(pvals.fnc(test29))
		#print(summary(test31))
		#print(pvals.fnc(test31))
		#print(anovatest)
		anovatest<-anova(test0,test1,test2,test3,test4,test5,test6,test7,test8,test9, test10,test11,test12,test13,test14,test15,test16,test17,test18,test19,test20,test21,test22,test23,test24,test25,test26,test27,test28,test29,test30,test31)
		anovatestord<-anovatest[order(rownames(anovatest)),]
		dev<-c(deviance(test0),deviance(test1),deviance(test10),deviance(test11),deviance(test12),deviance(test13),deviance(test14),deviance(test15),deviance(test16),deviance(test17),deviance(test18),deviance(test19),deviance(test2),deviance(test20),deviance(test21),deviance(test22),deviance(test23),deviance(test24),deviance(test25),deviance(test26),deviance(test27),deviance(test28),deviance(test29),deviance(test3),deviance(test30),deviance(test31),deviance(test4),deviance(test5),deviance(test6),deviance(test7),deviance(test8),deviance(test9))
		deltaAIC<-anovatestord[,2]-anovatestord[1,2]
		anovatest2<-cbind(anovatestord,deltaAIC,dev)
		rownum<-row(anovatest2)[which.min(anovatest2$AIC)]
		sitename<-stands[j]
		allsites<-cbind(sitename,allsites)
		otherbestmod<-anovatest2[which((anovatest2$AIC)<(2+anovatest2[rownum,2])),]#to see if any other simpler tests are within 2 AIC units
		simplebestmod<-otherbestmod[which.min(otherbestmod[,1]),]
		logliknull<-logLik(test0)
		nulldev<-deviance(test0)
		simplebestmod2<-cbind(rownames(simplebestmod),simplebestmod,logliknull,nulldev)
		colnames(simplebestmod2)<-c("bestmod",colnames	(simplebestmod),"logLiknull","nulldev")
		allsimplebestmod<-rbind(simplebestmod2,allsimplebestmod)
}
rownames(allsimplebestmod)<-allsites
(expdev<-round(100*((allsimplebestmod[,12]-allsimplebestmod[,10])/allsimplebestmod[,12]),digits=2))#to get the explained deviance of the best-fit models
allbestmod<-cbind(allsimplebestmod,expdev)
colnames(allbestmod)<-c(colnames(allsimplebestmod),"expdev")
print(allbestmod)
}


#Final models for each species (in "allsimplebestmod"):
#XANO: AM16=test18, AR07=test19,AE10=test24,PARA=test31
#TSME: AM16=test18, AR07=test31,AE10=test18,PARA=test31
#ABAM: TO04=test16,AG05=NULL,AV06=NULL,AM16=test27,AR07=test16,AE10=test1,PARA=test16
#TSHE: TO04=test16,TB13=test9,AG05=test5,AV06=test17,AX15=test20,AM16=test16
#THPL: TO04=test17,TB13=test16,AG05=NULL,AV06=test28,AX15=test29
#PSME: TO04=test6,TB13=test18,AG05=test6,AV06=test17,AX15=test31
#anovatest<-anova(test0,test1,test2,test3,test4,test5,test6,test7,test8,test9, test10,test11,test12,test13,test14,test15,test16,test17,test18,test19,test20,test21,test22,test23,test24,test25,test26,test27,test28,test29,test30,test31)


				
#Look at proportion of trees sensitive to each climate variable for the tests that were identified as best fit wit the linear mixed effects models
for (i in 1: length(species)){
	tmp<-dat.rwi.df[dat.rwi.df[,2]==species[i],1:dim(dat.rwi.df)[2]]
	propsig<-c()
	mnr2<-c()
	for(j in 1:length(stands)){
 		tmp2<-tmp[tmp[,1]==stands[j],]
 		if(dim(tmp2)[1]==0){next}
 		indtreesres<-c()
 			if(stands[j]=="PARA"){climst<-climstdPARA}
			if(stands[j]=="AE10"){climst<-climstdAE10}
			if(stands[j]=="AR07"){climst<-climstdAR07}
			if(stands[j]=="AV06"){climst<-climstdAV06}
			if(stands[j]=="AM16"){climst<-climstdAM16}
			if(stands[j]=="AX15"){climst<-climstdAX15}
			if(stands[j]=="AG05"){climst<-climstdAG05}
			if(stands[j]=="TB13"){climst<-climstdTB13}
			if(stands[j]=="TO04"){climst<-climstdTO04}

			yrs<-climst[,1]
			MAT<-climst[,2]
			PPT<-climst[,3]
			GPT<-climst[,4]
			DPT<-climst[,5]
			GST<-climst[,6]
			DST<-climst[,7]
			SWE<-climst[,8]
			SNDR<-climst[,9]
			GDD<-climst[,10]
  		for (k in 1:dim(tmp2)[[1]]){
 			indrwi<-as.numeric(tmp2[k,7:100])
 			temptest1<-lm(indrwi~MAT)
 			temptest5<-lm(indrwi~GST)
 			temptest6<-lm(indrwi~DST)
 			temptest9<-lm(indrwi~GST+DST)
 			temptest16<-lm(indrwi~SWE)
 			temptest17<-lm(indrwi~SNDR)
 			temptest18<-lm(indrwi~GDD)
 			temptest19<-lm(indrwi~SWE+GST)
 			temptest20<-lm(indrwi~SNDR+GST)
 			temptest24<-lm(indrwi~SWE+GST+GPT)
 			temptest27<-lm(indrwi~GDD+DST+GPT+DPT)
 			temptest28<-lm(indrwi~GDD+GPT+DPT)
 			temptest29<-lm(indrwi~GDD+GPT)
 			temptest31<-lm(indrwi~SWE+GDD)
			r1<-summary(temptest1)$r.squar
			r5<-summary(temptest5)$r.squar
			r6<-summary(temptest6)$r.squar
			r9<-summary(temptest9)$r.squar
			r16<-summary(temptest16)$r.squar
			r17<-summary(temptest17)$r.squar
			r18<-summary(temptest18)$r.squar
			r19<-summary(temptest19)$r.squar
			r20<-summary(temptest20)$r.squar
			r24<-summary(temptest24)$r.squar
			r27<-summary(temptest27)$r.squar
			r28<-summary(temptest28)$r.squar
			r29<-summary(temptest29)$r.squar
			r31<-summary(temptest31)$r.squar
			pvals1<-pf(summary(temptest1)$fstatistic[1],summary(temptest1)$fstatistic[2],summary(temptest1)$fstatistic[3],lower.tail=FALSE)
			pvals5<-pf(summary(temptest5)$fstatistic[1],summary(temptest5)$fstatistic[2],summary(temptest5)$fstatistic[3],lower.tail=FALSE)
			pvals6<-pf(summary(temptest6)$fstatistic[1],summary(temptest6)$fstatistic[2],summary(temptest6)$fstatistic[3],lower.tail=FALSE)
			pvals9<-pf(summary(temptest9)$fstatistic[1],summary(temptest9)$fstatistic[2],summary(temptest9)$fstatistic[3],lower.tail=FALSE)
			pvals16<-pf(summary(temptest16)$fstatistic[1],summary(temptest16)$fstatistic[2],summary(temptest16)$fstatistic[3],lower.tail=FALSE)
			pvals17<-pf(summary(temptest17)$fstatistic[1],summary(temptest17)$fstatistic[2],summary(temptest17)$fstatistic[3],lower.tail=FALSE)
			pvals18<-pf(summary(temptest18)$fstatistic[1],summary(temptest18)$fstatistic[2],summary(temptest18)$fstatistic[3],lower.tail=FALSE)
			pvals19<-pf(summary(temptest19)$fstatistic[1],summary(temptest19)$fstatistic[2],summary(temptest19)$fstatistic[3],lower.tail=FALSE)
			pvals20<-pf(summary(temptest20)$fstatistic[1],summary(temptest20)$fstatistic[2],summary(temptest20)$fstatistic[3],lower.tail=FALSE)
			pvals24<-pf(summary(temptest24)$fstatistic[1],summary(temptest24)$fstatistic[2],summary(temptest24)$fstatistic[3],lower.tail=FALSE)
			pvals27<-pf(summary(temptest27)$fstatistic[1],summary(temptest27)$fstatistic[2],summary(temptest27)$fstatistic[3],lower.tail=FALSE)
			pvals28<-pf(summary(temptest28)$fstatistic[1],summary(temptest28)$fstatistic[2],summary(temptest28)$fstatistic[3],lower.tail=FALSE)
			pvals29<-pf(summary(temptest29)$fstatistic[1],summary(temptest29)$fstatistic[2],summary(temptest29)$fstatistic[3],lower.tail=FALSE)
			pvals31<-pf(summary(temptest31)$fstatistic[1],summary(temptest31)$fstatistic[2],summary(temptest31)$fstatistic[3],lower.tail=FALSE)
			ifelse(pvals1<0.05,sensT<-1,sensT<-0)
			ifelse(pvals5<0.05,sensGT<-1,sensGT<-0)
			ifelse(pvals6<0.05,sensDT<-1,sensDT<-0)
			ifelse(pvals9<0.05,sensGTDT<-1,sensGTDT<-0)
			ifelse(pvals16<0.05,sensSWE<-1,sensSWE<-0)
			ifelse(pvals17<0.05,sensSNDR<-1,sensSNDR<-0)
			ifelse(pvals18<0.05,sensGDD<-1,sensGDD<-0)
			ifelse(pvals19<0.05,sensSWGT<-1,sensSWGT<-0)
			ifelse(pvals20<0.05,sensSDGT<-1,sensSDGT<-0)
			ifelse(pvals24<0.05,sensSWGTGP<-1,sensSWGTGP<-0)
			ifelse(pvals27<0.05,sensGDDTGPDP<-1,sensGDDTGPDP<-0)
			ifelse(pvals28<0.05,sensGDGPDP<-1,sensGDGPDP<-0)
			ifelse(pvals29<0.05,sensGDGP<-1,sensGDGP<-0)
			ifelse(pvals31<0.05,sensSWGD<-1,sensSWGD<-0)
			tmp1res<-c(round(r1,4),round(r5,4), round(r6,4), round(r9,4), round(r16,4), round(r17,4),round(r18,4),round(r19,4),round(r20,4), round(r24,4), round(r27,4), round(r28,4), round(r29,4),round(r31,4),sensT,sensGT,sensDT,sensGTDT,sensSWE,sensSNDR,sensGDD,sensSWGT,sensSDGT,sensSWGTGP,sensGDDTGPDP,sensGDGPDP,sensGDGP,sensSWGD)
			indtreesres<-rbind(indtreesres,tmp1res)
		}
		rownames(indtreesres)<-rownames(tmp2)
		colnames(indtreesres)<-c("Tr","GTr","DTr", "GTDTr", "SWEr","SNDRr","GDDr","SWGTr","SDGTr","SWGTGPr","GDDTGPDPr","GDGPDPr","GDGPr","SWGDr","Tsig","GTsig", "DTsig", "GTDTsig", "SWEsig","SNDRsig","GDDsig","SWGTsig","SDGTsig","SWGTGPsig","GDDTGPDPsig","GDGPDPsig","GDGPsig","SWGDsig")
propsigT<-round(sum(indtreesres[,15])/length(indtreesres[,8]),4)
propsigGT<-round(sum(indtreesres[,16])/length(indtreesres[,9]),4)
propsigDT<-round(sum(indtreesres[,17])/length(indtreesres[,10]),4)
propsigGTDT<-round(sum(indtreesres[,18])/length(indtreesres[,11]),4)
propsigSWE<-round(sum(indtreesres[,19])/length(indtreesres[,12]),4)
propsigSNDR<-round(sum(indtreesres[,20])/length(indtreesres[,13]),4)
propsigGDD<-round(sum(indtreesres[,21])/length(indtreesres[,14]),4)
propsigSWGT<-round(sum(indtreesres[,22])/length(indtreesres[,15]),4)
propsigSDGT<-round(sum(indtreesres[,23])/length(indtreesres[,16]),4)
propsigSWGTGP<-round(sum(indtreesres[,24])/length(indtreesres[,17]),4)
propsigGDDTGPDP<-round(sum(indtreesres[,25])/length(indtreesres[,18]),4)
propsigGDGPDP<-round(sum(indtreesres[,26])/length(indtreesres[,19]),4)
propsigGDGP<-round(sum(indtreesres[,27])/length(indtreesres[,20]),4)
propsigSWGD<-round(sum(indtreesres[,28])/length(indtreesres[,21]),4)
propsigALL<-c(as.character(stands[j]),propsigT,propsigGT,propsigDT,propsigGTDT,propsigSWE,propsigSNDR,propsigGDD,propsigSWGT,propsigSDGT,propsigSWGTGP,propsigGDDTGPDP,propsigGDGPDP,propsigGDGP,propsigSWGD)
propsig<-rbind(propsig,propsigALL)
mnr2T<-round(mean(indtreesres[,1]),3)
mnr2GT<-round(mean(indtreesres[,2]),3)
mnr2DT<-round(mean(indtreesres[,3]),3)
mnr2GTDT<-round(mean(indtreesres[,4]),3)
mnr2SWE<-round(mean(indtreesres[,5]),3)
mnr2SNDR<-round(mean(indtreesres[,6]),3)
mnr2GDD<-round(mean(indtreesres[,7]),3)
mnr2SWGT<-round(mean(indtreesres[,8]),3)
mnr2SDGT<-round(mean(indtreesres[,9]),3)
mnr2SWGTGP<-round(mean(indtreesres[,10]),3)
mnr2GDDTGPDP<-round(mean(indtreesres[,11]),3)
mnr2GDGPDP<-round(mean(indtreesres[,12]),3)
mnr2GDGP<-round(mean(indtreesres[,13]),3)
mnr2SWGD<-round(mean(indtreesres[,14]),3)
mnr2ALL<-c(as.character(stands[j]),mnr2T,mnr2GT,mnr2DT,mnr2GTDT,mnr2SWE,mnr2SNDR,mnr2GDD,mnr2SWGT,mnr2SDGT,mnr2SWGTGP,mnr2GDDTGPDP,mnr2GDGPDP,mnr2GDGP,mnr2SWGD)
mnr2<-rbind(mnr2,mnr2ALL)
}
colnames(propsig)<-c("stand","T","GT","DT", "GTDT", "SWE","SNDR","GDD","SWGT","SDGT","SWGTGP","GDDTGPDP","GDGPDP","GDGP","SWGD")
print(species[i])
print(propsig)
colnames(mnr2)<-c("stand","T","GT","DT", "GTDT", "SWE","SNDR","GDD","SWGT","SDGT","SWGTGP","GDDTGPDP","GDGPDP","GDGP","SWGD")
print(species[i])
print(mnr2)
}

############################################################
#For synchronized growth
stands<-c("TO04","TB13","AG05","AV06","AX15","AM16","AR07","AE10","PARA")
elev<-c(704,851,950,1064,1091,1197,1454,1460,1603)
species<-unique(dat.rwi[,2])
grwthcorall<-c()
for (i in 1: length(species)){
	tmp<-dat.rwi.df[dat.rwi.df[,2]==species[i],1:dim(dat.rwi.df)[2]]
	grwthcor<-c()
	col2<-c()
	for(j in 1:length(stands)){
 		tmp2<-tmp[tmp[,1]==stands[j],]
 		if(dim(tmp2)[1]==0){next}
		rwtmp<-tmp2[,7:105]
		rwtmp2<-t(rwtmp)
		r<- cor(rwtmp2, use = "na.or.complete")
		y<- which(lower.tri(r), TRUE)
		z <- data.frame(row = rownames(r)[y[, 1]],col = colnames(r)[y[, 2]],cor = r[y])
		tmpcor<-c(mean(z[,"cor"], na.rm=TRUE),sd(z[,"cor"], na.rm=TRUE), nrow(rwtmp))
		grwthcor<-rbind(grwthcor,tmpcor)
		#quartz()
		#hist(z[,"cor"], main = paste("",species[i],stands[j]))
		stdname<-stands[j]
		col2<-rbind(col2,stdname)
	}
col1<-rep(species[i], times=length(unique(tmp[,1])))
grwthcor2<-cbind(col1,col2,grwthcor)
grwthcorall<-rbind(grwthcorall,grwthcor2)
colnames(grwthcorall)<-c("species","stand","meancor","sdcor","numcors")
}
grwthcorall
stnderr<-as.numeric(grwthcorall[,3])/sqrt(as.numeric(grwthcorall[,5]))
allcors<-cbind(grwthcorall,stnderr)
allcors
########Figure2A for ABAM########################################
quartz(width=6, height=8)
par(mfrow = c(3,1),omi=c(1.5,0.1,0.1,0.1), mar=c(4,4.5,2,0.35))

#A
swecoefs<-tapply(as.numeric(test19allcoefs[,2]),test19allcoefs[,1],mean)
gstcoefs<-tapply(as.numeric(test19allcoefs[,3]),test19allcoefs[,1],mean)
tempcoefs<-c(NA,NA,NA,NA,0.022,0.048,0.079)
pptcoefs<-c(NA,NA,NA,NA,-0.035,-0.064,-0.060)
coefs<-matrix(c(tempcoefs,pptcoefs),2,7,byrow=TRUE)
barplot(coefs,beside=TRUE,col=c("red","blue"),ylim=c(-.1,.1),names.arg=c("","","","","","",""),xlab = "", ylab="Best-Fit Model Coefficient",cex.lab=1.3,cex.axis=1.2,axis.lty=1)
abline(h=0,lty=2)
mtext("A",side=3,line=.8,adj=0)

#2B. proportion of trees sensitive
abampropsig<-c(0.13,0,0,0.2,0.2,0.35,0.85)
abamelev<-c(704,950,1064,1197,1454,1460,1603)
barplot(abampropsig,xlab = "", ylab="Proportion of Trees Sensitive",names.arg=c("","","","","","",""),ylim=c(0,1.0),cex.lab=1.3,cex.axis=1.2,axis.lty=1,col="black")
mtext("B",side=3,line=.8,adj=0)

#2C mean correlations across elevations for ABAM
abamgrwthcor<-allcors[allcors[,1]=="Abam",]
abamgrwthcor2<-abamgrwthcor[-2,]#to remove TB13
abamelev<-c(704,950,1091,1197,1454,1460,1603)
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
stop("vectors must be same length")
arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
} 
y <- as.numeric(abamgrwthcor2[,3])
barx <- barplot(y, names.arg=abamelev,ylim=c(0,0.75), axis.lty=1, xlab="Elevation (m)", ylab="Mean correlation (r)",cex.lab=1.3,cex.axis=1.2,bty="l")
error.bar(barx,y, as.numeric(abamgrwthcor2[,6])) 
mtext("C",side=3,line=.8,adj=0)

########Figure2A for XANO########################################
quartz(width=6, height=7)
par(mfrow=c(3,1), mai=c(.57,.6,.2,0.1), omi=c(.31,.15,.15,.04), mgp=c(3.46,1,0))
#A
#use coefficients from test31 (GDD+SWE), arranged frmo low to high elevation
test31swecoefs<-c(-0.01437,-0.04991,-0.04843,-0.11066)
test31gddcoefs<-c(0.04476,0.01478,0.02937,0.03421)

test19swecoefs<-c(-0.02299,-0.05302,-0.05267,-0.11467 )
test19gstcoefs<-c(0.03693, 0.01059,0.02854,0.03166)

coefs<-matrix(c(test19gstcoefs,test19swecoefs),2,4,byrow=TRUE)
barplot(coefs,beside=TRUE,col=c("red","blue"),ylim=c(-.15,.08),names.arg=c("","","",""),xlab = "", ylab="Best-Fit Model Coefficient",cex.lab=1.3,cex.axis=1.3,axis.lty=1,las=1)
abline(h=0,lty=2)
mtext("A",side=3,line=.8,adj=0)

#2B. proportion of trees sensitive
xanopropsig<-c(0.35,0.45,0.65,0.95)
xanoelev<-c(1197,1454,1460,1603)
barplot(xanopropsig,xlab = "", ylab="Proportion Trees Sensitive",names.arg=c("","","",""),ylim=c(0,1.0),cex.lab=1.3,cex.axis=1.3,axis.lty=1,col="black",las=1)
mtext("B",side=3,line=.8,adj=0)

#2C mean correlations across elevations for ABAM
xanogrwthcor<-allcors[allcors[,1]=="Xano",]
xanoelev<-c(1197,1454,1460,1603)
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
stop("vectors must be same length")
arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
} 
y <- as.numeric(xanogrwthcor[,3])
barx <- barplot(y, names.arg=xanoelev,ylim=c(0,0.6), axis.lty=1, xlab="Elevation (m)", ylab="Mean correlation (r)",cex.lab=1.3,cex.axis=1.3,cex.names = 1.2,las=1)
error.bar(barx,y, as.numeric(xanogrwthcor[,6])) 
mtext("C",side=3,line=.8,adj=0)

######################New FIgure 2
#Select out a high elevation species (Xano) and a low-elevation sepcies (Tshe) and plot growth versus climate variables at upper and lower range margins.
stands<-c("TO04","AM16","PARA")
elev<-c(704,851,950,1064,1091,1197,1454,1460,1603)
r.elev<-elev[order(elev, decreasing=TRUE)]
#reorder stands by elevation
stands<-stands[order(elev, decreasing=TRUE)]
species<-unique(dat.rwi[,2])
#First, pull out xano and tshe data
	xanodat<-dat.rwi.df[dat.rwi[,2]=="Xano",1:dim(dat.rwi)[2]]
	tshedat<-dat.rwi.df[dat.rwi[,2]=="Tshe",1:dim(dat.rwi)[2]]
#now pull out upper and lower stands for Tshe and Xano
 		xanohighdat<-xanodat[xanodat[,1]=="PARA",]
 		xanohighmnrwi<-mean((xanohighdat[,7:100]),na.rm=TRUE) #pull out rwi for 1914-2007 and get mean
 		xanolodat<-xanodat[xanodat[,1]=="AM16",]
 		xanolomnrwi<-mean((xanolodat[,7:100]),na.rm=TRUE) 
 		tshehighdat<-tshedat[tshedat[,1]=="AM16",]
 		tshehighmnrwi<-mean((tshehighdat[,7:100]),na.rm=TRUE)
 		tshelodat<-tshedat[tshedat[,1]=="TO04",]
 		tshelomnrwi<-mean((tshelodat[,7:100]),na.rm=TRUE) 
 		swepara<-climstdPARA[,8]
 		sweam16<-climstdAM16[,8]
 		sweto04<-climstdTO04[,8]
 	quartz(width=7, height=7)
 	par(mfrow=c(2,2), mai=c(.75,.8,.25,0.1), omi=c(.3,.3,.3,.04), mgp=c(3.6,1,0))
	plot(swepara,xanohighmnrwi,pch=21,cex=1.2,bg="dark blue",ylab="Annual growth (mean RWI)",xlab="",bty="l",las=1,ylim=c(.4,1.6),xlim=c(-3,4), cex.lab=1.5, cex.axis=1.3)
 	abline(lm(xanohighmnrwi~swepara),lwd=2)
 	summary(lm(xanohighmnrwi~swepara),lwd=2)
 	mtext("R2=0.38", at=c(3), cex=1.2, line=-1.5)
 	mtext("p<0.001", at=c(3), cex=1.2, line=-3)
 	mtext("Upper Range Limit",side=3,line=.6,adj=0, cex=1.3)
 	
 	#mtext("A. Upper Range Limit",side=3,line=.6,adj=0)
 	#mtext("High Elevation Species (C. nootkatensis)",side=3,line=2,adj=0, cex=1.2)
 	plot(sweam16,xanolomnrwi,pch=21,cex=1.2,bg="dark blue",ylab="",xlab="",bty="l",las=1,ylim=c(.4,1.6),xlim=c(-3,4),cex.lab=1.5, cex.axis=1.3)
	abline(lm(xanolomnrwi~sweam16),lwd=2)
	summary(lm(xanolomnrwi~swepara),lwd=2)
	mtext("R2=0.06", at=c(3), cex=1.2, line=-1.5)
 	mtext("p=0.02", at=c(3), cex=1.2, line=-3)
 	mtext("Lower Range Limit",side=3,line=.6,adj=0, cex=1.3)

 	#mtext("B. Lower Range Limit",side=3,line=.6,adj=0)
 	plot(sweam16,tshehighmnrwi,pch=21,cex=1.2,bg="darkolivegreen3",ylab="Annual growth (mean RWI)",xlab="Snowpack (SWE, standardized)",bty="l",las=1,ylim=c(.4,1.6),xlim=c(-3,4),cex.lab=1.2)
 	abline(lm(tshehighmnrwi~sweam16),lwd=2)
 	mtext("C. Upper Range Limit",side=3,line=.6,adj=0)
 	mtext("Low Elevation Species (T. heterophylla)",side=3,line=2,adj=0, cex=1.2)
	plot(sweto04,tshelomnrwi,pch=21,cex=1.2,bg="darkolivegreen3",ylab="",xlab="Snowpack (SWE,standardized)",bty="l",las=1,ylim=c(.4,1.6),xlim=c(-3,4),cex.lab=1.2)
	abline(lm(tshelomnrwi~sweto04),lwd=2)
mtext("D. Mid-Range",side=3,line=.6,adj=0)
 	
 
 
 
 
 ####fig 2 version 2
 quartz(width=7, height=7)
 	par(mfrow=c(2,2), mai=c(.75,.75,.25,0.1), omi=c(.3,.3,.3,.04), mgp=c(3.6,1,0))
	plot(sweam16,tshehighmnrwi,pch=21,cex=1.2,bg="darkolivegreen3",ylab="Annual growth (mean RWI)",xlab="",bty="l",las=1,ylim=c(.4,1.6),xlim=c(-3,4),cex.lab=1.2)
 	abline(lm(tshehighmnrwi~sweam16),lwd=2)
 	mtext("A. Low Elevation Species (TSHE)",side=3,line=.6,adj=0) 
	mtext("Upper Range Limit",side=3,line=2,adj=0, cex=1.2)
	plot(swepara,xanohighmnrwi,pch=21,cex=1.2,bg="dark blue",ylab="",xlab="",bty="l",las=1,ylim=c(.4,1.6),xlim=c(-3,4))
 	abline(lm(xanohighmnrwi~swepara),lwd=2)
 	mtext("B. High Elevation Species (CANO)",side=3,line=.6,adj=0)
	plot(sweto04,tshelomnrwi,pch=21,cex=1.2,bg="darkolivegreen3",ylab="Annual growth (mean RWI)",xlab="Snowpack (SWE,standardized)",bty="l",las=1,ylim=c(.4,1.6),xlim=c(-3,4),cex.lab=1.2)
	abline(lm(tshelomnrwi~sweto04),lwd=2)
	mtext("Mid-Range and Lower Range Limit",side=3,line=2,adj=0, cex=1.2)
mtext("C. Low Elevation Species (TSHE)",side=3,line=.6,adj=0)
 	 plot(sweam16,xanolomnrwi,pch=21,cex=1.2,bg="dark blue",ylab="",xlab="Snowpack (SWE,standardized)",bty="l",las=1,ylim=c(.4,1.6),xlim=c(-3,4),cex.lab=1.2)
	abline(lm(xanolomnrwi~sweam16),lwd=2)
 	mtext("D. High Elevation Species (CANO)",side=3,line=.6,adj=0)
 
 ###Figure 3
lowspecies<-c("Thpl","Psme","Tshe")
highspecies<-c("Abam","Cano","Tsme")
lowelev<-c(704,851,950,1064,1091)
tsheelev<-c(704,851,950,1064,1091,1197)
abamelev<-c(704,950,1064,1197,1454,1460,1603)
highelev<-c(1197,1454,1460,1603)
psmeaic<-c(4,9,9,16,14)
thplaic<-c(43,8,0,17,7)
tsheaic<-c(7,10,3,9,23,4)
abamaic<-c(18,0,0,18,20,69,70)
canoaic<-c(15,39,32,67)
tsmeaic<-c(27,26,21,69)
psmeprop<-c(0.32,0.1,0.15,0.2,0.24)
thplprop<-c(0.21,0.25,0,0.5,0.27)
tsheprop<-c(0.1,0.2,0.1,0.15,0.6,0.15)
abamprop<-c(0.25,0,0,0.35,0.45,0.74,1.0)
canoprop<-c(0.35,0.45,0.65,0.95)
tsmeprop<-c(0.53,0.55,0.5,0.9)
psmesync<-c(0.25,.31,.34,.3,.41)
thplsync<-c(.26,.19,.31,.34,.3)
tshesync<-c(.19,.13,.14,.18,.38,.24)
abamsync<-c(.07,.23,.33,.34,.41,.38,.59)
canosync<-c(0.36,0.40,0.41,0.47)
tsmesync<-c(0.12,0.21,0.25,0.31)

quartz(width=7, height=6.5)
par(mfrow=c(3,2), mai=c(0.18,0.62,0.4,0.0001), omi=c(0.5,0.02,0.04,0.1), mgp=c(3.6,1,0))
plot(lowelev,psmeaic,type="l",ylab="Strength of climatic effect",xlab="",ylim=c(0,80),xlim=c(700,1200),col="dark orange",cex=1.4,bty="l",las=1,cex.lab=1.45,cex.axis=1.45, lwd=2,col.axis="white")
points(lowelev,psmeaic,bg="dark orange",pch=21,cex=1.5)
lines(lowelev,thplaic,col="red",lwd=2)
points(lowelev,thplaic,bg="red",pch=22,cex=1.5)
lines(tsheelev,tsheaic,col="darkolivegreen3",lwd=2)
points(tsheelev,tsheaic,bg="darkolivegreen3",pch=24,cex=1.5)
legend(700,80,lowspecies, cex=1.3,pt.lwd=1,pch=c(22,21,24),pt.bg=c("red","dark orange","darkolivegreen3"),pt.cex=1.5,bty="n")

mtext("(AICnull-AICbest)",side=2,line=2.7,adj=.5,cex=0.8)
mtext("A.",side=3,line=0.5,adj=0, cex=1.1)
mtext("Low Elevation Species",side=3,line=2,adj=0, cex=1.1)
axis(2,at=c(0,20,40,60,80),tick=TRUE, las=1,cex.axis=1.45)

plot(abamelev,abamaic,type="l",ylab="",xlab="",ylim=c(0,80),xlim=c(700,1650),col="dark green",cex=1.4,bty="l",las=1,cex.lab=1.45,cex.axis=1.45, lwd=2,col.axis="white")
points(abamelev,abamaic,bg="dark green",pch=22,cex=1.5)
lines(highelev,canoaic,col="dark blue",lwd=2)
points(highelev,canoaic,bg="dark blue",pch=21,cex=1.5)
lines(highelev,tsmeaic,col="purple",lwd=2)
points(highelev,tsmeaic,bg="purple",pch=24,cex=1.5)
mtext("B.",side=3,line=0.5,adj=0, cex=1.1)
legend(700,80,highspecies,cex=1.3,pt.lwd=1,pch=c(22,21,24),pt.bg=c("dark green","dark blue","purple"),pt.cex=1.5,bty="n")
mtext("High Elevation Species",side=3,line=2,adj=0, cex=1.1)
axis(2,at=c(0,20,40,60,80),tick=TRUE, las=1,cex.axis=1.45)

plot(lowelev,psmeprop,type="l",ylab="Proportion trees sensitive",xlab="",ylim=c(0,1.0),xlim=c(700,1200),col="dark orange",cex=1.4,bty="l",las=1,cex.lab=1.45,cex.axis=1.45, lwd=2,col.axis="white")
points(lowelev,psmeprop,bg="dark orange",pch=21,cex=1.5)
axis(2,at=c(0,0.20,0.40,0.60,0.80,1.0),tick=TRUE, las=1,cex.axis=1.45)
lines(lowelev,thplprop,col="red",lwd=2)
points(lowelev,thplprop,bg="red",pch=22,cex=1.5)
lines(tsheelev,tsheprop,col="darkolivegreen3",lwd=2)
points(tsheelev,tsheprop,bg="darkolivegreen3",pch=24,cex=1.5)
mtext("C.",side=3,line=0.5,adj=0, cex=1.1)

plot(abamelev,abamprop,type="l",ylab="",xlab="",ylim=c(0,1.0),xlim=c(700,1650),cex.lab=1.45,cex.axis=1.45,col="dark green",bty="l",las=1, lwd=2,col.axis="white")
points(abamelev,abamprop,bg="dark green",pch=22,cex=1.5)
lines(highelev,canoprop,col="dark blue",lwd=2)
points(highelev,canoprop,bg="dark blue",pch=21,cex=1.5)
lines(highelev,tsmeprop,col="purple",lwd=2)
points(highelev,tsmeprop,bg="purple",pch=24,cex=1.5)
mtext("D.",side=3,line=0.5,adj=0, cex=1.1)
axis(2,at=c(0,0.20,0.40,0.60,0.80,1.0),tick=TRUE, las=1,cex.axis=1.45)

plot(lowelev,psmesync,type="l",ylab="Growth synchrony",xlab="Elevation (m)",ylim=c(0,.6),xlim=c(700,1200),col="dark orange",cex=1.3,bty="l",las=1,cex.lab=1.45,cex.axis=1.45, lwd=2)
points(lowelev,psmesync,bg="dark orange",pch=21,cex=1.5)
lines(lowelev,thplsync,col="red",lwd=2)
points(lowelev,thplsync,bg="red",pch=22,cex=1.5)
lines(tsheelev,tshesync,col="darkolivegreen3",lwd=2)
points(tsheelev,tshesync,bg="darkolivegreen3",pch=24,cex=1.5)
mtext("(mean r)",side=2,line=2.7,adj=0.5,cex=0.8)
mtext("E.",side=3,line=0.5,adj=0, cex=1.1)
mtext("Elevation (m)",side=1,line=2.5,adj=.5, cex=1)

plot(abamelev,abamsync,type="l",ylab="",xlab="Elevation (m)",ylim=c(0,.6),xlim=c(700,1650),col="dark green",cex=1.3,bty="l",las=1,cex.lab=1.45,cex.axis=1.45, lwd=2)
points(abamelev,abamsync,bg="dark green",pch=22,cex=1.5)
lines(highelev,canosync,col="dark blue",lwd=2)
points(highelev,canosync,bg="dark blue",pch=21,cex=1.5)
lines(highelev,tsmesync,col="purple",lwd=2)
points(highelev,tsmesync,bg="purple",pch=24,cex=1.5)
mtext("F.",side=3,line=0.5,adj=0, cex=1.1)
#axis(1,at=c(800,1000,1200,1400,1600),tick=TRUE, las=1,cex.axis=1.45)
mtext("Elevation (m)",side=1,line=2.5,adj=.5, cex=1)

###alt figure 3 (with black lines and colored points)
quartz(width=7, height=6.5)
par(mfrow=c(3,2), mai=c(0.15,0.7,0.4,0.001), omi=c(0.5,0.02,0.04,0.2), mgp=c(3.6,1,0))
plot(lowelev,psmeaic,type="l",ylab="Strength of climatic effect",xlab="",ylim=c(0,80),xlim=c(700,1650),col="dark orange",cex=1.4,bty="l",las=1,lty=1,cex.lab=1.4,cex.axis=1.3, lwd=2,col.axis="white")
points(lowelev,psmeaic,bg="dark orange",pch=21,cex=1.5)
lines(lowelev,thplaic,lty=2,lwd=2)
points(lowelev,thplaic,bg="red",pch=22,cex=1.5)
lines(tsheelev,tsheaic,lty=3,lwd=2)
points(tsheelev,tsheaic,bg="darkolivegreen3",pch=24,cex=1.5)
legend(700,80,lowspecies, lty=c(1,2,3),pch=c(22,21,24),pt.bg=c("red","dark orange","darkolivegreen3"),pt.cex=1.5,bty="n")
mtext("(AICnull-AICbest)",side=2,line=2.7,adj=.5,cex=0.7)
mtext("A.",side=3,line=0.5,adj=0, cex=1.1)
mtext("Low Elevation Species",side=3,line=2,adj=0, cex=1.1)
axis(2,at=c(0,20,40,60,80),tick=TRUE, las=1,cex.axis=1.3)

plot(abamelev,abamaic,type="l",ylab="",xlab="",ylim=c(0,80),xlim=c(700,1650),cex=1.4,bty="l",las=1,cex.lab=1.4,cex.axis=1.1, lty=1,lwd=2,col.axis="white")
points(abamelev,abamaic,bg="dark green", pch=21,cex=1.5)
lines(highelev,canoaic,lwd=2,lty=2)
points(highelev,canoaic,bg="dark blue",pch=22,cex=1.5)
lines(highelev,tsmeaic,lty=3,lwd=2)
points(highelev,tsmeaic,bg="purple",pch=24,cex=1.5)
mtext("B.",side=3,line=0.5,adj=0, cex=1.1)
legend(700,80,highspecies,lty=c(1,2,3),pch=c(21,22,24),pt.bg=c("dark green","dark blue","purple"), pt.cex=1.5,bty="n")
mtext("High Elevation Species",side=3,line=2,adj=0, cex=1.1)

plot(lowelev,psmeprop,type="l",ylab="Proportion trees sensitive",xlab="",ylim=c(0,1.0),xlim=c(700,1650),lty=1,cex=1.4,bty="l",las=1,cex.lab=1.4,cex.axis=1.1, lwd=2,col.axis="white")
points(lowelev,psmeprop,bg="dark orange",pch=21,cex=1.5)
axis(2,at=c(0,0.20,0.40,0.60,0.80,1.0),tick=TRUE, las=1,cex.axis=1.3)
lines(lowelev,thplprop,lwd=2,lty=2)
points(lowelev,thplprop,bg="red",pch=22,cex=1.5)
lines(tsheelev,tsheprop,lwd=2,lty=3)
points(tsheelev,tsheprop,bg="darkolivegreen3",pch=24,cex=1.5)
mtext("C.",side=3,line=0.5,adj=0, cex=1.1)

plot(abamelev,abamprop,type="l",ylab="",xlab="",ylim=c(0,1.0),xlim=c(700,1650),lty=1,cex=1.3,bty="l",las=1,cex.lab=1.3,cex.axis=1.1, lwd=2,col.axis="white")
points(abamelev,abamprop,bg="dark green",pch=21,cex=1.5)
lines(highelev,canoprop,lty=2,lwd=2)
points(highelev,canoprop,bg="dark blue",pch=22,cex=1.5)
lines(highelev,tsmeprop,lty=3,lwd=2)
points(highelev,tsmeprop,bg="purple",pch=24,cex=1.5)
mtext("D.",side=3,line=0.5,adj=0, cex=1.1)

plot(lowelev,psmesync,type="l",ylab="Synchronicity",xlab="Elevation (m)",ylim=c(0,.6),xlim=c(700,1650),lty=1,cex=1.3,bty="l",las=1,cex.lab=1.4,cex.axis=1.3, lwd=2)
points(lowelev,psmesync,bg="dark orange",pch=21,cex=1.5)
lines(lowelev,thplsync,lty=2,lwd=2)
points(lowelev,thplsync,bg="red",pch=22,cex=1.5)
lines(tsheelev,tshesync,lty=3,lwd=2)
points(tsheelev,tshesync,bg="darkolivegreen3",pch=24,cex=1.5)
species<-c("Thpl","Psme","Tshe","Abam","Cano","Tsme")
mtext("(mean r)",side=2,line=2.7,adj=0.5,cex=0.7)
mtext("E.",side=3,line=0.5,adj=0, cex=1.1)
mtext("Elevation (m)",side=1,line=2.5,adj=.5, cex=1)

plot(abamelev,abamsync,type="l",ylab="",xlab="Elevation (m)",ylim=c(0,.6),xlim=c(700,1650),lty=1,cex=1.3,bty="l",las=1,cex.lab=1.4,cex.axis=1.1, lwd=2,col.axis="white")
points(abamelev,abamsync,bg="dark green",pch=21,cex=1.5)
lines(highelev,canosync,lty=2,lwd=2)
points(highelev,canosync,bg="dark blue",pch=22,cex=1.5)
lines(highelev,tsmesync,lty=3,lwd=2)
points(highelev,tsmesync,bg="purple",pch=24,cex=1.5)
mtext("F.",side=3,line=0.5,adj=0, cex=1.1)
axis(1,at=c(800,1000,1200,1400,1600),tick=TRUE, las=1,cex.axis=1.3)
mtext("Elevation (m)",side=1,line=2.5,adj=.5, cex=1)

###another alternative fig 3 (without points, and with markers at range limits)
quartz(width=7, height=6.5)
par(mfrow=c(3,2), mai=c(0.17,0.7,0.4,0.001), omi=c(0.5,0.02,0.04,0.2), mgp=c(3.6,1,0))
plot(lowelev,psmeaic,type="l",ylab="Strength of climatic effect",xlab="",ylim=c(0,80),xlim=c(700,1650),col="dark orange",cex=1.4,bty="l",las=1,cex.lab=1.45,cex.axis=1.45, lwd=2,col.axis="white")
points(lowelev[5],psmeaic[5],pch=124,cex=2,col="dark orange")
lines(lowelev,thplaic,col="red",lwd=2)
points(lowelev[5],thplaic[5],pch=124,cex=2,col="red")
lines(tsheelev,tsheaic,col="darkolivegreen3",lwd=2)
points(tsheelev[6],tsheaic[6],pch=124,cex=2,col="darkolivegreen3")
legend(700,80,lowspecies, lwd=1.5,col=c("red","dark orange","darkolivegreen3"),bty="n")
mtext("(AICnull-AICbest)",side=2,line=2.7,adj=.5,cex=0.8)
mtext("A.",side=3,line=0.5,adj=0, cex=1.1)
mtext("Low Elevation Species",side=3,line=2,adj=0, cex=1.1)
axis(2,at=c(0,20,40,60,80),tick=TRUE, las=1,cex.axis=1.3)

plot(abamelev,abamaic,type="l",ylab="",xlab="",ylim=c(0,80),xlim=c(700,1650),col="dark green",cex=1.4,bty="l",las=1,cex.lab=1.45,cex.axis=1.45, lwd=2,col.axis="white")
points(abamelev[8],abamaic[8],pch=124,cex=2,col="dark green")
points(abamelev[1],abamaic[1],pch=124,cex=2,col="dark green")
lines(highelev,canoaic,col="dark blue",lwd=2)
points(highelev[1],canoaic[1],pch=124,cex=2,col="dark blue")
points(highelev[4],canoaic[4],pch=124,cex=2,col="dark blue")
lines(highelev,tsmeaic,col="purple",lwd=2)
points(highelev[1],tsmeaic[1],pch=124,cex=2,col="purple")
points(highelev[4],tsmeaic[4],pch=124,cex=2,col="purple")
mtext("B.",side=3,line=0.5,adj=0, cex=1.1)
legend(700,80,highspecies,lwd=1.5,col=c("dark green","dark blue","purple"),bty="n")
mtext("High Elevation Species",side=3,line=2,adj=0, cex=1.1)

plot(lowelev,psmeprop,type="l",ylab="Proportion trees sensitive",xlab="",ylim=c(0,1.0),xlim=c(700,1650),col="dark orange",cex=1.4,bty="l",las=1,cex.lab=1.45,cex.axis=1.45, lwd=2,col.axis="white")
points(lowelev[5],psmeprop[5],pch=124,cex=2,col="dark orange")
axis(2,at=c(0,0.20,0.40,0.60,0.80,1.0),tick=TRUE, las=1,cex.axis=1.3)
lines(lowelev,thplprop,col="red",lwd=2)
points(lowelev[5],thplprop[5],col="red",pch=124,cex=2)
lines(tsheelev,tsheprop,col="darkolivegreen3",lwd=2)
points(tsheelev[6],tsheprop[6],col="darkolivegreen3",pch=124,cex=2)
mtext("C.",side=3,line=0.5,adj=0, cex=1.1)

plot(abamelev,abamprop,type="l",ylab="",xlab="",ylim=c(0,1.0),xlim=c(700,1650),col="dark green",cex=1.3,bty="l",las=1,cex.lab=1.45,cex.axis=1.45, lwd=2,col.axis="white")
points(abamelev[1],abamprop[1],col="dark green",pch=124,cex=2)
points(abamelev[8],abamprop[8],col="dark green",pch=124,cex=2)
lines(highelev,canoprop,col="dark blue",lwd=2)
points(highelev[1],canoprop[1],col="dark blue",pch=124,cex=2)
points(highelev[4],canoprop[4],col="dark blue",pch=124,cex=2)
lines(highelev,tsmeprop,col="purple",lwd=2)
points(highelev[1],tsmeprop[1],col="purple",pch=124,cex=2)
points(highelev[4],tsmeprop[4],col="purple",pch=124,cex=2)
mtext("D.",side=3,line=0.5,adj=0, cex=1.1)

plot(lowelev,psmesync,type="l",ylab="Synchronicity",xlab="Elevation (m)",ylim=c(0,.6),xlim=c(700,1650),col="dark orange",cex=1.3,bty="l",las=1,cex.lab=1.45,cex.axis=1.45, lwd=2)
points(lowelev[5],psmesync[5],col="dark orange",pch=124,cex=2)
lines(lowelev,thplsync,col="red",lwd=2)
points(lowelev[5],thplsync[5],col="red",pch=124,cex=2)
lines(tsheelev,tshesync,col="darkolivegreen3",lwd=2)
points(tsheelev[6],tshesync[6],col="darkolivegreen3",pch=124,cex=2)
mtext("(mean r)",side=2,line=2.7,adj=0.5,cex=0.8)
mtext("E.",side=3,line=0.5,adj=0, cex=1.1)
mtext("Elevation (m)",side=1,line=2.5,adj=.5, cex=1)

plot(abamelev,abamsync,type="l",ylab="",xlab="Elevation (m)",ylim=c(0,.6),xlim=c(700,1650),col="dark green",cex=1.3,bty="l",las=1,cex.lab=1.45,cex.axis=1.45, lwd=2,col.axis="white")
points(abamelev[1],abamsync[1],col="dark green",pch=124,cex=2)
points(abamelev[8],abamsync[8],col="dark green",pch=124,cex=2)
lines(highelev,canosync,col="dark blue",lwd=2)
points(highelev[1],canosync[1],col="dark blue",pch=124,cex=2)
points(highelev[4],canosync[4],col="dark blue",pch=124,cex=2)
lines(highelev,tsmesync,col="purple",lwd=2)
points(highelev[1],tsmesync[1],col="purple",pch=124,cex=2)
points(highelev[4],tsmesync[4],col="purple",pch=124,cex=2)
mtext("F.",side=3,line=0.5,adj=0, cex=1.1)
axis(1,at=c(800,1000,1200,1400,1600),tick=TRUE, las=1,cex.axis=1.45)
mtext("Elevation (m)",side=1,line=2.5,adj=.5, cex=1)



statsdat<-read.csv("statscor.csv", header=TRUE)
attach(statsdat)
head(statsdat)
aicsenstest<-lm(AIC~PropSens)
summary(aicsenstest)
aicsynchtest<-lm(AIC~Synch)
summary(aicsynchtest)
synchsenstest<-lm(Synch~PropSens)
summary(synchsenstest)
plot(AIC~PropSens, pch=16)
plot(AIC~Synch, pch=16)
plot(Synch~PropSens,pch=16)