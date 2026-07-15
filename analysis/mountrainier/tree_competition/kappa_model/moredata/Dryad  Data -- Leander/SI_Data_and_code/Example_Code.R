
#########################################################

###         LDL Anderegg and J Hille Ris Lambers        ###

#### Example code fore: Climatic and competitive constraints
#### on tree growth trade off at large scales but rarely 
#### explain local tree range boundaries

# additional data and full code available upon request from:
# Leander Anderegg (leanderegg@gmail.com)
#########################################################

# Note: this code shows worked examples to illustrate the analysis
# in Anderegg & HilleRisLambers (in submission). Due to the 
# size of the analysis, the goal of this code is to illustrate
# how critical analyses were performed, with a workable 
# example from Abies lasiocarpa from Colorado.
# These analyses were then applied to all species-replicates in the study.
# This example requires:
#  - Dataset S1: Mean growth and competitive environment data for CO A. lasiocarpa  
#  - Dataset S2: Raw tree ring widths for CO A. lasiocarpa
#  - Dataset S3: Summary data from mean growth, sensitivity to competition, growth synchrony, recruitment and survival analyses for all species-replicates
# provided as supplemental information with the submission
# 


# load required packages
library(betareg)
library(lmtest)
library(multcomp)
library(nlme)
library(lme4)
library(lmerTest)
library(car)
library(RColorBrewer)
library(dplyr)






#_____________________________________________________________________________________________
#_____________________________________________________________________________________________
########## BEGIN: Basal Area Growth ~ F(elev + comp)###########
#_____________________________________________________________________________________________
#_____________________________________________________________________________________________

# Goal: find the best model for predicting average BAI (2003-2012) as a function of Elevation band ($Band),
# tree size ($DBH), and the best single or combination of competitive metrics:
  # $BA_totLR = Basal area (from wedge prism) of live and recent dead trees around the focal tree
  # $N_Cr = # of crowns touching the focal tree
  # $in5_tot = # of trees w/in 5m of the focal tree
  # $ACFsc = active crown fraction of the focal tree
# plus their interactions with Band (if the effect of competition varies with elevation)

# These models allow estimation of size- and competition-standardized mean BAI
# as well as an estimate of the strength of competitive suppression for each elevation band.

# So I've made a suite of candidate models, and then automated model selection by AIC to get the best
# competition metric for each species. The family of functions takes a Treeinfo dataset and spits out
# the delta AIC table for all models <2 AIC from best model and returns the model object for the best model.
# I iteratively select the best model after determining the best variance structure 
# and whether a mixed model is necessary (with tree Pair as a random effect)
# This iterative approach is based on the suggested methods in Zuur et al. 2009.
# I then test the significance of an interaction term via a Likelihood Ratio Test




########### Create Model Selection Functions #############

# using lme with {nlme} for R side modeling with random effect, varname = BAI10yr
compselection.lme <- function(dataz, vf = varPower(form=~DBH), elevlevels = list(L="L", M="M", H="H"), aic.cutoff=4){
  levels(dataz$Band) <- elevlevels
  dataz$DBHsc <- scale(dataz$DBH) # scale predictors to avoid numerical instability
  dataz$BA_totLRsc <- scale(dataz$BA_totLR) # and to get BAI at species mean predictor level (intercept)
  dataz$N_Crsc <- scale(dataz$N_Cr)
  dataz$in5_totsc <- scale(dataz$in5_tot)
  dataz$ACFsc <- scale(dataz$ACF)
  dataz$BAI10yrsc <- scale(dataz$BAI10yr)
  na.flags <- apply(dataz[,c("DBH","BA_totLR", "N_Cr", "in5_tot", "ACF", "BAI10yr")], MARGIN=2, FUN= function(x){length(which(is.na(x)))})
  if(sum(na.flags>0)){print(paste("NAs in:", names(na.flags)[which(na.flags>0)]))}
  null0 <- lme(BAI10yr ~ Band , random = ~1|PairUn, weights = vf, data=dataz, method="ML")
  null1 <- lme(BAI10yr ~ Band + DBHsc, random = ~1|PairUn , weights = vf, data=dataz, method="ML")
  comp1 <- lme(BAI10yr ~ Band + DBHsc + N_Crsc, random = ~1|PairUn , weights = vf, data=dataz, method="ML")
  comp2 <- lme(BAI10yr ~ Band + DBHsc + BA_totLRsc, random = ~1|PairUn , weights = vf, data=dataz, method="ML")
  comp3 <- lme(BAI10yr ~ Band + DBHsc + in5_totsc, random = ~1|PairUn , weights = vf, data=dataz, method="ML")
  comp4 <- lme(BAI10yr ~ Band + DBHsc + ACFsc, random = ~1|PairUn , weights = vf, data=dataz, method="ML")
  comps1 <- lme(BAI10yr ~ Band + DBHsc + BA_sameLR + BA_otherLR, random = ~1|PairUn , weights = vf, data=dataz, method="ML")
  comps2 <- lme(BAI10yr ~ Band + DBHsc + in5_same + in5_other, random = ~1|PairUn , weights = vf, data=dataz, method="ML")
  compall1 <- lme(BAI10yr ~ Band + DBHsc + N_Crsc + BA_totLRsc , random = ~1|PairUn, weights = vf, data=dataz, method="ML")
  compall2 <- lme(BAI10yr ~ Band + DBHsc + N_Crsc + in5_totsc  , random = ~1|PairUn, weights = vf, data=dataz, method="ML")
  compall3 <- lme(BAI10yr ~ Band + DBHsc + N_Crsc + ACFsc , random = ~1|PairUn, weights = vf, data=dataz, method="ML")
  compall4 <- lme(BAI10yr ~ Band + DBHsc + ACFsc + BA_totLRsc , random = ~1|PairUn, weights = vf, data=dataz, method="ML")
  compall5 <- lme(BAI10yr ~ Band + DBHsc + ACFsc + in5_totsc , random = ~1|PairUn, weights = vf, data=dataz, method="ML")
  compall6 <- lme(BAI10yr ~ Band + DBHsc + BA_totLRsc + in5_totsc , random = ~1|PairUn, weights = vf, data=dataz, method="ML")
  compall7 <- lme(BAI10yr ~ Band + DBHsc + BA_totLRsc + in5_totsc + N_Crsc , random = ~1|PairUn, weights = vf, data=dataz, method="ML")
  compall8 <- lme(BAI10yr ~ Band + DBHsc + BA_totLRsc + in5_totsc + ACFsc , random = ~1|PairUn, weights = vf, data=dataz, method="ML")
  compall9 <- lme(BAI10yr ~ Band + DBHsc + BA_totLRsc + N_Crsc + ACFsc , random = ~1|PairUn, weights = vf, data=dataz, method="ML")
  compall10 <- lme(BAI10yr ~ Band + DBHsc + in5_totsc + N_Crsc + ACFsc , random = ~1|PairUn, weights = vf, data=dataz, method="ML")
  compall11 <- lme(BAI10yr ~ Band + DBHsc + BA_totLRsc + in5_totsc + N_Crsc + ACFsc , random = ~1|PairUn, weights = vf, data=dataz, method="ML")
  compint1 <- lme(BAI10yr ~ Band * N_Crsc + DBHsc , random = ~1|PairUn, weights = vf, data=dataz, method="ML")
  compint2 <- lme(BAI10yr ~ Band * BA_totLRsc + DBHsc , random = ~1|PairUn, weights = vf, data=dataz, method="ML")
  compint3 <- lme(BAI10yr ~ Band * in5_totsc +  DBHsc , random = ~1|PairUn, weights = vf, data=dataz, method="ML")
  compint4 <- lme(BAI10yr ~ Band * ACFsc + DBHsc , random = ~1|PairUn, weights = vf, data=dataz, method="ML")
  compint5 <- lme(BAI10yr ~ Band * N_Crsc + BA_totLRsc + Band:BA_totLRsc + DBHsc , random = ~1|PairUn, weights = vf, data=dataz, method="ML")
  compint6 <- lme(BAI10yr ~ Band * N_Crsc + in5_totsc + Band:in5_totsc + DBHsc , random = ~1|PairUn, weights = vf, data=dataz, method="ML")
  compint7 <- lme(BAI10yr ~ Band * N_Crsc + ACFsc + Band:ACFsc + DBHsc , random = ~1|PairUn, weights = vf, data=dataz, method="ML")
  compint8 <- lme(BAI10yr ~ Band * BA_totLRsc + in5_totsc + in5_totsc:Band +  DBHsc , random = ~1|PairUn, weights = vf, data=dataz, method="ML")
  compint9 <- lme(BAI10yr ~ Band * BA_totLRsc + ACFsc + ACFsc:Band +  DBHsc , random = ~1|PairUn, weights = vf, data=dataz, method="ML")
  compint10 <- lme(BAI10yr ~ Band * in5_totsc +  ACFsc + ACFsc:Band + DBHsc , random = ~1|PairUn, weights = vf, data=dataz, method="ML")
  compint11 <- lme(BAI10yr ~ Band + DBHsc + BA_totLRsc + in5_totsc + N_Crsc + ACFsc , random = ~1|PairUn, weights = vf, data=dataz, method="ML")
  
  ndcomp1 <- lme(BAI10yr ~ Band  + N_Crsc , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  ndcomp2 <- lme(BAI10yr ~ Band  + BA_totLRsc , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  ndcomp3 <- lme(BAI10yr ~ Band  + in5_totsc , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  ndcomp4 <- lme(BAI10yr ~ Band  + ACFsc , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  ndcomps1 <- lme(BAI10yr ~ Band  + BA_sameLR + BA_otherLR , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  ndcomps2 <- lme(BAI10yr ~ Band  + in5_same + in5_other , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  ndcompall1 <- lme(BAI10yr ~ Band  + N_Crsc + BA_totLRsc , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  ndcompall2 <- lme(BAI10yr ~ Band  + N_Crsc + in5_totsc , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  ndcompall3 <- lme(BAI10yr ~ Band  + N_Crsc + ACFsc , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  ndcompall4 <- lme(BAI10yr ~ Band  + ACFsc + BA_totLRsc , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  ndcompall5 <- lme(BAI10yr ~ Band  + ACFsc + in5_totsc , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  ndcompall6 <- lme(BAI10yr ~ Band  + BA_totLRsc + in5_totsc , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  ndcompall7 <- lme(BAI10yr ~ Band  + BA_totLRsc + in5_totsc + N_Crsc , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  ndcompall8 <- lme(BAI10yr ~ Band  + BA_totLRsc + in5_totsc + ACFsc , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  ndcompall9 <- lme(BAI10yr ~ Band  + BA_totLRsc + N_Crsc + ACFsc , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  ndcompall10 <- lme(BAI10yr ~ Band  + in5_totsc + N_Crsc + ACFsc , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  ndcompall11 <- lme(BAI10yr ~ Band  + BA_totLRsc + in5_totsc + N_Crsc + ACFsc , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  ndcompint1 <- lme(BAI10yr ~ Band * N_Crsc  , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  ndcompint2 <- lme(BAI10yr ~ Band * BA_totLRsc  , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  ndcompint3 <- lme(BAI10yr ~ Band * in5_totsc  , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  ndcompint4 <- lme(BAI10yr ~ Band * ACFsc  , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  ndcompint5 <- lme(BAI10yr ~ Band * N_Crsc + BA_totLRsc + Band:BA_totLRsc  , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  ndcompint6 <- lme(BAI10yr ~ Band * N_Crsc + in5_totsc + Band:in5_totsc  , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  ndcompint7 <- lme(BAI10yr ~ Band * N_Crsc + ACFsc + Band:ACFsc  , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  ndcompint8 <- lme(BAI10yr ~ Band * BA_totLRsc + in5_totsc + in5_totsc:Band , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  ndcompint9 <- lme(BAI10yr ~ Band * BA_totLRsc + ACFsc + ACFsc:Band , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  ndcompint10 <- lme(BAI10yr ~ Band * in5_totsc +  ACFsc + ACFsc:Band  , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  ## added 07.18.2016
  compint11 <- lme(BAI10yr ~ Band * BA_totLRsc + in5_totsc + Band:in5_totsc + N_Crsc + Band:N_Crsc + DBHsc , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  compint12 <- lme(BAI10yr ~ Band * BA_totLRsc + in5_totsc + Band:in5_totsc + ACFsc + Band:ACFsc + DBHsc , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  compint13 <- lme(BAI10yr ~ Band * BA_totLRsc + N_Crsc + Band:N_Crsc + ACFsc + Band:ACFsc + DBHsc , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  compint14 <- lme(BAI10yr ~ Band * in5_totsc + N_Crsc + Band:N_Crsc + ACFsc + Band:ACFsc + DBHsc , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  compint15 <- lme(BAI10yr ~ Band * BA_totLRsc + in5_totsc + Band:in5_totsc+ N_Crsc + Band:N_Crsc + ACFsc + Band:ACFsc + DBHsc , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  ndcompint11 <- lme(BAI10yr ~ Band * BA_totLRsc + in5_totsc + Band:in5_totsc + N_Crsc + Band:N_Crsc , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  ndcompint12 <- lme(BAI10yr ~ Band * BA_totLRsc + in5_totsc + Band:in5_totsc + ACFsc + Band:ACFsc , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  ndcompint13 <- lme(BAI10yr ~ Band * BA_totLRsc + N_Crsc + Band:N_Crsc + ACFsc + Band:ACFsc , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  ndcompint14 <- lme(BAI10yr ~ Band * in5_totsc + N_Crsc + Band:N_Crsc + ACFsc + Band:ACFsc , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  ndcompint15 <- lme(BAI10yr ~ Band * BA_totLRsc + in5_totsc + Band:in5_totsc+ N_Crsc + Band:N_Crsc + ACFsc + Band:ACFsc , random = ~1|PairUn, weights=vf, data=dataz, method="ML")
  ## added in 07/29/16 to test single interactions in models with two predictors
  compint5.1 <- lme(BAI10yr ~ Band * N_Crsc + BA_totLRsc + DBHsc , weights = vf, data=dataz, method="ML", random = ~1|PairUn)
  compint6.1 <- lme(BAI10yr ~ Band * N_Crsc + in5_totsc + DBHsc , weights = vf, data=dataz, method="ML", random = ~1|PairUn)
  compint7.1 <- lme(BAI10yr ~ Band * N_Crsc + ACFsc + DBHsc , weights = vf, data=dataz, method="ML", random = ~1|PairUn)
  compint8.1 <- lme(BAI10yr ~ Band * BA_totLRsc + in5_totsc +  DBHsc , weights = vf, data=dataz, method="ML", random = ~1|PairUn)
  compint9.1 <- lme(BAI10yr ~ Band * BA_totLRsc + ACFsc +  DBHsc , weights = vf, data=dataz, method="ML", random = ~1|PairUn)
  compint10.1 <- lme(BAI10yr ~ Band * in5_totsc +  ACFsc + DBHsc , weights = vf, data=dataz, method="ML", random = ~1|PairUn)
  compint5.2 <- lme(BAI10yr ~ Band + N_Crsc + BA_totLRsc + Band:BA_totLRsc + DBHsc , weights = vf, data=dataz, method="ML", random = ~1|PairUn)
  compint6.2 <- lme(BAI10yr ~ Band + N_Crsc + in5_totsc + Band:in5_totsc + DBHsc , weights = vf, data=dataz, method="ML", random = ~1|PairUn)
  compint7.2 <- lme(BAI10yr ~ Band + N_Crsc + ACFsc + Band:ACFsc + DBHsc , weights = vf, data=dataz, method="ML", random = ~1|PairUn)
  compint8.2 <- lme(BAI10yr ~ Band + BA_totLRsc + in5_totsc + Band:in5_totsc +  DBHsc , weights = vf, data=dataz, method="ML", random = ~1|PairUn)
  compint9.2 <- lme(BAI10yr ~ Band + BA_totLRsc + ACFsc + ACFsc:Band +  DBHsc , weights = vf, data=dataz, method="ML", random = ~1|PairUn)
  compint10.2 <- lme(BAI10yr ~ Band + in5_totsc +  ACFsc + ACFsc:Band + DBHsc , weights = vf, data=dataz, method="ML", random = ~1|PairUn)
  ndcompint5.1 <- lme(BAI10yr ~ Band * N_Crsc + BA_totLRsc  , weights = vf, data=dataz, method="ML", random = ~1|PairUn)
  ndcompint6.1 <- lme(BAI10yr ~ Band * N_Crsc + in5_totsc  , weights = vf, data=dataz, method="ML", random = ~1|PairUn)
  ndcompint7.1 <- lme(BAI10yr ~ Band * N_Crsc + ACFsc  , weights = vf, data=dataz, method="ML", random = ~1|PairUn)
  ndcompint8.1 <- lme(BAI10yr ~ Band * BA_totLRsc + in5_totsc  , weights = vf, data=dataz, method="ML", random = ~1|PairUn)
  ndcompint9.1 <- lme(BAI10yr ~ Band * BA_totLRsc + ACFsc  , weights = vf, data=dataz, method="ML", random = ~1|PairUn)
  ndcompint10.1 <- lme(BAI10yr ~ Band * in5_totsc +  ACFsc  , weights = vf, data=dataz, method="ML", random = ~1|PairUn)
  ndcompint5.2 <- lme(BAI10yr ~ Band + N_Crsc + BA_totLRsc + Band:BA_totLRsc  , weights = vf, data=dataz, method="ML", random = ~1|PairUn)
  ndcompint6.2 <- lme(BAI10yr ~ Band + N_Crsc + in5_totsc + Band:in5_totsc  , weights = vf, data=dataz, method="ML", random = ~1|PairUn)
  ndcompint7.2 <- lme(BAI10yr ~ Band + N_Crsc + ACFsc + Band:ACFsc  , weights = vf, data=dataz, method="ML", random = ~1|PairUn)
  ndcompint8.2 <- lme(BAI10yr ~ Band + BA_totLRsc + in5_totsc + in5_totsc:Band  , weights = vf, data=dataz, method="ML", random = ~1|PairUn)
  ndcompint9.2 <- lme(BAI10yr ~ Band + BA_totLRsc + ACFsc + ACFsc:Band  , weights = vf, data=dataz, method="ML", random = ~1|PairUn)
  ndcompint10.2 <- lme(BAI10yr ~ Band + in5_totsc +  ACFsc + ACFsc:Band  , weights = vf, data=dataz, method="ML", random = ~1|PairUn)
  
  
  aics <- AIC(null0, null1, comp1, comp2,comp3,comp4,comps1,comps2,compall1, compall2,compall3,compall4,compall5,compall6,compall7,compall8,compall9,compall10,compall11,compint1, compint2, compint3, compint4, compint5, compint6, compint7, compint8, compint9, compint10, compint11, compint12, compint13, compint14, compint15
              ,ndcomp1, ndcomp2,ndcomp3,ndcomp4,ndcomps1,ndcomps2,ndcompall1, ndcompall2,ndcompall3,ndcompall4,ndcompall5,ndcompall6,ndcompall7,ndcompall8,ndcompall9,ndcompall10,ndcompall11,ndcompint1, ndcompint2, ndcompint3, ndcompint4, ndcompint5, ndcompint6, ndcompint7, ndcompint8, ndcompint9, ndcompint10, ndcompint11, ndcompint12, ndcompint13, ndcompint14, ndcompint15
              ,compint5.1, compint6.1, compint7.1, compint8.1, compint9.1, compint10.1, compint5.2, compint6.2, compint7.2, compint8.2, compint9.2, compint10.2
              ,ndcompint5.1, ndcompint6.1, ndcompint7.1, ndcompint8.1, ndcompint9.1, ndcompint10.1, ndcompint5.2, ndcompint6.2, ndcompint7.2, ndcompint8.2, ndcompint9.2, ndcompint10.2)
  
  bics <- BIC(null0, null1, comp1, comp2,comp3,comp4,comps1,comps2,compall1, compall2,compall3,compall4,compall5,compall6,compall7,compall8,compall9,compall10,compall11,compint1, compint2, compint3, compint4, compint5, compint6, compint7, compint8, compint9, compint10, compint11, compint12, compint13, compint14, compint15
              ,ndcomp1, ndcomp2,ndcomp3,ndcomp4,ndcomps1,ndcomps2,ndcompall1, ndcompall2,ndcompall3,ndcompall4,ndcompall5,ndcompall6,ndcompall7,ndcompall8,ndcompall9,ndcompall10,ndcompall11,ndcompint1, ndcompint2, ndcompint3, ndcompint4, ndcompint5, ndcompint6, ndcompint7, ndcompint8, ndcompint9, ndcompint10, ndcompint11, ndcompint12, ndcompint13, ndcompint14, ndcompint15
              ,compint5.1, compint6.1, compint7.1, compint8.1, compint9.1, compint10.1, compint5.2, compint6.2, compint7.2, compint8.2, compint9.2, compint10.2
              ,ndcompint5.1, ndcompint6.1, ndcompint7.1, ndcompint8.1, ndcompint9.1, ndcompint10.1, ndcompint5.2, ndcompint6.2, ndcompint7.2, ndcompint8.2, ndcompint9.2, ndcompint10.2)
  
  
  
  aics$deltaAIC <- aics$AIC - min(aics$AIC)
  aics$BIC <- bics$BIC
  aics$deltaBIC <- aics$BIC - min(aics$BIC)
  print(aics[which(aics$deltaAIC <= aic.cutoff),])
  aics2 <- aics[which(aics$deltaAIC<2),]
  aics2$mod <- rownames(aics2)
  aicsbest <- arrange(aics2, df, deltaAIC )
  print(summary(get(aicsbest$mod[1])))
  
  return(get(aicsbest$mod[1]))
}

### compselection non-mixed model functions
# using gls with {nlme}, so can have R side modeling. 
compselection.gls <- function(dataz, vf = varPower(form=~DBH), elevlevels = list(L="L", M="M", H="H"), aic.cutoff=4, ...){
  levels(dataz$Band) <- elevlevels
  dataz$DBHsc <- scale(dataz$DBH) # scale predictors to avoid numerical instability
  dataz$BA_totLRsc <- scale(dataz$BA_totLR) # and to get BAI at species mean predictor level (intercept)
  dataz$N_Crsc <- scale(dataz$N_Cr)
  dataz$in5_totsc <- scale(dataz$in5_tot)
  dataz$ACFsc <- scale(dataz$ACF)
  dataz$BAI10yrsc <- scale(dataz$BAI10yr)
  na.flags <- apply(dataz[,c("DBH","BA_totLR", "N_Cr", "in5_tot", "ACF", "BAI10yr")], MARGIN=2, FUN= function(x){length(which(is.na(x)))})
  if(sum(na.flags>0)){print(paste("NAs in:", names(na.flags)[which(na.flags>0)]))}
  null0 <- gls(BAI10yr ~ Band , weights = vf, data=dataz, method = "ML",...)
  null1 <- gls(BAI10yr ~ Band + DBHsc , weights = vf, data=dataz, method="ML",...)
  comp1 <- gls(BAI10yr ~ Band + DBHsc + N_Crsc , weights = vf, data=dataz, method="ML",...)
  comp2 <- gls(BAI10yr ~ Band + DBHsc + BA_totLRsc , weights = vf, data=dataz, method="ML",...)
  comp3 <- gls(BAI10yr ~ Band + DBHsc + in5_totsc , weights = vf, data=dataz, method="ML",...)
  comp4 <- gls(BAI10yr ~ Band + DBHsc + ACFsc , weights = vf, data=dataz, method="ML",...)
  comps1 <- gls(BAI10yr ~ Band + DBHsc + BA_sameLR + BA_otherLR , weights = vf, data=dataz, method="ML",...)
  comps2 <- gls(BAI10yr ~ Band + DBHsc + in5_same + in5_other , weights = vf, data=dataz, method="ML",...)
  compall1 <- gls(BAI10yr ~ Band + DBHsc + N_Crsc + BA_totLRsc , weights = vf, data=dataz, method="ML",...)
  compall2 <- gls(BAI10yr ~ Band + DBHsc + N_Crsc + in5_totsc , weights = vf, data=dataz, method="ML",...)
  compall3 <- gls(BAI10yr ~ Band + DBHsc + N_Crsc + ACFsc , weights = vf, data=dataz, method="ML",...)
  compall4 <- gls(BAI10yr ~ Band + DBHsc + ACFsc + BA_totLRsc , weights = vf, data=dataz, method="ML",...)
  compall5 <- gls(BAI10yr ~ Band + DBHsc + ACFsc + in5_totsc , weights = vf, data=dataz, method="ML",...)
  compall6 <- gls(BAI10yr ~ Band + DBHsc + BA_totLRsc + in5_totsc , weights = vf, data=dataz, method="ML",...)
  compall7 <- gls(BAI10yr ~ Band + DBHsc + BA_totLRsc + in5_totsc + N_Crsc , weights = vf, data=dataz, method="ML",...)
  compall8 <- gls(BAI10yr ~ Band + DBHsc + BA_totLRsc + in5_totsc + ACFsc , weights = vf, data=dataz, method="ML",...)
  compall9 <- gls(BAI10yr ~ Band + DBHsc + BA_totLRsc + N_Crsc + ACFsc , weights = vf, data=dataz, method="ML",...)
  compall10 <- gls(BAI10yr ~ Band + DBHsc + in5_totsc + N_Crsc + ACFsc , weights = vf, data=dataz, method="ML",...)
  compall11 <- gls(BAI10yr ~ Band + DBHsc + BA_totLRsc + in5_totsc + N_Crsc + ACFsc , weights = vf, data=dataz, method="ML",...)
  compint1 <- gls(BAI10yr ~ Band * N_Crsc + DBHsc , weights = vf, data=dataz, method="ML",...)
  compint2 <- gls(BAI10yr ~ Band * BA_totLRsc + DBHsc , weights = vf, data=dataz, method="ML",...)
  compint3 <- gls(BAI10yr ~ Band * in5_totsc +  DBHsc , weights = vf, data=dataz, method="ML",...)
  compint4 <- gls(BAI10yr ~ Band * ACFsc + DBHsc , weights = vf, data=dataz, method="ML",...)
  compint5 <- gls(BAI10yr ~ Band * N_Crsc + BA_totLRsc + Band:BA_totLRsc + DBHsc , weights = vf, data=dataz, method="ML",...)
  compint6 <- gls(BAI10yr ~ Band * N_Crsc + in5_totsc + Band:in5_totsc + DBHsc , weights = vf, data=dataz, method="ML",...)
  compint7 <- gls(BAI10yr ~ Band * N_Crsc + ACFsc + Band:ACFsc + DBHsc , weights = vf, data=dataz, method="ML",...)
  compint8 <- gls(BAI10yr ~ Band * BA_totLRsc + in5_totsc + in5_totsc:Band +  DBHsc , weights = vf, data=dataz, method="ML",...)
  compint9 <- gls(BAI10yr ~ Band * BA_totLRsc + ACFsc + ACFsc:Band +  DBHsc , weights = vf, data=dataz, method="ML",...)
  compint10 <- gls(BAI10yr ~ Band * in5_totsc +  ACFsc + ACFsc:Band + DBHsc , weights = vf, data=dataz, method="ML",...)
  
  ndcomp1 <- gls(BAI10yr ~ Band  + N_Crsc , weights=vf, data=dataz, method="ML",...)
  ndcomp2 <- gls(BAI10yr ~ Band  + BA_totLRsc , weights=vf, data=dataz, method="ML",...)
  ndcomp3 <- gls(BAI10yr ~ Band  + in5_totsc , weights=vf, data=dataz, method="ML",...)
  ndcomp4 <- gls(BAI10yr ~ Band  + ACFsc , weights=vf, data=dataz, method="ML",...)
  ndcomps1 <- gls(BAI10yr ~ Band  + BA_sameLR + BA_otherLR , weights=vf, data=dataz, method="ML",...)
  ndcomps2 <- gls(BAI10yr ~ Band  + in5_same + in5_other , weights=vf, data=dataz, method="ML",...)
  ndcompall1 <- gls(BAI10yr ~ Band  + N_Crsc + BA_totLRsc , weights=vf, data=dataz, method="ML",...)
  ndcompall2 <- gls(BAI10yr ~ Band  + N_Crsc + in5_totsc , weights=vf, data=dataz, method="ML",...)
  ndcompall3 <- gls(BAI10yr ~ Band  + N_Crsc + ACFsc , weights=vf, data=dataz, method="ML",...)
  ndcompall4 <- gls(BAI10yr ~ Band  + ACFsc + BA_totLRsc , weights=vf, data=dataz, method="ML",...)
  ndcompall5 <- gls(BAI10yr ~ Band  + ACFsc + in5_totsc , weights=vf, data=dataz, method="ML",...)
  ndcompall6 <- gls(BAI10yr ~ Band  + BA_totLRsc + in5_totsc , weights=vf, data=dataz, method="ML",...)
  ndcompall7 <- gls(BAI10yr ~ Band  + BA_totLRsc + in5_totsc + N_Crsc , weights=vf, data=dataz, method="ML",...)
  ndcompall8 <- gls(BAI10yr ~ Band  + BA_totLRsc + in5_totsc + ACFsc , weights=vf, data=dataz, method="ML",...)
  ndcompall9 <- gls(BAI10yr ~ Band  + BA_totLRsc + N_Crsc + ACFsc , weights=vf, data=dataz, method="ML",...)
  ndcompall10 <- gls(BAI10yr ~ Band  + in5_totsc + N_Crsc + ACFsc , weights=vf, data=dataz, method="ML",...)
  ndcompall11 <- gls(BAI10yr ~ Band  + BA_totLRsc + in5_totsc + N_Crsc + ACFsc , weights=vf, data=dataz, method="ML",...)
  ndcompint1 <- gls(BAI10yr ~ Band * N_Crsc  , weights=vf, data=dataz, method="ML",...)
  ndcompint2 <- gls(BAI10yr ~ Band * BA_totLRsc  , weights=vf, data=dataz, method="ML",...)
  ndcompint3 <- gls(BAI10yr ~ Band * in5_totsc  , weights=vf, data=dataz, method="ML",...)
  ndcompint4 <- gls(BAI10yr ~ Band * ACFsc  , weights=vf, data=dataz, method="ML",...)
  ndcompint5 <- gls(BAI10yr ~ Band * N_Crsc + BA_totLRsc + Band:BA_totLRsc  , weights=vf, data=dataz, method="ML",...)
  ndcompint6 <- gls(BAI10yr ~ Band * N_Crsc + in5_totsc + Band:in5_totsc  , weights=vf, data=dataz, method="ML",...)
  ndcompint7 <- gls(BAI10yr ~ Band * N_Crsc + ACFsc + Band:ACFsc  , weights=vf, data=dataz, method="ML",...)
  ndcompint8 <- gls(BAI10yr ~ Band * BA_totLRsc + in5_totsc + in5_totsc:Band , weights=vf, data=dataz, method="ML",...)
  ndcompint9 <- gls(BAI10yr ~ Band * BA_totLRsc + ACFsc + ACFsc:Band , weights=vf, data=dataz, method="ML",...)
  ndcompint10 <- gls(BAI10yr ~ Band * in5_totsc +  ACFsc + ACFsc:Band  , weights=vf, data=dataz, method="ML",...)
  ## added 07.18.2016
  compint11 <- gls(BAI10yr ~ Band * BA_totLRsc + in5_totsc + Band:in5_totsc + N_Crsc + Band:N_Crsc + DBHsc , weights=vf, data=dataz, method="ML",...)
  compint12 <- gls(BAI10yr ~ Band * BA_totLRsc + in5_totsc + Band:in5_totsc + ACFsc + Band:ACFsc + DBHsc , weights=vf, data=dataz, method="ML",...)
  compint13 <- gls(BAI10yr ~ Band * BA_totLRsc + N_Crsc + Band:N_Crsc + ACFsc + Band:ACFsc + DBHsc , weights=vf, data=dataz, method="ML",...)
  compint14 <- gls(BAI10yr ~ Band * in5_totsc + N_Crsc + Band:N_Crsc + ACFsc + Band:ACFsc + DBHsc , weights=vf, data=dataz, method="ML",...)
  compint15 <- gls(BAI10yr ~ Band * BA_totLRsc + in5_totsc + Band:in5_totsc+ N_Crsc + Band:N_Crsc + ACFsc + Band:ACFsc + DBHsc , weights=vf, data=dataz, method="ML",...)
  ndcompint11 <- gls(BAI10yr ~ Band * BA_totLRsc + in5_totsc + Band:in5_totsc + N_Crsc + Band:N_Crsc , weights=vf, data=dataz, method="ML",...)
  ndcompint12 <- gls(BAI10yr ~ Band * BA_totLRsc + in5_totsc + Band:in5_totsc + ACFsc + Band:ACFsc , weights=vf, data=dataz, method="ML",...)
  ndcompint13 <- gls(BAI10yr ~ Band * BA_totLRsc + N_Crsc + Band:N_Crsc + ACFsc + Band:ACFsc , weights=vf, data=dataz, method="ML",...)
  ndcompint14 <- gls(BAI10yr ~ Band * in5_totsc + N_Crsc + Band:N_Crsc + ACFsc + Band:ACFsc , weights=vf, data=dataz, method="ML",...)
  ndcompint15 <- gls(BAI10yr ~ Band * BA_totLRsc + in5_totsc + Band:in5_totsc+ N_Crsc + Band:N_Crsc + ACFsc + Band:ACFsc , weights=vf, data=dataz, method="ML", control = list(maxIter=100,...),...)
  ## added in 07/29/16 to test single interactions in models with two predictors
  compint5.1 <- gls(BAI10yr ~ Band * N_Crsc + BA_totLRsc + DBHsc , weights = vf, data=dataz, method="ML",...)
  compint6.1 <- gls(BAI10yr ~ Band * N_Crsc + in5_totsc + DBHsc , weights = vf, data=dataz, method="ML",...)
  compint7.1 <- gls(BAI10yr ~ Band * N_Crsc + ACFsc + DBHsc , weights = vf, data=dataz, method="ML",...)
  compint8.1 <- gls(BAI10yr ~ Band * BA_totLRsc + in5_totsc +  DBHsc , weights = vf, data=dataz, method="ML",...)
  compint9.1 <- gls(BAI10yr ~ Band * BA_totLRsc + ACFsc +  DBHsc , weights = vf, data=dataz, method="ML",...)
  compint10.1 <- gls(BAI10yr ~ Band * in5_totsc +  ACFsc + DBHsc , weights = vf, data=dataz, method="ML",...)
  compint5.2 <- gls(BAI10yr ~ Band + N_Crsc + BA_totLRsc + Band:BA_totLRsc + DBHsc , weights = vf, data=dataz, method="ML",...)
  compint6.2 <- gls(BAI10yr ~ Band + N_Crsc + in5_totsc + Band:in5_totsc + DBHsc , weights = vf, data=dataz, method="ML",...)
  compint7.2 <- gls(BAI10yr ~ Band + N_Crsc + ACFsc + Band:ACFsc + DBHsc , weights = vf, data=dataz, method="ML",...)
  compint8.2 <- gls(BAI10yr ~ Band + BA_totLRsc + in5_totsc + Band:in5_totsc +  DBHsc , weights = vf, data=dataz, method="ML",...)
  compint9.2 <- gls(BAI10yr ~ Band + BA_totLRsc + ACFsc + ACFsc:Band +  DBHsc , weights = vf, data=dataz, method="ML",...)
  compint10.2 <- gls(BAI10yr ~ Band + in5_totsc +  ACFsc + ACFsc:Band + DBHsc , weights = vf, data=dataz, method="ML",...)
  ndcompint5.1 <- gls(BAI10yr ~ Band * N_Crsc + BA_totLRsc  , weights = vf, data=dataz, method="ML",...)
  ndcompint6.1 <- gls(BAI10yr ~ Band * N_Crsc + in5_totsc  , weights = vf, data=dataz, method="ML",...)
  ndcompint7.1 <- gls(BAI10yr ~ Band * N_Crsc + ACFsc  , weights = vf, data=dataz, method="ML",...)
  ndcompint8.1 <- gls(BAI10yr ~ Band * BA_totLRsc + in5_totsc  , weights = vf, data=dataz, method="ML",...)
  ndcompint9.1 <- gls(BAI10yr ~ Band * BA_totLRsc + ACFsc  , weights = vf, data=dataz, method="ML",...)
  ndcompint10.1 <- gls(BAI10yr ~ Band * in5_totsc +  ACFsc  , weights = vf, data=dataz, method="ML",...)
  ndcompint5.2 <- gls(BAI10yr ~ Band + N_Crsc + BA_totLRsc + Band:BA_totLRsc  , weights = vf, data=dataz, method="ML",...)
  ndcompint6.2 <- gls(BAI10yr ~ Band + N_Crsc + in5_totsc + Band:in5_totsc  , weights = vf, data=dataz, method="ML",...)
  ndcompint7.2 <- gls(BAI10yr ~ Band + N_Crsc + ACFsc + Band:ACFsc  , weights = vf, data=dataz, method="ML",...)
  ndcompint8.2 <- gls(BAI10yr ~ Band + BA_totLRsc + in5_totsc + in5_totsc:Band  , weights = vf, data=dataz, method="ML",...)
  ndcompint9.2 <- gls(BAI10yr ~ Band + BA_totLRsc + ACFsc + ACFsc:Band  , weights = vf, data=dataz, method="ML",...)
  ndcompint10.2 <- gls(BAI10yr ~ Band + in5_totsc +  ACFsc + ACFsc:Band  , weights = vf, data=dataz, method="ML",...)
  
  
  aics <- AIC(null0, null1, comp1, comp2,comp3,comp4,comps1,comps2,compall1, compall2,compall3,compall4,compall5,compall6,compall7,compall8,compall9,compall10,compall11,compint1, compint2, compint3, compint4, compint5, compint6, compint7, compint8, compint9, compint10, compint11, compint12, compint13, compint14, compint15
              ,ndcomp1, ndcomp2,ndcomp3,ndcomp4,ndcomps1,ndcomps2,ndcompall1, ndcompall2,ndcompall3,ndcompall4,ndcompall5,ndcompall6,ndcompall7,ndcompall8,ndcompall9,ndcompall10,ndcompall11,ndcompint1, ndcompint2, ndcompint3, ndcompint4, ndcompint5, ndcompint6, ndcompint7, ndcompint8, ndcompint9, ndcompint10, ndcompint11, ndcompint12, ndcompint13, ndcompint14, ndcompint15
              ,compint5.1, compint6.1, compint7.1, compint8.1, compint9.1, compint10.1, compint5.2, compint6.2, compint7.2, compint8.2, compint9.2, compint10.2
              ,ndcompint5.1, ndcompint6.1, ndcompint7.1, ndcompint8.1, ndcompint9.1, ndcompint10.1, ndcompint5.2, ndcompint6.2, ndcompint7.2, ndcompint8.2, ndcompint9.2, ndcompint10.2)
  
  bics <- BIC(null0, null1, comp1, comp2,comp3,comp4,comps1,comps2,compall1, compall2,compall3,compall4,compall5,compall6,compall7,compall8,compall9,compall10,compall11,compint1, compint2, compint3, compint4, compint5, compint6, compint7, compint8, compint9, compint10, compint11, compint12, compint13, compint14, compint15
              ,ndcomp1, ndcomp2,ndcomp3,ndcomp4,ndcomps1,ndcomps2,ndcompall1, ndcompall2,ndcompall3,ndcompall4,ndcompall5,ndcompall6,ndcompall7,ndcompall8,ndcompall9,ndcompall10,ndcompall11,ndcompint1, ndcompint2, ndcompint3, ndcompint4, ndcompint5, ndcompint6, ndcompint7, ndcompint8, ndcompint9, ndcompint10, ndcompint11, ndcompint12, ndcompint13, ndcompint14, ndcompint15
              ,compint5.1, compint6.1, compint7.1, compint8.1, compint9.1, compint10.1, compint5.2, compint6.2, compint7.2, compint8.2, compint9.2, compint10.2
              ,ndcompint5.1, ndcompint6.1, ndcompint7.1, ndcompint8.1, ndcompint9.1, ndcompint10.1, ndcompint5.2, ndcompint6.2, ndcompint7.2, ndcompint8.2, ndcompint9.2, ndcompint10.2)
  
  
  
  aics$deltaAIC <- aics$AIC - min(aics$AIC)
  aics$BIC <- bics$BIC
  aics$deltaBIC <- aics$BIC - min(aics$BIC)
  print(aics[which(aics$deltaAIC <= aic.cutoff),])
  aics2 <- aics[which(aics$deltaAIC<2),]
  aics2$mod <- rownames(aics2)
  aicsbest <- arrange(aics2, df, deltaAIC )
  print(summary(get(aicsbest$mod[1])))
  
  return(get(aicsbest$mod[1]))
}

# funciton sans the double interactions that sometimes don't converge with varPower() vfs
compselection.gls.simple <- function(dataz, vf = varPower(form=~DBH), elevlevels = list(L="L", M="M", H="H"), aic.cutoff=4){
  levels(dataz$Band) <- elevlevels
  dataz$DBHsc <- scale(dataz$DBH) # scale predictors to avoid numerical instability
  dataz$BA_totLRsc <- scale(dataz$BA_totLR) # and to get BAI at species mean predictor level (intercept)
  dataz$N_Crsc <- scale(dataz$N_Cr)
  dataz$in5_totsc <- scale(dataz$in5_tot)
  dataz$ACFsc <- scale(dataz$ACF)
  dataz$BAI10yrsc <- scale(dataz$BAI10yr)
  na.flags <- apply(dataz[,c("DBH","BA_totLR", "N_Cr", "in5_tot", "ACF", "BAI10yr")], MARGIN=2, FUN= function(x){length(which(is.na(x)))})
  if(sum(na.flags>0)){print(paste("NAs in:", names(na.flags)[which(na.flags>0)]))}
  null0 <- gls(BAI10yr ~ Band , weights = vf, data=dataz, method = "ML")
  null1 <- gls(BAI10yr ~ Band + DBHsc , weights = vf, data=dataz, method="ML")
  comp1 <- gls(BAI10yr ~ Band + DBHsc + N_Crsc , weights = vf, data=dataz, method="ML")
  comp2 <- gls(BAI10yr ~ Band + DBHsc + BA_totLRsc , weights = vf, data=dataz, method="ML")
  comp3 <- gls(BAI10yr ~ Band + DBHsc + in5_totsc , weights = vf, data=dataz, method="ML")
  comp4 <- gls(BAI10yr ~ Band + DBHsc + ACFsc , weights = vf, data=dataz, method="ML")
  comps1 <- gls(BAI10yr ~ Band + DBHsc + BA_sameLR + BA_otherLR , weights = vf, data=dataz, method="ML")
  comps2 <- gls(BAI10yr ~ Band + DBHsc + in5_same + in5_other , weights = vf, data=dataz, method="ML")
  compall1 <- gls(BAI10yr ~ Band + DBHsc + N_Crsc + BA_totLRsc , weights = vf, data=dataz, method="ML")
  compall2 <- gls(BAI10yr ~ Band + DBHsc + N_Crsc + in5_totsc , weights = vf, data=dataz, method="ML")
  compall3 <- gls(BAI10yr ~ Band + DBHsc + N_Crsc + ACFsc , weights = vf, data=dataz, method="ML")
  compall4 <- gls(BAI10yr ~ Band + DBHsc + ACFsc + BA_totLRsc , weights = vf, data=dataz, method="ML")
  compall5 <- gls(BAI10yr ~ Band + DBHsc + ACFsc + in5_totsc , weights = vf, data=dataz, method="ML")
  compall6 <- gls(BAI10yr ~ Band + DBHsc + BA_totLRsc + in5_totsc , weights = vf, data=dataz, method="ML")
  compall7 <- gls(BAI10yr ~ Band + DBHsc + BA_totLRsc + in5_totsc + N_Crsc , weights = vf, data=dataz, method="ML")
  compall8 <- gls(BAI10yr ~ Band + DBHsc + BA_totLRsc + in5_totsc + ACFsc , weights = vf, data=dataz, method="ML")
  compall9 <- gls(BAI10yr ~ Band + DBHsc + BA_totLRsc + N_Crsc + ACFsc , weights = vf, data=dataz, method="ML")
  compall10 <- gls(BAI10yr ~ Band + DBHsc + in5_totsc + N_Crsc + ACFsc , weights = vf, data=dataz, method="ML")
  compall11 <- gls(BAI10yr ~ Band + DBHsc + BA_totLRsc + in5_totsc + N_Crsc + ACFsc , weights = vf, data=dataz, method="ML")
  compint1 <- gls(BAI10yr ~ Band * N_Crsc + DBHsc , weights = vf, data=dataz, method="ML")
  compint2 <- gls(BAI10yr ~ Band * BA_totLRsc + DBHsc , weights = vf, data=dataz, method="ML")
  compint3 <- gls(BAI10yr ~ Band * in5_totsc +  DBHsc , weights = vf, data=dataz, method="ML")
  compint4 <- gls(BAI10yr ~ Band * ACFsc + DBHsc , weights = vf, data=dataz, method="ML")
  compint5 <- gls(BAI10yr ~ Band * N_Crsc + BA_totLRsc + Band:BA_totLRsc + DBHsc , weights = vf, data=dataz, method="ML")
  compint6 <- gls(BAI10yr ~ Band * N_Crsc + in5_totsc + Band:in5_totsc + DBHsc , weights = vf, data=dataz, method="ML")
  compint7 <- gls(BAI10yr ~ Band * N_Crsc + ACFsc + Band:ACFsc + DBHsc , weights = vf, data=dataz, method="ML")
  compint8 <- gls(BAI10yr ~ Band * BA_totLRsc + in5_totsc + in5_totsc:Band +  DBHsc , weights = vf, data=dataz, method="ML")
  compint9 <- gls(BAI10yr ~ Band * BA_totLRsc + ACFsc + ACFsc:Band +  DBHsc , weights = vf, data=dataz, method="ML")
  compint10 <- gls(BAI10yr ~ Band * in5_totsc +  ACFsc + ACFsc:Band + DBHsc , weights = vf, data=dataz, method="ML")
  
  ndcomp1 <- gls(BAI10yr ~ Band  + N_Crsc , weights=vf, data=dataz, method="ML")
  ndcomp2 <- gls(BAI10yr ~ Band  + BA_totLRsc , weights=vf, data=dataz, method="ML")
  ndcomp3 <- gls(BAI10yr ~ Band  + in5_totsc , weights=vf, data=dataz, method="ML")
  ndcomp4 <- gls(BAI10yr ~ Band  + ACFsc , weights=vf, data=dataz, method="ML")
  ndcomps1 <- gls(BAI10yr ~ Band  + BA_sameLR + BA_otherLR , weights=vf, data=dataz, method="ML")
  ndcomps2 <- gls(BAI10yr ~ Band  + in5_same + in5_other , weights=vf, data=dataz, method="ML")
  ndcompall1 <- gls(BAI10yr ~ Band  + N_Crsc + BA_totLRsc , weights=vf, data=dataz, method="ML")
  ndcompall2 <- gls(BAI10yr ~ Band  + N_Crsc + in5_totsc , weights=vf, data=dataz, method="ML")
  ndcompall3 <- gls(BAI10yr ~ Band  + N_Crsc + ACFsc , weights=vf, data=dataz, method="ML")
  ndcompall4 <- gls(BAI10yr ~ Band  + ACFsc + BA_totLRsc , weights=vf, data=dataz, method="ML")
  ndcompall5 <- gls(BAI10yr ~ Band  + ACFsc + in5_totsc , weights=vf, data=dataz, method="ML")
  ndcompall6 <- gls(BAI10yr ~ Band  + BA_totLRsc + in5_totsc , weights=vf, data=dataz, method="ML")
  ndcompall7 <- gls(BAI10yr ~ Band  + BA_totLRsc + in5_totsc + N_Crsc , weights=vf, data=dataz, method="ML")
  ndcompall8 <- gls(BAI10yr ~ Band  + BA_totLRsc + in5_totsc + ACFsc , weights=vf, data=dataz, method="ML")
  ndcompall9 <- gls(BAI10yr ~ Band  + BA_totLRsc + N_Crsc + ACFsc , weights=vf, data=dataz, method="ML")
  ndcompall10 <- gls(BAI10yr ~ Band  + in5_totsc + N_Crsc + ACFsc , weights=vf, data=dataz, method="ML")
  ndcompall11 <- gls(BAI10yr ~ Band  + BA_totLRsc + in5_totsc + N_Crsc + ACFsc , weights=vf, data=dataz, method="ML")
  ndcompint1 <- gls(BAI10yr ~ Band * N_Crsc  , weights=vf, data=dataz, method="ML")
  ndcompint2 <- gls(BAI10yr ~ Band * BA_totLRsc  , weights=vf, data=dataz, method="ML")
  ndcompint3 <- gls(BAI10yr ~ Band * in5_totsc  , weights=vf, data=dataz, method="ML")
  ndcompint4 <- gls(BAI10yr ~ Band * ACFsc  , weights=vf, data=dataz, method="ML")
  ndcompint5 <- gls(BAI10yr ~ Band * N_Crsc + BA_totLRsc + Band:BA_totLRsc  , weights=vf, data=dataz, method="ML")
  ndcompint6 <- gls(BAI10yr ~ Band * N_Crsc + in5_totsc + Band:in5_totsc  , weights=vf, data=dataz, method="ML")
  ndcompint7 <- gls(BAI10yr ~ Band * N_Crsc + ACFsc + Band:ACFsc  , weights=vf, data=dataz, method="ML")
  ndcompint8 <- gls(BAI10yr ~ Band * BA_totLRsc + in5_totsc + in5_totsc:Band , weights=vf, data=dataz, method="ML")
  ndcompint9 <- gls(BAI10yr ~ Band * BA_totLRsc + ACFsc + ACFsc:Band , weights=vf, data=dataz, method="ML")
  ndcompint10 <- gls(BAI10yr ~ Band * in5_totsc +  ACFsc + ACFsc:Band  , weights=vf, data=dataz, method="ML")
  ## added 07.18.2016
  # compint11 <- gls(BAI10yr ~ Band * BA_totLRsc + in5_totsc + Band:in5_totsc + N_Crsc + Band:N_Crsc + DBHsc , weights=vf, data=dataz, method="ML")
  # compint12 <- gls(BAI10yr ~ Band * BA_totLRsc + in5_totsc + Band:in5_totsc + ACFsc + Band:ACFsc + DBHsc , weights=vf, data=dataz, method="ML")
  # compint13 <- gls(BAI10yr ~ Band * BA_totLRsc + N_Crsc + Band:N_Crsc + ACFsc + Band:ACFsc + DBHsc , weights=vf, data=dataz, method="ML")
  # compint14 <- gls(BAI10yr ~ Band * in5_totsc + N_Crsc + Band:N_Crsc + ACFsc + Band:ACFsc + DBHsc , weights=vf, data=dataz, method="ML")
  # compint15 <- gls(BAI10yr ~ Band * BA_totLRsc + in5_totsc + Band:in5_totsc+ N_Crsc + Band:N_Crsc + ACFsc + Band:ACFsc + DBHsc , weights=vf, data=dataz, method="ML")
  # ndcompint11 <- gls(BAI10yr ~ Band * BA_totLRsc + in5_totsc + Band:in5_totsc + N_Crsc + Band:N_Crsc , weights=vf, data=dataz, method="ML")
  # ndcompint12 <- gls(BAI10yr ~ Band * BA_totLRsc + in5_totsc + Band:in5_totsc + ACFsc + Band:ACFsc , weights=vf, data=dataz, method="ML")
  # ndcompint13 <- gls(BAI10yr ~ Band * BA_totLRsc + N_Crsc + Band:N_Crsc + ACFsc + Band:ACFsc , weights=vf, data=dataz, method="ML")
  # ndcompint14 <- gls(BAI10yr ~ Band * in5_totsc + N_Crsc + Band:N_Crsc + ACFsc + Band:ACFsc , weights=vf, data=dataz, method="ML")
  # ndcompint15 <- gls(BAI10yr ~ Band * BA_totLRsc + in5_totsc + Band:in5_totsc+ N_Crsc + Band:N_Crsc + ACFsc + Band:ACFsc , weights=vf, data=dataz, method="ML", control = list(maxIter=100))
  ## added in 07/29/16 to test single interactions in models with two predictors
  # compint5.1 <- gls(BAI10yr ~ Band * N_Crsc + BA_totLRsc + DBHsc , weights = vf, data=dataz, method="ML")
  # compint6.1 <- gls(BAI10yr ~ Band * N_Crsc + in5_totsc + DBHsc , weights = vf, data=dataz, method="ML")
  # compint7.1 <- gls(BAI10yr ~ Band * N_Crsc + ACFsc + DBHsc , weights = vf, data=dataz, method="ML")
  # compint8.1 <- gls(BAI10yr ~ Band * BA_totLRsc + in5_totsc +  DBHsc , weights = vf, data=dataz, method="ML")
  # compint9.1 <- gls(BAI10yr ~ Band * BA_totLRsc + ACFsc +  DBHsc , weights = vf, data=dataz, method="ML")
  # compint10.1 <- gls(BAI10yr ~ Band * in5_totsc +  ACFsc + DBHsc , weights = vf, data=dataz, method="ML")
  # compint5.2 <- gls(BAI10yr ~ Band + N_Crsc + BA_totLRsc + Band:BA_totLRsc + DBHsc , weights = vf, data=dataz, method="ML")
  # compint6.2 <- gls(BAI10yr ~ Band + N_Crsc + in5_totsc + Band:in5_totsc + DBHsc , weights = vf, data=dataz, method="ML")
  # compint7.2 <- gls(BAI10yr ~ Band + N_Crsc + ACFsc + Band:ACFsc + DBHsc , weights = vf, data=dataz, method="ML")
  # compint8.2 <- gls(BAI10yr ~ Band + BA_totLRsc + in5_totsc + Band:in5_totsc +  DBHsc , weights = vf, data=dataz, method="ML")
  # compint9.2 <- gls(BAI10yr ~ Band + BA_totLRsc + ACFsc + ACFsc:Band +  DBHsc , weights = vf, data=dataz, method="ML")
  # compint10.2 <- gls(BAI10yr ~ Band + in5_totsc +  ACFsc + ACFsc:Band + DBHsc , weights = vf, data=dataz, method="ML")
  # ndcompint5.1 <- gls(BAI10yr ~ Band * N_Crsc + BA_totLRsc  , weights = vf, data=dataz, method="ML")
  # ndcompint6.1 <- gls(BAI10yr ~ Band * N_Crsc + in5_totsc  , weights = vf, data=dataz, method="ML")
  # ndcompint7.1 <- gls(BAI10yr ~ Band * N_Crsc + ACFsc  , weights = vf, data=dataz, method="ML")
  # ndcompint8.1 <- gls(BAI10yr ~ Band * BA_totLRsc + in5_totsc  , weights = vf, data=dataz, method="ML")
  # ndcompint9.1 <- gls(BAI10yr ~ Band * BA_totLRsc + ACFsc  , weights = vf, data=dataz, method="ML")
  # ndcompint10.1 <- gls(BAI10yr ~ Band * in5_totsc +  ACFsc  , weights = vf, data=dataz, method="ML")
  # ndcompint5.2 <- gls(BAI10yr ~ Band + N_Crsc + BA_totLRsc + Band:BA_totLRsc  , weights = vf, data=dataz, method="ML")
  # ndcompint6.2 <- gls(BAI10yr ~ Band + N_Crsc + in5_totsc + Band:in5_totsc  , weights = vf, data=dataz, method="ML")
  # ndcompint7.2 <- gls(BAI10yr ~ Band + N_Crsc + ACFsc + Band:ACFsc  , weights = vf, data=dataz, method="ML")
  # ndcompint8.2 <- gls(BAI10yr ~ Band + BA_totLRsc + in5_totsc + in5_totsc:Band  , weights = vf, data=dataz, method="ML")
  # ndcompint9.2 <- gls(BAI10yr ~ Band + BA_totLRsc + ACFsc + ACFsc:Band  , weights = vf, data=dataz, method="ML")
  # ndcompint10.2 <- gls(BAI10yr ~ Band + in5_totsc +  ACFsc + ACFsc:Band  , weights = vf, data=dataz, method="ML")
  # 
  
  aics <- AIC(null0, null1, comp1, comp2,comp3,comp4,comps1,comps2,compall1, compall2,compall3,compall4,compall5,compall6,compall7,compall8,compall9,compall10,compall11,compint1, compint2, compint3, compint4 #, compint5, compint6, compint7, compint8, compint9, compint10, compint11, compint12, compint13, compint14, compint15
              ,ndcomp1, ndcomp2,ndcomp3,ndcomp4,ndcomps1,ndcomps2,ndcompall1, ndcompall2,ndcompall3,ndcompall4,ndcompall5,ndcompall6,ndcompall7,ndcompall8,ndcompall9,ndcompall10,ndcompall11,ndcompint1, ndcompint2, ndcompint3, ndcompint4 #, ndcompint5, ndcompint6, ndcompint7, ndcompint8, ndcompint9, ndcompint10, ndcompint11, ndcompint12, ndcompint13, ndcompint14, ndcompint15
  )#,compint5.1, compint6.1, compint7.1, compint8.1, compint9.1, compint10.1, compint5.2, compint6.2, compint7.2, compint8.2, compint9.2, compint10.2
  #,ndcompint5.1, ndcompint6.1, ndcompint7.1, ndcompint8.1, ndcompint9.1, ndcompint10.1, ndcompint5.2, ndcompint6.2, ndcompint7.2, ndcompint8.2, ndcompint9.2, ndcompint10.2)
  
  bics <- BIC(null0, null1, comp1, comp2,comp3,comp4,comps1,comps2,compall1, compall2,compall3,compall4,compall5,compall6,compall7,compall8,compall9,compall10,compall11,compint1, compint2, compint3, compint4 #, compint5, compint6, compint7, compint8, compint9, compint10, compint11, compint12, compint13, compint14, compint15
              ,ndcomp1, ndcomp2,ndcomp3,ndcomp4,ndcomps1,ndcomps2,ndcompall1, ndcompall2,ndcompall3,ndcompall4,ndcompall5,ndcompall6,ndcompall7,ndcompall8,ndcompall9,ndcompall10,ndcompall11,ndcompint1, ndcompint2, ndcompint3, ndcompint4 #, ndcompint5, ndcompint6, ndcompint7, ndcompint8, ndcompint9, ndcompint10, ndcompint11, ndcompint12, ndcompint13, ndcompint14, ndcompint15
  )#,compint5.1, compint6.1, compint7.1, compint8.1, compint9.1, compint10.1, compint5.2, compint6.2, compint7.2, compint8.2, compint9.2, compint10.2
  #,ndcompint5.1, ndcompint6.1, ndcompint7.1, ndcompint8.1, ndcompint9.1, ndcompint10.1, ndcompint5.2, ndcompint6.2, ndcompint7.2, ndcompint8.2, ndcompint9.2, ndcompint10.2)
  
  
  
  aics$deltaAIC <- aics$AIC - min(aics$AIC)
  aics$BIC <- bics$BIC
  aics$deltaBIC <- aics$BIC - min(aics$BIC)
  print(aics[which(aics$deltaAIC <= aic.cutoff),])
  aics2 <- aics[which(aics$deltaAIC<2),]
  aics2$mod <- rownames(aics2)
  aicsbest <- arrange(aics2, df, deltaAIC )
  print(summary(get(aicsbest$mod[1])))
  
  return(get(aicsbest$mod[1]))
}



# useful model diagnostics plots. Not all that are suggested by icebreakR (https://cran.r-project.org/doc/contrib/Robinson-icebreaker.pdf), but a reasonable sample
residplots <- function(model, dataz, res.type="normalized", Title = "", mixed = F){
  quartz(width=7, height=6, title=paste(dataz$State[1], dataz$Species[1], Title))
  par(mfrow=c(3,2), mar=c(4,4,1,1), oma=c(0,0,2,0))
  if(mixed==T){
    ref.group <- ranef(model)[[1]]
    ref.var.group <- tapply(residuals(model, type="pearson", level=1),
                            dataz$PairUn, var)
    qqnorm(ref.group, main="Q-Q Normal - group Random Effects")
    qqline(ref.group, col="red")
  }
  else{
    qqnorm(residuals(model), main = "raw"); qqline(residuals(model))
  }
  qqnorm(residuals(model, type=res.type), main="normalized"); qqline(residuals(model, type=res.type))
  mod.resids <- residuals(model, type=res.type)
  scatter.smooth(mod.resids~dataz$DBH); abline(h=0, col="red")
  scatter.smooth(mod.resids~dataz$in5_tot); abline(h=0, col="red")
  scatter.smooth(mod.resids~dataz$BAI10yr); abline(h=0, col="red")
  boxplot(mod.resids~dataz$Band, notch=T, varwidth=T); abline(h=0)
  
}





#____________________________
### EXAMPLE: iterative model selection technique to find best random, fixed, and variance structure #########
#____________________________
# Strategy: 
# - test whether random effect for tree pair is needed
# - select the best variance structure, searching across all fixed effects structures for each variance structure
# - validate model (model criticism plots, testing for non-linear responses, etc.)
# - refit best model with REML and draw statistical inferences from t-tests on parameters and LRT on Elevation X Competition interaction


### Example for Colorado Abies lasiocarpa (data in CO_ABLAinfo)

## read in Growth and Competitive Environment data from Abies lasiocarpa from Colorado (SI Data S1)
CO_ABLAinfo <- read.csv("SIData_S1_CO_ABLA_GrowthCompetition_20180428.csv", header=T)[,-c(1:2)]

## remove na's from dataset for feeding to lme
CO_ABLAinfo.na <- CO_ABLAinfo[!is.na(CO_ABLAinfo$BAI10yr),]

# test whether pair random effect is needed with 'over the top' fixed model
prelimCO_ABLA1 <- gls(BAI10yr ~ Band * BA_totLRsc + ACFsc * Band + DBHsc, data=CO_ABLAinfo.na)
prelimCO_ABLA2 <- lme(BAI10yr ~ Band * BA_totLRsc + ACFsc * Band + DBHsc, random= ~1|PairUn, data=CO_ABLAinfo.na)
AIC(prelimCO_ABLA1, prelimCO_ABLA2); anova(prelimCO_ABLA1, prelimCO_ABLA2) # no need for PairUn


### use compselection.gls() function 
# to find the best fixed effects structure for each possible variance
# structure, and then compare the AIC of the best models
# note: this proceedure is needed because the optimal fixed structure sometimes
#       changes depending on the variance structure applied. This searches all combos of fixed and variance structure space
CO_ABLAbaigls0 <- compselection.gls(CO_ABLAinfo.na, vf= NULL, aic.cutoff = 4)
CO_ABLAbaigls1 <- compselection.gls(CO_ABLAinfo.na, vf= varPower(form=~DBH), aic.cutoff = 4)
CO_ABLAbaigls2 <- compselection.gls(CO_ABLAinfo.na, vf= varExp(form=~DBH), aic.cutoff = 4)
CO_ABLAbaigls3 <- compselection.gls(CO_ABLAinfo.na, vf= varFixed(~DBH), aic.cutoff = 4)
#CO_ABLAbaigls4 <- compselection.gls(CO_ABLAinfo.na, vf= varPower(), aic.cutoff = 4) # didn't converge...
#CO_ABLAbaigls4 <- compselection.gls.simple(CO_ABLAinfo.na, vf= varPower(), aic.cutoff = 4) # didn't converge even without complex interactions
CO_ABLAbaigls5 <- compselection.gls(CO_ABLAinfo.na, vf= varPower(form=~DBH|Band))
## need to fit gls4 with update because some models don't converge even with compselection.gls.simple
dataz <- CO_ABLAinfo.na
CO_ABLAbaigls4 <- update(CO_ABLAbaigls5, weights=varPower()) # still better even if using suboptimal model (i.e gls3 with no int rather than gls5 with int...)

AIC(CO_ABLAbaigls0, CO_ABLAbaigls1, CO_ABLAbaigls2, CO_ABLAbaigls3, CO_ABLAbaigls4, CO_ABLAbaigls5) 
# varPower() is the best variance structure

# model criticism plots
residplots(CO_ABLAbaigls4, CO_ABLAinfo.na, Title="BAI")
plot(fitted(CO_ABLAbaigls4, level=0)~BAI10yr, CO_ABLAinfo.na, col=Band)
abline(0,1)
# doesn't violate assumptions, but could explain more variance
# so we'll test to make sure whe're not missing anything
dataz <- CO_ABLAinfo.na
vf <- varPower()
ablatest1 <- update(CO_ABLAbaigls4, .~. + I(DBHsc^2))
ablatest2 <- update(CO_ABLAbaigls4, .~. + I(DBHsc^2) + I(BA_totLRsc ^2))
ablatest3 <- update(CO_ABLAbaigls4, .~. + DBHsc : Band)
plot(BAI10yr~BA_totLRsc, CO_ABLAinfo.na, col=Band)
points(fitted(CO_ABLAbaigls4)~BA_totLRsc, CO_ABLAinfo.na, col=Band, pch=17)
points(fitted(ablatest1)~DBHsc, CO_ABLAinfo.na, col=Band, pch=16, cex=.8)
AIC(ablatest1, ablatest2, ablatest3, CO_ABLAbaigls4) # nothing improves the model

#### Best BAI model 
levels(CO_ABLAinfo.na$Band) <- list(M="M",L="L",H="H") # set mid elevation range center as intercept for statistical inference
CO_ABLAbaibest <- gls(BAI10yr ~ Band * BA_totLRsc + DBHsc, weights= varPower(), CO_ABLAinfo.na, method="REML")
### use LRT test to assess the significance of the Band * Competition interaction
CO_ABLAbaibest.0 <- gls(BAI10yr ~ Band + BA_totLRsc + DBHsc, weights= varPower(), CO_ABLAinfo.na, method="ML")
CO_ABLAbaibest.1 <- gls(BAI10yr ~ Band * BA_totLRsc + DBHsc, weights= varPower(), CO_ABLAinfo.na, method="ML")
anova(CO_ABLAbaibest.0, CO_ABLAbaibest.1) # p=0.0146
# calculate BAI change from a 1sd increase in competition
CO_ABLApreddata <-  data.frame(Band=c("L","M","H"), DBHsc=c(0,0,0), BA_totLRsc = c(1,1,1))

summary(CO_ABLAbaibest)
# Coefficients reported in Table S3
#                      Value Std.Error   t-value p-value
# (Intercept)      1475.5527 120.21584 12.274196  0.0000
# BandL              79.2796 170.58202  0.464760  0.6434
# BandH            -385.1298 150.55303 -2.558101  0.0124
# BA_totLRsc       -289.5943  91.87255 -3.152130  0.0023
# DBHsc             371.3288  61.57888  6.030132  0.0000
# BandL:BA_totLRsc  -69.4021 133.57564 -0.519572  0.6048
# BandH:BA_totLRsc  227.1680 120.57180  1.884089  0.0631








#___________________________________________________________________________
####### Absolute best Growth and Competition models ################
#___________________________________________________________________________
# (for statistical inferences reported in Table S2, Table S3)
  # these are the best models selected using the iterative model selection technique (first selecting random effects, then variance structures, than fixed effects)
  # suggested in Zuur et al. 2007 and demonstrated above

##CO
CO_PIPOcompmod <-gls(BAI10yr ~ Band * N_Crsc + ACFsc + DBHsc , CO_PIPOinfo, weights=varPower(), method="REML")#update(CO_PIPObai0, method="REML")
CO_PIPOpreddata <- data.frame(Band = c("L","M","H"), N_Crsc=c(1,1,1),ACFsc=c(-1,-1,-1), DBHsc=c(0,0,0))
CO_PIPOpreds <- predictSE(CO_PIPOcompmod, newdata=CO_PIPOpreddata, se.fit=T)

CO_POTRcompmod <- lme(BAI10yr ~ Band + DBHsc + BA_totLRsc + ACFsc, random= ~1|PairUn, CO_POTRinfo, weights=varExp(form=~DBH), method="REML")
CO_POTRpreddata <- data.frame(Band=c("L","M","H"), DBHsc=c(0,0,0), BA_totLRsc = c(1,1,1), in5_totsc = c(1,1,1), ACFsc = c(-1,-1,-1))
CO_POTRpreds <- predictSE(CO_POTRcompmod, newdata = CO_POTRpreddata, level=0, se.fit=T )


CO_ABLAinfo.na <- CO_ABLAinfo[!is.na(CO_ABLAinfo$BAI10yr),] #remove NAs
CO_ABLAcompmod <- gls(BAI10yr ~ Band * BA_totLRsc + DBHsc, weights= varPower(), CO_ABLAinfo.na, method="REML")
CO_ABLApreddata <-  data.frame(Band=c("L","M","H"), DBHsc=c(0,0,0), BA_totLRsc = c(1,1,1))
CO_ABLApreds <- predictSE(CO_ABLAcompmod, newdata = CO_ABLApreddata, se.fit=T)

## MT
MT_TSHEcompmod <- gls(BAI10yr ~ Band * BA_totLRsc + DBHsc, weights=varPower(form=~DBHsc|Band), MT_TSHEinfo, method="REML")
MT_TSHEpreddata <-  data.frame(Band=c("L","M","H"), DBHsc=c(0,0,0), BA_totLRsc = c(1,1,1))
MT_TSHEpreds <- predictSE(MT_TSHEcompmod, newdata = MT_TSHEpreddata, se.fit=T)

MT_PSMEinfo.cl <- MT_PSMEinfo[-76,] # remove one outlier tree 
MT_PSMEcompmod <- gls(BAI10yr~ Band * ACFsc, MT_PSMEinfo.cl)
MT_PSMEpreddata <- data.frame(Band=c("L","M","H"), ACFsc=c(-1,-1,-1))
MT_PSMEpreds <- predictSE(MT_PSMEcompmod, newdata=MT_PSMEpreddata, se.fit=T) # this is identical if you use lm() and predict()

MT_ABLAcompmod  <- lme(BAI10yr ~ Band + DBHsc + ACFsc + Band:DBHsc , random= ~1|PairUn, data=MT_ABLAinfo, weights=varPower(form=~DBH), method="REML")
MT_ABLApreddata <- data.frame(Band=c("L","M","H"), ACFsc=c(-1,-1,-1), DBHsc = c(0,0,0))
MT_ABLApreds <- predictSE(MT_ABLAcompmod, newdata = MT_ABLApreddata, level=0, se.fit=T )

### WA
WA_TSHEinfo.cl <- WA_TSHEinfo[-c(20,26,43),] # remove three outlier trees
WA_TSHEcompmod <- gls(BAI10yr ~ Band + DBHsc + ACFsc,weights=varExp(form=~DBH), WA_TSHEinfo.cl, method="REML")
WA_TSHEpreddata <- data.frame(Band=c("M","H"), ACFsc=c(-1,-1), DBHsc = c(0,0))
WA_TSHEpreds <- predictSE(WA_TSHEcompmod, newdata = WA_TSHEpreddata, se.fit=T )

WA_PSMEinfo.cl <- WA_PSMEinfo[-1,] # remove one oultier tree
WA_PSMEcompmod <- lme(BAI10yr~ Band + DBHsc + BA_totLRsc, random= ~ 1|PairUn, weights=varPower(form=~DBH), WA_PSMEinfo.cl, method="REML")
WA_PSMEpreddata <- data.frame(Band=c("L","M","H"), BA_totLRsc=c(1,1,1), DBHsc = c(0,0,0))
WA_PSMEpreds <- predictSE(WA_PSMEcompmod, newdata = WA_PSMEpreddata, level=0, se.fit=T )

WA_ABLAcompmod <- gls(BAI10yr ~ Band * BA_totLRsc + DBHsc, weights=varPower(), WA_ABLAinfo, control=list(maxIter=100, msMaxIter=100), method="REML")
WA_ABLApreddata  <- data.frame(Band = c("L","M","H"), BA_totLRsc = c(1,1,1), DBHsc=c(0,0,0))
WA_ABLApreds <- predictSE(WA_ABLAcompmod, newdata=WA_ABLApreddata, se.fit=T)


#___________________________________________________________________________
##### full model for best parameter estimates ##############################
#___________________________________________________________________________
# (used for Figures 2,4 and S2)
  # note: these full models include Band and Band*Competition terms even where they weren't siginificant
  # this provides the best estimates of mean growth and competitive effects for all species, all elevation bands
  # even when there was no statistically significant difference between elevation bands. 

##CO
CO_PIPOcompmod.int <-gls(BAI10yr ~ Band * N_Crsc + ACFsc + DBHsc , CO_PIPOinfo, weights=varPower(), method="REML")#update(CO_PIPObai0, method="REML")
CO_PIPOpreddata <- data.frame(Band = c("L","M","H"), N_Crsc=c(1,1,1),ACFsc=c(-1,-1,-1), DBHsc=c(0,0,0))
CO_PIPOpreds.int <- predictSE(CO_PIPOcompmod.int, newdata=CO_PIPOpreddata, se.fit=T)

#included ns interactions
CO_POTRcompmod.int <- lme(BAI10yr ~ Band  * BA_totLRsc + ACFsc + Band:ACFsc + DBHsc, random= ~1|PairUn, CO_POTRinfo, weights=varExp(form=~DBH), method="REML")
CO_POTRpreddata <- data.frame(Band=c("L","M","H"), DBHsc=c(0,0,0), BA_totLRsc = c(1,1,1), in5_totsc = c(1,1,1), ACFsc = c(-1,-1,-1))
CO_POTRpreds.int <- predictSE(CO_POTRcompmod.int, newdata = CO_POTRpreddata, level=0, se.fit=T )

CO_ABLAinfo.na <- CO_ABLAinfo[!is.na(CO_ABLAinfo$BAI10yr),]
CO_ABLAcompmod.int <- gls(BAI10yr ~ Band * BA_totLRsc + DBHsc, weights= varPower(), CO_ABLAinfo.na, method="REML")
CO_ABLApreddata <-  data.frame(Band=c("L","M","H"), DBHsc=c(0,0,0), BA_totLRsc = c(1,1,1))
CO_ABLApreds.int <- predictSE(CO_ABLAcompmod.int, newdata = CO_ABLApreddata, se.fit=T)

## MT
MT_TSHEcompmod.int <- gls(BAI10yr ~ Band * BA_totLRsc + DBHsc, weights=varPower(form=~DBHsc|Band), MT_TSHEinfo, method="REML")
MT_TSHEpreddata <-  data.frame(Band=c("L","M","H"), DBHsc=c(0,0,0), BA_totLRsc = c(1,1,1))
MT_TSHEpreds.int <- predictSE(MT_TSHEcompmod.int, newdata = MT_TSHEpreddata, se.fit=T)

MT_PSMEinfo.cl <- MT_PSMEinfo[-76,]
MT_PSMEcompmod.int <- gls(BAI10yr~ Band * ACFsc, MT_PSMEinfo.cl)
MT_PSMEpreddata <- data.frame(Band=c("L","M","H"), ACFsc=c(-1,-1,-1))
MT_PSMEpreds.int <- predictSE(MT_PSMEcompmod.int, newdata=MT_PSMEpreddata, se.fit=T) # this is identical if you use lm() and predict()
#included ns interactions
MT_ABLAcompmod.int  <- lme(BAI10yr ~ Band * ACFsc + DBHsc + Band:DBHsc , random= ~1|PairUn, data=MT_ABLAinfo, weights=varPower(form=~DBH), method="REML")
MT_ABLApreddata <- data.frame(Band=c("L","M","H"), ACFsc=c(-1,-1,-1), DBHsc = c(0,0,0))
MT_ABLApreds.int <- predictSE(MT_ABLAcompmod.int, newdata = MT_ABLApreddata, level=0, se.fit=T )

### WA
# included ns interaction
WA_TSHEinfo.cl <- WA_TSHEinfo[-c(20,26,43),]
WA_TSHEcompmod.int <- gls(BAI10yr ~ Band * ACFsc+ DBHsc ,weights=varExp(form=~DBH), WA_TSHEinfo.cl, method="REML")
WA_TSHEpreddata <- data.frame(Band=c("M","H"), ACFsc=c(-1,-1), DBHsc = c(0,0))
WA_TSHEpreds.int <- predictSE(WA_TSHEcompmod.int, newdata = WA_TSHEpreddata, se.fit=T )
# included ns interaction
WA_PSMEinfo.cl <- WA_PSMEinfo[-1,]
WA_PSMEcompmod.int <- lme(BAI10yr~ Band * BA_totLRsc + DBHsc, random= ~ 1|PairUn, weights=varPower(form=~DBH), WA_PSMEinfo.cl, method="REML")
WA_PSMEpreddata <- data.frame(Band=c("L","M","H"), BA_totLRsc=c(1,1,1), DBHsc = c(0,0,0))
WA_PSMEpreds.int <- predictSE(WA_PSMEcompmod.int, newdata = WA_PSMEpreddata, level=0, se.fit=T )

WA_ABLAcompmod.int <- gls(BAI10yr ~ Band * BA_totLRsc + DBHsc, weights=varPower(), WA_ABLAinfo, control=list(maxIter=100, msMaxIter=100), method="REML")
WA_ABLApreddata  <- data.frame(Band = c("L","M","H"), BA_totLRsc = c(1,1,1), DBHsc=c(0,0,0))
WA_ABLApreds.int <- predictSE(WA_ABLAcompmod.int, newdata=WA_ABLApreddata, se.fit=T)

# combine the predicted growth with 1sd increased competition from all species into one df
growthint <- c(CO_PIPOpreds.int[[1]],CO_POTRpreds.int[[1]],CO_ABLApreds.int[[1]],MT_TSHEpreds.int[[1]],MT_PSMEpreds.int[[1]],MT_ABLApreds.int[[1]],NA, WA_TSHEpreds.int[[1]],WA_PSMEpreds.int[[1]],WA_ABLApreds.int[[1]])
# calculate the growth sensitvity to competition as the difference between growth at mean competition and growth at +1sd competition.
Compsens.mean.new <- growthint - Sumstatsall$BAI.mod.mean 
  # note: the sign of sensitivity is negative, more sensitive populations = more negative numbers




####### END: Mean Growth and Sensitivity to Competition modeling ############################
#_____________________________________________________________________________________________







#_____________________________________________________________________________________________
#_____________________________________________________________________________________________
############## BEGIN: Growth Synchrony/Climate Sensitivity Analysis ##########################################
#_____________________________________________________________________________________________
#_____________________________________________________________________________________________



# example of data processing for Colorado species-replciates (CO_PIPO, CO_POTR, CO_ABLA), data provided for CO_ABLA: 
# 1) import raw ringwidth data, crossdated and averaged per tree
# 2) create truncated dataframes cut to the length of the shortest tree chronology (to avoid biases earlier in the record when trees start dropping out),
#    and detrend these truncated tree chronologies
# 3) calculate all the pairwise correlations between trees in an elevation band
# 4) transform correlations [-1,1] to [0,1] and perform beta regressions to test
#    for significant differences between range margin and range center growth synchrony
#    The results of this regression then have to be transformed back into correlation scale (from logit scale)





########### . Required Functions ########################################


#--- Function to identify first NA to appear in a chronology --
firstNA <- function(x){
  x.na <- is.na(x)
  ifelse(length(which(x.na==T))==0, return(length(x)), return(min(which(x.na==T))))
}



#-------- core_trunc() ---

# Use: takes a data set (years as columns, trees as rows)
# and turns cores into columns and years into rows,
# with Elev-Tree is colnames
# also truncates to desired length of time
# this format is need to melt this into a long form df and to calculate growth synchrony

# Inupts: workingdata= data frame produced by averaging ring widths from two cores per tree
# first = year to start time series, as character
# last = year to end time series, as character
# detmethod = method (default "Spline") with which to detrend
# detrend = logical, whether to return detrended time series
#       or undetrended ring widths, default = TRUE
# drop = number of the oldest rings to drop (if many cores hit pith the first few years of growth are sometimes weird)
# yrcored = the year the core was collected (last year in the rw dataset)
# prewhiten = logical, whether to subtract out an AR1 model to pull out autocorrelation

#______________________________________
core_trunc <- function (workingdata, first=start, last=end, detmethod="Spline", detrend=TRUE, drop=3, yrcored=2013, prewhiten=F, ...){
  require(dplR)
  library(gmp)
  # Use: takes a data set produced by core_average()
  # and turns cores into columns and years into rows,
  # with Elev-Tree is colnames
  # also truncates to desired length of time
  # if colnames are in X2013 format, this renames them to be pure numerics
  if(is.na(suppressWarnings(as.numeric(colnames(workingdata)[ncol(workingdata)])))){
    colnames(workingdata)[6:ncol(workingdata)]<- seq(yrcored,yrcored-(ncol(workingdata)-6), by=-1)
  }
  # grabbing desired rows + drop from dataframe, unless there aren't enough rows to begin with
  First <- as.character(max((as.numeric(first) - drop), as.numeric(colnames(workingdata)[ncol(workingdata)])))
  if (as.numeric(First)>as.numeric(first)) print ("Time series shorter than desired")  
  # establish vector of column numbers that correspond to the years that we want to pull
  yrs <- seq(from=which(colnames(workingdata)==last),to=which(colnames(workingdata)==First), by=1)
  
  # subsetting working data to just desired range, holding on to the first 4 columns of ID variables
  treedata2011 <- workingdata[,c(1:4,yrs)]
  
  # create a column with a full Tree ID, and move it to be with the other 4 ID columns at beginning of dataframe
  treedata2011$ID <- paste(treedata2011$State,treedata2011$Species,treedata2011$Elev, treedata2011$Tree, sep = "-")
  idcol <- as.numeric(ncol(treedata2011))
  treedata<- treedata2011[,c(1:4,idcol,5:(idcol-1))]  # good to go!
  #________________________________________________________________
  #------------------------------------------------------------
  ## Creating a detrended rw data frame 
  lastyr <- which(colnames(treedata)==last)
  spstdat.rwl<-t(treedata[,lastyr:ncol(treedata)]) ### 
  rownames(spstdat.rwl)<-colnames(treedata[lastyr:ncol(treedata)])
  colnames(spstdat.rwl)<-treedata$ID
  spstdat.rwl.df<-as.data.frame(spstdat.rwl)
  
  if (drop > 0){
    for (i in 1:ncol(spstdat.rwl.df)){
      tmp <- spstdat.rwl.df[,i]
      if (is.na(spstdat.rwl.df[length(tmp),i])){
        l <- min(which(is.na(tmp)))
        drops <- seq(from=l-drop, to=l-1)
        spstdat.rwl.df[drops,i] <- NA
      }
    }
    spstdat.drop <- spstdat.rwl.df[1:which(rownames(spstdat.rwl.df)==first),]
  }
  else {spstdat.drop <- spstdat.rwl.df}
  #__________ Detrend each individual tree chronology (if detrend==T)________
  rw.det <- detrend(rwl=spstdat.drop, method=detmethod)#, ...) # final object if detrending
  rw.nodet <- spstdat.drop # final object if not detrending
  #-----------------------------------------------------------
  #==========================================================
  if(prewhiten==TRUE & detrend==TRUE){
    years <- rownames(rw.det)
    rw.det <- data.frame(apply(rw.det, MARGIN=2, FUN = det.ar))#,...)
    rownames(rw.det) <- years
    colnames(rw.det) <- str_replace_all(colnames(rw.det),pattern = "\\.",replacement = "-")
  }
  if(prewhiten==TRUE & detrend==FALSE){
    years <- rownames(rw.nodet)
    rw.nodet <- data.frame(apply(rw.nodet, MARGIN=2, FUN = det.ar))#,...)
    rownames(rw.nodet) <- years
    colnames(rw.nodet) <- str_replace_all(colnames(rw.nodet),pattern = "\\.",replacement = "-")
  }
  if (detrend==TRUE)
    return(rw.det)
  else
    return(rw.nodet)
}


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




##_________Function to calculate pairwise correlations__________________
# takes data output from core_trunc(), breaks it into elevation bands based on the Tree ID, and returns either all the pairwise correlations or a summary per elevation band
paircors <- function(tmpdata, use.method="pairwise.complete.obs", summaries.only=F){
  
  ## Mid Elevation
  tmpdata.M <- tmpdata[,grep("-M-",colnames(tmpdata))]
  r.M<- cor(tmpdata.M, use = use.method) # get corr matrix
  y.M<- which(lower.tri(r.M), TRUE) # get values to subset to lower triangular
  z.M <- data.frame(Tree1 = rownames(r.M)[y.M[, 1]],Tree2 = colnames(r.M)[y.M[, 2]],cor = r.M[y.M], Elev = rep("M", times=length(r.M[y.M]))) # make a dataframe of each pairwise correlation
  # make a row with the average Pearson Corr coef, the sd, and the number of trees in the stand
  summary.M <-c(mean(z.M[,"cor"], na.rm=TRUE),sd(z.M[,"cor"], na.rm=TRUE), ncol(tmpdata.M))
  
  ## High Elevation
  tmpdata.H <- tmpdata[,grep("-H-",colnames(tmpdata))]
  r.H<- cor(tmpdata.H, use = use.method) # get corr matrix
  y.H<- which(lower.tri(r.H), TRUE) # get values to subset to lower triangular
  z.H <- data.frame(Tree1 = rownames(r.H)[y.H[, 1]],Tree2 = colnames(r.H)[y.H[, 2]],cor = r.H[y.H], Elev = rep("H", times=length(r.H[y.H]))) # make a dataframe of each pairwise correlation
  # make a row with the average Pearson Corr coef, the sd, and the number of trees in the stand
  summary.H <-c(mean(z.H[,"cor"], na.rm=TRUE),sd(z.H[,"cor"], na.rm=TRUE), ncol(tmpdata.H))
  ## Low Elevation
  if(length(grep("-L", colnames(tmpdata))>0)){
    tmpdata.L <- tmpdata[,grep("-L-",colnames(tmpdata))]
    r.L<- cor(tmpdata.L, use = use.method) # get corr matrix
    y.L<- which(lower.tri(r.L), TRUE) # get values to subset to lower triangular
    z.L <- data.frame(Tree1 = rownames(r.L)[y.L[, 1]],Tree2 = colnames(r.L)[y.L[, 2]],cor = r.L[y.L], Elev = rep("L", times=length(r.L[y.L]))) # make a dataframe of each pairwise correlation
    # make a row with the average Pearson Corr coef, the sd, and the number of trees in the stand
    summary.L <-c(mean(z.L[,"cor"], na.rm=TRUE),sd(z.L[,"cor"], na.rm=TRUE), ncol(tmpdata.L))
    sumstats <- cbind(c("L","M","H"),data.frame(rbind(summary.L,summary.M,summary.H)))
    colnames(sumstats) <- c("Elevation", "Meancorr", "SDcorr","Ntrees")
    indcors <- rbind(z.L,z.M,z.H)
  }
  else{
    sumstats <- cbind(c("M","H"),data.frame(rbind(summary.M,summary.H)))
    colnames(sumstats) <- c("Elevation", "Meancorr", "SDcorr","Ntrees")
    indcors <- rbind(z.M,z.H)
  }
  
  if(summaries.only==T){
    return(sumstats)
  }
  else{
    return(indcors)
  }
}




#++++++++++++++ . End import needed functions++++++++++++++++++++++++







#+++++++++++++ 1) Import raw RW dataframes ++++++++++++++++++++++++++++
# Pinus ponderosa (low elev)
# CO_PIPOrw <- read.csv("PIPO_avgRingWidth_2-19-15.csv", header=T)[,-1]
# Populus tremuloides (mid elev)
# CO_POTRrw <- read.csv("POTR_avgRingWidth_2-19-15.csv", header=T)[,-1]
# Abies lasiocarpa (subalpine/high elev)
CO_ABLArw <- read.csv("SIData_S2_RingWidths_CO_ABLA_20180428.csv", header=T)[,-1]



#++++++++++++++ 2) Creat Detrended and truncated RWIs +++++++++++++++

#### Parameters for core_trunc
end <- "2012"     # last date of time series
method <- "Spline" # detrending method desired. c("Spline", "ModNegExp", "Mean")
dropyrs <- 3
####

# note: earliest year ('first') is the earliest year to grap, which we get by finding the column of the first NA for a species-replicate 
#       and then calculating the year based on the column number (2013 + 6 ID columns + 1 because firstNA gives the column before the NA - the column number = the year of the first NA)
# dat <- CO_PIPOrw
# CO_PIPOcut <- core_trunc(workingdata=dat, first=2013+7-min(apply(X = dat,MARGIN = 1,FUN = firstNA)), last=end, detmethod=method, detrend=T, drop=0)
# dat <- CO_POTRrw
# CO_POTRcut <- core_trunc(workingdata=dat, first=2013+7-min(apply(X = dat,MARGIN = 1,FUN = firstNA)), last=end, detmethod=method, detrend=T, drop=0)
dat <- CO_ABLArw
CO_ABLAcut <- core_trunc(workingdata=dat, first=2013+7-min(apply(X = dat,MARGIN = 1,FUN = firstNA)), last=end, detmethod=method, detrend=T, drop=0)




#++++++++++++++ 3) Calcualte all pairwise correlations between trees +++++++++++++++

# 
# CO_PIPOcorscut <- paircors(CO_PIPOcut)
# CO_PIPOsumstats <- paircors(CO_PIPOcut,summaries.only = T)
# 
# CO_POTRcorscut <- paircors(CO_POTRcut)
# CO_POTRsumstats <- paircors(CO_POTRcut,summaries.only = T)

CO_ABLAcorscut <- paircors(CO_ABLAcut)
#CO_ABLAsumstats <- paircors(CO_ABLAcut,summaries.only = T)



#++++++++++++++ 4) perform beta regressions to test for significant difference +++++++++++++++

# 
# #+++++++++++++++++ CO-PIPO beta regression +++++++++++++++++++++
# levels(CO_PIPOcorscut$Elev) <- list(L="L", M="M", H="H")
# CO_PIPOcorscut$cor_trans <- (CO_PIPOcorscut$cor + 1 )/2
# 
# copipo_logit <- betareg(cor_trans ~ Elev - 1 , CO_PIPOcorscut) # default uses logit link
# copipo_logit2 <- betareg(cor_trans ~ Elev - 1 | Elev, CO_PIPOcorscut)
# lrtest(copipo_logit, copipo_logit2) # looks like link function doesn't much matter, but different precisions per elev are important
# summary(copipo_logit2)
# 
# # transforming coefficients from logit scale back to correlation coefficeint scale
# bestmod <- copipo_logit2
# copipo.out <- summary(bestmod) # get output from best model
# betas.logit <- data.frame(copipo.out$coefficients$mean[,1:2],confint(bestmod)[1:3,])
# betas.logit$use <- betas.logit[,1] + betas.logit[,2]
# betas.logit$lse <- betas.logit[,1] - betas.logit[,2]
# betas.prop <- exp(betas.logit)/(1 + exp(betas.logit)) # get back to proportion scale
# betas.corr <- betas.prop * 2 - 1 # then back transform to pearson R scale
# CO_PIPOsumstats$corr.mod.mean <- betas.corr[,1]
# CO_PIPOsumstats$corr.mod.lse <- betas.corr[,6]
# CO_PIPOsumstats$corr.mod.use <- betas.corr[,5]
# CO_PIPOsumstats$corr.mod.ci2.5 <- betas.corr[,3]
# CO_PIPOsumstats$corr.mod.ci97.5 <- betas.corr[,4]
# ## Statistics to report
# levels(CO_PIPOcorscut$Elev) <- list( M="M", L="L", H="H")
# copipo_logitb <- betareg(cor_trans ~ Elev | Elev, CO_PIPOcorscut)
# summary(copipo_logitb) # L diff p= 2.7e-07, H diff p < 2e-16
# # compared to 
# summary(glht(copipo_logit2, linfct = c("ElevH-ElevM=0", "ElevL-ElevM=0", "ElevH-ElevL=0"))) #was based off of http://stats.stackexchange.com/questions/167159/post-hoc-test-for-betareg-model-r
# # HvM p < 1e-06, LvM p<1e-06, LvH p<1e-06
# 
# 
# #++++++++++++++++++ CO-POTR beta regression +++++++++++++++++++
# levels(CO_POTRcorscut$Elev) <- list(L="L", M="M", H="H")
# CO_POTRcorscut$cor_trans <- (CO_POTRcorscut$cor + 1 )/2
# 
# copotr_logit <- betareg(cor_trans ~ Elev - 1 , CO_POTRcorscut) # default uses logit link
# copotr_logit2 <- betareg(cor_trans ~ Elev - 1 | Elev, CO_POTRcorscut)
# lrtest(copotr_logit, copotr_logit2) # looks like link function doesn't much matter, but different precisions per elev are important
# # qqnorm(resid(copotr_logit2)); qqline(resid(copotr_logit2))
# # plot(copotr_logit2, which=1:5, type="pearson") # Whoa shoot! There's a major outlier
# # which(resid(copotr_logit2)>10) # 1204 is the outlier
# copotr_logit2n <- betareg(cor_trans ~ Elev - 1 | Elev, CO_POTRcorscut, subset= -1204)
# # MUCH better
# 
# # transforming coefficients from logit scale back to correlation coefficeint scale
# bestmod <- copotr_logit2n
# copotr.out <- summary(bestmod) # get output from best model
# betas.logit <- data.frame(copotr.out$coefficients$mean[,1:2],confint(bestmod)[1:3,])
# betas.logit$use <- betas.logit[,1] + betas.logit[,2]
# betas.logit$lse <- betas.logit[,1] - betas.logit[,2]
# betas.prop <- exp(betas.logit)/(1 + exp(betas.logit)) # get back to proportion scale
# betas.corr <- betas.prop * 2 - 1 # then back transform to pearson R scale
# CO_POTRsumstats$corr.mod.mean <- betas.corr[,1]
# CO_POTRsumstats$corr.mod.lse <- betas.corr[,6]
# CO_POTRsumstats$corr.mod.use <- betas.corr[,5]
# CO_POTRsumstats$corr.mod.ci2.5 <- betas.corr[,3]
# CO_POTRsumstats$corr.mod.ci97.5 <- betas.corr[,4]
# ## Statistics to report
# levels(CO_POTRcorscut$Elev) <- list( M="M", L="L", H="H")
# copotr_logitb <- betareg(cor_trans ~ Elev | Elev, CO_POTRcorscut)
# summary(copotr_logitb) # L diff p= 0.0816, H diff p < 2e-16
# # compared to 
# summary(glht(copotr_logit2, linfct = c("ElevH-ElevM=0", "ElevL-ElevM=0", "ElevH-ElevL=0"))) #was based off of http://stats.stackexchange.com/questions/167159/post-hoc-test-for-betareg-model-r
# # HvM p < 1e-04, LvM p=0.19, LvH p<1e-04
# 


#+++++++++++++++++++++ CO-ABLA beta regression +++++++++++++++++
levels(CO_ABLAcorscut$Elev) <- list(L="L", M="M", H="H")
# transform correlations to be [0,1]
CO_ABLAcorscut$cor_trans <- (CO_ABLAcorscut$cor + 1 )/2

coabla_logit <- betareg(cor_trans ~ Elev - 1 , CO_ABLAcorscut) # default uses logit link
coabla_logit2 <- betareg(cor_trans ~ Elev -1 | Elev, CO_ABLAcorscut)
lrtest(coabla_logit, coabla_logit2) # different precisions per elev are important
# qqnorm(resid(coabla_logit2)); qqline(resid(coabla_logit2)) # looks ok
# plot(coabla_logit2, which=1:5, type="pearson") #actually looks pretty good. better than POTR even
coabla_betareg <- update(coabla_logit2, ~ +1)

# transforming coefficients from logit scale back to correlation coefficeint scale
bestmod <- coabla_logit2
coabla.out <- summary(bestmod) # get output from best model
betas.logit <- data.frame(coabla.out$coefficients$mean[,1:2],confint(bestmod)[1:3,])
betas.logit$use <- betas.logit[,1] + betas.logit[,2]
betas.logit$lse <- betas.logit[,1] - betas.logit[,2]
betas.prop <- exp(betas.logit)/(1 + exp(betas.logit)) # get back to proportion scale
betas.corr <- betas.prop * 2 - 1 # then back transform to pearson R scale
## adding results to species-replicate summary
# CO_ABLAsumstats$corr.mod.mean <- betas.corr[,1]
# CO_ABLAsumstats$corr.mod.lse <- betas.corr[,6]
# CO_ABLAsumstats$corr.mod.use <- betas.corr[,5]
# CO_ABLAsumstats$corr.mod.ci2.5 <- betas.corr[,3]
# CO_ABLAsumstats$corr.mod.ci97.5 <- betas.corr[,4]
## Statistics to report
levels(CO_ABLAcorscut$Elev) <- list( M="M", L="L", H="H")
coabla_logitb <- betareg(cor_trans ~ Elev | Elev, CO_ABLAcorscut)
summary(coabla_logitb) # L diff p= 1.63e-14, H diff p =.251




####### END: Growth Synchrony ############################
#_____________________________________________________________________________________________













#_____________________________________________________________________________________________
#_____________________________________________________________________________________________
############## BEGIN: Survival Modeling Example ##########################################
#_____________________________________________________________________________________________
#_____________________________________________________________________________________________
# for the recruitment and survival data used in this section will be archived on Dryad upon publication.


## Goal: analyze the survival rates of conpsecific trees captured in the variable radius plots
# centered on each cored tree. Trees falling inside the plots were id'ed to species and characterized
# as Live, red-dead (still had red needles in crown, usually 0-2 years post mortality),
# gray-dead (still had all bark and small branches, usually <5 years post mortality), and snags.
# For this analysis, we combined red-dead and gray-dead trees to asses survival rates over the previous ~5 years
# (decomposition rates may vary by species and/or by mountain, but probably not between the range margins of a species-replicate)

# Approach: 1) Turn data from variable radius plots (# of live and dead conspecific trees)
#             into binary data (0= live tree, 1= dead tree)
#           2) model survival as a function of Elevation band + a plot random effect (termed 'Tag') using
#             a binomial error distribution and logit link
#           3) acquire confidence intervals in back transformed response (survival probability) space for model predictions
#             based only on the fixed effects using the predictInterval() function from {merTools}



## turn # of trees captured BAF prism variable radius plot into binary data frame
  # first combine red-dead (_R) and gray-dead (_G) trees into total dead (_d)
PIPOba <- CO_PIPOinfo %>% mutate(PIPO_d = PIPO_R + PIPO_G) %>% select(Tag, Band, PIPO_L, PIPO_d)
## live trees
PIPObin.baL <- PIPOba[rep(1:nrow(PIPOba), PIPOba$PIPO_L),c("Tag","Band")]
PIPObin.baL$L_D <- rep(1, nrow(PIPObin.baL))
## dead trees
PIPObin.baD <- PIPOba[rep(1:nrow(PIPOba), PIPOba$PIPO_d),c("Tag","Band")]
PIPObin.baD$L_D <- rep(0, nrow(PIPObin.baD))
## combine all
PIPObin.ba <- rbind(PIPObin.baL, PIPObin.baD) %>% arrange(Tag)
PIPObin.ba$Tag <- as.numeric(as.factor(PIPObin.ba$Tag))


  # set mid elevation range center as intercept
levels(PIPObin.ba$Band) <- list(M="M",L="L",H="H")
# model survival as glmm using {lme4} with a plot random effect, bionial error and logit link
PIPOmort <- glmer(L_D ~ Band + (1|Tag), PIPObin.ba, family=binomial(link='logit'))
# get confidence intervals in terms of survival probability
CO_PIPOmortprd <- predictInterval(PIPOmort, newdata=data.frame(Band=c("L","M","H"), Tag=c(1,1,1)),which="fixed", type="probability")



##### POTR
## turn BAF into binary df
CO_POTRba <- CO_POTRinfo %>% mutate(POTR_d = POTR_R + POTR_G) %>% select(Tag, Band, POTR_L, POTR_d)
## live trees
CO_POTRbin.baL <- CO_POTRba[rep(1:nrow(CO_POTRba), CO_POTRba$POTR_L),c("Tag","Band")]
CO_POTRbin.baL$L_D <- rep(1, nrow(CO_POTRbin.baL))
## dead trees
CO_POTRbin.baD <- CO_POTRba[rep(1:nrow(CO_POTRba), CO_POTRba$POTR_d),c("Tag","Band")]
CO_POTRbin.baD$L_D <- rep(0, nrow(CO_POTRbin.baD))
## combine all
CO_POTRbin.ba <- rbind(CO_POTRbin.baL, CO_POTRbin.baD) %>% arrange(Tag)
CO_POTRbin.ba$Tag <- as.numeric(as.factor(CO_POTRbin.ba$Tag))

CO_POTRbin <- CO_POTRbin.ba 
levels(CO_POTRbin$Band) <- list(M="M",L="L",H="H")
CO_POTRmort <- glmer(L_D ~ Band + (1|Tag), CO_POTRbin, family=binomial(link='logit'))
CO_POTRmortprd <- predictInterval(CO_POTRmort, newdata=data.frame(Band=c("L","M","H"), Tag=c(1,1,1)),which="fixed", type="probability")



##### CO-ABLA
## turn BAF into binary df
CO_ABLAba <- CO_ABLAinfo %>% mutate(ABLA_d = ABLA_R + ABLA_G) %>% select(Tag, Band, ABLA_L, ABLA_d)
## live trees
CO_ABLAbin.baL <- CO_ABLAba[rep(1:nrow(CO_ABLAba), CO_ABLAba$ABLA_L),c("Tag","Band")]
CO_ABLAbin.baL$L_D <- rep(1, nrow(CO_ABLAbin.baL))
## dead trees
CO_ABLAbin.baD <- CO_ABLAba[rep(1:nrow(CO_ABLAba), CO_ABLAba$ABLA_d),c("Tag","Band")]
CO_ABLAbin.baD$L_D <- rep(0, nrow(CO_ABLAbin.baD))
## combine all
CO_ABLAbin.ba <- rbind(CO_ABLAbin.baL, CO_ABLAbin.baD) %>% arrange(Tag)
CO_ABLAbin.ba$Tag <- as.numeric(as.factor(CO_ABLAbin.ba$Tag))

CO_ABLAbin <- CO_ABLAbin.ba 
levels(CO_ABLAbin$Band) <- list(M="M",L="L",H="H")
CO_ABLAmort <- glmer(L_D ~ Band + (1|Tag), CO_ABLAbin, family=binomial(link='logit'))
CO_ABLAmortprd <- predictInterval(CO_ABLAmort, newdata=data.frame(Band=c("L","M","H"), Tag=c(1,1,1)),which="fixed", type="probability")



########## Montana Mortality ############
##### MT-TSHE
## turn BAF into binary df
MT_TSHEba <- MT_TSHEinfo %>% mutate(TSHE_d = TSHE_R + TSHE_G) %>% select(Tag, Band, TSHE_L, TSHE_d)
## live trees
MT_TSHEbin.baL <- MT_TSHEba[rep(1:nrow(MT_TSHEba), MT_TSHEba$TSHE_L),c("Tag","Band")]
MT_TSHEbin.baL$L_D <- rep(1, nrow(MT_TSHEbin.baL))
## dead trees
MT_TSHEbin.baD <- MT_TSHEba[rep(1:nrow(MT_TSHEba), MT_TSHEba$TSHE_d),c("Tag","Band")]
MT_TSHEbin.baD$L_D <- rep(0, nrow(MT_TSHEbin.baD))
## combine all
MT_TSHEbin.ba <- rbind(MT_TSHEbin.baL, MT_TSHEbin.baD) %>% arrange(Tag)
MT_TSHEbin.ba$Tag <- as.numeric(as.factor(MT_TSHEbin.ba$Tag))

MT_TSHEbin <- MT_TSHEbin.ba 
levels(MT_TSHEbin$Band) <- list(M="M",L="L",H="H")
MT_TSHEmort <- glmer(L_D ~ Band + (1|Tag), MT_TSHEbin, family=binomial(link='logit'))
MT_TSHEmortprd <- predictInterval(MT_TSHEmort, newdata=data.frame(Band=c("L","M","H"), Tag=c(1,1,1)),which="fixed", type="probability")



##### MT-PSME
## turn BAF into binary df
MT_PSMEba <- MT_PSMEinfo %>% mutate(PSME_d = PSME_R + PSME_G) %>% select(Tag, Band, PSME_L, PSME_d)
## live trees
MT_PSMEbin.baL <- MT_PSMEba[rep(1:nrow(MT_PSMEba), MT_PSMEba$PSME_L),c("Tag","Band")]
MT_PSMEbin.baL$L_D <- rep(1, nrow(MT_PSMEbin.baL))
## dead trees
MT_PSMEbin.baD <- MT_PSMEba[rep(1:nrow(MT_PSMEba), MT_PSMEba$PSME_d),c("Tag","Band")]
MT_PSMEbin.baD$L_D <- rep(0, nrow(MT_PSMEbin.baD))
## combine all
MT_PSMEbin.ba <- rbind(MT_PSMEbin.baL, MT_PSMEbin.baD) %>% arrange(Tag)
MT_PSMEbin.ba$Tag <- as.numeric(as.factor(MT_PSMEbin.ba$Tag))

MT_PSMEbin <- MT_PSMEbin.ba 
levels(MT_PSMEbin$Band) <- list(M="M",L="L",H="H")
MT_PSMEmort <- glmer(L_D ~ Band + (1|Tag), MT_PSMEbin, family=binomial(link='logit'))
MT_PSMEmortprd <- predictInterval(MT_PSMEmort, newdata=data.frame(Band=c("L","M","H"), Tag=c(1,1,1)),which="fixed", type="probability")


##### MT-ABLA
## turn BAF into binary df
MT_ABLAba <- MT_ABLAinfo %>% mutate(ABLA_d = ABLA_R + ABLA_G) %>% select(Tag, Band, ABLA_L, ABLA_d)
## live trees
MT_ABLAbin.baL <- MT_ABLAba[rep(1:nrow(MT_ABLAba), MT_ABLAba$ABLA_L),c("Tag","Band")]
MT_ABLAbin.baL$L_D <- rep(1, nrow(MT_ABLAbin.baL))
## dead trees
MT_ABLAbin.baD <- MT_ABLAba[rep(1:nrow(MT_ABLAba), MT_ABLAba$ABLA_d),c("Tag","Band")]
MT_ABLAbin.baD$L_D <- rep(0, nrow(MT_ABLAbin.baD))
## combine all
MT_ABLAbin.ba <- rbind(MT_ABLAbin.baL, MT_ABLAbin.baD) %>% arrange(Tag)
MT_ABLAbin.ba$Tag <- as.numeric(as.factor(MT_ABLAbin.ba$Tag))

MT_ABLAbin <- MT_ABLAbin.ba 
levels(MT_ABLAbin$Band) <- list(M="M",L="L",H="H")
MT_ABLAmort <- glmer(L_D ~ Band + (1|Tag), MT_ABLAbin, family=binomial(link='logit'))
MT_ABLAmortprd <- predictInterval(MT_ABLAmort, newdata=data.frame(Band=c("L","M","H"), Tag=c(1,1,1)),which="fixed", type="probability")




########## WA Mortality ############
##### WA-TSHE
## turn BAF into binary df
WA_TSHEba <- WA_TSHEinfo %>% mutate(TSHE_d = TSHE_R + TSHE_G) %>% select(Tag, Band, TSHE_L, TSHE_d)
## live trees
WA_TSHEbin.baL <- WA_TSHEba[rep(1:nrow(WA_TSHEba), WA_TSHEba$TSHE_L),c("Tag","Band")]
WA_TSHEbin.baL$L_D <- rep(1, nrow(WA_TSHEbin.baL))
## dead trees
WA_TSHEbin.baD <- WA_TSHEba[rep(1:nrow(WA_TSHEba), WA_TSHEba$TSHE_d),c("Tag","Band")]
WA_TSHEbin.baD$L_D <- rep(0, nrow(WA_TSHEbin.baD))
## combine all
WA_TSHEbin.ba <- rbind(WA_TSHEbin.baL, WA_TSHEbin.baD) %>% arrange(Tag)
WA_TSHEbin.ba$Tag <- as.numeric(as.factor(WA_TSHEbin.ba$Tag))

WA_TSHEbin <- WA_TSHEbin.ba 
levels(WA_TSHEbin$Band) <- list(M="M",H="H")
WA_TSHEmort <- glmer(L_D ~ Band + (1|Tag), WA_TSHEbin, family=binomial(link='logit'))
WA_TSHEmortprd <- rbind(c(NA,NA,NA),predictInterval(WA_TSHEmort, newdata=data.frame(Band=c("M","H"), Tag=c(1,1)),which="fixed", type="probability"))



##### WA-PSME
## turn BAF into binary df
WA_PSMEba <- WA_PSMEinfo %>% mutate(PSME_d = PSME_R + PSME_G) %>% select(Tag, Band, PSME_L, PSME_d)
## live trees
WA_PSMEbin.baL <- WA_PSMEba[rep(1:nrow(WA_PSMEba), WA_PSMEba$PSME_L),c("Tag","Band")]
WA_PSMEbin.baL$L_D <- rep(1, nrow(WA_PSMEbin.baL))
## dead trees
WA_PSMEbin.baD <- WA_PSMEba[rep(1:nrow(WA_PSMEba), WA_PSMEba$PSME_d),c("Tag","Band")]
WA_PSMEbin.baD$L_D <- rep(0, nrow(WA_PSMEbin.baD))
## combine all
WA_PSMEbin.ba <- rbind(WA_PSMEbin.baL, WA_PSMEbin.baD) %>% arrange(Tag)
WA_PSMEbin.ba$Tag <- as.numeric(as.factor(WA_PSMEbin.ba$Tag))

WA_PSMEbin <- WA_PSMEbin.ba 
levels(WA_PSMEbin$Band) <- list(M="M",L="L",H="H")
WA_PSMEmort <- glmer(L_D ~ Band + (1|Tag), WA_PSMEbin, family=binomial(link='logit'))
WA_PSMEmortprd <- predictInterval(WA_PSMEmort, newdata=data.frame(Band=c("L","M","H"), Tag=c(1,1,1)),which="fixed", type="probability")


##### WA-ABLA
## turn BAF into binary df
WA_ABLAba <- WA_ABLAinfo %>% mutate(ABLA_d = ABLA_R + ABLA_G) %>% select(Tag, Band, ABLA_L, ABLA_d)
## live trees
WA_ABLAbin.baL <- WA_ABLAba[rep(1:nrow(WA_ABLAba), WA_ABLAba$ABLA_L),c("Tag","Band")]
WA_ABLAbin.baL$L_D <- rep(1, nrow(WA_ABLAbin.baL))
## dead trees
WA_ABLAbin.baD <- WA_ABLAba[rep(1:nrow(WA_ABLAba), WA_ABLAba$ABLA_d),c("Tag","Band")]
WA_ABLAbin.baD$L_D <- rep(0, nrow(WA_ABLAbin.baD))
## combine all
WA_ABLAbin.ba <- rbind(WA_ABLAbin.baL, WA_ABLAbin.baD) %>% arrange(Tag)
WA_ABLAbin.ba$Tag <- as.numeric(as.factor(WA_ABLAbin.ba$Tag))

WA_ABLAbin <- WA_ABLAbin.ba 
levels(WA_ABLAbin$Band) <- list(M="M",L="L",H="H")
WA_ABLAmort <- glmer(L_D ~ Band + (1|Tag), WA_ABLAbin, family=binomial(link='logit'))
WA_ABLAmortprd <- predictInterval(WA_ABLAmort, newdata=data.frame(Band=c("L","M","H"), Tag=c(1,1,1)),which="fixed", type="probability")


#### combine all the predicted mortality probabilities #######
morts <- rbind(CO_PIPOmortprd, CO_POTRmortprd, CO_ABLAmortprd, MT_TSHEmortprd,MT_PSMEmortprd, MT_ABLAmortprd, WA_TSHEmortprd,WA_PSMEmortprd,WA_ABLAmortprd)



####### END: Survival Modeling ############################
#_____________________________________________________________________________________________













#_____________________________________________________________________________________________
#_____________________________________________________________________________________________
############## BEGIN: Regeneration Density Modeling ##########################################
#_____________________________________________________________________________________________
#_____________________________________________________________________________________________
# data used in this analysis will be archived on Dryad upon publication.

## goal: use the number of total saplings/seedlings within 5m of all cored trees to estimate
# the regeneration density of conspecifics at each elevation band.
# We will model $rtot 
# (total conspecific regeneration, originally divided into size classes but summed for analysis)
# as a function of Elevation band ($Band), some metric of competitive environment
#  - $BA_tot = total basal area around focal tree asses with wedge prism
#  - $in5_tot = number of trees w/in 5m of the focal tree
#  - $BA_same_all = stand basal area of conspecific trees only from wedge prism plots

# note: for extracting final model estimates for plotting in Figure 5b, competitive environment
#       metrics were mean centered and sd scaled (XXXXsc variables) so that predicting regen densities
#       at a competitive metric of 0 = densities at mean competitive density.
# also, MT-TSHE showed a positive relationship with BA_same and negative relationship with in5_tot.
# This was interpreted as both competitive suppression and seed limitation, and only competitive suppression
# is shown in Figure S7 (i.e. BA_same is fixed at the mean for prediction lines drawn in figure)


################ v Colorado v #########################


#### CO_PIPO
  # note: not enough regen was detected near focal trees, so we had to augment with regeneration transects
# read in regeneration transect results
#PIPOregen <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/TreeCoresRanges_LDLA/CO_RingWidths/PIPO_regendsurveys_090816.csv", header=T)
PIPOregen <- read.csv("PIPO_regendsurveys_090816.csv", header=T)
PIPOregen$TransUn <- paste(PIPOregen$Elev,PIPOregen$Transect)
PIPOtotregen <- PIPOregen %>% group_by (Elev, Transect) %>% summarise(RegenAll=sum(Count))

levels(PIPOtotregen$Elev) <- list(L="L", M="M", H="H")
regmod <- glm(RegenAll~Elev, data=PIPOtotregen, family = poisson) # significant diffs for both margins
predregmean <- predict(regmod, newdata = list( Elev = c("L","M","H")), se.fit=TRUE, type="response")
predregci <- exp(confint(update(regmod, .~.-1))[1:3,])
cs.reg <- matrix(NA, nrow=3, ncol=2)
compsens.reg <- c(cs.reg[1,1], cs.reg[1,1] + cs.reg[2,1], cs.reg[1,1] + cs.reg[3,1])
predsCO_PIPO <- data.frame(reg.mod.mean=predregmean$fit, reg.mod.lci = predregci[,1], reg.mod.uci = predregci[,2], reg.comspens = compsens.reg, reg.compsens.se = cs.reg[,2])




#### CO_POTR
### simplifying it to just model total regen
levels(CO_POTRinfo$Band) <- list(M="M",L="L", H="H")
CO_POTRinfo$rtot <- with(CO_POTRinfo, POTR_r1 + POTR_r2 + POTR_r3 + POTR_r4) # combined regeneration of all sizes

# model selection to determine which index of competition is the best
regtotCOPOTR0 <- glm(rtot~Band , data=CO_POTRinfo, family = poisson)
regtotCOPOTR1 <- glm(rtot~Band + BA_tot, data=CO_POTRinfo, family = poisson)
regtotCOPOTR1.1 <- glm(rtot~Band * BA_tot, data=CO_POTRinfo, family = poisson)
regtotCOPOTR2 <- glm(rtot~Band + in5_tot, data=CO_POTRinfo, family = poisson)
regtotCOPOTR2.1 <- glm(rtot~Band * in5_tot, data=CO_POTRinfo, family = poisson)
regtotCOPOTR3 <- glm(rtot~Band + BA_same_all, data=CO_POTRinfo, family = poisson)
regtotCOPOTR3.1 <- glm(rtot~Band * BA_same_all, data=CO_POTRinfo, family = poisson)
# BIC(regtotCOPOTR0, regtotCOPOTR1,regtotCOPOTR1.1, regtotCOPOTR2, regtotCOPOTR2.1, regtotCOPOTR3, regtotCOPOTR3.1)
AIC(regtotCOPOTR0, regtotCOPOTR1,regtotCOPOTR1.1, regtotCOPOTR2, regtotCOPOTR2.1, regtotCOPOTR3, regtotCOPOTR3.1)
CO_POTRinfo$reghat<- predict(regtotCOPOTR2.1, type="response")


### extracting model parameters
levels(CO_POTRinfo$Band) <- list(L="L", M="M", H="H")
regmod <- glm(rtot~Band * in5_totsc, data=CO_POTRinfo, family = poisson)
predregmean <- predict(regmod, newdata = list(in5_totsc = rep(0, times = 3), Band = c("L","M","H")), se.fit=TRUE, type="response")
predregci <- exp(confint(update(regmod, .~.-1))[1:3,])
cs.reg <- coef(summary(regmod))[4:6,1:2]
compsens.reg <- c(cs.reg[1,1], cs.reg[1,1] + cs.reg[2,1], cs.reg[1,1] + cs.reg[3,1])
predsCO_POTR <- data.frame(reg.mod.mean=predregmean$fit, reg.mod.lci = predregci[,1], reg.mod.uci = predregci[,2], reg.comspens = compsens.reg, reg.compsens.se = cs.reg[,2])



#### CO_ABLA
CO_ABLAinfo$rtot <- with(CO_ABLAinfo, ABLA_r1 + ABLA_r2 + ABLA_r3 + ABLA_r4)
levels(CO_ABLAinfo$Band) <- list(M="M",L="L", H="H")
regtotCOABLA0 <- glm(rtot~Band , data=CO_ABLAinfo, family = poisson)
regtotCOABLA1 <- glm(rtot~Band + BA_tot, data=CO_ABLAinfo, family = poisson)
regtotCOABLA1.1 <- glm(rtot~Band * BA_tot, data=CO_ABLAinfo, family = poisson)
regtotCOABLA2 <- glm(rtot~Band + in5_tot, data=CO_ABLAinfo, family = poisson)
regtotCOABLA2.1 <- glm(rtot~Band * in5_tot, data=CO_ABLAinfo, family = poisson)
regtotCOABLA3 <- glm(rtot~Band + BA_same_all, data=CO_ABLAinfo, family = poisson)
regtotCOABLA3.1 <- glm(rtot~Band * BA_same_all, data=CO_ABLAinfo, family = poisson)
# BIC(regtotCOABLA0, regtotCOABLA1,regtotCOABLA1.1, regtotCOABLA2, regtotCOABLA2.1, regtotCOABLA3, regtotCOABLA3.1)
AIC(regtotCOABLA0, regtotCOABLA1,regtotCOABLA1.1, regtotCOABLA2, regtotCOABLA2.1, regtotCOABLA3, regtotCOABLA3.1)
  # either 2 (BIC) or 2.1 (AIC) is the best, depending on how heavily you penalize complexity
CO_ABLAinfo$reghat <- NA
CO_ABLAinfo$reghat[-which(is.na(CO_ABLAinfo$rtot))]<- predict(regtotCOABLA2.1, type="response")
CO_ABLAreg <- CO_ABLAinfo[-which(is.na(CO_ABLAinfo$rtot)),]
CO_ABLAreg$reghat<- predict(regtotCOABLA2.1, type="response")

### extracting model parameters
levels(CO_ABLAreg$Band) <- list(L="L", M="M", H="H")
# full model from which model parameters are drawn for Figure 5b (essentially identical to above)
regmod <- glm(rtot~Band * in5_totsc, data=CO_ABLAinfo, family = poisson)
predregmean <- predict(regmod, newdata = list(in5_totsc = rep(0, times = 3), Band = c("L","M","H")), se.fit=TRUE, type="response")
predregci <- exp(confint(update(regmod, .~.-1))[1:3,])
cs.reg <- coef(summary(regmod))[4:6,1:2] # no intercept, so could either fit an intercept model (like I did for BA), or make them all identical
compsens.reg <- c(cs.reg[1,1], cs.reg[1,1] + cs.reg[2,1], cs.reg[1,1] + cs.reg[3,1])
predsCO_ABLA <- data.frame(reg.mod.mean=predregmean$fit, reg.mod.lci = predregci[,1], reg.mod.uci = predregci[,2], reg.comspens = compsens.reg, reg.compsens.se = cs.reg[,2])


#--------- ^ Colorado ^ ----------------------------------------------



################ v Montana v #########################

###   MT_TSHE
### simplifying it to just model total regen
MT_TSHEinfo$rtot <- with(MT_TSHEinfo, TSHE_r1 + TSHE_r2 + TSHE_r3 + TSHE_r4)
levels(MT_TSHEinfo$Band) <- list(M="M",L="L", H="H")
# model selection
regtotMTTSHE0 <- glm(rtot~Band , data=MT_TSHEinfo, family = poisson)
regtotMTTSHE1 <- glm(rtot~Band + BA_tot, data=MT_TSHEinfo, family = poisson)
regtotMTTSHE1.1 <- glm(rtot~Band * BA_tot, data=MT_TSHEinfo, family = poisson)
regtotMTTSHE2 <- glm(rtot~Band + in5_tot, data=MT_TSHEinfo, family = poisson)
regtotMTTSHE2.1 <- glm(rtot~Band * in5_tot, data=MT_TSHEinfo, family = poisson)
regtotMTTSHE3 <- glm(rtot~Band + BA_same_all, data=MT_TSHEinfo, family = poisson)
regtotMTTSHE3.1 <- glm(rtot~Band * BA_same_all, data=MT_TSHEinfo, family = poisson)
# BIC(regtotMTTSHE0, regtotMTTSHE1,regtotMTTSHE1.1, regtotMTTSHE2, regtotMTTSHE2.1, regtotMTTSHE3, regtotMTTSHE3.1)
AIC(regtotMTTSHE0, regtotMTTSHE1,regtotMTTSHE1.1, regtotMTTSHE2, regtotMTTSHE2.1, regtotMTTSHE3, regtotMTTSHE3.1)
# fitted values
MT_TSHEinfo$reghat<- predict(regtotMTTSHE3.1, type="response")


### extracting model estimates
levels(MT_TSHEinfo$Band) <- list(L="L", M="M", H="H")
MT_TSHEinfo$BA_same_allsc <- as.numeric(scale(MT_TSHEinfo$BA_same_all))
regmod <- glm(rtot~Band * BA_same_allsc, data=MT_TSHEinfo, family = poisson)
predregmean <- predict(regmod, newdata = list(BA_same_allsc = rep(0, times = 3), Band = c("L","M","H")), se.fit=TRUE, type="response")
predregci <- exp(confint(update(regmod, .~.-1))[1:3,])
cs.reg <- coef(summary(regmod))[4:6,1:2] # no intercept, so could either fit an intercept model (like I did for BA), or make them all identical
compsens.reg <- c(cs.reg[1,1], cs.reg[1,1] + cs.reg[2,1], cs.reg[1,1] + cs.reg[3,1])
predsMT_TSHE <- data.frame(reg.mod.mean=predregmean$fit, reg.mod.lci = predregci[,1], reg.mod.uci = predregci[,2], reg.comspens = compsens.reg, reg.compsens.se = cs.reg[,2])



#### MT_PSME
### simplifying it to just model total regen
MT_PSMEinfo$rtot <- with(MT_PSMEinfo, PSME_r1 + PSME_r2 + PSME_r3 + PSME_r4)
levels(MT_PSMEinfo$Band) <- list(M="M",L="L", H="H")
regtotMTPSME0 <- glm(rtot~Band , data=MT_PSMEinfo, family = poisson)
regtotMTPSME1 <- glm(rtot~Band + BA_tot, data=MT_PSMEinfo, family = poisson)
regtotMTPSME1.1 <- glm(rtot~Band * BA_tot, data=MT_PSMEinfo, family = poisson)
regtotMTPSME2 <- glm(rtot~Band + in5_tot, data=MT_PSMEinfo, family = poisson)
regtotMTPSME2.1 <- glm(rtot~Band * in5_tot, data=MT_PSMEinfo, family = poisson)
regtotMTPSME3 <- glm(rtot~Band + BA_same_all, data=MT_PSMEinfo, family = poisson)
regtotMTPSME3.1 <- glm(rtot~Band * BA_same_all, data=MT_PSMEinfo, family = poisson)
# BIC(regtotMTPSME0, regtotMTPSME1,regtotMTPSME1.1, regtotMTPSME2, regtotMTPSME2.1, regtotMTPSME3, regtotMTPSME3.1)
AIC(regtotMTPSME0, regtotMTPSME1,regtotMTPSME1.1, regtotMTPSME2, regtotMTPSME2.1, regtotMTPSME3, regtotMTPSME3.1)

MT_PSMEinfo$reghat<- predict(regtotMTPSME3.1, type="response")

# extracing model estimates
levels(MT_PSMEinfo$Band) <- list(L="L", M="M", H="H")
MT_PSMEinfo$BA_same_allsc <- as.numeric(scale(MT_PSMEinfo$BA_same_all))
regmod <- glm(rtot~Band * BA_same_allsc, data=MT_PSMEinfo, family = poisson)
predregmean <- predict(regmod, newdata = list(BA_same_allsc = rep(0, times = 3), Band = c("L","M","H")), se.fit=TRUE, type="response")
predregci <- exp(confint(update(regmod, .~.-1))[1:3,])
cs.reg <- coef(summary(regmod))[4:6,1:2] # no intercept, so could either fit an intercept model (like I did for BA), or make them all identical
compsens.reg <- c(cs.reg[1,1], cs.reg[1,1] + cs.reg[2,1], cs.reg[1,1] + cs.reg[3,1])
predsMT_PSME <- data.frame(reg.mod.mean=predregmean$fit, reg.mod.lci = predregci[,1], reg.mod.uci = predregci[,2], reg.comspens = compsens.reg, reg.compsens.se = cs.reg[,2])



#### MT_ABLA
### simplifying it to just model total regen
MT_ABLAinfo$rtot <- with(MT_ABLAinfo, ABLA_r1 + ABLA_r2 + ABLA_r3 + ABLA_r4)
levels(MT_ABLAinfo$Band) <- list(M="M", L="L",H="H")
regtotMTABLA0 <- glm(rtot~Band , data=MT_ABLAinfo, family = poisson)
regtotMTABLA1 <- glm(rtot~Band + BA_tot, data=MT_ABLAinfo, family = poisson)
regtotMTABLA1.1 <- glm(rtot~Band * BA_tot, data=MT_ABLAinfo, family = poisson)
regtotMTABLA2 <- glm(rtot~Band + in5_tot, data=MT_ABLAinfo, family = poisson)
regtotMTABLA2.1 <- glm(rtot~Band * in5_tot, data=MT_ABLAinfo, family = poisson)
regtotMTABLA3 <- glm(rtot~Band + BA_same_all, data=MT_ABLAinfo, family = poisson)
regtotMTABLA3.1 <- glm(rtot~Band * BA_same_all, data=MT_ABLAinfo, family = poisson)
# BIC(regtotMTABLA0, regtotMTABLA1,regtotMTABLA1.1, regtotMTABLA2, regtotMTABLA2.1, regtotMTABLA3, regtotMTABLA3.1)
AIC(regtotMTABLA0, regtotMTABLA1,regtotMTABLA1.1, regtotMTABLA2, regtotMTABLA2.1, regtotMTABLA3, regtotMTABLA3.1)
MT_ABLAinfo$reghat<- predict(regtotMTABLA2.1, type="response")

# extract model estimates
levels(MT_ABLAinfo$Band) <- list(L="L", M="M", H="H")
regmod <- glm(rtot~Band * in5_totsc, data=MT_ABLAinfo, family = poisson)
predregmean <- predict(regmod, newdata = list(in5_totsc = rep(0, times = 3), Band = c("L","M","H")), se.fit=TRUE, type="response")
predregci <- exp(confint(update(regmod, .~.-1))[1:3,])
cs.reg <- coef(summary(regmod))[4:6,1:2] # no intercept, so could either fit an intercept model (like I did for BA), or make them all identical
compsens.reg <- c(cs.reg[1,1], cs.reg[1,1] + cs.reg[2,1], cs.reg[1,1] + cs.reg[3,1])
predsMT_ABLA <- data.frame(reg.mod.mean=predregmean$fit, reg.mod.lci = predregci[,1], reg.mod.uci = predregci[,2], reg.comspens = compsens.reg, reg.compsens.se = cs.reg[,2])

# note: outliers are present, but removing them doesn't qualitatively influence the inference

#--------- ^ Montanta ^ ----------------------------------------------



################ v Washington v #########################

#### WA_TSHE
### simplifying it to just model total regen
WA_TSHEinfo$rtot <- with(WA_TSHEinfo, TSHE_r1 + TSHE_r2 + TSHE_r3 + TSHE_r4)
levels(WA_TSHEinfo$Band) <- list(M="M", H="H")
#WA_TSHEinfo.short <- WA_TSHEinfo[-which(is.na(WA_TSHEinfo$rtot)),]
regtotWATSHE0 <- glm(rtot~Band , data=WA_TSHEinfo, family = poisson)
regtotWATSHE1 <- glm(rtot~Band + BA_tot, data=WA_TSHEinfo, family = poisson)
regtotWATSHE1.1 <- glm(rtot~Band * BA_tot, data=WA_TSHEinfo, family = poisson)
regtotWATSHE2 <- glm(rtot~Band + in5_tot, data=WA_TSHEinfo, family = poisson)
regtotWATSHE2.1 <- glm(rtot~Band * in5_tot, data=WA_TSHEinfo, family = poisson)
regtotWATSHE3 <- glm(rtot~Band + BA_same_all, data=WA_TSHEinfo, family = poisson)
regtotWATSHE3.1 <- glm(rtot~Band * BA_same_all, data=WA_TSHEinfo, family = poisson)
regtotWATSHE4 <- glm(rtot~Band + in5_tot + BA_same_all, data=WA_TSHEinfo, family=poisson)
#BIC(regtotWATSHE0, regtotWATSHE1,regtotWATSHE1.1, regtotWATSHE2, regtotWATSHE2.1, regtotWATSHE3, regtotWATSHE3.1)
AIC(regtotWATSHE0, regtotWATSHE1,regtotWATSHE1.1, regtotWATSHE2, regtotWATSHE2.1, regtotWATSHE3, regtotWATSHE3.1, regtotWATSHE4)
#WA_TSHEinfo$reghat <- NA
WA_TSHEinfo$reghat<- predict(regtotWATSHE4, type="response") # 3 and 2 are almost indestinguishable, but show different signs.
  # so added 4 w/ both and it is by far and away the best model
# for Figure S7, will show effect of trees in 5 at a constant BA
watshe <- data.frame(in5_tot= rep(seq(0,8, by=1), times=2), BA_same_all = rep(mean(WA_TSHEinfo$BA_same_all), times = 18), Band = rep(c("M","H"), each=9))
reghat <- predict(regtotWATSHE4, type="response", newdata = watshe)
watshe$reghat <- reghat

## extracting model estimates
levels(WA_TSHEinfo$Band) <- list( M="M", H="H")
WA_TSHEinfo$BA_same_allsc <- as.numeric(scale(WA_TSHEinfo$BA_same_all))
regmod <- glm(rtot~Band + BA_same_allsc + in5_totsc, data=WA_TSHEinfo, family = poisson)
predregmean <- predict(regmod, newdata = list(in5_totsc= c(0,0), BA_same_allsc = rep(0, times = 2), Band = c("M","H")), se.fit=TRUE, type="response")
predregci <- exp(confint(update(regmod, .~.-1))[1:2,])
#cs.reg <- coef(summary(regmod))[3:4,1:2] # no intercept, so could either fit an intercept model (like I did for BA), or make them all identical
#compsens.reg <- c(cs.reg[1,1], cs.reg[1,1] + cs.reg[2,1] )#, cs.reg[1,1] + cs.reg[3,1])
predsWA_TSHE <- data.frame(reg.mod.mean=c(NA,predregmean$fit), reg.mod.lci = c(NA,predregci[,1]), reg.mod.uci = c(NA,predregci[,2]), reg.comspens = c(NA,NA,NA), reg.compsens.se = c(NA,NA,NA))



#### WA_PSME
### simplifying it to just model total regen
WA_PSMEinfo$rtot <- with(WA_PSMEinfo, PSME_r1 + PSME_r2 + PSME_r3 + PSME_r4)
levels(WA_PSMEinfo$Band) <- list(M="M",L="L", H="H")
## NOTE: point 45 is a major outlier that results in an unrealistic parameter estimates.
# I chose to remove it
regtotWAPSME0 <- glm(rtot~Band , data=WA_PSMEinfo[-48,], family = poisson)
regtotWAPSME1 <- glm(rtot~Band + BA_tot, data=WA_PSMEinfo[-48,], family = poisson)
regtotWAPSME1.1 <- glm(rtot~Band * BA_tot, data=WA_PSMEinfo[-48,], family = poisson)
regtotWAPSME2 <- glm(rtot~Band + in5_tot, data=WA_PSMEinfo[-48,], family = poisson)
regtotWAPSME2.1 <- glm(rtot~Band * in5_tot, data=WA_PSMEinfo[-48,], family = poisson)
regtotWAPSME3 <- glm(rtot~Band + BA_same_all, data=WA_PSMEinfo[-48,], family = poisson)
regtotWAPSME3.1 <- glm(rtot~Band * BA_same_all, data=WA_PSMEinfo[-48,], family = poisson)
#BIC(regtotWAPSME0, regtotWAPSME1,regtotWAPSME1.1, regtotWAPSME2, regtotWAPSME2.1, regtotWAPSME3, regtotWAPSME3.1)
AIC(regtotWAPSME0, regtotWAPSME1,regtotWAPSME1.1, regtotWAPSME2, regtotWAPSME2.1, regtotWAPSME3, regtotWAPSME3.1)

WA_PSMEinfo$reghat <- NA
WA_PSMEinfo$reghat[-48]<- predict(regtotWAPSME3.1, type="response")
# 48 is a MAJOR outlier

levels(WA_PSMEinfo$Band) <- list(L="L", M="M", H="H")
WA_PSMEinfo$BA_same_allsc <- as.numeric(scale(WA_PSMEinfo$BA_same_all)) # create a mean center and sd scaled predictor variable
regmod <- glm(rtot~Band * BA_same_allsc, data=WA_PSMEinfo, family = poisson)
predregmean <- predict(regmod, newdata = list(BA_same_allsc = rep(0, times = 3), Band = c("L","M","H")), se.fit=TRUE, type="response")
# best model via BIC (less complex). Results qualitatively identical, except for high comp M elev.
# WA_PSMEinfo$BA_totsc <- as.numeric(scale(WA_PSMEinfo$BA_tot)) # create a mean center and sd scaled predictor variable
# regmod <- glm(rtot~Band * BA_totsc, data=WA_PSMEinfo[-48], family = poisson)
# predregmean <- predict(regmod, newdata = list(BA_totsc = rep(0, times = 3), Band = c("L","M","H")), se.fit=TRUE, type="response")

predregci <- exp(confint(update(regmod, .~.-1))[1:3,])
cs.reg <- coef(summary(regmod))[4:6,1:2] # no intercept, so could either fit an intercept model (like I did for BA), or make them all identical
compsens.reg <- c(cs.reg[1,1], cs.reg[1,1] + cs.reg[2,1], cs.reg[1,1] + cs.reg[3,1])
predsWA_PSME <- data.frame(reg.mod.mean=predregmean$fit, reg.mod.lci = predregci[,1], reg.mod.uci = predregci[,2], reg.comspens = compsens.reg, reg.compsens.se = cs.reg[,2])



#### WA_ABLA
### simplifying it to just model total regen
WA_ABLAinfo$rtot <- with(WA_ABLAinfo, ABLA_r1 + ABLA_r2 + ABLA_r3 + ABLA_r4)
levels(WA_ABLAinfo$Band) <- list(M="M",L="L", H="H")
#WA_ABLAinfo.short <- WA_ABLAinfo[-which(is.na(WA_ABLAinfo$rtot)),]
regtotWAABLA0 <- glm(rtot~Band , data=WA_ABLAinfo, family = poisson)
regtotWAABLA1 <- glm(rtot~Band + BA_tot, data=WA_ABLAinfo, family = poisson)
regtotWAABLA1.1 <- glm(rtot~Band * BA_tot, data=WA_ABLAinfo, family = poisson)
regtotWAABLA2 <- glm(rtot~Band + in5_tot, data=WA_ABLAinfo, family = poisson)
regtotWAABLA2.1 <- glm(rtot~Band * in5_tot, data=WA_ABLAinfo, family = poisson)
regtotWAABLA3 <- glm(rtot~Band + BA_same_all, data=WA_ABLAinfo, family = poisson)
regtotWAABLA3.1 <- glm(rtot~Band * BA_same_all, data=WA_ABLAinfo, family = poisson)
#BIC(regtotWAABLA0, regtotWAABLA1,regtotWAABLA1.1, regtotWAABLA2, regtotWAABLA2.1, regtotWAABLA3, regtotWAABLA3.1)
AIC(regtotWAABLA0, regtotWAABLA1,regtotWAABLA1.1, regtotWAABLA2, regtotWAABLA2.1, regtotWAABLA3, regtotWAABLA3.1)
#WA_ABLAinfo$reghat <- NA
WA_ABLAregen <- WA_ABLAinfo[-which(is.na(WA_ABLAinfo$rtot)),]
WA_ABLAregen$reghat<- predict(regtotWAABLA3.1, type="response")
# 19 might be a high leverage point


levels(WA_ABLAregen$Band) <- list(L="L", M="M", H="H")
WA_ABLAinfo$BA_same_allsc <- as.numeric(scale(WA_ABLAinfo$BA_same_all))
regmod <- glm(rtot~Band * BA_same_allsc, data=WA_ABLAinfo, family = poisson)
predregmean <- predict(regmod, newdata = list(BA_same_allsc = rep(0, times = 3), Band = c("L","M","H")), se.fit=TRUE, type="response")
predregci <- exp(confint(update(regmod, .~.-1))[1:3,])
cs.reg <- coef(summary(regmod))[4:6,1:2] # no intercept, so could either fit an intercept model (like I did for BA), or make them all identical
compsens.reg <- c(cs.reg[1,1], cs.reg[1,1] + cs.reg[2,1], cs.reg[1,1] + cs.reg[3,1])
predsWA_ABLA <- data.frame(reg.mod.mean=predregmean$fit, reg.mod.lci = predregci[,1], reg.mod.uci = predregci[,2], reg.comspens = compsens.reg, reg.compsens.se = cs.reg[,2])




####### combining all Regen model output #####
regen <- rbind(predsCO_PIPO, predsCO_POTR, predsCO_ABLA, predsMT_TSHE, predsMT_PSME, predsMT_ABLA, predsWA_TSHE, predsWA_PSME, predsWA_ABLA)





######### Figure S7: Plotting Regen as f(Elev, Comp) #########

## boxplot of CO-PIPO regeneration densities from regeneration transects
pCOPIPO <- ggplot(PIPOtotregen, aes(y=RegenAll,x=Elev, colour=Elev)) + 
  geom_boxplot() +
  ggtitle("CO-PIPO")

pCOPOTR<- ggplot(CO_POTRinfo, aes(x = in5_tot, y = reghat, colour = Band)) +
  geom_point(aes(y = rtot), alpha=.5) +
  geom_line(size = 1) +
  labs( y = "Expected Regen") +
  ggtitle("CO-POTR")

pCOABLA <- ggplot(CO_ABLAreg, aes(x = in5_tot, y = reghat, colour = Band)) +
  geom_point(aes(y = rtot), alpha=.5) +
  geom_line(size = 1) +
  labs( y = "Expected Regen") +
  ggtitle("CO-ABLA")

pMTTSHE<- ggplot(MT_TSHEinfo, aes(x = BA_same_all, y = reghat, colour = Band)) +
  geom_point(aes(y = rtot), alpha=.5) +
  geom_line(size = 1) +
  labs( y = "Expected Regen") +
  ggtitle("MT-TSHE")

pMTPSME<- ggplot(MT_PSMEinfo, aes(x = BA_same_all, y = reghat, colour = Band)) +
  geom_point(aes(y = rtot), alpha=.5) +
  geom_line(size = 1) +
  labs( y = "Expected Regen") +
  ggtitle("MT-PSME")

pMTABLA <- ggplot(MT_ABLAinfo, aes(x = in5_tot, y = reghat, colour = Band)) +
  geom_point(aes(y = rtot), alpha=.5) +
  geom_line(size = 1) +
  labs( y = "Expected Regen") +
  ggtitle("MT-ABLA")

pWATSHE<- ggplot(WA_TSHEinfo, aes(x = in5_tot, y = reghat, colour = Band)) +
  geom_point(aes(y = rtot), alpha=.5) +
  geom_line(data=watshe, size = 1) +
  labs( y = "Expected Regen") +
  ggtitle("WA-TSHE")

pWAPSME<- ggplot(WA_PSMEinfo[-48,], aes(x = BA_same_all, y = reghat, colour = Band)) +
  geom_point(aes(y = rtot), alpha=.5) +
  geom_line(size = 1) +
  labs( y = "Expected Regen") +
  ggtitle("WA-PSME")

pWAABLA <- ggplot(WA_ABLAregen, aes(x = BA_same_all, y = reghat, colour = Band)) +
  geom_point(aes(y = rtot), alpha=.5) +
  geom_line(size = 1) +
  labs( y = "Expected Regen") +
  ggtitle("WA-ABLA")


quartz(width=6, height=5)
multiplot(plotlist =list(pCOPIPO,pMTTSHE,pWATSHE, pCOPOTR, pMTPSME, pWAPSME, pCOABLA, pMTABLA,pWAABLA), cols=3) #list(pCOPIPO, pCOPOTR, pCOABLA,pMTTSHE,pMTPSME, pMTABLA,pWATSHE,pWAPSME,pWAABLA), cols=3)


####### END: Regeneration Density Modeling ############################
#_____________________________________________________________________________________________














#_____________________________________________________________________________________________
#_____________________________________________________________________________________________
########## BEGIN: Large-scale Trade-off Analysis ###########
#_____________________________________________________________________________________________
#_____________________________________________________________________________________________

# This code takes the compiled output of the Basal Area Growth analysis 
# and Growth Synchrony Analysis for all species replicates, performs statistical tests
# for large-scale relationships between climatic and competitive growth constraints
# and mean growth rate, and plots the results.
# The Supplemental Dataset S3: Growth_Summaries_20180428.csv provides the full output derived from performing the example analyses
# above for all species-replicates.

Sumstatsall <- read.csv("SIData_S3_All_Summaries_20180428.csv")

#Sumstatsall <- read.csv("/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/TreeCoresRanges_LDLA/Manuscript Prep/Example_Data_forMS/SIData_S3_All_Summaries_20180428.csv")

######## Mean Growth v Synchrony Stats ##################
# since the SEs of synchrony estimates aren't symmetrical 
# (because of backtransforming from a beta regression)
# I'll use an inverse variance-weighting based on the lower SE only.

## raw BAI and bai reduction
# invCorrSD <- 1/Sumstatsall$SDcorr^2 / sum(1/Sumstatsall$SDcorr^2, na.rm=T) # inverse Variance weighting
invCorr.se <- 1/(Sumstatsall$corr.mod.mean-Sumstatsall$corr.mod.lse) / sum(1/(Sumstatsall$corr.mod.mean-Sumstatsall$corr.mod.lse), na.rm=T) # inverse SE weighting


baisynchglobal.int <- lmer(corr.mod.mean~scale(BAI.mod.mean.int) + (scale(BAI.mod.mean.int)|St_Sp), Sumstatsall)
baisynchglobalweighted.int <- lmer(corr.mod.mean~scale(BAI.mod.mean.int) + (scale(BAI.mod.mean.int)|St_Sp), Sumstatsall, weights=invCorr.se)
#baisynchglobalweighted.int <- lmer(corr.mod.mean~scale(BAI.mod.mean.int) + (scale(BAI.mod.mean.int)|St_Sp), Sumstatsall, weights=invCorrSD)

# examine residuals
qqp(resid(baisynchglobalweighted.int))
rfs <- ranef(baisynchglobalweighted.int)[[1]]
qqp(rfs[,2])

summary(baisynchglobal.int) # p = 0.066
summary(baisynchglobalweighted.int) # p=0.0446




######## Mean Growth v CompSens Stats ##################
## RAW BAI and RAW comp sens (bai reduction)
# create an inverse se metric with which to weight mixed effects model 
invCompSens.se.int <- 1/(Sumstatsall$CompSens.se.int) / sum(1/Sumstatsall$CompSens.se.int, na.rm=T)

baicompglobal.int <- lmer(CompSens.mean.int~BAI.mod.mean.sc.int + (BAI.mod.mean.sc.int|St_Sp), Sumstatsall)
baicompglobalweighted.int <- lmer(CompSens.mean.int~BAI.mod.mean.sc.int + (BAI.mod.mean.sc.int|St_Sp), Sumstatsall, weights = invCompSens.se.int)
#lmer(CompSens.mean.int~scale(BAI.mod.mean.int) + (scale(BAI.mod.mean.int)|St_Sp), Sumstatsall, weights= invCompSens.se.int)

# model criticism plots
qqp(resid(baicompglobal.int)) # innermost residuals look normal?
rfs <- ranef(baicompglobal.int)[[1]]
qqp(rfs[,2]) # random effects look normal?


summary(baicompglobal.int) # p=0.000434
summary(baicompglobalweighted.int) #p=0.00138 




#### Mean growth vs CompSens, everything as a proportion of max BAI ###
baicompglobal.prop.int <- lmer(CompSens.prop.int~BAI.prop.int + (BAI.prop.int|St_Sp), Sumstatsall)
invCompSens.se <- 1/(Sumstatsall$CompSens.se.int) / sum(1/Sumstatsall$CompSens.se.int, na.rm=T)
baicompglobal.propweighted.int <- lmer(CompSens.prop.int~BAI.prop.int + (BAI.prop.int|St_Sp), Sumstatsall , weights=invCompSens.se)

qqp(resid(baicompglobal.propweighted.int))
rfs <- ranef(baicompglobal.propweighted.int)[[1]]
qqp(rfs[,1]) # not the world's most normal

summary(baicompglobal.prop.int) # p =0.155
summary(baicompglobal.propweighted.int) # p=0.351




#======================================================
######## ** FIGURE 3: Large Scale Trade-off ################################
#======================================================
#. light St_Sp trend lines, emphasize points

# set color palet to be paired light and dark blue and green
prescols <- brewer.pal(n=4, "Paired")[c(2,4,1,3)]

# set lwd, lt for species-replicate lines
stsp_lwd <- 1
stsp_lty <- 2

quartz(width=6, height=4)
par(mar=c(4,4.4,2,0.2), oma=c(0,0,0,1), mfrow=c(1,2), mgp=c(2.5,.8,0))
#palette(paste0(colors.ord,"AA"))

plot(corr.mod.mean~BAI.mod.mean.int,Sumstatsall, col=prescols[1], pch=16, cex=1.3,
     ylab="Growth Synchrony",
     xlab=expression(paste("Mean Growth ("*mm^2/yr*")")),
     ylim=c(0, .9), xlim=c(200,2800), cex.lab=1.2, type="n")
# error_bars(xvar="BAI.mod.mean.int", yvar = "corr.mod.mean", upper = "SDcorr", errordata = Sumstatsall,length = 0, color = prescols[3])
#arrows(x0=Sumstatsall$BAI.mod.mean.int, y0=Sumstatsall$corr.mod.mean, y1=Sumstatsall$corr.mod.ci97.5,length=0, lwd=2, col=prescols[3] )
#arrows(x0=Sumstatsall$BAI.mod.mean.int, y0=Sumstatsall$corr.mod.mean, y1=Sumstatsall$corr.mod.ci2.5,length=0, lwd=2, col=prescols[3] )
 arrows(x0=Sumstatsall$BAI.mod.mean.int, y0=Sumstatsall$corr.mod.mean, y1=Sumstatsall$corr.mod.use,length=0, lwd=2, col=prescols[3] )
 arrows(x0=Sumstatsall$BAI.mod.mean.int, y0=Sumstatsall$corr.mod.mean, y1=Sumstatsall$corr.mod.lse,length=0, lwd=2, col=prescols[3] )



points(corr.mod.mean~BAI.mod.mean.int,Sumstatsall, col=prescols[1], pch=16, cex=1.3)
# St_Sp trend lines
for(i in unique(Sumstatsall$St_Sp)){
  d <- Sumstatsall[which(Sumstatsall$St_Sp==i),]
  lines(fitted(lm(corr.mod.mean~BAI.mod.mean.int, d))~na.omit(d$BAI.mod.mean.int), col=prescols[1],lty=stsp_lty,lwd=stsp_lwd) #lty=2)
}
## global trend line
bc <- 1298.616
bs <- 595.8935#615.95
bcoef <- coef(summary(baisynchglobal.int))[,1]
bcoefa <- (bcoef[1]*bs - bcoef[2]*bc)/bs
bcoefb <- bcoef[2]/bs
abline(b=bcoefb, a=bcoefa, lwd=2)
mtext(text="p=0.045", side = 3, line = -1, adj=0.1)

#legend('bottomright', legend=levels(Sumstatsall$State), pch=c(16,17,18), ncol=1, bty="n", cex=.8)

plot(CompSens.mean.int.rev~BAI.mod.mean.int,Sumstatsall, col=prescols[2], pch=16, cex=1.3,
     ylab="Sens to Comp (growth reduction)", #\n(mm2 BAI reduction)",
     xlab=expression(paste("Mean Growth ("*mm^2/yr*")")),
     xlim=c(200,2800), cex.lab=1.2 ,ylim=c(-350, 750))
error_bars(xvar="BAI.mod.mean.int", yvar = "CompSens.mean.int.rev", upper = "CompSens.se.int", errordata = Sumstatsall,length = 0, color = prescols[4])
points(CompSens.mean.int.rev~BAI.mod.mean.int,Sumstatsall, col=prescols[2], pch=16, cex=1.3)
# individual St_Sp trendlines
for(i in unique(Sumstatsall$St_Sp)){
  d <- Sumstatsall[which(Sumstatsall$St_Sp==i),]
  lines(fitted(lm(CompSens.mean.int.rev~BAI.mod.mean.int, d))~na.omit(d$BAI.mod.mean.int), col=prescols[2], lty=stsp_lty, lwd=stsp_lwd)
}
# adding global trend line
bc <- 1298.616 # BAI centering component
bs <- 595.8935#615.95 # BAI scaling component
baicompglobal.rev <- lmer(CompSens.mean.int.rev~BAI.mod.mean.sc.int + (BAI.mod.mean.sc.int|St_Sp), Sumstatsall)
bcoef <- coef(summary(baicompglobal.rev))[,1]
bcoefa <- (bcoef[1]*bs - bcoef[2]*bc)/bs
bcoefb <- bcoef[2]/bs
abline(b=bcoefb, a=bcoefa, lwd=2)
# abline(lm(CompSens.mean~BAI.mod.mean,Sumstatsall), lwd=2, col="blue")
mtext(text="p=0.001", side = 3, line=-1, adj=.1)
# legend in bottom-right corner
#legend("bottomright", legend=levels(Sumstatsall$Species), col=colors.ord, pch=16, xpd=T, cex=.8, ncol=2, bty="n")




####### END: Large Scale Trade-off analysis ############################
#_____________________________________________________________________________________________





#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
####### ** Fig 4: Survival and Regeneration #################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# new horizontal for 1 col width figure (cribbed from Fig 4)
quartz(width=6, height=3)
par(mar=c(3,4.4,1,0.2), oma=c(0,0,0,1), mfrow=c(1,2), mgp=c(2.5,.8,0))


################## A) Modeled Survival vs Range position
regenhyps <- c(1,1,1,2,2,2,4,3,3) # 1 = conforms w. growth, 2= doesn't, 3= partial, 4= no pattern
survhyps <- c(3,1,1,4,3,2,4,2,1)
colsnegative <- brewer.pal(n=11, "RdGy")[c(8,3,8)] # colors to highlight the opposite of expected patterns
colsnegative <- c("#AAAAAA", "#E41A1C", "#984EA3", "#AAAAAA") # pink = "#D6705D"
plot(Survival.mod.mean~xval, Sumstatsall, type="n", ylab="Survival Rate", xlim=c(0,3.5),xaxt="n", xlab="", ylim=c(0.6,1.01)) # previous xlim=c(0,3.33)
axis(1, at=Sumstatsall$xval[1:3],labels=c("Harsh","", "Benign"))# labels=c("","", ""))
mtext("a)", side=3, adj=0.05,line=-1)
for (j in 1:9){
  i <- unique(Sumstatsall$St_Sp)[j]
  if(Sumstatsall$BAI.prop.int[which(Sumstatsall$St_Sp==i)][3]==min(Sumstatsall$BAI.prop.int[which(Sumstatsall$St_Sp==i)], na.rm=T)){
    points(Survival.mod.mean~I(xval*-1 +.6666 + 2.6666), Sumstatsall, subset=St_Sp==i, type="l", pch=16, col=colsnegative[survhyps[j]], lwd=2
    ) #, col=as.numeric(Sumstatsall$State[which(Sumstatsall$St_Sp==i)])[1]) # , col=synch.agree[which(Sumstatsall$St_Sp==i)], pch=synch.sig[which(Sumstatsall$St_Sp==i)],
  }
  else{
    points(Survival.mod.mean~xval, Sumstatsall, subset=St_Sp==i, type="l", pch=16, col=colsnegative[survhyps[j]], lwd=2
    )#, col=as.numeric(Sumstatsall$State[which(Sumstatsall$St_Sp==i)])[1]) #col=synch.agree[which(Sumstatsall$St_Sp==i)], pch=synch.sig[which(Sumstatsall$St_Sp==i)]
    
  }
}
lines(x = Sumstatsall$xval[c(1,3)], y=c(.77,1), lwd=3)
legend('bottomright', legend = c("Prediction","Follows/ns","Partially","Opposite"), col = c("black",colsnegative[c(1,3,2)]), lwd=c(3,2,2,2), bty="n", cex=.8)
text(x = 3.15, y=Sumstatsall$Survival.mod.mean[c(16,24)]+.02,labels = c("MT-ABLA","WA-PSME"), cex=.7, srt=30)


######################## B) Modeled Regen vs range position
plot(reg.mod.prop~xval, Sumstatsall, type="n", ylab="Prop. max Regeneration", xlim=c(0,3.5),xaxt="n", xlab="", ylim=c(0,1))
axis(1, at=Sumstatsall$xval[1:3], labels=c("Harsh","", "Benign"))
mtext("b)", side=3, adj=0.05,line=-1)
for (j in 1:9){
  i <- unique(Sumstatsall$St_Sp)[j]
  if(Sumstatsall$BAI.prop.int[which(Sumstatsall$St_Sp==i)][3]==min(Sumstatsall$BAI.prop.int[which(Sumstatsall$St_Sp==i)], na.rm=T)){
    points(reg.mod.prop~I(xval*-1 +.6666 + 2.6666), Sumstatsall, subset=St_Sp==i, type="l", pch=16, col=colsnegative[regenhyps[j]],lwd=2
    )#, col=as.numeric(Sumstatsall$State[which(Sumstatsall$St_Sp==i)])[1])# , col=synch.agree[which(Sumstatsall$St_Sp==i)], pch=synch.sig[which(Sumstatsall$St_Sp==i)],
  }
  else{
    points(reg.mod.prop~xval, Sumstatsall, subset=St_Sp==i, type="l", pch=16, col=colsnegative[regenhyps[j]], lwd=2
    )#, col=as.numeric(Sumstatsall$State[which(Sumstatsall$St_Sp==i)])[1]) #col=synch.agree[which(Sumstatsall$St_Sp==i)], pch=synch.sig[which(Sumstatsall$St_Sp==i)]
    
  }
}

lines(x = Sumstatsall$xval[c(1,3)], y=c(.3,1), lwd=3)
text(x = 3.15, y=Sumstatsall$reg.mod.prop[c(12,13,16)]+.07,labels = c("MT-TSHE","MT-PSME","MT-ABLA"), cex=.7, srt=30)






