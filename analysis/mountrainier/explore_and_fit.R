
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

wd <- "~/projects/climategrowthshifts/analysis/mountrainier"

# read and prepare data
tv01002 <- read.csv(file.path(wd, "data/trees/TV01002_v17.csv"))
mtr_data <- tv01002[tv01002$PSP_STUDYID == "MRRS",]
length(unique(mtr_data$TREEID))

mtr_data <- na.omit(mtr_data[c("STANDID", "TREEID", "DBH", "YEAR")]) # remove NA (death)
mtr_data <- mtr_data[ave(seq_len(nrow(mtr_data)), mtr_data$TREEID, FUN = length) >= 4, ] # keep only if at least 4 obs.

# ------------------------------------ #
# First test with M. Betancourt subset #
# ------------------------------------ #
mbetsubids <-  c("AV06000100001", "AV06000100002", "AV06000100003", "AV06000100004", "AV06000100005", "AV06000100006", "AV06000100007", "AV06000100008", "AV06000100009", "AV06000100010", "AV06000100011",
             "AV06000100012", "AV06000100013", "AV06000100014", "AV06000100015", "AV06000100016", "AV06000100017", "AV06000100018", "AV06000100019", "AV06000100020", "AV06000100021", "AV06000100022")
mtr_data_mbet <- mtr_data[mtr_data$TREEID %in% mbetsubids,]
# mtr_data <- mtr_data[mtr_data$TREEID %in% sample(mtr_data$TREEID, 100),]

mtr_data_mbet$NUMID <- as.integer(factor(mtr_data_mbet$TREEID))
mtr_data_mbet$NUMOBS <- ave(seq_along(mtr_data_mbet$TREEID), mtr_data_mbet$TREEID, FUN = seq_along)
Ntrees <- max(mtr_data_mbet$NUMID)
Nobs_pertree <- aggregate(NUMOBS ~ NUMID, mtr_data_mbet, max)$NUMOBS
Sobs <- c(0,cumsum(Nobs_pertree))

gompertz <- stan_model(file.path(wd, "stan/heterogeneous_gompertz.stan"))
 
mdl.data <- list(N = nrow(mtr_data_mbet),
                 N_trees =  max(mtr_data_mbet$NUMID),
                 Sobs = Sobs,
                 years = mtr_data_mbet$YEAR,
                 dbhs = mtr_data_mbet$DBH)

fit <- sampling(gompertz, mdl.data, 
                chains = 1, iter = 4000,
                verbose = TRUE)
# it's working!

# -------------------------------------- #
# Second test with random sampled subset #
# -------------------------------------- #
mtr_data_samp <- mtr_data[mtr_data$TREEID %in% sample(mtr_data$TREEID, 100),]

mtr_data_samp$NUMID <- as.integer(factor(mtr_data_samp$TREEID))
mtr_data_samp$NUMOBS <- ave(seq_along(mtr_data_samp$TREEID), mtr_data_samp$TREEID, FUN = seq_along)
Ntrees <- max(mtr_data_samp$NUMID)
Nobs_pertree <- aggregate(NUMOBS ~ NUMID, mtr_data_samp, max)$NUMOBS
Sobs <- c(0,cumsum(Nobs_pertree))

gompertz <- stan_model(file.path(wd, "stan/heterogeneous_gompertz.stan"))

mdl.data <- list(N = nrow(mtr_data_samp),
                 N_trees =  max(mtr_data_samp$NUMID),
                 Sobs = Sobs,
                 years = mtr_data_samp$YEAR,
                 dbhs = mtr_data_samp$DBH)

fit <- sampling(gompertz, mdl.data, 
                chains = 1, iter = 4000,
                verbose = TRUE)
# not working:
# Chain 1: Initialization between (-2, 2) failed after 100 attempts. 
# Chain 1:  Try specifying initial values, reducing ranges of constrained values, or reparameterizing the model.
# [1] "Error : Initialization failed."




