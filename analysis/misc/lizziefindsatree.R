## Started (finally) 25 September 2026 ##
## Lizzie made the file by copying from ...
## https://github.com/vvandermeersch/climategrowthshifts/issues/160#issuecomment-5528581597 ##

## Looking for a long-term tree to try out a GP ##

# housekeeping
rm(list=ls())
options(stringsAsFactors = FALSE)

set.seed(11)

wd <- "/Users/lizzie/Documents/git/projects/grephon/climategrowthshifts/analysis"

library(dplR)

datasets_summary <- readRDS(file.path(wd, "/pnwandmore/input/itrdb", "datasets_summary_usonly.rds"))
raw_dir <- paste0(wd, "/pnwandmore/data/itrdb/usa")

# selection <- datasets_summary[datasets_summary$species_code == 'PSME' & datasets_summary$first_year < 1200, 'dataset']
selection <- datasets_summary[datasets_summary$species_code == 'PIED' & datasets_summary$first_year < 1200, 'dataset']

random_dataset <- sample(na.omit(selection), 1)

f <- paste0(random_dataset, '.rwl')
rwdat <- dplR::read.rwl(file.path(raw_dir, f), format = 'tucson', verbose = FALSE)

idxs <- dplR::autoread.ids(rwdat, fix.typos = FALSE)
idxs$original <- colnames(rwdat)[1:nrow(idxs)] 
colnames(rwdat)[1:nrow(idxs)] <- paste0('tree', idxs$tree, '_core', idxs$core)

rwdat$year <- rownames(rwdat)

rwdat_long <- reshape(rwdat, varying = names(rwdat)[names(rwdat) != "year"],
                      v.names = "rw_mm", timevar = "full_id",
                      times = names(rwdat)[names(rwdat) != "year"],
                      direction = "long")
rwdat_long$tree_id <- sub(".*tree(\\d+)_core.*", "\\1", rwdat_long$full_id)
rwdat_long$core_id <- sub(".*core(\\d+).*", "\\1", rwdat_long$full_id)

years_bytree <- aggregate(year ~ tree_id, data = rwdat_long[!is.na(rwdat_long$rw_mm), ], 
          function(x) c(min = min(x), max = max(x)))

years_bytree$yearn <- as.numeric(years_bytree$year[,2])-as.numeric(years_bytree$year[,1])

years_bytree[with(years_bytree, order(-yearn)), ]

tree30 <- rwdat_long[rwdat_long$tree_id == "30",] 
tree29 <- rwdat_long[rwdat_long$tree_id == "29",] 

plot(rw_mm~year, tree30)
plot(rw_mm~year, tree29)
write.csv(tree30, paste0(wd, "/misc/output/piedtree30.csv"), row.names=FALSE)
