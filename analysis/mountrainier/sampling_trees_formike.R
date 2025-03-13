
rm(list=ls())
wd <- "~/projects/climategrowthshifts/analysis/mountrainier"

# read Ailene data
treerings <- read.csv(file.path(wd, "data", "treerings_ailene", "SouthSidecoresX.csv"), header=TRUE)
treerings$uniqueID <- paste(treerings$Tag, treerings$Core, sep = "_")                      
meas_cols <-  names(treerings)[grep("X", names(treerings))]
treerings <- reshape(treerings, dir = "long", 
                     varying = meas_cols, 
                     times = as.numeric(gsub("X", "", meas_cols)),
                     v.names = "rw", timevar = "year")
treerings <- treerings[c("Stand", "Species", "Tag", "Core", "uniqueID", "DBH", "RingCount", "year", "rw")]
treerings <- na.omit(treerings)
treerings <- treerings[order(treerings$year),]

# sample one tree
treerings <- treerings[treerings$Stand == "AV06",] # same plot as in the workshop
sampled_tree <- treerings[treerings$Tag == 'E60', ] # nice looking tree, starting in 1791 (French first constitution)

sampled_core1 <- sampled_tree[sampled_tree$Core == 1, c('year', 'rw')]
sampled_core2 <- sampled_tree[sampled_tree$Core == 2, c('year', 'rw')]

sampled_cores <- merge(sampled_core1, sampled_core2, , by = 'year')
names(sampled_cores) <- c('year', 'rw_core1 (mm)', 'rw_core2 (mm)')

write.csv(sampled_cores, file = file.path(wd, 'input/treerings/tag_E60.csv'))

# sample trees from the same species
treerings <- treerings[treerings$Stand == "AV06",] # same plot as in the workshop
sampled_tree <- treerings[treerings$Species == 'Abam', ]

sampled_core1 <- sampled_tree[sampled_tree$Core == 1, c('Tag', 'year', 'rw')]
sampled_core2 <- sampled_tree[sampled_tree$Core == 2, c('Tag', 'year', 'rw')]

sampled_cores <- merge(sampled_core1, sampled_core2, , by = c('Tag', 'year'), all.x = TRUE)
names(sampled_cores) <- c('tree_id', 'year', 'rw_core1 (mm)', 'rw_core2 (mm)')
length(unique(sampled_cores$tree_id))

write.csv(sampled_cores, file = file.path(wd, 'input/treerings/abam_AV06.csv'))
