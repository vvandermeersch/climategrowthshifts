
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

# sample one tree, one plot
treerings <- treerings[treerings$Stand == "AV06",] # same plot as in the workshop
sampled_tree <- treerings[treerings$Tag == 'E60', ] # nice looking tree, starting in 1791 (French first constitution)

sampled_core1 <- sampled_tree[sampled_tree$Core == 1, c('year', 'rw')]
sampled_core2 <- sampled_tree[sampled_tree$Core == 2, c('year', 'rw')]

sampled_cores <- merge(sampled_core1, sampled_core2, , by = 'year')
names(sampled_cores) <- c('year', 'rw_core1 (mm)', 'rw_core2 (mm)')

write.csv(sampled_cores, file = file.path(wd, 'input/treerings/tag_E60.csv'))

# sample trees from the same species, one plot
treerings <- treerings[treerings$Stand == "AV06",] # same plot as in the workshop
sampled_tree <- treerings[treerings$Species == 'Abam', ]

sampled_core1 <- sampled_tree[sampled_tree$Core == 1, c('Tag', 'year', 'rw')]
sampled_core2 <- sampled_tree[sampled_tree$Core == 2, c('Tag', 'year', 'rw')]

sampled_cores <- merge(sampled_core1, sampled_core2, , by = c('Tag', 'year'), all.x = TRUE)
names(sampled_cores) <- c('tree_id', 'year', 'rw_core1 (mm)', 'rw_core2 (mm)')
length(unique(sampled_cores$tree_id))

write.csv(sampled_cores, file = file.path(wd, 'input/treerings/abam_AV06.csv'))


# sample trees from the same species, all plots
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

sampled_tree <- treerings[treerings$Species == 'Abam', ]
sampled_tree$tree_id <- paste0(sampled_tree$Stand, '_', sampled_tree$Tag)

sampled_core1 <- sampled_tree[sampled_tree$Core == 1, c('Stand', 'tree_id', 'year', 'rw')]
sampled_core2 <- sampled_tree[sampled_tree$Core == 2, c('Stand', 'tree_id', 'year', 'rw')]

sampled_cores <- merge(sampled_core1, sampled_core2, , by = c('Stand', 'tree_id', 'year'), all = TRUE)
names(sampled_cores) <- c('stand', 'tree_id', 'year', 'rw_core1 (mm)', 'rw_core2 (mm)')
length(unique(sampled_cores$tree_id))


any(aggregate(`rw_core1 (mm)` ~ tree_id, data = sampled_cores, FUN = function(x) all(is.na(x))) == TRUE)

write.csv(sampled_cores, file = file.path(wd, 'input/treerings/abam_allplots.csv'))


# sample trees from all species, all plots
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

sampled_tree <- treerings
sampled_tree$tree_id <- paste0(sampled_tree$Stand, '_', sampled_tree$Tag)

sampled_core1 <- sampled_tree[sampled_tree$Core == 1, c('Stand', 'Species', 'tree_id', 'year', 'rw')]
sampled_core2 <- sampled_tree[sampled_tree$Core == 2, c('Stand', 'Species', 'tree_id', 'year', 'rw')]

sampled_cores <- merge(sampled_core1, sampled_core2, , by = c('Stand', 'Species', 'tree_id', 'year'), all = TRUE)
names(sampled_cores) <- c('stand', 'species', 'tree_id', 'year', 'rw_core1 (mm)', 'rw_core2 (mm)')
length(unique(sampled_cores$tree_id))

write.csv(sampled_cores, file = file.path(wd, 'input/treerings/allspecies_allplots.csv'))

span <- data.frame()
newsampledcores <- data.frame()
for(t in unique(sampled_cores$tree_id)){
  
  df <- sampled_cores[sampled_cores$tree_id == t , ]
  
  span1 = length(na.omit(df$`rw_core1 (mm)`))
  span2 = length(na.omit(df$`rw_core2 (mm)`))
  
  if(span1 < span2){
    
    save <- df$`rw_core1 (mm)`
    df$`rw_core1 (mm)` <-  df$`rw_core2 (mm)` 
    df$`rw_core2 (mm)` <-  save
    
  }
  
  newsampledcores <- rbind(newsampledcores, df)
  
  span1 = length(na.omit(df$`rw_core1 (mm)`))
  span2 = length(na.omit(df$`rw_core2 (mm)`))
  
  span <- rbind(
    span,
    data.frame(tree_id = t, species = unique(df$species), span1, span2))
  
}

newspan <- span[span$span1 > 50,]
newsampledcores <- newsampledcores[newsampledcores$tree_id %in% newspan$tree_id, ]
length(unique(newsampledcores$tree_id))

write.csv(newsampledcores, file = file.path(wd, 'input/treerings/mosttrees_allspecies_allplots.csv'))


# sample trees from the same species, one plot
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

treerings <- treerings[treerings$Stand == "PARA",] # same plot as in the workshop
sampled_tree <- treerings[treerings$Species == 'Tsme', ]

sampled_core1 <- sampled_tree[sampled_tree$Core == 1, c('Tag', 'year', 'rw')]
sampled_core2 <- sampled_tree[sampled_tree$Core == 2, c('Tag', 'year', 'rw')]

sampled_cores <- merge(sampled_core1, sampled_core2, , by = c('Tag', 'year'), all.x = TRUE)
names(sampled_cores) <- c('tree_id', 'year', 'rw_core1 (mm)', 'rw_core2 (mm)')
length(unique(sampled_cores$tree_id))

write.csv(sampled_cores, file = file.path(wd, 'input/treerings/tsme_PARA.csv'))
