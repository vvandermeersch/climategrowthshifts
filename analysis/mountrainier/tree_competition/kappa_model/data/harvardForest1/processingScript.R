# started on Aug 23, 2026 by Devina
library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(tidyverse)
library(dbscan)

setwd("D:/ubc_study/udergrad_research/temporal ecology lab/climategrowthshifts/analysis/mountrainier/tree_competition/kappa_model/data/harvardForest1")

#alive stems-> treeid, stemid match -> dbh difference/2
# tree2014 <- read_csv("hf253-05-stems-2014.csv",show_col_types = FALSE)
# tree2019 <- read_csv("hf253-06-stems-2019.csv",show_col_types = FALSE)

#treeid, stemid, sp, gx, gy, dbh
# alive status
stemdataName <- list("hf253-05-stems-2014.csv","hf253-06-stems-2019.csv")
stem_list <- map(stemdataName,function(f){
  dat <- read_csv(f,show_col_types = FALSE)
  dat%>%
    filter(status=='A')%>%
    select(tree.id,stem.id,sp,gx,gy,dbh)
})

growthchange <- stem_list[[1]]%>%
  inner_join(stem_list[[2]],
             by = c("tree.id","stem.id"),
             suffix = c("_2014","_2019"))%>%
  filter(sp_2014==sp_2019)%>%
  mutate(growth = (dbh_2019-dbh_2014)/2,
         locationchange = sqrt((gx_2014-gx_2019)^2+(gy_2014-gy_2019)^2))%>%
  filter(locationchange==0)%>%
  select(tree.id,stem.id,sp_2014,gx_2014,gy_2014,growth)%>%
  rename(sp=sp_2014,gx=gx_2014,gy=gy_2014)

# neighborhood_2014 <- stem_list[[1]] %>%
#   select(focal_treeid = tree.id,
#          focal_stemid = stem.id,
#          focal_x = gx,
#          focal_y = gy) %>%
#   crossing(
#     stem_list[[1]] %>%
#       select(neighbor_treeid = tree.id,
#              neighbor_stemid = stem.id,
#              neighbor_sp = sp,
#              neighbor_x = gx,
#              neighbor_y = gy,
#              neighbor_dbh = dbh)
#   ) %>%
#   mutate(dist = sqrt((neighbor_x-focal_x)^2+(neighbor_y-focal_y)^2))%>%
#   filter(dist > 0 & dist <= 10)
neighborhood_list <- map(stem_list,function(d) {
  dat <- d %>% filter(!is.na(gx),!is.na(gy))
  coords <- as.matrix(dat[,c("gx","gy")])
  neighbors <- frNN(coords,eps=10,sort=TRUE)
  
  neighborhood <- map_dfr(seq_len(nrow(dat)),function(i){
    neighbor_idx <- neighbors$id[[i]]
    neighbor_idx <- neighbor_idx[neighbor_idx!=i]
    
    if (length(neighbor_idx)==0) {
      return (NULL)
    }
    
    tibble(
      focal_treeid = dat$tree.id[i],
      focal_stemid = dat$stem.id[i],
      
      neighbor_treeid = dat$tree.id[neighbor_idx],
      neighbor_stemid = dat$stem.id[neighbor_idx],
      neighbor_sp = dat$sp[neighbor_idx],
      neighbor_dbh = dat$dbh[neighbor_idx],
      
      dist = sqrt((dat$gx[neighbor_idx]-dat$gx[i])^2+
                    (dat$gy[neighbor_idx]-dat$gy[i])^2)
    )
    }) %>%
    filter(focal_treeid!=neighbor_treeid)
  
  return(neighborhood)
})

write.csv(neighborhood_list[[1]],"processed/neighborhood_2014.csv")
write.csv(neighborhood_list[[2]],"processed/neighborhood_2019.csv")
write.csv(growthchange,"processed/tree_growth.csv")
